package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;
import com.duckandcover.html.Transforms;
import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.swing.IValueToColor;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.FastIndex;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * The UncertainValuesBase class serves to store a collection of related
 * variable values along with the associated covariance/uncertainty matrix.
 * </p>
 * <p>
 * UncertainValuesBase instances can be the input into calculations class or can
 * be the results of calculations performed with the
 * {@link UncertaintyPropagator} class.
 * </p>
 * <p>
 * To facilitate bookkeeping, the variables within the object are labeled with
 * Object derived labels. Make sure that all the classes implementing labels
 * implement hashCode() and equals() as the labels are stored in hash maps. The
 * entries in the values vector and covariance matrix can be indexed by label or
 * by integer index.
 * </p>
 * <p>
 * Covariance matrices have certain properties (square, symmetric, positive
 * definite, and the covariance elements must be related to the variances via a
 * corrolation coefficient which can take on values between -1 and 1.) These
 * properties are checked in the constructor.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
abstract public class UncertainValuesBase<H> //
		implements IToHTML, IUncertainValues<H> {

	private static class CombinedUncertainValues<J> extends UncertainValuesBase<J> {

		private final List<Pair<UncertainValuesBase<? extends J>, Integer>> mIndices;

		private final SimplyLazy<RealVector> mValues = new SimplyLazy<RealVector>() {

			@Override
			protected RealVector initialize() {
				final RealVector res = new ArrayRealVector(mIndices.size());
				for (int i = 0; i < res.getDimension(); ++i) {
					final Pair<UncertainValuesBase<? extends J>, Integer> pr = mIndices.get(i);
					final int idx = pr.getSecond().intValue();
					res.setEntry(i, pr.getFirst().getEntry(idx));
				}
				return res;
			}
		};

		private final SimplyLazy<RealMatrix> mCovariances = new SimplyLazy<RealMatrix>() {

			@Override
			protected RealMatrix initialize() {
				final RealMatrix res = MatrixUtils.createRealMatrix(mIndices.size(), mIndices.size());
				for (int r = 0; r < mIndices.size(); ++r)
					for (int c = 0; c < mIndices.size(); ++c) {
						final Pair<UncertainValuesBase<? extends J>, Integer> rPr = mIndices.get(r);
						final Pair<UncertainValuesBase<? extends J>, Integer> cPr = mIndices.get(c);
						if (rPr.getFirst() == cPr.getFirst()) {
							final RealMatrix baseRes = rPr.getFirst().getCovariances();
							final int rIdx = rPr.getSecond().intValue();
							final int cIdx = cPr.getSecond().intValue();
							res.setEntry(r, c, baseRes.getEntry(rIdx, cIdx));
						}
					}
				return res;
			}

		};

		protected CombinedUncertainValues( //
				final List<J> labels, //
				final List<? extends UncertainValuesBase<? extends J>> bases, //
				final boolean first //
		) throws ArgumentException {
			super(labels);
			mIndices = new ArrayList<>();
			for (final J lbl : labels) {
				int cx = 0;
				for (final UncertainValuesBase<? extends J> uvb : bases) {
					final int idx = uvb.indexOf(lbl);
					if (idx >= 0) {
						mIndices.add(Pair.create(uvb, idx));
						cx++;
						if (first)
							break;
					}
				}
				if (cx == 0)
					throw new ArgumentException("The argument " + lbl + " has not been defined.");
				if (cx > 1)
					throw new ArgumentException("The argument " + lbl + " has been multiply defined.");

			}
		}

		@Override
		public RealMatrix getCovariances() {
			return mCovariances.get();
		}

		@Override
		public RealVector getValues() {
			return mValues.get();
		}

	};

	private static class ReorderedUncertainValues<L> //
			extends UncertainValuesBase<L> {

		private final UncertainValuesBase<? super L> mBase;
		private final int[] mIndexes;

		private final SimplyLazy<RealVector> mValues = new SimplyLazy<RealVector>() {

			@Override
			protected RealVector initialize() {
				final RealVector baseVals = mBase.getValues();
				final RealVector res = new ArrayRealVector(mIndexes.length);
				for (int i = 0; i < mIndexes.length; ++i)
					res.setEntry(i, baseVals.getEntry(mIndexes[i]));
				return res;
			}
		};

		private final SimplyLazy<RealMatrix> mCovariances = new SimplyLazy<RealMatrix>() {

			@Override
			protected RealMatrix initialize() {
				final RealMatrix res = MatrixUtils.createRealMatrix(mIndexes.length, mIndexes.length);
				final RealMatrix baseRes = mBase.getCovariances();
				for (int r = 0; r < mIndexes.length; ++r)
					for (int c = 0; c < mIndexes.length; ++c)
						res.setEntry(r, c, baseRes.getEntry(mIndexes[r], mIndexes[c]));
				return res;
			}
		};

		protected ReorderedUncertainValues(final List<? extends L> labels, final UncertainValuesBase<? super L> base)//
				throws ArgumentException {
			super(new ArrayList<L>(labels));
			mIndexes = new int[labels.size()];
			mBase = base;
			List<L> missing = new ArrayList<>();
			for (int i = 0; i < mIndexes.length; ++i) {
				final L label = labels.get(i);
				final int idx = mBase.indexOf(label);
				if (idx < 0)
					missing.add(label);
				mIndexes[i] = idx;
			}
			if (missing.size() > 0)
				throw new ArgumentException("The argument labels " + missing + " are missing.");
		}

		@Override
		public RealMatrix getCovariances() {
			return mCovariances.get();
		}

		@Override
		public RealVector getValues() {
			return mValues.get();
		}

	};

	private static class UVSOutOfRangeException extends OutOfRangeException {

		private static final long serialVersionUID = -8886319004986992747L;

		private final String mExtra;

		public UVSOutOfRangeException(final Number wrong, final Number low, final Number high, final String extra) {
			super(wrong, low, high);
			mExtra = extra;
		}

		@Override
		public String toString() {
			return super.toString() + " - " + mExtra;
		}

	}

	protected static final double MAX_CORR = 1.001;

	/**
	 * Builds an {@link UncertainValuesBase} object representing the specified
	 * labeled quantities as extracted from the list of {@link UncertainValuesBase}
	 * objects.
	 *
	 * @param labels
	 * @param uvs
	 * @return {@link UncertainValuesBase}
	 * @throws ArgumentException
	 */
	public static <J> UncertainValues<J> build(//
			final List<J> labels, //
			final List<UncertainValuesBase<J>> uvs //
	) throws ArgumentException {
		// Test that each requested label is defined once and only once.
		for (final J label : labels) {
			int count = 0;
			for (final UncertainValuesBase<? extends J> uv : uvs)
				if (uv.hasEntry(label))
					count++;
			if (count < 1)
				throw new ArgumentException(
						"The label " + label + " is not defined in one of the input UncertainValues.");
			if (count > 1)
				throw new ArgumentException(
						"The label " + label + " is muliply defined in one of the input UncertainValues.");
		}
		final UncertainValues<J> res = new UncertainValues<>(labels);
		for (final UncertainValuesBase<J> uv : uvs)
			for (final J label1 : labels)
				if (uv.hasEntry(label1)) {
					res.set(label1, uv.getEntry(label1), uv.getVariance(label1));
					for (final J label2 : uv.getLabels())
						if (res.hasEntry(label2)) {
							final double cv = uv.getCovariance(label1, label2);
							if (cv != 0.0)
								res.setCovariance(label1, label2, cv);
						}
				}
		return res;
	}

	/**
	 * Combines a disjoint set of {@link UncertainValues} into a single one.
	 * (Disjoint meaning not sharing a common label.)
	 *
	 * @param uvs List&lt;UncertainValues&gt;
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public static <J> UncertainValuesBase<J> combine( //
			final List<? extends UncertainValuesBase<? extends J>> uvs, //
			final boolean takeFirst//
	) throws ArgumentException {
		// Test that each requested label is defined once and only once.
		final List<J> labels = new ArrayList<>();
		for (final UncertainValuesBase<? extends J> uv : uvs)
			labels.addAll(uv.getLabels());
		return new CombinedUncertainValues<>(labels, uvs, takeFirst);
	}

	/**
	 * Creates a bitmap that represents the difference between uncertainties
	 * associated with these two sets of UncertainValuesBase. The difference between
	 * the covariances is plotted.
	 *
	 * @param uvs1
	 * @param uvs2
	 * @param corr
	 * @param pixDim
	 * @return BufferedImage
	 * @throws ArgumentException
	 */
	public static <J> BufferedImage compareAsBitmap(//
			final UncertainValuesBase<J> uvs1, //
			final UncertainValuesBase<J> uvs2, //
			final IValueToColor corr, //
			final int pixDim //
	) throws ArgumentException {
		if (uvs1.getDimension() != uvs2.getDimension())
			throw new DimensionMismatchException(uvs2.getDimension(), uvs1.getDimension());
		final int dim = uvs1.getDimension();
		final boolean[] disp = new boolean[dim];
		final UncertainValuesBase<J> uvs2r = UncertainValues.extract(uvs1.getLabels(), uvs2);
		final BufferedImage bi = new BufferedImage(pixDim * (dim + 2), pixDim * dim, BufferedImage.TYPE_3BYTE_BGR);
		final Graphics2D g2 = bi.createGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, bi.getWidth(), bi.getHeight());
		for (int rr = 0; rr < dim; ++rr) {
			final double v1 = uvs1.getEntry(rr), v2 = uvs2r.getEntry(rr);
			// Red if values are substantially different
			if ((Math.abs(v1 - v2) > 0.01 * Math.max(v1, v2)) && (Math.max(v1, v2) > 1.0e-10))
				g2.setColor(Color.RED);
			else
				g2.setColor(Color.white);
			g2.fillRect(0, rr * pixDim, pixDim, pixDim);
			final double dv1 = uvs1.getUncertainty(rr), dv2 = uvs2r.getUncertainty(rr);
			if ((dv1 < 1.0e-6 * Math.abs(v1)) && (dv2 < 1.0e-6 * Math.abs(v2))) {
				g2.setColor(Color.WHITE);
				disp[rr] = false;
			} else {
				disp[rr] = true;
				final double delta = 0.5 * Math.abs((dv1 - dv2) / (dv1 + dv2));
				g2.setColor(corr.map(delta));
			}
			g2.fillRect((2 + rr) * pixDim, rr * pixDim, pixDim, pixDim);
		}
		for (int rr = 0; rr < dim; ++rr) {
			if (disp[rr]) {
				for (int cc = 0; cc < dim; ++cc) {
					if (disp[cc]) {
						final double cc1 = uvs1.getCorrelationCoefficient(rr, cc);
						final double cc2 = uvs2r.getCorrelationCoefficient(rr, cc);
						final double delta = 0.5 * Math.abs(cc1 - cc2);
						g2.setColor(corr.map(delta));
						g2.fillRect((2 + rr) * pixDim, rr * pixDim, pixDim, pixDim);
					}
				}
			}
		}
		return bi;
	}

	public static <J> UncertainValues<J> forceMinCovariance(final UncertainValuesBase<J> vals,
			final RealVector minCov) {
		assert vals.getDimension() == minCov.getDimension();
		final RealMatrix cov = vals.getCovariances().copy();
		for (int rc = 0; rc < vals.getDimension(); ++rc)
			if (cov.getEntry(rc, rc) < minCov.getEntry(rc))
				cov.setEntry(rc, rc, minCov.getEntry(rc));
		return new UncertainValues<J>(vals.getLabels(), vals.getValues(), cov);
	}

	/**
	 * Returns the UncertainValues that result from applying the function/Jacobian
	 * in <code>nmjf</code> to the input values/variances in <code>input</code>.
	 *
	 *
	 * @param nmjf
	 * @param input
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static <J> UncertainValuesCalculator<J> propagate( //
			final LabeledMultivariateJacobianFunction<? extends J, ? extends J> nmjf, //
			final UncertainValuesBase<J> input //
	) throws ArgumentException {
		return new UncertainValuesCalculator<J>(nmjf, input);
	}

	/**
	 * Returns the UncertainValues that result from applying the function/Jacobian
	 * in <code>nmjf</code> to the input values/variances in <code>input</code>.
	 *
	 *
	 * @param nmjf
	 * @param input
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static <J> UncertainValuesCalculator<J> propagateFiniteDifference( //
			final LabeledMultivariateJacobianFunction<? extends J, ? extends J> nmjf, //
			final UncertainValuesBase<J> input, //
			final RealVector dinp //
	) throws ArgumentException {
		UncertainValuesCalculator<J> res = new UncertainValuesCalculator<J>(nmjf, input);
		res.setCalculator(res.new FiniteDifference(dinp));
		return res;
	}

	/**
	 * Returns the UncertainValues that result from applying the function/Jacobian
	 * in <code>nmjf</code> to the input values/variances in <code>input</code>.
	 *
	 *
	 * @param nmjf
	 * @param input
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static <J> UncertainValuesCalculator<J> propagateMonteCarlo( //
			final LabeledMultivariateJacobianFunction<? extends J, ? extends J> nmjf, //
			final UncertainValuesBase<J> input, //
			final int nEvals //
	) throws ArgumentException {
		UncertainValuesCalculator<J> res = new UncertainValuesCalculator<J>(nmjf, input);
		res.setCalculator(res.new MonteCarlo(nEvals));
		return res;
	}

	public static <J> boolean testEquality(final UncertainValuesBase<J> uvs1, final UncertainValuesBase<J> uvs2) {
		if (uvs1.getDimension() != uvs2.getDimension())
			return false;
		for (final J label1 : uvs1.getLabels()) {
			if (uvs1.getEntry(label1) != uvs2.getEntry(label1))
				return false;
			for (final J label2 : uvs1.getLabels())
				if (uvs1.getCovariance(label1, label2) != uvs2.getCovariance(label1, label2))
					return false;
		}
		return true;
	}

	/**
	 * Return an {@link UncertainValues} object with the same dimension and values
	 * as input except all the covariances except those associated with label are
	 * zeroed.
	 *
	 * @param label
	 * @param input
	 * @return {@link UncertainValues}
	 */
	public static <J> UncertainValues<J> zeroBut(final J label, final UncertainValuesBase<J> input) {
		final UncertainValues<J> res = new UncertainValues<J>(input.getLabels(), input.getValues(), 0.0);
		final int idx = input.indexOf(label);
		res.setCovariance(idx, idx, input.getCovariance(idx, idx));
		return res;
	}

	/**
	 * Checks that the covariance matrix is m x m, the variances are non-negative,
	 * the covariances are bounded correlation coefficients between -1<r<1.
	 *
	 * @param m
	 * @param covar
	 */
	private static void validateCovariance(final int m, final RealMatrix covar, final double maxVal) {
		if (covar.getRowDimension() != m)
			throw new DimensionMismatchException(covar.getRowDimension(), m);
		if (covar.getColumnDimension() != m)
			throw new DimensionMismatchException(covar.getColumnDimension(), m);
		final double EPS = 1.0e-10, SREPS = Math.sqrt(EPS);
		for (int r = 0; r < covar.getRowDimension(); ++r) {
			final double entryRR = covar.getEntry(r, r);
			if (entryRR < 0.0) {
				// This is necessary because of subtle rounding errors
				if (entryRR > -EPS * maxVal)
					covar.setEntry(r, r, 0.0);
				else
					throw new UVSOutOfRangeException(entryRR, 0.0, Double.MAX_VALUE, "Row=" + r);
			}
		}
		for (int r = 0; r < covar.getRowDimension(); ++r)
			for (int c = r + 1; c < covar.getColumnDimension(); ++c) {
				final double entryRC = covar.getEntry(r, c);
				if (Math.abs(entryRC - covar.getEntry(c, r)) > SREPS * Math.abs(entryRC) + EPS)
					throw new OutOfRangeException(entryRC, covar.getEntry(c, r) - EPS, covar.getEntry(c, r) + EPS);
				final double max = Math.sqrt(covar.getEntry(c, c) * covar.getEntry(r, r)) + EPS;
				final double rr = entryRC / max;
				if ((rr > MAX_CORR) || (rr < -MAX_CORR)) {
					if (Math.abs(entryRC) > EPS)
						throw new UVSOutOfRangeException(entryRC, -max, max, "Row=" + r + ", Col=" + c);
					else {
						covar.setEntry(r, c, Math.signum(rr) * max);
						covar.setEntry(c, r, Math.signum(rr) * max);
					}
				}
			}
	}

	private final List<H> mLabels;

	private int mHashCode = 0;

	/**
	 * Constructs a UncertainValuesBase object based on the specified labels with
	 * zero values and NaN covariances.
	 *
	 * @param labels List&lt;Object&gt; A list of objects implementing hashCode()
	 *               and equals().
	 */
	protected UncertainValuesBase(//
			final List<H> labels //
	) {
		mLabels = new FastIndex<>(labels);
	}

	public BufferedImage asCovarianceBitmap(final int dim, final IValueToColor sigma, final IValueToColor corr) {
		final RealMatrix sc = NullableRealMatrix.build(getDimension(), getDimension());
		final RealVector values = getValues();
		for (int r = 0; r < getDimension(); ++r) {
			final double crr = getCovariance(r, r);
			final double rVal = values.getEntry(r);
			sc.setEntry(r, r, Math.sqrt(crr) / rVal);
			if (Math.sqrt(crr) > 1.0e-8 * Math.abs(rVal)) {
				for (int c = r + 1; c < getDimension(); ++c) {
					final double ccc = getCovariance(c, c);
					if (Math.sqrt(ccc) > 1.0e-8 * Math.abs(values.getEntry(c))) {
						final double rr = getCovariance(r, c) / (Math.sqrt(crr * ccc));
						assert rr >= -MAX_CORR : rr + " at " + getLabel(r) + ", " + getCovariance(r, c) + ", " + crr
								+ ", " + ccc;
						assert rr <= MAX_CORR : rr + " at " + getLabel(r) + ", " + getCovariance(r, c) + ", " + crr
								+ ", " + ccc;
						sc.setEntry(r, c, rr);
						sc.setEntry(c, r, rr);
					}
				}
			}
		}
		final BufferedImage bi = new BufferedImage(dim * getDimension(), dim * getDimension(),
				BufferedImage.TYPE_3BYTE_BGR);
		final Graphics2D g2 = bi.createGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, bi.getWidth(), bi.getHeight());
		for (int r = 0; r < getDimension(); ++r) {
			g2.setColor(sigma.map(Math.sqrt(sc.getEntry(r, r))));
			g2.fillRect(r * dim, r * dim, dim, dim);
			for (int c = r + 1; c < getDimension(); ++c) {
				final double entry = sc.getEntry(r, c);
				if (!Double.isNaN(entry))
					g2.setColor(corr.map(entry));
				else
					g2.setColor(Color.yellow);
				g2.fillRect(c * dim, r * dim, dim, dim);
				g2.fillRect(r * dim, c * dim, dim, dim);
			}
		}
		return bi;
	}

	/**
	 * Check that the indices are valid and not repeated. Uses assert rather than an
	 * Exception.
	 *
	 * @param indices
	 * @return true
	 */
	public boolean assertIndices(final int[] indices) {
		for (int i = 0; i < indices.length; ++i) {
			assert indices[i] >= 0 : "Index[" + i + "] is less than zero.";
			assert indices[i] < mLabels.size() : "Index[" + i + "] is larger than the number of labels.";
			for (int j = i + 1; j < indices.length; ++j)
				assert indices[i] != indices[j] : "Duplicated index: Index[" + i + "] equals Index[" + j + "]";
		}
		return true;
	}

	/**
	 * Takes the input UncertainValues and creates a new, reordered UncertainValues
	 * object in which the correlated rows/columns are grouped together. THe
	 * non-zero covariances are moved away from the edges and upper-right and
	 * lower-left corners towards the diagonal.
	 *
	 *
	 * @return UncertainValues
	 */

	public UncertainValuesBase<H> blockDiagnonalize() {
		final int[] count = new int[getDimension()];
		for (int r = 0; r < count.length; ++r)
			for (int c = 0; c < count.length; ++c)
				if ((r != c) && (getCovariance(r, c) != 0.0))
					++count[r];
		final int[] idx = new int[getDimension()];
		final boolean[] done = new boolean[getDimension()];
		Arrays.fill(done, false);
		int next = 0;
		for (int covs = 0; (next < count.length) && (covs < count.length); ++covs) {
			for (int r = 0; (next < count.length) && (r < count.length); ++r) {
				if ((!done[r]) && (count[r] == covs)) {
					final List<Integer> cols = new ArrayList<>();
					idx[next] = r;
					++next;
					done[r] = true;
					cols.add(r);
					while (!cols.isEmpty()) {
						final int col = cols.remove(0);
						done[col] = true;
						for (int covs2 = covs; covs2 < count.length; ++covs2) {
							for (int rr = 0; rr < done.length; ++rr)
								if ((!done[rr]) && (count[rr] == covs2) && (getCovariance(rr, col) != 0.0)) {
									idx[next] = rr;
									++next;
									done[rr] = true;
									cols.add(rr);
								}
						}
					}
				}
			}
		}
		assert next == count.length;
		try {
			return extract(idx);
		} catch (final ArgumentException e) {
			System.err.println("Should never happen!");
			e.printStackTrace();
		}
		return null;
	}

	final public RealMatrix delta(final UncertainValuesBase<H> zeroed) {
		return getCovariances().subtract(zeroed.getCovariances());
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final IUncertainValues<?> other = (IUncertainValues<?>) obj;
		return Objects.equals(mLabels, other.getLabels());
	}

	/**
	 * Tests the equality of the values, variances and correlation coefficients to
	 * within the fractional tolerance.
	 * 
	 * @param uv2
	 * @param tol
	 * @return boolean
	 */
	public boolean equals(final UncertainValuesBase<H> uv2, final double tol) {
		for (int r = 0; r < this.getDimension(); ++r) {
			final double tmp = Math.abs((getEntry(r) - uv2.getEntry(mLabels.get(r))) / getEntry(r));
			if (Double.isFinite(tmp))
				if (tmp > tol)
					return false;
			final double var = Math.abs((getVariance(r) - getVariance(mLabels.get(r))) / getEntry(r));
			if (Double.isFinite(var))
				if (var > tol)
					return false;
			for (int c = r; c < this.getDimension(); ++c)
				if (Math.abs(getCorrelationCoefficient(r, c)
						- uv2.getCorrelationCoefficient(mLabels.get(r), mLabels.get(c))) > tol)
					return false;
		}
		return true;
	}

	/**
	 * Creates a new {@link UncertainValues} object representing only those
	 * rows/columns whose indexes are in idx. Can be used to create a sub-set of
	 * this {@link UncertainValues} or reorder this {@link UncertainValues}.
	 *
	 * @param idx
	 * @return {@link UncertainValues} A new object
	 * @throws ArgumentException
	 */
	final public UncertainValuesBase<H> extract(final int[] indices) //
			throws ArgumentException {
		return extract(getLabels(indices));
	}
	
	/**
	 * Creates a new {@link UncertainValues} object representing only those
	 * rows/columns whose indexes are in idx. Can be used to create a sub-set of
	 * this {@link UncertainValues} or reorder this {@link UncertainValues}.
	 *
	 * @param idx
	 * @return {@link UncertainValues} A new object
	 * @throws ArgumentException
	 */
	final public UncertainValuesBase<H> extract(final List<? extends H> labels) //
			throws ArgumentException {
		return reorder(labels);
	}


	/**
	 * Extract an array of values from this UncertainValue object for the specified
	 * list of labels in the order specified by the label list.
	 *
	 * @param labels List&lt;H&gt;
	 * @return RealVector
	 */
	final public RealMatrix extractCovariances(final List<H> labels) {
		final RealMatrix res = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		for (int r = 0; r < labels.size(); ++r) {
			final int ridx = indexOf(labels.get(r));
			res.setEntry(r, r, getCovariance(r, r));
			for (int c = r + 1; c < labels.size(); ++c) {
				final int cidx = indexOf(labels.get(c));
				final double cv = getCovariance(ridx, cidx);
				res.setEntry(r, c, cv);
				res.setEntry(c, r, cv);
			}
		}
		return res;
	}

	/**
	 * Extract an array of values from this UncertainValue object for the specified
	 * list of labels in the order specified by the label list.
	 *
	 * @param labels List&lt;H&gt;
	 * @return RealVector
	 */
	final public RealVector extractValues(final List<? extends H> labels) {
		final RealVector res = new ArrayRealVector(labels.size());
		int i = 0;
		for (final H label : labels) {
			assert indexOf(label) != -1 : label + " is missing in extractValues(...)";
			res.setEntry(i, getEntry(label));
			++i;
		}
		return res;
	}

	public List<H> findMissing(final List<H> inputLabels) {
		final List<H> res = new ArrayList<>();
		for (final H lbl : inputLabels)
			if (indexOf(lbl) == -1)
				res.add(lbl);
		return res;
	}

	/**
	 * Returns the covariance associated with the specific labels.
	 *
	 * @param p1 H
	 * @param p2 H
	 * @return double The correlation coefficient associated with indices p1 and p2
	 */
	final public double getCorrelationCoefficient(final H p1, final H p2) {
		return getCorrelationCoefficient(indexOf(p1), indexOf(p2));
	}

	/**
	 * Returns the covariance associated with the specific integer indices.
	 *
	 * @param p1 int
	 * @param p2 int
	 * @return double The correlation coefficient associated with indices p1 and p2
	 */
	final public double getCorrelationCoefficient(final int p1, final int p2) {
		final RealMatrix cov = getCovariances();
		final double c12 = cov.getEntry(p1, p2);
		if (c12 == 0.0)
			return 0.0;
		else {
			final double c11 = cov.getEntry(p1, p1);
			final double c22 = cov.getEntry(p2, p2);
			if (c11 * c22 > 0) {
				final double rho = c12 / Math.sqrt(c11 * c22);
				assert (rho > -1.00000001) && (rho < 1.00000001) : rho;
				return Math.min(1.0, Math.max(-1.0, rho));
			} else
				return Double.NaN;
		}
	}

	final public RealMatrix getCorrelationCoefficients() {
		final int dim = getDimension();
		final RealMatrix rm = MatrixUtils.createRealMatrix(dim, dim);
		for (int r = 0; r < dim; ++r)
			for (int c = r + 1; c < dim; ++c) {
				final double cc = getCorrelationCoefficient(r, c);
				rm.setEntry(r, c, cc);
				rm.setEntry(c, r, cc);
			}
		return rm;
	}

	/**
	 * Returns the covariance associated with the specific labels
	 *
	 * @param label1 H
	 * @param label2 H
	 * @return double The covariance associated with labels label1 and label2 or 0
	 *         if one or both labels unknown
	 */
	final public double getCovariance(final H label1, final H label2) {
		final int p1 = indexOf(label1), p2 = indexOf(label2);
		return (p1 != -1) && (p2 != -1) ? getCovariance(p1, p2) : 0.0;
	}

	/**
	 * Returns the covariance associated with the specific integer indices.
	 *
	 * @param p1    int
	 * @param p2int
	 * @return double The covariance associated with indices p1 and p2
	 */
	final public double getCovariance(final int p1, final int p2) {
		return getCovariances().getEntry(p1, p2);
	}

	/**
	 * Returns a map of the variance and covariances relative to the specified
	 * label.
	 *
	 * @param label
	 * @return Map&lt;H, Double&gt;
	 */
	final public Map<H, Double> getCovariances(final H label) {
		final Map<H, Double> res = new HashMap<>();
		final int p = indexOf(label);
		if (p != -1) {
			for (int i = 0; i < mLabels.size(); ++i) {
				final RealMatrix cov = getCovariances();
				if (cov.getEntry(p, i) != 0.0)
					res.put(mLabels.get(i), cov.getEntry(p, i));
			}
		}
		return res;
	}

	/**
	 * Returns a RealMatrix of the variance and covariances relative to the
	 * specified list of labels.
	 *
	 * @param label
	 * @return Map&lt;H, Double&gt;
	 */
	final public RealMatrix getCovariances(final List<H> labels) {
		final int[] idx = new int[labels.size()];
		for (int i = 0; i < idx.length; ++i)
			idx[i] = indexOf(labels.get(i));
		final RealMatrix res = MatrixUtils.createRealMatrix(idx.length, idx.length);
		final RealMatrix cov = getCovariances();
		for (int i = 0; i < idx.length; ++i)
			for (int j = 0; j < idx.length; ++j) {
				res.setEntry(i, j, cov.getEntry(idx[i], idx[j]));
			}
		return res;
	}

	/**
	 * Returns the number of labels which is equivalent to the number of values and
	 * row/columns in the covariance matrix.
	 *
	 * @return int
	 */
	final public int getDimension() {
		return getLabels().size();
	}

	/**
	 * Returns the value associated with the entry associated with label.
	 *
	 * @param label Object
	 * @return double
	 */
	final public double getEntry(final H label) {
		final int p = indexOf(label);
		assert p >= 0 : "Label " + label.toString() + " missing.";
		return getValues().getEntry(p);
	}

	/**
	 * Returns the function value associated with the label specified by index.
	 *
	 * @param p
	 * @return double
	 */
	final public double getEntry(final int p) {
		return getValues().getEntry(p);
	}

	/**
	 * Returns the value associated with the entry associated with label if the
	 * label is associated with a value or returns defVal otherwise. Useful for a
	 * quantity that may or may not have an associated uncertainty.
	 *
	 * @param label
	 * @param defVal
	 * @return double
	 */
	final public double getEntryWithDefault(final H label, final double defVal) {
		final int p = indexOf(label);
		return p >= 0 ? getEntry(p) : defVal;
	}

	/**
	 * Returns the label at the p-th entry in the values vector and in the p-th row
	 * and column in the covariance matrix.
	 *
	 * @param p
	 * @return Object
	 */
	final public H getLabel(final int p) {
		return getLabels().get(p);
	}

	/**
	 * Returns an ordered list of the Object labels associated with the values and
	 * covariances in this object.
	 *
	 * @return List&lt;Object&gt;
	 */
	@Override
	final public List<H> getLabels() {
		return Collections.unmodifiableList(mLabels);
	}

	/**
	 * Extracts all labels assignable as cls
	 *
	 * @param cls<T> The class type
	 * @return List&lt;T&gt;
	 */
	final public <T> List<T> getLabels(final Class<T> cls) {
		final List<T> res = new ArrayList<>();
		for (final Object tag : getLabels())
			if (cls.isInstance(tag))
				res.add(cls.cast(tag));
		return Collections.unmodifiableList(res);
	}

	/**
	 * Returns the uncertainty associated with the specific label
	 *
	 * @param label H
	 * @return double The uncertainty associated with label or 0.0 if label unknown
	 */
	final public double getUncertainty(final H label) {
		final int p = indexOf(label);
		return p != -1 ? Math.sqrt(getCovariance(p, p)) : 0.0;
	}

	/**
	 * Returns the uncertainty associated with the specific index
	 *
	 * @param p int
	 * @return double The uncertainty associated with index p
	 */
	final public double getUncertainty(final int p) {
		return Math.sqrt(getCovariance(p, p));
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified
	 * label. The value is the same as getEntry(...) and the uncertainty is the
	 * on-diagonal covariance matrix entry associated with that entry.
	 *
	 * @param label
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValue(final H label) {
		final int p = indexOf(label);
		return new UncertainValue(getEntry(p), label, Math.sqrt(getCovariance(p, p)));
	}

	final public Map<H, UncertainValue> getUncertainValueMap() {
		final Map<H, UncertainValue> res = new HashMap<>();
		for (final H label : getLabels())
			res.put(label, getUncertainValue(label));
		return res;
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified
	 * label (if defined). The value is the same as getEntry(...) and the
	 * uncertainty is the on-diagonal covariance matrix entry associated with that
	 * entry. If the label isn't defined, then the value defVal will be returned.
	 *
	 * @param label
	 * @param defVal
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValueWithDefault(final H label, final UncertainValue defVal) {
		final int p = indexOf(label);
		return p >= 0 ? new UncertainValue(getEntry(p), label, getCovariance(p, p)) : defVal;

	}

	/**
	 * Returns the value and the associated uncertainty as an {@link UncertainValue}
	 * object.
	 *
	 * @param label
	 * @return {@link UncertainValue}
	 */
	final public UncertainValue getValue(final H label) {
		return new UncertainValue(getEntry(label), label, Math.sqrt(getVariance(label)));
	}

	final public Map<H, Double> getValueMap() {
		return getValueMap(getLabels());
	}

	final public <T> Map<T, Double> getValueMap(final Class<T> cls) {
		final Map<T, Double> res = new HashMap<>();
		for (final H label : getLabels())
			if (cls.isInstance(label))
				res.put(cls.cast(label), getEntry(label));
		return res;
	}

	final public Map<H, Double> getValueMap(final List<H> labels) {
		final Map<H, Double> res = new HashMap<>();
		for (final H label : labels)
			res.put(label, getEntry(label));
		return res;
	}

	/**
	 * Returns a {@link RealVector} containing the values associated with this
	 * object in the order specified by the List labels.
	 *
	 * @param labels
	 * @return {@link RealVector}
	 */
	final public RealVector getValues(final List<? extends H> labels) {
		final RealVector res = new ArrayRealVector(labels.size());
		for (int i = 0; i < res.getDimension(); ++i)
			res.setEntry(i, getValues().getEntry(indexOf(labels.get(i))));
		return res;
	}

	/**
	 * Returns the variance associated with the specific label
	 *
	 * @param label H
	 * @return double The variance associated with label
	 */
	final public double getVariance(final H label) {
		final int p = indexOf(label);
		return getCovariance(p, p);
	}

	/**
	 * Returns the variance associated with the specific index
	 *
	 * @param idx int
	 * @return double The variance associated with label
	 */
	final public double getVariance(final int idx) {
		return getCovariance(idx, idx);
	}

	/**
	 * Is there a value and covariances associated with the specified label?
	 *
	 * @param label
	 * @return boolean
	 */
	@Override
	final public boolean hasEntry(final Object label) {
		return getLabels().indexOf(label) != -1;
	}

	@Override
	public int hashCode() {
		if (mHashCode == 0)
			mHashCode = Objects.hash(mLabels);
		return mHashCode;
	}

	/**
	 * Returns the index associated with the specified label or -1 if not found.
	 *
	 * @param p
	 * @return Object
	 */
	final public int indexOf(final Object label) {
		return getLabels().indexOf(label);
	}

	/**
	 * Note: Returns an index of -1 if the label is missing
	 *
	 * @param labels
	 * @return Returns an array of integer indices for the specified labels in order
	 */
	final public int[] indices(final List<?> labels) {
		final int[] res = new int[labels.size()];
		for (int i = 0; i < res.length; ++i)
			res[i] = indexOf(labels.get(i));
		return res;
	}

	/**
	 * Checks the values and covariances to determine whether any are equivalent to
	 * NaN.
	 *
	 * @return boolean true if one value is NaN, false otherwise.
	 */
	final public boolean isNaN() {
		if (getValues().isNaN())
			return true;
		final RealMatrix cov = getCovariances();
		for (int r = 0; r < cov.getRowDimension(); ++r)
			for (int c = 0; c < cov.getColumnDimension(); ++c)
				if (Double.isNaN(getCovariance(r, c)))
					return true;
		return false;
	}

	/**
	 * Returns an UncertainValuesBase object with the labels in the order specified
	 * by the argument. If this is already in this order, this is returned;
	 * otherwise a new {@link UncertainValues} object is created.
	 *
	 * @param list
	 * @return {@link UncertainValuesBase}
	 * @throws ArgumentException
	 */
	@SuppressWarnings("unchecked")
	public <L extends H> UncertainValuesBase<L> reorder(final List<? extends L> list) throws ArgumentException {
		if (mLabels.size() == list.size()) {
			// Check if already in correct order...
			boolean eq = true;
			for (int i = 0; (i < list.size()) && eq; ++i)
				if (!mLabels.get(i).equals(list.get(i))) {
					eq = false;
					break;
				}
			if (eq)
				return (UncertainValuesBase<L>) this;
		}
		return new ReorderedUncertainValues<L>(list, this);
	}

	/**
	 * Build a new UncertainValues from this one in which the labels have been
	 * sorted into alphabetical order by label.toString().
	 *
	 * @return {@link UncertainValues} A new instance with the same data reordered.
	 */
	final public UncertainValuesBase<H> sort() {
		return sort(new Comparator<H>() {

			@Override
			public int compare(final H o1, final H o2) {
				return o1.toString().compareTo(o2.toString());
			}
		});
	}

	/**
	 * Build a new UncertainValues from this one in which the labels have been
	 * sorted into an order determined by the specified {@link Comparator}.
	 *
	 * @param compare
	 * @return {@link UncertainValues} A new instance with the same data reordered.
	 */
	final public UncertainValuesBase<H> sort(final Comparator<H> compare) {
		final List<H> labels = new ArrayList<>(getLabels());
		labels.sort(compare);
		try {
			return reorder(labels);
		} catch (final ArgumentException e) {
			System.err.println("This should never happen.");
			e.printStackTrace();
		}
		return null;
	}

	final public String toCSV() {
		final StringBuffer sb = new StringBuffer(4096);
		sb.append("\"Name\",\"Value\"");
		final RealMatrix cov = getCovariances();
		final List<H> labels = getLabels();
		for (int c = 0; c < cov.getColumnDimension(); ++c) {
			sb.append(",\"");
			sb.append(labels.get(c));
			sb.append("\"");
		}
		sb.append("\n");
		final RealVector vals = getValues();
		for (int r = 0; r < vals.getDimension(); ++r) {
			sb.append("\"");
			sb.append(labels.get(r));
			sb.append("\",");
			sb.append(vals.getEntry(r));
			for (int c = 0; c < cov.getColumnDimension(); ++c) {
				sb.append(",");
				sb.append(cov.getEntry(r, c));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	/**
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mode, new BasicNumberFormat("0.00E0"));
	}

	/**
	 * Provides a mechanism to convert this {@link UncertainValuesBase} object to
	 * HTML with a little more control over number formating.
	 *
	 * @param mode
	 * @param nf
	 * @return String
	 */
	final public String toHTML(final Mode mode, final BasicNumberFormat nf) {
		switch (mode) {
		case TERSE: {
			final Table table = new Table();
			final List<Item> header = new ArrayList<>();
			final List<Item> vals = new ArrayList<>();
			for (final H rowLabel : getLabels()) {
				header.add(Table.th(HTML.toHTML(rowLabel, Mode.TERSE)));
				vals.add(Table.tdc(
						nf.formatHTML(getEntry(rowLabel)) + "&pm;" + nf.formatHTML(Math.sqrt(getVariance(rowLabel)))));
			}
			table.addRow(header);
			table.addRow(vals);
			return table.toHTML(Mode.NORMAL);
		}
		case NORMAL: {
			final Map<H, UncertainValue> tmp = getUncertainValueMap();
			final Table t = new Table();
			final Map<String, H> tagMap = new TreeMap<>();
			for (final H tag : tmp.keySet())
				tagMap.put(HTML.toHTML(tag, Mode.TERSE), tag);
			t.addRow(Table.th("Label"), Table.th("Value"), Table.th("Uncertainty"), Table.th("Fractional"));
			final HalfUpFormat df = new HalfUpFormat("0.0%");
			for (final Map.Entry<String, H> me : tagMap.entrySet()) {
				final UncertainValue uv = tmp.get(me.getValue());
				t.addRow(Table.td(me.getKey()), Table.td(uv.doubleValue()), Table.td(uv.uncertainty()),
						Table.td(df.format(uv.fractionalUncertainty())));
			}
			return t.toHTML(Mode.NORMAL);
		}
		case VERBOSE:
		default: {
			final Table vals = new Table();
			final List<Item> all = new ArrayList<>();
			all.add(Table.td("Name"));
			all.add(Table.tdc("Quantity"));
			all.add(Table.td());
			final List<H> labels = getLabels();
			for (final H colLabel : labels)
				all.add(Table.tdc(HTML.toHTML(colLabel, Mode.TERSE)));
			vals.addRow(all);
			for (int r = 0; r < labels.size(); ++r) {
				final H rowLabel = labels.get(r);
				final List<Item> row = new ArrayList<>();
				row.add(Table.td(HTML.toHTML(rowLabel, Mode.TERSE)));
				row.add(Table.tdc(nf.formatHTML(getEntry(rowLabel))));
				if (r == (labels.size() - 1) / 2)
					row.add(Table.tdc("&nbsp;&nbsp;&plusmn;&nbsp;&nbsp;"));
				else
					row.add(Table.td());
				for (int c = 0; c < mLabels.size(); ++c)
					row.add(Table.tdc(toHTML_Covariance(r, c)));
				vals.addRow(row);
			}
			return vals.toHTML(Mode.NORMAL);
		}
		}
	}

	/**
	 * Convert the covariance at r,c into HTML in a human-friendly manner. Variances
	 * are converted into "(v)^2" and covariances into the correlation coefficient
	 * times sR sC.
	 *
	 * @param r
	 * @param c
	 * @return String
	 */
	final public String toHTML_Covariance(final int r, final int c) {
		final double val = getCovariance(r, c);
		if (r == c) {
			final BasicNumberFormat nf = new BasicNumberFormat("0.00E0");
			final String html = "(" + nf.formatHTML(Math.sqrt(Math.abs(val))) + ")<sup>2</sup>";
			if (c >= 0)
				return html;
			else // This is a problem! Highlight it....
				return HTML.fontColor(Transforms.NON_BREAKING_DASH + html, Color.RED);
		} else {
			final double vr = getCovariance(r, r), vc = getCovariance(c, c);
			if ((vr != 0.0) && (vc != 0.0)) {
				final BasicNumberFormat nf = new BasicNumberFormat("0.000");
				final double cc = val / Math.sqrt(vr * vc);
				if (Math.abs(cc) < 1.0e-5)
					return "&nbsp;";
				else {
					final String html = nf.formatHTML(cc) + "&middot;&sigma;<sub>R</sub>&sigma;<sub>C</sub>";
					if ((cc >= -MAX_CORR) || (cc <= MAX_CORR))
						return html;
					else
						return "&nbsp;";
				}
			} else {
				if (val == 0.0)
					return "&nbsp";
				else {
					final BasicNumberFormat nf = new BasicNumberFormat("0.00E0");
					return HTML.fontColor(nf.format(val), Color.RED);
				}
			}
		}
	}

	final public String toSimpleHTML(final BasicNumberFormat bnf) {
		final Table t0 = new Table();
		{
			final List<Table.Item> row = new ArrayList<>();
			row.add(Table.th("Label"));
			row.add(Table.thc("Value"));
			row.add(Table.thc("&nbsp;"));
			for (int c = 0; c < getDimension(); ++c)
				row.add(Table.thc(HTML.toHTML(getLabel(c), Mode.NORMAL)));
			t0.addRow(row);
		}
		for (int r = 0; r < getDimension(); ++r) {
			final List<Table.Item> row = new ArrayList<>();
			row.add(Table.thc(HTML.toHTML(getLabel(r), Mode.NORMAL)));
			row.add(Table.tdc(bnf.format(getEntry(r))));
			row.add(Table.tdc(r == getDimension() / 2 ? "&#177;" : "&nbsp;"));
			for (int c = 0; c < getDimension(); ++c)
				row.add(Table.tdc(bnf.format(getCovariance(r, c))));
			t0.addRow(row);
		}
		return t0.toHTML(Mode.NORMAL);
	}

	@Override
	public String toString() {
		final List<H> labels = getLabels();
		final List<H> labelSub = new ArrayList<>();
		for (int i = 0; (i < labels.size()) && (i < 5); ++i)
			labelSub.add(labels.get(i));
		final String lblStr = labelSub.toString();
		return "UVS[" + lblStr.substring(1, lblStr.length() - 1)
				+ (labelSub.size() < labels.size() ? "+" + (labels.size() - labelSub.size()) + " more" : "") + "]";
	}

	final public void validateCovariance() throws ArgumentException {
		final List<H> labels = getLabels();
		final Set<H> dupLabels = new HashSet<>(labels);
		if (dupLabels.size() != labels.size())
			for (final H label : labels)
				if (!dupLabels.contains(label))
					throw new ArgumentException("The label " + label + " is repeated.");
		validateCovariance(labels.size(), getCovariances(),
				Math.max(getValues().getMaxValue(), -getValues().getMinValue()));
	}

	/**
	 * Extracts the labels specified by the indices
	 *
	 * @param indices An integer list of indices
	 * @return List&lt;T&gt;
	 */
	final private List<H> getLabels(final int[] indices) {
		final List<H> res = new ArrayList<>();
		for (int i = 0; i < indices.length; ++i)
			res.add(getLabel(i));
		return Collections.unmodifiableList(res);
	}

}
