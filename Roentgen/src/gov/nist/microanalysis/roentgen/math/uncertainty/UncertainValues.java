package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.swing.IValueToColor;

/**
 * <p>
 * The UncertainValues class serves to store a collection of related variable
 * values along with the associated covariance/uncertainty matrix.
 * </p>
 * <p>
 * UncertainValues instances can be the input into calculations class or can be
 * the results of calculations performed with the {@link UncertaintyPropagator}
 * class.
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
public class UncertainValues //
		extends UncertainValuesBase {

	private final RealVector mValues;
	private final RealMatrix mCovariance;
	private int mHashCode;

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

	public static final UncertainValues NULL = new UncertainValues(Collections.emptyList());

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
		final double EPS = 1.0e-9, SREPS = Math.sqrt(EPS);
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
				final double max = Math.sqrt(covar.getEntry(c, c) * covar.getEntry(r, r));
				final double rr = entryRC / max;
				if (Double.isFinite(rr)) {
					if (((rr < -MAX_CORR) || (rr > MAX_CORR)) && (Math.abs(max) > 1.0e-6))
						throw new UVSOutOfRangeException(entryRC, -max, max, "Row=" + r + ", Col=" + c);
				} else
					covar.setEntry(r, c, 0.0);
			}
	}

	private static final RealMatrix buildCovariances(final RealVector variances, final RealMatrix corrCoeffs) {
		final int dim = variances.getDimension();
		final RealMatrix res = MatrixUtils.createRealMatrix(dim, dim);
		for (int r = 0; r < dim; ++r) {
			final double vr = variances.getEntry(r);
			assert vr >= 0.0;
			res.setEntry(r, r, vr);
			for (int c = r + 1; c < dim; ++c) {
				final double vc = variances.getEntry(c);
				final double cc = corrCoeffs.getEntry(r, c);
				assert cc >= -1.0 : "Correlation coefficient less than negative unity. " + cc;
				assert cc <= 1.0 : "Correlation coefficient greater than unity. " + cc;
				final double cov = vr * vc > 0.0 ? Math.sqrt(vr * vc) * cc : 0.0;
				res.setEntry(r, c, cov);
				res.setEntry(c, r, cov);
			}
		}
		return res;
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * covariance matrix.
	 *
	 * @param labels List&lt;Object&gt; A list of objects implementing hashCode()
	 *               and equals().
	 * @param vals   {@link RealVector}
	 * @param covar  {@link RealMatrix}
	 */
	public UncertainValues( //
			final List<? extends Object> labels, //
			final RealVector vals, //
			final RealMatrix covar //
	) {
		super(labels);
		if (vals.getDimension() != labels.size())
			throw new DimensionMismatchException(vals.getDimension(), labels.size());
		if (covar.getRowDimension() != labels.size())
			throw new DimensionMismatchException(covar.getRowDimension(), labels.size());
		if (covar.getColumnDimension() != labels.size())
			throw new DimensionMismatchException(covar.getColumnDimension(), labels.size());
		assert noReplicateLabels(labels);
		assert noCollectionLabels(labels);
		mValues = vals;
		validateCovariance(getLabels().size(), covar, Math.max(mValues.getMaxValue(), -mValues.getMinValue()));
		mCovariance = covar;
	}

	public UncertainValues(final Object[] labels, final double[] vals, final double[] cov) {
		this(Arrays.asList(labels), new ArrayRealVector(vals), new ArrayRealVector(cov));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels with zero
	 * values and NaN covariances.
	 *
	 * @param labels List&lt;Object&gt; A list of objects implementing hashCode()
	 *               and equals().
	 */
	public UncertainValues(//
			final List<? extends Object> labels //
	) {
		this(labels, new ArrayRealVector(labels.size()), NullableRealMatrix.build(labels.size(), labels.size()));
		mValues.set(Double.NaN);
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param labels    List&lt;Object&gt; A list of objects implementing hashCode()
	 *                  and equals().
	 * @param vals      {@link RealVector}
	 * @param variances {@link RealVector} covariances are zero.
	 */
	public UncertainValues( //
			final List<? extends Object> labels, //
			final RealVector vals, //
			final RealVector variances //
	) {
		this(labels, vals, MatrixUtils.createRealDiagonalMatrix(variances.toArray()));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param labels    List&lt;Object&gt; A list of objects implementing hashCode()
	 *                  and equals().
	 * @param vals      {@link RealVector}
	 * @param variances {@link RealVector} covariances are zero.
	 */
	public UncertainValues( //
			final List<? extends Object> labels, //
			final RealVector vals, //
			final RealVector variances, //
			final RealMatrix corrCoeffs) {
		this(labels, vals, buildCovariances(variances, corrCoeffs));
	}

	private static final List<Pair<? extends Object, ? extends Number>> convert(
			final Map<? extends Object, ? extends Number> vals) {
		final List<Pair<? extends Object, ? extends Number>> res = new ArrayList<>();
		for (final Map.Entry<? extends Object, ? extends Number> me : vals.entrySet())
			res.add(Pair.create(me.getKey(), me.getValue()));
		return res;
	}

	private static final ArrayList<? extends Object> extractLabels(
			final List<Pair<? extends Object, ? extends Number>> vals) {
		final ArrayList<Object> res = new ArrayList<>();
		for (final Pair<? extends Object, ? extends Number> pr : vals)
			res.add(pr.getFirst());
		return res;
	}

	private static final RealVector extractVals(final List<Pair<? extends Object, ? extends Number>> vals,
			final boolean val) {
		final RealVector res = new ArrayRealVector(vals.size());
		for (int i = 0; i < vals.size(); ++i) {
			final Pair<? extends Object, ? extends Number> pr = vals.get(i);
			final Number n = pr.getSecond();
			if (val)
				res.setEntry(i, n.doubleValue());
			else {
				if (n instanceof UncertainValue) {
					final UncertainValue uv = (UncertainValue) n;
					res.setEntry(i, uv.variance());
				}
			}
		}
		return res;
	}

	private UncertainValues(final List<Pair<? extends Object, ? extends Number>> vals, final boolean extra) {
		this(extractLabels(vals), extractVals(vals, true), extractVals(vals, false));
	}

	public UncertainValues(//
			final Map<? extends Object, ? extends Number> vals //
	) {
		this(convert(vals), true);
	}

	public static UncertainValues force(//
			final UncertainValuesBase base//
	) {
		if (base instanceof UncertainValues)
			return (UncertainValues) base;
		else
			return new UncertainValues(base.getLabels(), base.getValues(), base.getCovariances());
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances. "Homoskedastic-case"
	 *
	 * @param labels   List&lt;Object&gt; A list of objects implementing hashCode()
	 *                 and equals().
	 * @param vals     {@link RealVector}
	 * @param variance double (common to all vals)
	 */
	public UncertainValues(final List<? extends Object> labels, final RealVector vals, final double variance) {
		this(labels, vals, MatrixUtils.createRealIdentityMatrix(labels.size()).scalarMultiply(variance));
	}

	private static boolean noReplicateLabels(final List<? extends Object> labels) {
		final Set<Object> rep = new HashSet<>();
		for (final Object label : labels) {
			if (rep.contains(label)) {
				System.err.println("The label " + label + " is replicated.");
				return false;
			}
			rep.add(label);
		}
		return true;
	}

	private static boolean noCollectionLabels(final List<? extends Object> labels) {
		for (final Object label : labels) {
			if (label instanceof Collection) {
				System.err.println("The label " + label + " is a collection.");
				return false;
			}
		}
		return true;
	}

	/**
	 * Builds an {@link UncertainValues} object representing the specified labeled
	 * quantities as extracted from the list of {@link UncertainValues} objects.
	 *
	 * @param labels
	 * @param uvs
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public static UncertainValues build(final List<? extends Object> labels, final UncertainValues... uvs)
			throws ArgumentException {
		// Test that each requested label is defined once and only once.
		for (final Object label : labels) {
			int count = 0;
			for (final UncertainValues uv : uvs)
				if (uv.hasEntry(label))
					count++;
			if (count < 1)
				throw new ArgumentException(
						"The label " + label + " is not defined in one of the input UncertainValues.");
			if (count > 1)
				throw new ArgumentException(
						"The label " + label + " is muliply defined in one of the input UncertainValues.");
		}
		final UncertainValues res = new UncertainValues(labels);
		for (final UncertainValues uv : uvs)
			for (final Object label1 : labels)
				if (uv.hasEntry(label1)) {
					res.set(label1, uv.getEntry(label1), uv.getVariance(label1));
					for (final Object label2 : uv.getLabels())
						if (res.hasEntry(label2)) {
							final double cv = uv.getCovariance(label1, label2);
							if (cv != 0.0)
								res.setCovariance(label1, label2, cv);
						}
				}
		return res;
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
	public static UncertainValues zeroBut(final Object label, final UncertainValues input) {
		final UncertainValues res = new UncertainValues(input.getLabels(), input.getValues(), 0.0);
		final int idx = input.indexOf(label);
		res.setCovariance(idx, idx, input.getCovariance(idx, idx));
		return res;
	}

	/**
	 * Return an {@link UncertainValues} object with the same dimension and values
	 * as input except all the covariances except those associated with the labels
	 * are zeroed.
	 *
	 * @param labels
	 * @param input
	 * @return {@link UncertainValues}
	 */
	public static UncertainValues zeroBut(final Collection<? extends Object> labels, final UncertainValues input) {
		final UncertainValues res = new UncertainValues(input.getLabels(), input.getValues(), 0.0);
		final List<Integer> idxs = new ArrayList<>();
		for (final Object label : labels)
			idxs.add(input.indexOf(label));
		for (final int ridx : idxs)
			for (final int cidx : idxs)
				res.mCovariance.setEntry(ridx, cidx, input.getCovariance(ridx, cidx));
		return res;
	}

	/**
	 * Check that the indices are valid and not repeated. Uses assert rather than an
	 * Exception.
	 *
	 * @param indices
	 * @return true
	 */
	@Override
	public boolean assertIndices(final int[] indices) {
		for (int i = 0; i < indices.length; ++i) {
			assert indices[i] >= 0 : "Index[" + i + "] is less than zero.";
			final List<Object> labels = getLabels();
			assert indices[i] < labels.size() : "Index[" + i + "] is larger than the number of labels.";
			for (int j = i + 1; j < indices.length; ++j)
				assert indices[i] != indices[j] : "Duplicated index: Index[" + i + "] equals Index[" + j + "]";
		}
		return true;
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
	public static UncertainValuesBase propagate( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues input //
	) throws ArgumentException {
		return new UncertainValuesCalculator(nmjf, input);
	}

	/**
	 * Returns the UncertainValues that result from applying the function/Jacobian
	 * in <code>nmjf</code> to the input values/variances in <code>input</code>.
	 *
	 *
	 * @param nmjf
	 * @param inputs
	 * @param dinp   Delta in the input values in the same order as inputs...
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static UncertainValues propagateFiniteDifference( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues inputs, //
			final RealVector dinp //
	) throws ArgumentException {
		UncertainValuesCalculator uvc = new UncertainValuesCalculator(nmjf, inputs);
		final List<Object> labels = uvc.getInputs().getLabels();
		RealVector reordered = new ArrayRealVector(labels.size());
		for (int i = 0; i < reordered.getDimension(); ++i)
			reordered.setEntry(i, dinp.getEntry(inputs.indexOf(labels.get(i))));
		uvc.setCalculator(new UncertainValuesCalculator.FiniteDifference(reordered));
		return UncertainValues.force(uvc);
	}

	public static UncertainValues forceMinCovariance(final UncertainValues vals, final RealVector minCov) {
		assert vals.getDimension() == minCov.getDimension();
		final RealMatrix cov = vals.getCovariances().copy();
		for (int rc = 0; rc < vals.getDimension(); ++rc)
			if (cov.getEntry(rc, rc) < minCov.getEntry(rc))
				cov.setEntry(rc, rc, minCov.getEntry(rc));
		return new UncertainValues(vals.getLabels(), vals.getValues(), cov);
	}

	/**
	 * Similar to <code>propagate(...)</code> except uses a Monte Carlo-style
	 * evaluation rather than the Jacobian to propagate the uncertainties.
	 *
	 * @param nmjf
	 * @param input
	 * @param nEvals
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static UncertainValues propagateMonteCarlo(//
			final UncertainValuesCalculator uvc, //
			final int nEvals //
	) throws ArgumentException {
		return propagateMC(uvc.getFunction(), uvc.getInputs(), nEvals);
	}

	/**
	 * Similar to <code>propagate(...)</code> except uses a Monte Carlo-style
	 * evaluation rather than the Jacobian to propagate the uncertainties.
	 *
	 * @param nmjf
	 * @param inputs
	 * @param nEvals
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static UncertainValues propagateMC(//
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValuesBase inputs, //
			final int nEvals //
	) throws ArgumentException {
		UncertainValuesCalculator uvc = new UncertainValuesCalculator(nmjf, inputs);
		final MCPropagator mcp = new MCPropagator(uvc);
		return mcp.compute(nEvals);
	}

	/**
	 * Returns a {@link RealVector} containing the values associated with this
	 * object.
	 *
	 * @return {@link RealVector}
	 */
	@Override
	final public RealVector getValues() {
		return mValues;
	}

	/**
	 * Returns the matrix containing the covariances associate with this uncertainty
	 * calculation.
	 *
	 * @return RealMatrix
	 */
	@Override
	final public RealMatrix getCovariances() {
		return mCovariance;
	}

	/**
	 * Extracts a subset of the {@link UncertainValues} uvs associated with the
	 * specified labels into a new {@link UncertainValues} object. All labels must
	 * exists in uvs.
	 *
	 * @param labels
	 * @param uvs
	 * @return UncertainValues
	 */
	public static UncertainValues extract(final List<? extends Object> labels, final UncertainValuesBase uvs) {
		final RealVector vals = new ArrayRealVector(labels.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		final int[] idx = uvs.indices(labels);
		for (int ri = 0; ri < labels.size(); ++ri) {
			final int r = idx[ri];
			assert r >= 0 : //
			labels.get(r) + " is unavailable in UncertainValues.extract(...)";
			vals.setEntry(ri, uvs.getEntry(r));
			cov.setEntry(ri, ri, uvs.getCovariance(r, r));
			for (int ci = 0; ci < ri; ++ci) {
				final int c = idx[ci];
				final double cc = uvs.getCovariance(r, c);
				cov.setEntry(ri, ci, cc);
				cov.setEntry(ci, ri, cc);
			}
		}
		return new UncertainValues(labels, vals, cov);
	}

	/**
	 * Extracts a subset of the {@link UncertainValues} uvs associated with the
	 * specified labels into a new {@link UncertainValues} object. All labels must
	 * exists in auvs - if an item exists in more than one {@link UncertainValues}
	 * input then only those values and covariances from the first matching are
	 * extracted.
	 *
	 * @param labels
	 * @param inputs A list containing an arbitrary number of
	 *               {@link UncertainValues} objects
	 * @return UncertainValues
	 */
	public static UncertainValues extract(final List<? extends Object> labels, final UncertainValuesBase... inputs) {
		final RealVector vals = new ArrayRealVector(labels.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		for (int ri = 0; ri < labels.size(); ++ri) {
			boolean found = false;
			for (final UncertainValuesBase uvs : inputs) {
				final int r = uvs.indexOf(labels.get(ri));
				if (r >= 0) {
					found = true;
					vals.setEntry(ri, uvs.getEntry(r));
					cov.setEntry(ri, ri, uvs.getCovariance(r, r));
					for (int ci = r + 1; ci < labels.size(); ++ci) {
						final int c = uvs.indexOf(labels.get(ci));
						if (c >= 0) {
							final double cc = uvs.getCovariance(r, c);
							cov.setEntry(ri, ci, cc);
							cov.setEntry(ci, ri, cc);
						}
					}
				}
			}
			assert found : labels.get(ri) + " is unavailable in any input to UncertainValues.extract(...)";
		}
		return new UncertainValues(labels, vals, cov);
	}

	/**
	 * Creates a bitmap that represents the difference between uncertainties
	 * associated with these two sets of UncertainValues. The difference between the
	 * covariances is plotted.
	 *
	 * @param uvs1
	 * @param uvs2
	 * @param corr
	 * @param pixDim
	 * @return BufferedImage
	 */
	public static BufferedImage compareAsBitmap(final UncertainValuesBase uvs1, final UncertainValuesBase uvs2,
			final IValueToColor corr, final int pixDim) {
		if (uvs1.getDimension() != uvs2.getDimension())
			throw new DimensionMismatchException(uvs2.getDimension(), uvs1.getDimension());
		final int dim = uvs1.getDimension();
		final boolean[] disp = new boolean[dim];
		final UncertainValues uvs2r = UncertainValues.extract(uvs1.getLabels(), uvs2);
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

	/**
	 * <p>
	 * Copies those entries labeled in <i>from</i> that are also labeled in
	 * <i>to</i>. Useful for replicating the variances and covariances in
	 * <i>from</i> in <i>to</i>.
	 * </p>
	 * <p>
	 * Will overwrite values in <i>to</i>
	 * </p>
	 *
	 * @param from
	 * @param to
	 */
	public static void copy(final UncertainValuesBase from, final UncertainValues to) {
		final List<? extends Object> fromLabels = from.getLabels();
		final Map<Object, Integer> toMap = new HashMap<>();
		final Map<Object, Integer> fromMap = new HashMap<>();
		for (int i = 0; i < fromLabels.size(); ++i) {
			final Object label = fromLabels.get(i);
			final int p = to.indexOf(label);
			if (p >= 0) {
				toMap.put(label, Integer.valueOf(p));
				fromMap.put(label, i);
			}
		}
		final RealVector toValues = to.mValues;
		final RealMatrix toCovs = to.mCovariance;
		for (final Map.Entry<Object, Integer> me1 : toMap.entrySet()) {
			final int p1 = me1.getValue().intValue();
			toValues.setEntry(p1, from.getEntry(me1.getKey()));
			for (final Map.Entry<Object, Integer> me2 : toMap.entrySet()) {
				final int p2 = me2.getValue().intValue();
				toCovs.setEntry(p1, p2, from.getCovariance(me1.getKey(), me2.getKey()));
			}
		}
	}

	/**
	 * Initializes the value and variance associated with the specified label.
	 *
	 * @param label
	 * @param val
	 * @param var
	 */
	final public void set(final Object label, final double val, final double var) {
		final int p = indexOf(label);
		mValues.setEntry(p, val);
		mCovariance.setEntry(p, p, var);
	}

	/**
	 * Sets the covariance value associated with the specified values (both i,j and
	 * j,i)
	 *
	 * @param label1
	 * @param label2
	 * @param cov
	 */
	final public void setCovariance(final Object label1, final Object label2, final double cov) {
		final int p1 = indexOf(label1), p2 = indexOf(label2);
		setCovariance(p1, p2, cov);
	}

	/**
	 * Sets the covariance value associated with the specified values (both i,j and
	 * j,i)
	 *
	 * @param p1
	 * @param p2
	 * @param cov
	 */
	final public void setCovariance(final int p1, final int p2, final double cov) {
		mCovariance.setEntry(p1, p2, cov);
		mCovariance.setEntry(p2, p1, cov);
	}

	/**
	 * Initializes the specified value and variance for the specified label.
	 *
	 * @param label
	 * @param uv    (uses the doubleValue() and variance())
	 */
	final public void set(final Object label, final UncertainValue uv) {
		final int p = indexOf(label);
		mValues.setEntry(p, uv.doubleValue());
		mCovariance.setEntry(p, p, uv.variance());
	}

	@Override
	public int hashCode() {
		if (mHashCode == 0)
			mHashCode = Objects.hash(super.hashCode(), mValues, mCovariance);
		return mHashCode;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValues other = (UncertainValues) obj;
		return super.equals(obj) && //
				Objects.equals(mCovariance, other.mCovariance) && //
				Objects.equals(mValues, other.mValues);
	}

	public UncertainValues copy() {
		return new UncertainValues(getLabels(), mValues.copy(), mCovariance.copy());
	}

	public static boolean testEquality(final UncertainValues uvs1, final UncertainValues uvs2) {
		if (uvs1.getDimension() != uvs2.getDimension())
			return false;
		for (final Object label1 : uvs1.getLabels()) {
			if (uvs1.getEntry(label1) != uvs2.getEntry(label1))
				return false;
			for (final Object label2 : uvs1.getLabels())
				if (uvs1.getCovariance(label1, label2) != uvs2.getCovariance(label1, label2))
					return false;
		}
		return true;
	}

	@Override
	public String toString() {
		final List<? extends Object> labels = getLabels();
		final List<Object> labelSub = new ArrayList<>();
		for (int i = 0; (i < labels.size()) && (i < 5); ++i)
			labelSub.add(labels.get(i));
		final String lblStr = labelSub.toString();
		return "UVS[" + lblStr.substring(1, lblStr.length() - 1)
				+ (labelSub.size() < labels.size() ? "+" + (labels.size() - labelSub.size()) + " more" : "") + "]";
	}
}
