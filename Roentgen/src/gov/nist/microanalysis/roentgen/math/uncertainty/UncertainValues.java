package gov.nist.microanalysis.roentgen.math.uncertainty;

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
public class UncertainValues<H> //
		extends UncertainValuesBase<H> {

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
	/**
	 * Builds an {@link UncertainValues} object representing the specified labeled
	 * quantities as extracted from the list of {@link UncertainValues} objects.
	 *
	 * @param labels
	 * @param uvs
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public static <J> UncertainValues<J> build( //
			final List<J> labels, //
			@SuppressWarnings("unchecked") final UncertainValues<J>... uvs //
	) throws ArgumentException {
		// Test that each requested label is defined once and only once.
		for (final J label : labels) {
			int count = 0;
			for (final UncertainValues<J> uv : uvs)
				if (uv.hasEntry(label))
					count++;
			if (count < 1)
				throw new ArgumentException(
						"The label " + label + " is not defined in one of the input UncertainValues.");
			if (count > 1)
				throw new ArgumentException(
						"The label " + label + " is muliply defined in one of the input UncertainValues.");
		}
		final UncertainValues<J> res = new UncertainValues<J>(labels);
		for (final UncertainValues<J> uv : uvs)
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
	public static <H> void copy(final UncertainValuesBase<H> from, final UncertainValues<H> to) {
		final List<H> fromLabels = from.getLabels();
		final Map<H, Integer> toMap = new HashMap<>();
		final Map<H, Integer> fromMap = new HashMap<>();
		for (int i = 0; i < fromLabels.size(); ++i) {
			final H label = fromLabels.get(i);
			final int p = to.indexOf(label);
			if (p >= 0) {
				toMap.put(label, Integer.valueOf(p));
				fromMap.put(label, i);
			}
		}
		final RealVector toValues = to.mValues;
		final RealMatrix toCovs = to.mCovariance;
		for (final Map.Entry<H, Integer> me1 : toMap.entrySet()) {
			final int p1 = me1.getValue().intValue();
			toValues.setEntry(p1, from.getEntry(me1.getKey()));
			for (final Map.Entry<H, Integer> me2 : toMap.entrySet()) {
				final int p2 = me2.getValue().intValue();
				toCovs.setEntry(p1, p2, from.getCovariance(me1.getKey(), me2.getKey()));
			}
		}
	}

	/**
	 * Extracts a subset of the {@link UncertainValues} uvs associated with the
	 * specified labels into a new {@link UncertainValues} object. All labels must
	 * exists in uvs.
	 *
	 * @param labels
	 * @param luvs
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static <J> UncertainValues<J> extract(final List<J> labels, final List<UncertainValuesBase<?>> luvs)
			throws ArgumentException {
		final Set<J> remaining = new HashSet<J>(labels);
		final RealVector vals = new ArrayRealVector(labels.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		for (final UncertainValuesBase<?> uvs : luvs) {
			final int[] idx = uvs.indices(labels);
			for (int ri = 0; ri < idx.length; ++ri) {
				if (idx[ri] != -1) {
					remaining.remove(labels.get(ri));
					final int r = idx[ri];
					vals.setEntry(ri, uvs.getEntry(r));
					cov.setEntry(ri, ri, uvs.getCovariance(r, r));
					for (int ci = 0; ci < ri; ++ci) {
						if (idx[ci] != -1) {
							final int c = idx[ci];
							final double cc = uvs.getCovariance(r, c);
							cov.setEntry(ri, ci, cc);
							cov.setEntry(ci, ri, cc);
						}
					}
				}
			}
		}
		if (remaining.size() > 0)
			throw new ArgumentException("The labels " + remaining + " are missing in extract(...)");
		return new UncertainValues<J>(labels, vals, cov);
	}

	/**
	 * Extracts a subset of the {@link UncertainValues} uvs associated with the
	 * specified labels into a new {@link UncertainValues} object. All labels must
	 * exists in uvs.
	 *
	 * @param labels
	 * @param luvs
	 * @return UncertainValues
	 * @throws ArgumentException
	 */
	public static <J> UncertainValues<J> extract(final List<J> labels, final UncertainValuesBase<?> uvs) //
			throws ArgumentException {
		return extract(labels, Collections.singletonList(uvs));
	}

	@SuppressWarnings("unchecked")
	public static <J> UncertainValues<J> force(//
			final UncertainValuesBase<? extends J> base//
	) {
		if (base instanceof UncertainValues)
			return (UncertainValues<J>) base;
		else {
			// Compute getCovariances() first since this forces getValues() to be calculated
			// simultaneously.
			final RealMatrix covs = base.getCovariances();
			return new UncertainValues<J>(base.getLabels(), base.getValues(), covs);

		}
	}

	public static <J> UncertainValues<J> forceMinCovariance(final UncertainValues<J> vals, final RealVector minCov) {
		assert vals.getDimension() == minCov.getDimension();
		final RealMatrix cov = vals.getCovariances().copy();
		for (int rc = 0; rc < vals.getDimension(); ++rc)
			if (cov.getEntry(rc, rc) < minCov.getEntry(rc))
				cov.setEntry(rc, rc, minCov.getEntry(rc));
		return new UncertainValues<J>(vals.getLabels(), vals.getValues(), cov);
	}

	public static <H> boolean testEquality(final UncertainValues<H> uvs1, final UncertainValues<H> uvs2) {
		if (uvs1.getDimension() != uvs2.getDimension())
			return false;
		for (final H label1 : uvs1.getLabels()) {
			if (uvs1.getEntry(label1) != uvs2.getEntry(label1))
				return false;
			for (final H label2 : uvs1.getLabels())
				if (uvs1.getCovariance(label1, label2) != uvs2.getCovariance(label1, label2))
					return false;
		}
		return true;
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
	public static <J> UncertainValues<J> zeroBut(final Collection<? extends J> labels, final UncertainValues<J> input) {
		final UncertainValues<J> res = new UncertainValues<J>(input.getLabels(), input.getValues(), 0.0);
		final List<Integer> idxs = new ArrayList<>();
		for (final J label : labels)
			idxs.add(input.indexOf(label));
		for (final int ridx : idxs)
			for (final int cidx : idxs)
				res.mCovariance.setEntry(ridx, cidx, input.getCovariance(ridx, cidx));
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
	public static <J> UncertainValues<J> zeroBut(final J label, final UncertainValues<J> input) {
		final UncertainValues<J> res = new UncertainValues<J>(input.getLabels(), input.getValues(), 0.0);
		final int idx = input.indexOf(label);
		res.setCovariance(idx, idx, input.getCovariance(idx, idx));
		return res;
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

	private static final <J> List<Pair<J, ? extends Number>> convert(final Map<J, ? extends Number> vals) {
		final List<Pair<J, ? extends Number>> res = new ArrayList<>();
		for (final Map.Entry<J, ? extends Number> me : vals.entrySet())
			res.add(Pair.create(me.getKey(), me.getValue()));
		return res;
	}

	private static final <J> ArrayList<J> extractLabels(final List<Pair<J, ? extends Number>> vals) {
		final ArrayList<J> res = new ArrayList<>();
		for (final Pair<J, ? extends Number> pr : vals)
			res.add(pr.getFirst());
		return res;
	}

	private static final <J> RealVector extractVals(final List<Pair<J, ? extends Number>> vals, final boolean val) {
		final RealVector res = new ArrayRealVector(vals.size());
		for (int i = 0; i < vals.size(); ++i) {
			final Pair<J, ? extends Number> pr = vals.get(i);
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

	private static <J> boolean noCollectionLabels(final List<J> labels) {
		for (final J label : labels) {
			if (label instanceof Collection) {
				System.err.println("The label " + label + " is a collection.");
				return false;
			}
		}
		return true;
	}

	private static <J> boolean noReplicateLabels(final List<J> labels) {
		final Set<J> rep = new HashSet<>();
		for (final J label : labels) {
			if (rep.contains(label)) {
				System.err.println("The label " + label + " is replicated.");
				return false;
			}
			rep.add(label);
		}
		return true;
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
		final double EPS = 1.0e8, SREPS = Math.sqrt(EPS);
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

	private final RealVector mValues;

	private final RealMatrix mCovariance;

	private int mHashCode;

	public UncertainValues(final H[] labels, final double[] vals, final double[] cov) {
		this(Arrays.asList(labels), new ArrayRealVector(vals), new ArrayRealVector(cov));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * covariance matrix.
	 *
	 * @param labels List&lt;? extends H&gt; A list of labels of size L.
	 * @param vals   {@link RealVector} A vector of dimension L
	 * @param covar  {@link RealMatrix} A matrix of dimension L x L
	 */
	public UncertainValues( //
			final List<? extends H> labels, //
			final RealVector vals, //
			final RealMatrix covar //
	) {
		super(new ArrayList<>(labels));
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

	/**
	 * Constructs a UncertainValues object based on the specified labels with zero
	 * values and NaN covariances.
	 *
	 * @param labels List&lt;H&gt; A list of objects implementing hashCode() and
	 *               equals().
	 */
	public UncertainValues(//
			final List<H> labels //
	) {
		this(labels, new ArrayRealVector(labels.size()), NullableRealMatrix.build(labels.size(), labels.size()));
		mValues.set(Double.NaN);
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances. "Homoskedastic-case"
	 *
	 * @param labels   List&lt;H&gt; A list of objects implementing hashCode() and
	 *                 equals().
	 * @param vals     {@link RealVector}
	 * @param variance double (common to all vals)
	 */
	public UncertainValues(final List<H> labels, final RealVector vals, final double variance) {
		this(labels, vals, MatrixUtils.createRealIdentityMatrix(labels.size()).scalarMultiply(variance));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param labels    List&lt;H&gt; A list of objects implementing hashCode() and
	 *                  equals().
	 * @param vals      {@link RealVector}
	 * @param variances {@link RealVector} covariances are zero.
	 */
	public UncertainValues( //
			final List<H> labels, //
			final RealVector vals, //
			final RealVector variances //
	) {
		this(labels, vals, MatrixUtils.createRealDiagonalMatrix(variances.toArray()));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param labels    List&lt;H&gt; A list of objects implementing hashCode() and
	 *                  equals().
	 * @param vals      {@link RealVector}
	 * @param variances {@link RealVector} covariances are zero.
	 */
	public UncertainValues( //
			final List<H> labels, //
			final RealVector vals, //
			final RealVector variances, //
			final RealMatrix corrCoeffs) {
		this(labels, vals, buildCovariances(variances, corrCoeffs));
	}

	public UncertainValues(//
			final Map<H, ? extends Number> vals //
	) {
		this(convert(vals), true);
	}

	private UncertainValues(final List<Pair<H, ? extends Number>> vals, final boolean extra) {
		this(extractLabels(vals), extractVals(vals, true), extractVals(vals, false));
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
			final List<H> labels = getLabels();
			assert indices[i] < labels.size() : "Index[" + i + "] is larger than the number of labels.";
			for (int j = i + 1; j < indices.length; ++j)
				assert indices[i] != indices[j] : "Duplicated index: Index[" + i + "] equals Index[" + j + "]";
		}
		return true;
	}

	public UncertainValues<H> copy() {
		return new UncertainValues<H>(getLabels(), mValues.copy(), mCovariance.copy());
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValues<?> other = (UncertainValues<?>) obj;
		return super.equals(obj) && //
				Objects.equals(mCovariance, other.mCovariance) && //
				Objects.equals(mValues, other.mValues);
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
	 * Returns a {@link RealVector} containing the values.
	 *
	 * @return {@link RealVector}
	 */
	@Override
	final public RealVector getValues() {
		return mValues;
	}

	@Override
	public int hashCode() {
		if (mHashCode == 0)
			mHashCode = Objects.hash(super.hashCode(), mValues, mCovariance);
		return mHashCode;
	}

	/**
	 * Initializes the value and variance associated with the specified label.
	 *
	 * @param label
	 * @param val
	 * @param var
	 */
	final public void set(final H label, final double val, final double var) {
		final int p = indexOf(label);
		mValues.setEntry(p, val);
		mCovariance.setEntry(p, p, var);
	}

	/**
	 * Initializes the specified value and variance for the specified label.
	 *
	 * @param label
	 * @param uv    (uses the doubleValue() and variance())
	 */
	final public void set(final H label, final UncertainValue uv) {
		final int p = indexOf(label);
		mValues.setEntry(p, uv.doubleValue());
		mCovariance.setEntry(p, p, uv.variance());
	}

	/**
	 * Sets the covariance value associated with the specified values (both i,j and
	 * j,i)
	 *
	 * @param label1
	 * @param label2
	 * @param cov
	 */
	final public void setCovariance(final H label1, final H label2, final double cov) {
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
}
