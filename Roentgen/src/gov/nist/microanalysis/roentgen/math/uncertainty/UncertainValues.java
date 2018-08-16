package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Transforms;
import com.duckandcover.html.Table.Item;
import com.google.common.base.Objects;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.swing.IValueToColor;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

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
 * To facilicate bookkeeping, the variables within the object are labeled with
 * Object derived tags. Make sure that all the classes implementing labels
 * implement hashCode() and equals() as the labels are stored in hash maps. The
 * entries in the values vector and covariance matrix can be indexed by tag or
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
public class UncertainValues implements IToHTML {

	private static final double MAX_CORR = 1.00000001;
	private final RealVector mValues;
	private final List<Object> mTags;
	private final RealMatrix mCovariance;

	private final Map<Object, Integer> mIndex;

	/**
	 * Checks that the covariance matrix is m x m, the variances are non-negative,
	 * the covariances are bounded correlation coefficients between -1<r<1.
	 *
	 * @param m
	 * @param covar
	 */
	public static void validateCovariance(final int m, final RealMatrix covar, final double maxVal) {
		if (covar.getRowDimension() != m)
			throw new DimensionMismatchException(covar.getRowDimension(), m);
		if (covar.getColumnDimension() != m)
			throw new DimensionMismatchException(covar.getColumnDimension(), m);
		final double EPS = 1.0e-15;
		for (int r = 0; r < covar.getRowDimension(); ++r) {
			final double entryRR = covar.getEntry(r, r);
			if (entryRR < 0.0) {
				// This is necessary because of subtle rounding errors
				if (entryRR > -EPS * maxVal)
					covar.setEntry(r, r, 0.0);
				else
					throw new OutOfRangeException(entryRR, 0.0, Double.MAX_VALUE);
			}
		}
		for (int r = 0; r < covar.getRowDimension(); ++r)
			for (int c = r + 1; c < covar.getColumnDimension(); ++c) {
				final double entryRC = covar.getEntry(r, c);
				if (Math.abs(entryRC - covar.getEntry(c, r)) > 1.0e-6 * Math.abs(entryRC))
					throw new OutOfRangeException(entryRC, covar.getEntry(c, r), covar.getEntry(c, r));
				final double max = Math.sqrt(covar.getEntry(c, c) * covar.getEntry(r, r));
				final double rr = entryRC / max;
				if ((Math.abs(entryRC) > 1.0e-15) && ((rr > MAX_CORR) || (rr < -MAX_CORR)))
					throw new OutOfRangeException(entryRC, -max, max);
			}
	}

	/**
	 * Constructs a UncertainValues object based on the specified tags, values and
	 * covariance matrix.
	 *
	 * @param tags  List&lt;Object&gt; A list of objects implementing hashCode() and
	 *              equals().
	 * @param vals  {@link RealVector}
	 * @param covar {@link RealMatrix}
	 */
	public UncertainValues(final List<? extends Object> tags, final RealVector vals, final RealMatrix covar) {
		if (vals.getDimension() != tags.size())
			throw new DimensionMismatchException(covar.getRowDimension(), tags.size());
		if (covar.getRowDimension() != tags.size())
			throw new DimensionMismatchException(covar.getRowDimension(), tags.size());
		if (covar.getColumnDimension() != tags.size())
			throw new DimensionMismatchException(covar.getColumnDimension(), tags.size());
		mTags = Collections.unmodifiableList(tags);
		final HashMap<Object, Integer> index = new HashMap<>();
		for (int i = 0; i < mTags.size(); ++i)
			index.put(mTags.get(i), Integer.valueOf(i));
		mIndex = Collections.unmodifiableMap(index);
		mValues = vals;
		validateCovariance(mTags.size(), covar, Math.max(mValues.getMaxValue(), -mValues.getMinValue()));
		mCovariance = covar;
	}

	/**
	 * Constructs a UncertainValues object based on the specified tags with zero
	 * values and NaN covariances.
	 *
	 * @param tags List&lt;Object&gt; A list of objects implementing hashCode() and
	 *             equals().
	 */
	public UncertainValues(final List<? extends Object> tags) {
		this(tags, new ArrayRealVector(tags.size()), new Array2DRowRealMatrix(tags.size(), tags.size()));
		mValues.set(Double.NaN);
	}

	/**
	 * Constructs a UncertainValues object based on the specified tags, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param tags      List&lt;Object&gt; A list of objects implementing hashCode()
	 *                  and equals().
	 * @param vals      {@link RealVector}
	 * @param variances {@link RealVector} covariances are zero.
	 */
	public UncertainValues(final List<? extends Object> tags, final RealVector vals, final RealVector variances) {
		this(tags, vals, MatrixUtils.createRealDiagonalMatrix(variances.toArray()));
	}

	private static final RealVector extractValues(final Map<? extends Object, UncertainValue> vals) {
		final double[] d = new double[vals.size()];
		int i = 0;
		for (final Object key : vals.keySet()) {
			d[i] = vals.get(key).doubleValue();
			++i;
		}
		return new ArrayRealVector(d);
	}

	private static final RealVector extractVariances(final Map<? extends Object, UncertainValue> vals) {
		final double[] d = new double[vals.size()];
		int i = 0;
		for (final Object key : vals.keySet()) {
			d[i] = vals.get(key).variance();
			++i;
		}
		return new ArrayRealVector(d);
	}

	public UncertainValues(final Map<? extends Object, UncertainValue> vals) {
		this(Arrays.asList(vals.keySet().toArray()), extractValues(vals), extractVariances(vals));
	}

	/**
	 * Constructs a UncertainValues object based on the specified tags, values and
	 * variances. "Homoskedastic-case"
	 *
	 * @param tags     List&lt;Object&gt; A list of objects implementing hashCode()
	 *                 and equals().
	 * @param vals     {@link RealVector}
	 * @param variance double (common to all vals)
	 */
	public UncertainValues(final List<? extends Object> tags, final RealVector vals, final double variance) {
		this(tags, vals, MatrixUtils.createRealIdentityMatrix(tags.size()).scalarMultiply(variance));
	}

	public static UncertainValues build(final List<? extends Object> tags, final UncertainValues... uvs)
			throws ArgumentException {
		// Test that each requested tag is defined once and only once.
		for (final Object tag : tags) {
			int count = 0;
			for (final UncertainValues uv : uvs)
				if (uv.hasEntry(tag))
					count++;
			if (count < 1)
				throw new ArgumentException("The tag " + tag + " is not defined in one of the input UncertainValues.");
			if (count > 1)
				throw new ArgumentException(
						"The tag " + tag + " is muliply defined in one of the input UncertainValues.");
		}
		final UncertainValues res = new UncertainValues(tags);
		for (final UncertainValues uv : uvs)
			for (final Object tag1 : tags)
				if (uv.hasEntry(tag1)) {
					res.set(tag1, uv.getEntry(tag1), uv.getVariance(tag1));
					for (final Object tag2 : uv.getTags())
						if (res.hasEntry(tag2)) {
							final double cv = uv.getCovariance(tag1, tag2);
							if (cv != 0.0)
								res.setCovariance(tag1, tag2, cv);
						}
				}
		return res;
	}

	/**
	 * Return an {@link UncertainValues} object with the same dimension and values
	 * as input except all the covariances except those associated with tag are
	 * zeroed.
	 *
	 * @param tag
	 * @param input
	 * @return {@link UncertainValues}
	 */
	public static UncertainValues zeroBut(final Object tag, final UncertainValues input) {
		final UncertainValues res = new UncertainValues(input.getTags(), input.getValues(), 0.0);
		final int idx = input.indexOf(tag);
		for (int i = 0; i < res.getDimension(); ++i) {
			res.mCovariance.setEntry(i, idx, input.getCovariance(i, idx));
			res.mCovariance.setEntry(idx, i, input.getCovariance(idx, i));
		}
		return res;
	}

	public static UncertainValues propagate(final NamedMultivariateJacobianFunction nmjf, final UncertainValues vals) {
		assert nmjf.getInputTags().equals(vals.getTags());
		final Pair<RealVector, RealMatrix> eval = nmjf.evaluate(vals.getValues());
		final RealMatrix jac = eval.getSecond();
		return new UncertainValues(nmjf.getOutputTags(), eval.getFirst(),
				jac.multiply(vals.getCovariances().multiply(jac.transpose())));
	}

	public static UncertainValues forceMinCovariance(final UncertainValues vals, final RealVector minCov) {
		assert vals.getDimension() == minCov.getDimension();
		final RealMatrix cov = vals.getCovariances().copy();
		for (int rc = 0; rc < vals.getDimension(); ++rc)
			if (cov.getEntry(rc, rc) < minCov.getEntry(rc))
				cov.setEntry(rc, rc, minCov.getEntry(rc));
		return new UncertainValues(vals.getTags(), vals.getValues(), cov);
	}

	public static UncertainValues propagateMC(final NamedMultivariateJacobianFunction nmjf, final UncertainValues vals,
			final int nEvals) {
		final MultivariateNormalDistribution mnd = new MultivariateNormalDistribution(vals.mValues.toArray(),
				vals.mCovariance.getData());
		final List<RealVector> evals = new ArrayList<>();
		for (int eval = 0; eval < nEvals; ++eval) {
			final RealVector pt = MatrixUtils.createRealVector(mnd.sample());
			evals.add(nmjf.evaluate(pt).getFirst());
		}
		final int len = nmjf.getOutputDimension();
		final RealVector sum = new ArrayRealVector(len);
		for (final RealVector eval : evals)
			sum.combineToSelf(1.0, 1.0, eval);
		final RealVector avg = sum.mapDivide(evals.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(len, len);
		for (int r = 0; r < len; r++)
			for (int c = r; c < len; c++) {
				double tmp = 0.0;
				for (final RealVector eval : evals)
					tmp += (eval.getEntry(r) - avg.getEntry(r)) * (eval.getEntry(c) - avg.getEntry(c));
				tmp /= evals.size();
				cov.setEntry(r, c, tmp);
				cov.setEntry(c, r, tmp);
			}
		return new UncertainValues(nmjf.getOutputTags(), avg, cov);
	}
	
	
	/**
	 * Extract an array of values from this UncertainValue object for the 
	 * specified list of tags
	 * 
	 * @param tags List&lt;? extends Object&gt;
	 * @return RealVector
	 */
	public RealVector extractValues(List<? extends Object> tags) {
		RealVector res = new ArrayRealVector(tags.size());
		int i=0;
		for(Object tag : tags) {
			res.setEntry(i, getEntry(tag));
			++i;
		}
		return res;
	}

	public static UncertainValues implicit(final NamedMultivariateJacobian jY,
			final NamedMultivariateJacobian jX, UncertainValues uv) {
		final Pair<RealVector, RealMatrix> jXe = jX.evaluate(uv.extractValues(jX.getInputTags()));
		final Pair<RealVector, RealMatrix> jYe = jY.evaluate(uv.extractValues(jY.getInputTags()));
		// From page 64 of NPL Report DEM-ES-011
		// Jy.Vy.JyT = Jx.Vx.JxT solve for Vy.
		final boolean NAIVE = true;
		if (NAIVE) {
			// Order the 
			UncertainValues ouvs = UncertainValues.extract(jX.getInputTags(), uv);
			final RealMatrix right = jXe.getSecond().multiply(ouvs.mCovariance.multiply(jXe.getSecond().transpose()));
			final RealMatrix ijYe=MatrixUtils.inverse(jYe.getSecond());
			return new UncertainValues(jY.getOutputTags(), jYe.getFirst(), ijYe.multiply(right.multiply(ijYe.transpose())));
		} else {
			// 1. Form the Cholesky factor Rx of Vx, i.e., the upper triangular matrix such
			// that RxT.Rx = Vx.
			final CholeskyDecomposition chyd = new CholeskyDecomposition(uv.getCovariances());
			final RealMatrix Rx = chyd.getL();
			// 2. Factor Jx as the product Jx = Qx.Ux, where Qx is an orthogonal matrix and
			// Ux is upper triangular.
			final QRDecomposition qrd = new QRDecomposition(jXe.getSecond());
			final RealMatrix Qx = qrd.getQ(), Ux = qrd.getR();
			// 3. Factor Jy as the product Jy = Ly.Uy, where Ly is lower triangular and Uy
			// is upper triangular.
			final LUDecomposition lud = new LUDecomposition(jYe.getSecond());
			// ????What about pivoting?????
			final RealMatrix Ly = lud.getL(), Uy = lud.getU();
			// 4. Solve the matrix equation UyT.M1 = I for M1.
			final RealMatrix M1 = MatrixUtils.inverse(Uy.transpose());
			// 5. Solve LyT.M2 = M1 for M2,
			final RealMatrix M2 = MatrixUtils.inverse(Ly.transpose()).multiply(M1);
			// 6. Form M3 = QxT.M2.
			final RealMatrix M3 = Qx.transpose().multiply(M2);
			// 7. Form M4 = UxT.M3.
			final RealMatrix M4 = Ux.transpose().multiply(M3);
			// 8. Form M = Rx.M4.
			final RealMatrix M = Rx.multiply(M4);
			// 9. Orthogonally triangularize M to give the upper triangular matrix R.
			final RealMatrix R = new QRDecomposition(M).getR();
			// 10. Form Vy = RT.R.
			return new UncertainValues(jY.getOutputTags(), jYe.getFirst(), (R.transpose()).multiply(R));
		}
	}

	/**
	 * Returns a {@link RealVector} containing the values associated with this
	 * object.
	 *
	 * @return {@link RealVector}
	 */
	final public RealVector getValues() {
		return mValues;
	}

	final public RealMatrix getCovariances() {
		return mCovariance;
	}

	/**
	 * Extracts a subset of the {@link UncertainValues} uvs associated with the
	 * specified tags into a new {@link UncertainValues} object. All tags must
	 * exists in uvs.
	 *
	 * @param tags
	 * @param uvs
	 * @return UncertainValues
	 */
	public static UncertainValues extract(final List<? extends Object> tags, final UncertainValues uvs) {
		final RealVector vals = new ArrayRealVector(tags.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(tags.size(), tags.size());
		for (int ri = 0; ri < tags.size(); ++ri) {
			final int r = uvs.indexOf(tags.get(ri));
			vals.setEntry(ri, uvs.getEntry(r));
			cov.setEntry(ri, ri, uvs.getCovariance(r, r));
			for (int ci = r + 1; ci < tags.size(); ++ci) {
				final int c = uvs.indexOf(tags.get(ci));
				final double cc = uvs.getCovariance(r, c);
				cov.setEntry(ri, ci, cc);
				cov.setEntry(ci, ri, cc);
			}
		}
		return new UncertainValues(tags, vals, cov);

	}

	/**
	 * Returns an ordered list of the Object tags associated with the values and
	 * covariances in this object.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<Object> getTags() {
		return mTags;
	}

	/**
	 * Is there a value and covariances associated with the specified tag?
	 *
	 * @param tag
	 * @return boolean
	 */
	final public boolean hasEntry(final Object tag) {
		return mTags.contains(tag);
	}

	/**
	 * Returns the tag at the p-th entry in the values vector and in the p-th row
	 * and column in the covariance matrix.
	 *
	 * @param p
	 * @return Object
	 */
	final public Object getTag(final int p) {
		return mTags.get(p);
	}

	/**
	 * Returns the index associated with the specified tag or -1 if not found.
	 *
	 * @param p
	 * @return Object
	 */
	final public int indexOf(final Object tag) {
		return mIndex.getOrDefault(tag, Integer.valueOf(-1));
	}

	/**
	 * Returns the value associated with the entry associated with tag.
	 *
	 * @param tag Object
	 * @return double
	 */
	final public double getEntry(final Object tag) {
		final int p = indexOf(tag);
		assert p >= 0 : "Tag " + tag.toString() + " missing.";
		return mValues.getEntry(p);
	}

	/**
	 * Returns the value associated with the entry associated with tag if the tag is
	 * associated with a value or returns defVal otherwise. Useful for a quantity
	 * that may or may not have an associated uncertainty.
	 *
	 * @param tag
	 * @param defVal
	 * @return double
	 */
	final public double getEntryWithDefault(final Object tag, final double defVal) {
		final int p = indexOf(tag);
		return p >= 0 ? mValues.getEntry(p) : defVal;
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified tag.
	 * The value is the same as getEntry(...) and the uncertainty is the on-diagonal
	 * covariance matrix entry associated with that entry.
	 *
	 * @param tag
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValue(final Object tag) {
		final int p = indexOf(tag);
		return new UncertainValue(getEntry(p), tag, Math.sqrt(mCovariance.getEntry(p, p)));
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified tag
	 * (if defined). The value is the same as getEntry(...) and the uncertainty is
	 * the on-diagonal covariance matrix entry associated with that entry. If the
	 * tag isn't defined, then the value defVal will be returned.
	 *
	 * @param tag
	 * @param defVal
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValueWithDefault(final Object tag, final UncertainValue defVal) {
		final int p = indexOf(tag);
		return p >= 0 ? new UncertainValue(getEntry(p), tag, mCovariance.getEntry(p, p)) : defVal;

	}

	/**
	 * Returns the number of tags which is equivalent to the number of values and
	 * row/columns in the covariance matrix.
	 *
	 * @return int
	 */
	final public int getDimension() {
		return mTags.size();
	}

	final public double getEntry(final int p) {
		return mValues.getEntry(p);
	}

	/**
	 * Returns the variance associated with the specific tag
	 *
	 * @param tag Object
	 * @return double The variance associated with tag
	 */
	public double getVariance(final Object tag) {
		final int p = indexOf(tag);
		return mCovariance.getEntry(p, p);
	}

	/**
	 * Returns the variance associated with the specific index
	 *
	 * @param idx int
	 * @return double The variance associated with tag
	 */
	public double getVariance(final int idx) {
		return mCovariance.getEntry(idx, idx);
	}

	/**
	 * Returns the uncertainty associated with the specific tag
	 *
	 * @param tag Object
	 * @return double The uncertainty associated with tag or 0.0 if tag unknown
	 */
	public double getUncertainty(final Object tag) {
		final int p = indexOf(tag);
		return p != -1 ? Math.sqrt(mCovariance.getEntry(p, p)) : 0.0;
	}

	/**
	 * Returns the uncertainty associated with the specific index
	 *
	 * @param p int
	 * @return double The uncertainty associated with index p
	 */
	public double getUncertainty(final int p) {
		return Math.sqrt(mCovariance.getEntry(p, p));
	}

	/**
	 * Returns the covariance associated with the specific tags
	 *
	 * @param tag1 Object
	 * @param tag2 Object
	 * @return double The covariance associated with tags tag1 and tag2 or 0 if one
	 *         or both tags unknown
	 */
	public double getCovariance(final Object tag1, final Object tag2) {
		final int p1 = indexOf(tag1), p2 = indexOf(tag2);
		// assert mCovariance.getEntry(p1, p2) == mCovariance.getEntry(p2, p1) :
		// mCovariance.getEntry(p1, p2) + "!="
		// + mCovariance.getEntry(p2, p1);
		return (p1 != -1) && (p2 != -1) ? mCovariance.getEntry(p1, p2) : 0.0;
	}

	/**
	 * Returns the covariance associated with the specific integer indices.
	 *
	 * @param p1    int
	 * @param p2int
	 * @return double The covariance associated with indices p1 and p2
	 */
	public double getCovariance(final int p1, final int p2) {
		return mCovariance.getEntry(p1, p2);
	}

	/**
	 * Returns the covariance associated with the specific integer indices.
	 *
	 * @param p1 int
	 * @param p2 int
	 * @return double The correlation coefficient associated with indices p1 and p2
	 */
	public double getCorrelationCoefficient(final int p1, final int p2) {
		return mCovariance.getEntry(p1, p2) / Math.sqrt(mCovariance.getEntry(p1, p1) * mCovariance.getEntry(p2, p2));
	}

	/**
	 * Returns a map of the variance and covariances relative to the specified tag.
	 *
	 * @param tag
	 * @return Map&lt;Object, Double&gt;
	 */
	public Map<Object, Double> getCovariances(final Object tag) {
		final Map<Object, Double> res = new HashMap<>();
		final int p = indexOf(tag);
		if (p != -1) {
			for (int i = 0; i < mTags.size(); ++i)
				if (mCovariance.getEntry(p, i) != 0.0)
					res.put(mTags.get(i), mCovariance.getEntry(p, i));
		}
		return res;
	}

	/**
	 * Returns the value and the associated uncertainty as an {@link UncertainValue}
	 * object.
	 *
	 * @param tag
	 * @return {@link UncertainValue}
	 */
	public UncertainValue getValue(final Object tag) {
		return new UncertainValue(getEntry(tag), tag, Math.sqrt(getVariance(tag)));
	}

	public String toCSV() {
		final StringBuffer sb = new StringBuffer(4096);
		sb.append("\"Name\",\"Value\"");
		for (int c = 0; c < mCovariance.getColumnDimension(); ++c) {
			sb.append(",\"");
			sb.append(mTags.get(c));
			sb.append("\"");
		}
		sb.append("\n");
		for (int r = 0; r < mValues.getDimension(); ++r) {
			sb.append("\"");
			sb.append(mTags.get(r));
			sb.append("\",");
			sb.append(mValues.getEntry(r));
			for (int c = 0; c < mCovariance.getColumnDimension(); ++c) {
				sb.append(",");
				sb.append(mCovariance.getEntry(r, c));
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
	 * Provides a mechanism to convert this {@link UncertainValues} object to HTML
	 * with a little more control over number formating.
	 *
	 * @param mode
	 * @param nf
	 * @return String
	 */
	public String toHTML(final Mode mode, final BasicNumberFormat nf) {
		switch (mode) {
		case TERSE: {
			final Table table = new Table();
			final List<Item> header = new ArrayList<>();
			final List<Item> vals = new ArrayList<>();
			for (final Object rowTag : mTags) {
				header.add(Table.th(HTML.toHTML(rowTag, Mode.TERSE)));
				vals.add(Table
						.tdc(nf.formatHTML(getEntry(rowTag)) + "&pm;" + nf.formatHTML(Math.sqrt(getVariance(rowTag)))));
			}
			table.addRow(header);
			table.addRow(vals);
			return table.toHTML(Mode.NORMAL);
		}
		case NORMAL:
		case VERBOSE:
		default: {
			final Table vals = new Table();
			final List<Item> all = new ArrayList<>();
			all.add(Table.td("Name"));
			all.add(Table.tdc("Quantity"));
			all.add(Table.td());
			for (final Object colTag : mTags)
				all.add(Table.tdc(HTML.toHTML(colTag, Mode.TERSE)));
			vals.addRow(all);
			final BasicNumberFormat bnf2 = new BasicNumberFormat("0.00");
			for (int r = 0; r < mTags.size(); ++r) {
				final Object rowTag = mTags.get(r);
				final List<Item> row = new ArrayList<>();
				row.add(Table.td(HTML.toHTML(rowTag, Mode.TERSE)));
				row.add(Table.tdc(nf.formatHTML(getEntry(rowTag))));
				if (r == (mTags.size() - 1) / 2)
					row.add(Table.tdc("&nbsp;&nbsp;&plusmn;&nbsp;&nbsp;"));
				else
					row.add(Table.td());
				for (final Object colTag : mTags) {
					final double c = getCovariance(rowTag, colTag);
					if ((rowTag == colTag) || (mode != Mode.VERBOSE))
						row.add(Table.tdc((c > 0.0 ? "" : Transforms.NON_BREAKING_DASH) + "("
								+ nf.formatHTML(Math.sqrt(Math.abs(c))) + ")<sup>2</sup>"));
					else {
						final double cv = Math.sqrt(getVariance(rowTag) * getVariance(colTag));
						final double cc = cv > 0.0 ? c / cv : 0.0;
						// Compute the corrolatation coefficient
						assert (cc >= -MAX_CORR) && (cc <= MAX_CORR) : cc;
						if (cc == 0.0)
							row.add(Table.td());
						else
							row.add(Table.tdc((cc < 0.0 ? Transforms.NON_BREAKING_DASH : "")
									+ bnf2.formatHTML(Math.sqrt(Math.abs(cc)))
									+ "&middot;&sigma;<sub>R</sub>&sigma;<sub>C</sub>"));
					}
				}
				vals.addRow(row);
			}
			return vals.toHTML(Mode.NORMAL);
		}
		}
	}

	public BufferedImage asCovarianceBitmap(final int dim, final IValueToColor sigma, final IValueToColor corr) {
		final Array2DRowRealMatrix sc = new Array2DRowRealMatrix(getDimension(), getDimension());
		for (int r = 0; r < getDimension(); ++r) {
			sc.setEntry(r, r, Math.sqrt(mCovariance.getEntry(r, r)) / mValues.getEntry(r));
			for (int c = r + 1; c < getDimension(); ++c) {
				final double rr = mCovariance.getEntry(r, c)
						/ (Math.sqrt(mCovariance.getEntry(r, r) * mCovariance.getEntry(c, c)));
				if (!Double.isNaN(rr)) {
					assert rr >= -MAX_CORR : rr;
					assert rr <= MAX_CORR : rr;
					sc.setEntry(r, c, rr);
					sc.setEntry(c, r, rr);
				}
			}
		}
		final BufferedImage bi = new BufferedImage(dim * getDimension(), dim * getDimension(),
				BufferedImage.TYPE_3BYTE_BGR);
		final Graphics2D g2 = bi.createGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, bi.getWidth(), bi.getHeight());
		for (int r = 0; r < getDimension(); ++r) {
			g2.setColor(sigma.map(sc.getEntry(r, r)));
			g2.fillRect(r * dim, r * dim, dim, dim);
			for (int c = r + 1; c < getDimension(); ++c) {
				g2.setColor(corr.map(sc.getEntry(r, c)));
				g2.fillRect(c * dim, r * dim, dim, dim);
				g2.fillRect(r * dim, c * dim, dim, dim);
			}
		}
		return bi;
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
	public static BufferedImage compareAsBitmap(UncertainValues uvs1, UncertainValues uvs2, final IValueToColor corr,
			int pixDim) {
		if (uvs1.getDimension() != uvs2.getDimension())
			throw new DimensionMismatchException(uvs2.getDimension(), uvs1.getDimension());
		int dim = uvs1.getDimension();
		for (int i = 0; i < dim; ++i)
			assert uvs1.getTag(i) == uvs2.getTag(i) : uvs1.getTag(i) + "!=" + uvs2.getTag(i);
		final Array2DRowRealMatrix sc = new Array2DRowRealMatrix(dim, dim);
		for (int r = 0; r < dim; ++r) {
			for (int c = 0; c < dim; ++c) {
				final double rr = 0.5 * (uvs1.getCovariance(r, c) - uvs2.getCovariance(r, c))
						/ (uvs1.getCovariance(r, c) + uvs2.getCovariance(r, c));
				if (!Double.isNaN(rr)) {
					sc.setEntry(r, c, rr);
					sc.setEntry(c, r, rr);
				}
			}
		}
		final BufferedImage bi = new BufferedImage(pixDim * dim, pixDim * dim, BufferedImage.TYPE_3BYTE_BGR);
		final Graphics2D g2 = bi.createGraphics();
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, bi.getWidth(), bi.getHeight());
		for (int r = 0; r < dim; ++r) {
			g2.setColor(corr.map(sc.getEntry(r, r)));
			g2.fillRect(r * dim, r * dim, dim, dim);
			for (int c = r + 1; c < dim; ++c) {
				g2.setColor(corr.map(sc.getEntry(r, c)));
				g2.fillRect(c * pixDim, r * pixDim, pixDim, pixDim);
				g2.fillRect(r * pixDim, c * pixDim, pixDim, pixDim);
			}
		}
		return bi;
	}

	public boolean equals(final UncertainValues uv2, final double tol) {
		for (int r = 0; r < this.getDimension(); ++r) {
			if (Math.abs(getEntry(r) - uv2.getEntry(mTags.get(r))) > tol)
				return false;
			for (int c = r; c < this.getDimension(); ++c)
				if (Math.abs(getCovariance(r, c) - uv2.getCovariance(mTags.get(r), mTags.get(c))) > tol)
					return false;
		}
		return true;
	}

	/**
	 * <p>
	 * Copies those entries tagged in <i>from</i> that are also tagged in <i>to</i>.
	 * Useful for replicating the variances and covariances in <i>from</i> in
	 * <i>to</i>.
	 * </p>
	 * <p>
	 * Will overwrite values in <i>to</i>
	 * </p>
	 *
	 * @param from
	 * @param to
	 */
	public static void copy(final UncertainValues from, final UncertainValues to) {
		final List<? extends Object> fromTags = from.getTags();
		final Map<Object, Integer> toMap = new HashMap<>();
		final Map<Object, Integer> fromMap = new HashMap<>();
		for (int i = 0; i < fromTags.size(); ++i) {
			final Object tag = fromTags.get(i);
			final int p = to.indexOf(tag);
			if (p >= 0) {
				toMap.put(tag, Integer.valueOf(p));
				fromMap.put(tag, i);
			}
		}
		for (final Map.Entry<Object, Integer> me1 : toMap.entrySet()) {
			final int p1 = me1.getValue().intValue();
			to.mValues.setEntry(p1, from.getEntry(me1.getKey()));
			for (final Map.Entry<Object, Integer> me2 : toMap.entrySet()) {
				final int p2 = me2.getValue().intValue();
				to.mCovariance.setEntry(p1, p2, from.getCovariance(me1.getKey(), me2.getKey()));
			}
		}
	}

	/**
	 * Initializes the value and variance associated with the specified tag.
	 *
	 * @param tag
	 * @param val
	 * @param var
	 */
	public void set(final Object tag, final double val, final double var) {
		final int p = indexOf(tag);
		mValues.setEntry(p, val);
		mCovariance.setEntry(p, p, var);
	}

	/**
	 * Sets the covariance value associated with the specified values (both i,j and
	 * j,i)
	 *
	 * @param tag1
	 * @param tag2
	 * @param cov
	 */
	public void setCovariance(final Object tag1, final Object tag2, final double cov) {
		final int p1 = indexOf(tag1), p2 = indexOf(tag2);
		mCovariance.setEntry(p1, p2, cov);
		mCovariance.setEntry(p2, p1, cov);
	}

	/**
	 * Initializes the specified value and variance for the specified tag.
	 *
	 * @param tag
	 * @param uv  (uses the doubleValue() and variance())
	 */
	public void set(final Object tag, final UncertainValue uv) {
		final int p = indexOf(tag);
		mValues.setEntry(p, uv.doubleValue());
		mCovariance.setEntry(p, p, uv.variance());
	}

	/**
	 * Checks the values and covariances to determine whether any are equivalent to
	 * NaN.
	 *
	 * @return boolean true if one value is NaN, false otherwise.
	 */
	public boolean isNaN() {
		if (mValues.isNaN())
			return true;
		for (int r = 0; r < mCovariance.getRowDimension(); ++r)
			for (int c = 0; c < mCovariance.getColumnDimension(); ++c)
				if (Double.isNaN(mCovariance.getEntry(r, c)))
					return true;
		return false;
	}

	public RealMatrix delta(final UncertainValues zeroed) {
		return mCovariance.subtract(zeroed.mCovariance);
	}

	private int rmHashCode(final RealMatrix rm) {
		return rm.getData().hashCode();
	}

	private boolean rmEquals(final RealMatrix rm1, final RealMatrix rm2) {
		for (int i = 0; i < rm1.getRowDimension(); ++i)
			if (!Arrays.equals(rm1.getRow(i), rm2.getRow(i)))
				return false;
		return true;
	}

	@Override
	public int hashCode() {
		return 37 * Objects.hashCode(mValues, mTags) + rmHashCode(mCovariance);
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
		return rmEquals(mCovariance, other.mCovariance) && Objects.equal(mTags, other.mTags)
				&& Objects.equal(mValues, other.mValues);
	}

	/**
	 * Build a new UncertainValues from this one in which the tags have been sorted
	 * into an order determined by the specified {@link Comparator}.
	 *
	 * @param compare
	 * @return {@link UncertainValues} A new instance with the same data reordered.
	 */
	public UncertainValues sort(final Comparator<Object> compare) {
		final List<Object> tags = new ArrayList<>(mTags);
		tags.sort(compare);
		final UncertainValues res = new UncertainValues(tags);
		for (int r = 0; r < mTags.size(); ++r) {
			res.mValues.setEntry(res.indexOf(mTags.get(r)), mValues.getEntry(r));
			for (int c = 0; c < mTags.size(); ++c)
				res.setCovariance(mTags.get(r), mTags.get(c), mCovariance.getEntry(r, c));
		}
		return res;
	}

	/**
	 * Build a new UncertainValues from this one in which the tags have been sorted
	 * into alphabetical order by tag.toString().
	 *
	 * @return {@link UncertainValues} A new instance with the same data reordered.
	 */
	public UncertainValues sort() {
		return sort(new Comparator<Object>() {

			@Override
			public int compare(final Object o1, final Object o2) {
				return o1.toString().compareTo(o2.toString());
			}
		});
	}

}
