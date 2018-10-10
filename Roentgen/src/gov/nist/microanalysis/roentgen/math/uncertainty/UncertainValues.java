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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;
import com.duckandcover.html.Transforms;
import com.google.common.base.Objects;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
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
		implements IToHTML {

	private static final double MAX_CORR = 1.00000001;
	private final RealVector mValues;
	private final List<Object> mLabels;
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
		final double EPS = 1.0e-10, SREPS = Math.sqrt(EPS);
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
				if (Math.abs(entryRC - covar.getEntry(c, r)) > SREPS * Math.abs(entryRC))
					throw new OutOfRangeException(entryRC, covar.getEntry(c, r), covar.getEntry(c, r));
				final double max = Math.sqrt(covar.getEntry(c, c) * covar.getEntry(r, r));
				final double rr = entryRC / max;
				if ((rr > MAX_CORR) || (rr < -MAX_CORR)) {
					if (Math.abs(entryRC) > EPS)
						throw new OutOfRangeException(entryRC, -max, max);
					else {
						covar.setEntry(r, c, Math.signum(rr) * max);
						covar.setEntry(c, r, Math.signum(rr) * max);
					}
				}
			}
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * covariance matrix.
	 *
	 * @param labels  List&lt;Object&gt; A list of objects implementing hashCode() and
	 *              equals().
	 * @param vals  {@link RealVector}
	 * @param covar {@link RealMatrix}
	 */
	public UncertainValues( //
			final List<? extends Object> labels, //
			final RealVector vals, //
			final RealMatrix covar //
	) {
		if (vals.getDimension() != labels.size())
			throw new DimensionMismatchException(covar.getRowDimension(), labels.size());
		if (covar.getRowDimension() != labels.size())
			throw new DimensionMismatchException(covar.getRowDimension(), labels.size());
		if (covar.getColumnDimension() != labels.size())
			throw new DimensionMismatchException(covar.getColumnDimension(), labels.size());
		mLabels = Collections.unmodifiableList(labels);
		final HashMap<Object, Integer> index = new HashMap<>();
		for (int i = 0; i < mLabels.size(); ++i)
			index.put(mLabels.get(i), Integer.valueOf(i));
		mIndex = Collections.unmodifiableMap(index);
		mValues = vals;
		validateCovariance(mLabels.size(), covar, Math.max(mValues.getMaxValue(), -mValues.getMinValue()));
		mCovariance = covar;
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels with zero
	 * values and NaN covariances.
	 *
	 * @param labels List&lt;Object&gt; A list of objects implementing hashCode() and
	 *             equals().
	 */
	public UncertainValues(//
			final List<? extends Object> labels //
	) {
		this(labels, new ArrayRealVector(labels.size()), new Array2DRowRealMatrix(labels.size(), labels.size()));
		mValues.set(Double.NaN);
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances array. "Heteroskedastic-case"
	 *
	 * @param labels      List&lt;Object&gt; A list of objects implementing hashCode()
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

	public UncertainValues(//
			final Map<? extends Object, UncertainValue> vals //
	) {
		this(Arrays.asList(vals.keySet().toArray()), extractValues(vals), extractVariances(vals));
	}

	/**
	 * Constructs a UncertainValues object based on the specified labels, values and
	 * variances. "Homoskedastic-case"
	 *
	 * @param labels     List&lt;Object&gt; A list of objects implementing hashCode()
	 *                 and equals().
	 * @param vals     {@link RealVector}
	 * @param variance double (common to all vals)
	 */
	public UncertainValues(final List<? extends Object> labels, final RealVector vals, final double variance) {
		this(labels, vals, MatrixUtils.createRealIdentityMatrix(labels.size()).scalarMultiply(variance));
	}

	/**
	 * Builds an {@link UncertainValues} object representing the specified labeled quantities
	 * as extracted from the list of {@link UncertainValues} objects.
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
				throw new ArgumentException("The label " + label + " is not defined in one of the input UncertainValues.");
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
	 * Combines a disjoint set of {@link UncertainValues} into a single one.
	 * (Disjoint meaning not sharing a common label.)
	 * 
	 * @param uvs
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public static UncertainValues combine(final UncertainValues... uvs) throws ArgumentException {
		// Test that each requested label is defined once and only once.
		final List<Object> labels = new ArrayList<>();
		for (final UncertainValues uv : uvs) {
			for (final Object label : uv.getLabels())
				if (labels.contains(label))
					throw new ArgumentException("The label " + label + " is duplicated in UncertainValues.combine(...)");
				else
					labels.add(label);
		}
		final UncertainValues res = new UncertainValues(labels);
		int r0 = 0;
		for (final UncertainValues uv : uvs) {
			int dr = 0;
			for (; dr < uv.getDimension(); ++dr) {
				res.mValues.setEntry(r0 + dr, uv.getEntry(dr));
				for (int dc = 0; dc < uv.getDimension(); ++dc)
					res.mCovariance.setEntry(r0 + dr, r0 + dc, uv.mCovariance.getEntry(dr, dc));
			}
			r0 += dr;
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
		for (int i = 0; i < res.getDimension(); ++i) {
			res.mCovariance.setEntry(i, idx, input.getCovariance(i, idx));
			res.mCovariance.setEntry(idx, i, input.getCovariance(idx, i));
		}
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
	 */
	public static UncertainValues propagate( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues input //
	) {
		assert nmjf.getInputLabels().equals(input.getLabels());
		final Pair<RealVector, RealMatrix> eval = nmjf.evaluate(input.getValues());
		final RealMatrix jac = eval.getSecond();
		return new UncertainValues(nmjf.getOutputLabels(), eval.getFirst(),
				jac.multiply(input.getCovariances().multiply(jac.transpose())));
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
	 * @param vals
	 * @param nEvals
	 * @return UncertainValues
	 */
	public static UncertainValues propagateMC(//
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues vals, //
			final int nEvals //
	) {
		MCPropagator mcp = new MCPropagator(nmjf, vals);
		return mcp.compute(nEvals);	}

	/**
	 * Extract an array of values from this UncertainValue object for the specified
	 * list of labels in the order specified by the label list.
	 *
	 * @param labels List&lt;? extends Object&gt;
	 * @return RealVector
	 */
	public RealVector extractValues(final List<? extends Object> labels) {
		final RealVector res = new ArrayRealVector(labels.size());
		int i = 0;
		for (final Object label : labels) {
			res.setEntry(i, getEntry(label));
			++i;
		}
		return res;
	}

	/**
	 * Extract an array of values from this UncertainValue object for the specified
	 * list of labels in the order specified by the label list.
	 *
	 * @param labels List&lt;? extends Object&gt;
	 * @return RealVector
	 */
	public RealMatrix extractCovariances(final List<? extends Object> labels) {
		final RealMatrix res = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		for (int r = 0; r < labels.size(); ++r) {
			int ridx = indexOf(labels.get(r));
			res.setEntry(r, r, getCovariance(r, r));
			for (int c = r+1; c < labels.size(); ++c) {
				int cidx = indexOf(labels.get(c));
				final double cv = getCovariance(ridx, cidx);
				res.setEntry(r, c, cv);
				res.setEntry(c, r, cv);
			}
		}
		return res;
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

	/**
	 * Returns a {@link RealVector} containing the values associated with this
	 * object in the order specified by the List labels.
	 *
	 * @param labels
	 * @return {@link RealVector}
	 */
	final public RealVector getValues(final List<? extends Object> labels) {
		final RealVector res = new ArrayRealVector(labels.size());
		for (int i = 0; i < res.getDimension(); ++i)
			res.setEntry(i, mValues.getEntry(indexOf(labels.get(i))));
		return res;
	}

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
	public static UncertainValues extract(final List<? extends Object> labels, final UncertainValues uvs) {
		final RealVector vals = new ArrayRealVector(labels.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(labels.size(), labels.size());
		for (int ri = 0; ri < labels.size(); ++ri) {
			final int r = uvs.indexOf(labels.get(ri));
			vals.setEntry(ri, uvs.getEntry(r));
			cov.setEntry(ri, ri, uvs.getCovariance(r, r));
			for (int ci = r + 1; ci < labels.size(); ++ci) {
				final int c = uvs.indexOf(labels.get(ci));
				final double cc = uvs.getCovariance(r, c);
				cov.setEntry(ri, ci, cc);
				cov.setEntry(ci, ri, cc);
			}
		}
		return new UncertainValues(labels, vals, cov);

	}

	/**
	 * Returns an ordered list of the Object labels associated with the values and
	 * covariances in this object.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<Object> getLabels() {
		return mLabels;
	}

	/**
	 * Is there a value and covariances associated with the specified label?
	 *
	 * @param label
	 * @return boolean
	 */
	final public boolean hasEntry(final Object label) {
		return mLabels.contains(label);
	}

	/**
	 * Returns the label at the p-th entry in the values vector and in the p-th row
	 * and column in the covariance matrix.
	 *
	 * @param p
	 * @return Object
	 */
	final public Object getLabel(final int p) {
		return mLabels.get(p);
	}

	/**
	 * Returns the index associated with the specified label or -1 if not found.
	 *
	 * @param p
	 * @return Object
	 */
	final public int indexOf(final Object label) {
		return mIndex.getOrDefault(label, Integer.valueOf(-1));
	}

	/**
	 * Returns the value associated with the entry associated with label.
	 *
	 * @param label Object
	 * @return double
	 */
	final public double getEntry(final Object label) {
		final int p = indexOf(label);
		assert p >= 0 : "Label " + label.toString() + " missing.";
		return mValues.getEntry(p);
	}

	/**
	 * Returns the value associated with the entry associated with label if the label is
	 * associated with a value or returns defVal otherwise. Useful for a quantity
	 * that may or may not have an associated uncertainty.
	 *
	 * @param label
	 * @param defVal
	 * @return double
	 */
	final public double getEntryWithDefault(final Object label, final double defVal) {
		final int p = indexOf(label);
		return p >= 0 ? mValues.getEntry(p) : defVal;
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified label.
	 * The value is the same as getEntry(...) and the uncertainty is the on-diagonal
	 * covariance matrix entry associated with that entry.
	 *
	 * @param label
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValue(final Object label) {
		final int p = indexOf(label);
		return new UncertainValue(getEntry(p), label, Math.sqrt(mCovariance.getEntry(p, p)));
	}

	/**
	 * Returns an UncertainValue for the quantity assoicated with the specified label
	 * (if defined). The value is the same as getEntry(...) and the uncertainty is
	 * the on-diagonal covariance matrix entry associated with that entry. If the
	 * label isn't defined, then the value defVal will be returned.
	 *
	 * @param label
	 * @param defVal
	 * @return UncertanValue
	 */
	final public UncertainValue getUncertainValueWithDefault(final Object label, final UncertainValue defVal) {
		final int p = indexOf(label);
		return p >= 0 ? new UncertainValue(getEntry(p), label, mCovariance.getEntry(p, p)) : defVal;

	}

	/**
	 * Returns the number of labels which is equivalent to the number of values and
	 * row/columns in the covariance matrix.
	 *
	 * @return int
	 */
	final public int getDimension() {
		return mLabels.size();
	}

	final public double getEntry(final int p) {
		return mValues.getEntry(p);
	}

	/**
	 * Returns the variance associated with the specific label
	 *
	 * @param label Object
	 * @return double The variance associated with label
	 */
	public double getVariance(final Object label) {
		final int p = indexOf(label);
		return mCovariance.getEntry(p, p);
	}

	/**
	 * Returns the variance associated with the specific index
	 *
	 * @param idx int
	 * @return double The variance associated with label
	 */
	public double getVariance(final int idx) {
		return mCovariance.getEntry(idx, idx);
	}

	/**
	 * Returns the uncertainty associated with the specific label
	 *
	 * @param label Object
	 * @return double The uncertainty associated with label or 0.0 if label unknown
	 */
	public double getUncertainty(final Object label) {
		final int p = indexOf(label);
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
	 * Returns the covariance associated with the specific labels
	 *
	 * @param label1 Object
	 * @param label2 Object
	 * @return double The covariance associated with labels label1 and label2 or 0 if one
	 *         or both labels unknown
	 */
	public double getCovariance(final Object label1, final Object label2) {
		final int p1 = indexOf(label1), p2 = indexOf(label2);
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
	 * Returns a map of the variance and covariances relative to the specified label.
	 *
	 * @param label
	 * @return Map&lt;Object, Double&gt;
	 */
	public Map<Object, Double> getCovariances(final Object label) {
		final Map<Object, Double> res = new HashMap<>();
		final int p = indexOf(label);
		if (p != -1) {
			for (int i = 0; i < mLabels.size(); ++i)
				if (mCovariance.getEntry(p, i) != 0.0)
					res.put(mLabels.get(i), mCovariance.getEntry(p, i));
		}
		return res;
	}

	public Map<Object, Double> getValueMap() {
		final Map<Object, Double> res = new HashMap<>();
		for (final Object label : getLabels())
			res.put(label, getEntry(label));
		return res;
	}

	/**
	 * Returns the value and the associated uncertainty as an {@link UncertainValue}
	 * object.
	 *
	 * @param label
	 * @return {@link UncertainValue}
	 */
	public UncertainValue getValue(final Object label) {
		return new UncertainValue(getEntry(label), label, Math.sqrt(getVariance(label)));
	}

	public String toCSV() {
		final StringBuffer sb = new StringBuffer(4096);
		sb.append("\"Name\",\"Value\"");
		for (int c = 0; c < mCovariance.getColumnDimension(); ++c) {
			sb.append(",\"");
			sb.append(mLabels.get(c));
			sb.append("\"");
		}
		sb.append("\n");
		for (int r = 0; r < mValues.getDimension(); ++r) {
			sb.append("\"");
			sb.append(mLabels.get(r));
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
			for (final Object rowLabel : mLabels) {
				header.add(Table.th(HTML.toHTML(rowLabel, Mode.TERSE)));
				vals.add(Table
						.tdc(nf.formatHTML(getEntry(rowLabel)) + "&pm;" + nf.formatHTML(Math.sqrt(getVariance(rowLabel)))));
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
			for (final Object colLabel : mLabels)
				all.add(Table.tdc(HTML.toHTML(colLabel, Mode.TERSE)));
			vals.addRow(all);
			final BasicNumberFormat bnf2 = new BasicNumberFormat("0.00");
			for (int r = 0; r < mLabels.size(); ++r) {
				final Object rowLabel = mLabels.get(r);
				final List<Item> row = new ArrayList<>();
				row.add(Table.td(HTML.toHTML(rowLabel, Mode.TERSE)));
				row.add(Table.tdc(nf.formatHTML(getEntry(rowLabel))));
				if (r == (mLabels.size() - 1) / 2)
					row.add(Table.tdc("&nbsp;&nbsp;&plusmn;&nbsp;&nbsp;"));
				else
					row.add(Table.td());
				for (final Object colLabel : mLabels) {
					final double c = getCovariance(rowLabel, colLabel);
					if ((rowLabel == colLabel) || (mode != Mode.VERBOSE))
						row.add(Table.tdc((c > 0.0 ? "" : Transforms.NON_BREAKING_DASH) + "("
								+ nf.formatHTML(Math.sqrt(Math.abs(c))) + ")<sup>2</sup>"));
					else {
						final double cv = Math.sqrt(getVariance(rowLabel) * getVariance(colLabel));
						final double cc = cv > 0.0 ? c / cv : 0.0;
						// Compute the corrolatation coefficient
						// assert (cc >= -MAX_CORR) && (cc <= MAX_CORR) : cc;
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
		final RealMatrix sc = new NullableRealMatrix(getDimension(), getDimension());
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
	public static BufferedImage compareAsBitmap(final UncertainValues uvs1, final UncertainValues uvs2,
			final IValueToColor corr, final int pixDim) {
		if (uvs1.getDimension() != uvs2.getDimension())
			throw new DimensionMismatchException(uvs2.getDimension(), uvs1.getDimension());
		final int dim = uvs1.getDimension();
		for (int i = 0; i < dim; ++i)
			assert uvs1.getLabel(i) == uvs2.getLabel(i) : uvs1.getLabel(i) + "!=" + uvs2.getLabel(i);
		final Array2DRowRealMatrix sc = new Array2DRowRealMatrix(dim, dim);
		for (int r = 0; r < dim; ++r) {
			final double rr = (uvs1.getCovariance(r, r) - uvs2.getCovariance(r, r))
					/ Math.max(1.0 - 100, (uvs1.getCovariance(r, r) + uvs2.getCovariance(r, r)));
			sc.setEntry(r, r, Math.max(-1, Math.min(1.0, rr)));
			for (int c = r + 1; c < dim; ++c) {
				final double rc = (uvs1.getCorrelationCoefficient(r, c) - uvs2.getCorrelationCoefficient(r, c)) / 2.0;
				if (!Double.isNaN(rc)) {
					sc.setEntry(r, c, rc);
					sc.setEntry(c, r, rc);
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
				if (!Double.isNaN(sc.getEntry(r, c)))
					g2.setColor(corr.map(sc.getEntry(r, c)));
				else
					g2.setColor(Color.yellow);
				g2.fillRect(c * pixDim, r * pixDim, pixDim, pixDim);
				g2.fillRect(r * pixDim, c * pixDim, pixDim, pixDim);
			}
		}
		return bi;
	}

	public boolean equals(final UncertainValues uv2, final double tol) {
		for (int r = 0; r < this.getDimension(); ++r) {
			if (Math.abs(getEntry(r) - uv2.getEntry(mLabels.get(r))) > tol)
				return false;
			for (int c = r; c < this.getDimension(); ++c)
				if (Math.abs(getCovariance(r, c) - uv2.getCovariance(mLabels.get(r), mLabels.get(c))) > tol)
					return false;
		}
		return true;
	}

	/**
	 * <p>
	 * Copies those entries labeled in <i>from</i> that are also labeled in <i>to</i>.
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
	 * Initializes the value and variance associated with the specified label.
	 *
	 * @param label
	 * @param val
	 * @param var
	 */
	public void set(final Object label, final double val, final double var) {
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
	public void setCovariance(final Object label1, final Object label2, final double cov) {
		final int p1 = indexOf(label1), p2 = indexOf(label2);
		mCovariance.setEntry(p1, p2, cov);
		mCovariance.setEntry(p2, p1, cov);
	}

	/**
	 * Initializes the specified value and variance for the specified label.
	 *
	 * @param label
	 * @param uv  (uses the doubleValue() and variance())
	 */
	public void set(final Object label, final UncertainValue uv) {
		final int p = indexOf(label);
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

	@Override
	public int hashCode() {
		return Objects.hashCode(mValues, mLabels, mCovariance);
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
		return Objects.equal(mCovariance, other.mCovariance) && //
				Objects.equal(mLabels, other.mLabels) && //
				Objects.equal(mValues, other.mValues);
	}

	/**
	 * Build a new UncertainValues from this one in which the labels have been sorted
	 * into an order determined by the specified {@link Comparator}.
	 *
	 * @param compare
	 * @return {@link UncertainValues} A new instance with the same data reordered.
	 */
	public UncertainValues sort(final Comparator<Object> compare) {
		final List<Object> labels = new ArrayList<>(mLabels);
		labels.sort(compare);
		final UncertainValues res = new UncertainValues(labels);
		for (int r = 0; r < mLabels.size(); ++r) {
			res.mValues.setEntry(res.indexOf(mLabels.get(r)), mValues.getEntry(r));
			for (int c = 0; c < mLabels.size(); ++c)
				res.setCovariance(mLabels.get(r), mLabels.get(c), mCovariance.getEntry(r, c));
		}
		return res;
	}

	/**
	 * Build a new UncertainValues from this one in which the labels have been sorted
	 * into alphabetical order by label.toString().
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
