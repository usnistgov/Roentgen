package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.utility.FastIndex;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * The Jacobian is a matrix consisting of n partial derivatives associated with
 * m functions. This implementation of the Jacobian is designed to work with the
 * CovarianceMatrix class as these two classes are necessary to implement a
 * linear algebra-based implementation
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
abstract public class LabeledMultivariateJacobianFunction //
		implements MultivariateJacobianFunction, IToHTML //
{

	/***
	 * Unique object labels identifying each of the random variable arguments to the
	 * functions
	 */
	private final List<? extends Object> mInputLabels;
	/***
	 * Unique object labels identifying each of the functions (in order)
	 */
	private final List<? extends Object> mOutputLabels;

	public static PrintStream sDump = null;

	protected void dumpArguments(RealVector point, LabeledMultivariateJacobianFunction parent) {
		if (sDump != null) {
			StringBuffer sb = new StringBuffer();
			NumberFormat nf = new HalfUpFormat("0.00E0");
			sb.append(toString());
			sb.append("[");
			for (int i = 0; i < getInputDimension(); ++i) {
				if (i != 0)
					sb.append(",");
				final Object lbl = getInputLabel(i);
				sb.append(lbl);
				sb.append("=");
				sb.append(nf.format(point.getEntry(i)));
			}
			sb.append("] in ");
			sb.append(parent);
			sb.append("\n");
			sDump.append(sb);
		}
	}

	/**
	 * Check labels are only used once.
	 *
	 * @param labels
	 */
	private void validateLabels(final Collection<? extends Object> labels) {
		final Set<? extends Object> set = new HashSet<>(labels);
		for (Object lbl : labels)
			if (!set.remove(lbl))
				throw new RuntimeException("The label " + lbl + " is duplicated.");
	}

	public LabeledMultivariateJacobianFunction( //
			final List<? extends Object> inputLabels, //
			final List<? extends Object> outputs) {
		validateLabels(inputLabels);
		// assert inputLabels != outputLabels;
		mInputLabels = new FastIndex<>(inputLabels);
		validateLabels(outputs);
		mOutputLabels = new FastIndex<>(outputs);
	}

	/**
	 * Returns an array consisting of the labels associated with the n input random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<? extends Object> getInputLabels() {
		return mInputLabels;
	}

	/**
	 * Returns an array consisting of the labels associated with the m output random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<? extends Object> getOutputLabels() {
		return mOutputLabels;
	}

	final public Object getOutputLabel(final int idx) {
		return mOutputLabels.get(idx);
	}

	final public Object getInputLabel(final int idx) {
		return mInputLabels.get(idx);
	}

	/**
	 * Number of input random variables expected.
	 *
	 * @return int
	 */
	final public int getInputDimension() {
		return mInputLabels.size();
	}

	/**
	 * Number of output random variable produced.
	 *
	 * @return int
	 */
	final public int getOutputDimension() {
		return mOutputLabels.size();
	}

	/**
	 * Returns the index of the input variable identified by the specified Object
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	final public int inputIndex(final Object label) {
		return mInputLabels.indexOf(label);
	}

	/**
	 * Returns the index of the output variable identified by the specified Object
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	final public int outputIndex(final Object label) {
		return mOutputLabels.indexOf(label);
	}

	/**
	 * Only writes <code>value</code> to the Jacobian matrix <code>jacob</code> if
	 * <code>label</code> is an input variable (not a constant.) This is useful for
	 * implementing compute(...) to handle situations in which an input may either
	 * be considered a constant or a random variable.
	 *
	 *
	 * @param row   Row index
	 * @param label Column label
	 * @param value Quantity to place at row, inputIndex(label) in
	 *              <code>jacob</code>
	 * @param jacob The result Jacobian (modified)
	 */
	public void writeJacobian(final int row, final Object label, final double value, final RealMatrix jacob) {
		final int p = inputIndex(label);
		if (p != -1)
			jacob.setEntry(row, p, value);
	}

	/**
	 * Extracts from <code>point</code> representing an argument to
	 * <code>nmvj.compute(point)</code> the input RealArray that is suitable as the
	 * argument to <code>this.comptue(...)</code>.
	 *
	 * @param nmvj  The outer {@link LabeledMultivariateJacobianFunction}
	 * @param point A {@link RealVector} of length nmvj.getInputDimension()
	 * @return {@link RealVector} of length this.getInputDimension() containing the
	 *         appropriate input values from point.
	 */
	public RealVector extractArgument(final LabeledMultivariateJacobianFunction nmvj, final RealVector point) {
		assert point.getDimension() == nmvj.getInputDimension();
		final int dim = getInputDimension();
		final RealVector res = new ArrayRealVector(dim);
		final List<? extends Object> labels = getInputLabels();
		for (int i = 0; i < dim; ++i) {
			final int idx = nmvj.inputIndex(labels.get(i));
			assert idx != -1 : "Can't find " + labels.get(i) + " in the arguments to " + nmvj.toString();
			res.setEntry(i, point.getEntry(idx));
		}
		return res;
	}

	/**
	 * <p>
	 * A safer version of <code>value(x)</code>.
	 * </p>
	 * <p>
	 * Checks the length of the input argument x, calls <code>value(x)</code> and
	 * then checks the output {@link RealVector} and {@link RealMatrix} to ensure
	 * that they are all the correct dimensions.
	 * </p>
	 *
	 * @param x
	 * @return Pair&lt;{@link RealVector}, {@link RealMatrix}&gt; As from a call to
	 *         <code>value(x)</code>.
	 */
	final public Pair<RealVector, RealMatrix> evaluate(final RealVector x) {
		if (x.getDimension() != getInputDimension())
			throw new DimensionMismatchException(x.getDimension(), getInputDimension());
		final Pair<RealVector, RealMatrix> res = value(x);
		if (res.getFirst().getDimension() != getOutputDimension())
			throw new DimensionMismatchException(res.getFirst().getDimension(), getOutputDimension());
		if (res.getSecond().getRowDimension() != getOutputDimension())
			throw new DimensionMismatchException(res.getSecond().getRowDimension(), getOutputDimension());
		if (res.getSecond().getColumnDimension() != getInputDimension())
			throw new DimensionMismatchException(res.getSecond().getColumnDimension(), getInputDimension());
		return res;
	}

	/**
	 * Computes the result {@link RealVector} using the most efficient available
	 * mechanism depending upon whether the instantiated class implements
	 * {@link ILabeledMultivariateFunction} or not. If
	 * {@link ILabeledMultivariateFunction} is available then this is called.
	 * Otherwise value(...) is called and only the {@link RealVector} value portion
	 * is returned.
	 *
	 *
	 * @param inp Evaluation point
	 * @return {@link RealVector} The result values
	 */
	final public RealVector compute(final RealVector inp) {
		if (this instanceof ILabeledMultivariateFunction)
			return ((ILabeledMultivariateFunction) this).optimized(inp);
		else
			return value(inp).getFirst();
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects. Equivalent to
	 * <code>getOutputValues(uvs, 1.0e-6)</code>
	 *
	 * @param uvs Input {@link UncertainValues}
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues(final UncertainValuesBase uvs)
			throws ArgumentException {
		return getOutputValues(uvs, 1.0e-6);
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects.
	 *
	 * @param uvs Input {@link UncertainValues}
	 * @param tol Tolerance relative to value
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues( //
			final UncertainValuesBase uvs, //
			final double tol//
	) throws ArgumentException {
		UncertainValuesCalculator uvc = new UncertainValuesCalculator(this, uvs, true);
		return uvc.getOutputValues(tol);
	}

	/**
	 * <p>
	 * Get the resulting UncertainValue for each output quantity with sources
	 * grouped and named according to the map of names to a collection of labels.
	 * </p>
	 * </p>
	 * Returns a map from output value to an UncertainValue with discretely broken
	 * out uncertainty components according to the labels map.
	 * </p>
	 *
	 * @param uvs
	 * @param labels Group according to this mapping
	 * @param tol    Minimum size uncertainty to include
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues( //
			final UncertainValues uvs, //
			final Map<String, Collection<? extends Object>> labels, //
			final double tol //
	) throws ArgumentException {
		UncertainValuesCalculator uvc = new UncertainValuesCalculator(this, uvs);
		return uvc.getOutputValues(labels, tol);
	}

	/**
	 * <p>
	 * First checks to see if <code>label</code> is an input variable in which case
	 * it gets the associated value from <code>point</code>. Otherwise, it returns
	 * the constant value associate with <code>label</code>. This is useful for
	 * implementing compute(...) in situations in which an input is sometimes a
	 * constant and other times a random variable.
	 * </p>
	 * <p>
	 * getValue(...) is typically used for input parameters and not used for
	 * intermediate parameters in multi-step calculations.
	 * </p>
	 *
	 * @param label
	 * @param point
	 * @return double
	 */
	final public double getValue(final Object label, final RealVector point) {
		return point.getEntry(inputIndex(label));
	}

	/**
	 * A value either with or without uncertainty has been defined for this label.
	 *
	 * @param label
	 * @return true if a value has been defined.
	 */
	final public boolean hasValue(final Object label) {
		return (inputIndex(label) != -1);
	}

	@Override
	public int hashCode() {
		return Objects.hash(mInputLabels, mOutputLabels);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final LabeledMultivariateJacobianFunction other = (LabeledMultivariateJacobianFunction) obj;
		return Objects.equals(mInputLabels, other.mInputLabels) && //
				Objects.equals(mOutputLabels, other.mOutputLabels);
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE: {
			return HTML
					.escape("V[" + getOutputLabels().size() + " values]=F(" + getInputLabels().size() + " arguments)");
		}
		case NORMAL:
			return HTML.escape(toString());
		default:
		case VERBOSE: {
			final Table t = new Table();
			t.addRow(Table.th(HTML.escape(toString()), 2));
			{
				final StringBuffer sb = new StringBuffer();
				for (final Object label : getOutputLabels()) {
					if (sb.length() != 0)
						sb.append("<br/>");
					sb.append(HTML.toHTML(label, Mode.TERSE));
				}
				t.addRow(Table.td("Ouputs"), Table.td(sb.toString()));
			}
			{
				final StringBuffer sb = new StringBuffer();
				for (final Object label : getInputLabels()) {
					if (sb.length() != 0)
						sb.append("<br/>");
					sb.append(HTML.toHTML(label, Mode.TERSE));
				}
				t.addRow(Table.td("Inputs"), Table.td(sb.toString()));
			}
			return t.toHTML(Mode.NORMAL);
		}
		}
	}
}
