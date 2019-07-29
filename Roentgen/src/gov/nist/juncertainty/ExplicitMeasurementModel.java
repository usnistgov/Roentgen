package gov.nist.juncertainty;

import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Set;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ValueAndJacobianFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.juncertainty.utility.FastIndex;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.Constraint;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * {@link ExplicitMeasurementModel} is an abstract class to represent
 * calculational steps in an explicit measurement model. Instances of this class
 * may represent a full measurement model or a step in a larger measurement
 * model.
 * </p>
 *
 * <p>
 * {@link ExplicitMeasurementModel} extends {@link MultivariateJacobianFunction}
 * by adding support for input and output variable labels. The
 * {@link MultivariateJacobianFunction} interface defines a single method that
 * calculates a {@link Pair}&lt;{@link RealVector}, {@link RealMatrix}&gt;. The
 * {@link RealVector} represents the values associated with a vector function
 * and the {@link RealMatrix} represents the Jacobian for the vector function
 * with respect to the input variables. The {@link ExplicitMeasurementModel}
 * class makes it easier to keep track of the input and output variables
 * particularly in circumstances in which the number and type of the input and
 * output variables may change dependent on the nature of the calculation.
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
abstract public class ExplicitMeasurementModel<G, H> //
		implements MultivariateJacobianFunction, ValueAndJacobianFunction, IToHTML //
{

	public static PrintStream sDump = null;

	/**
	 * Check labels are only used once.
	 *
	 * @param inLabels
	 * @param outLabels
	 * @throws ArgumentException
	 */
	private static <G, H> void validateLabels(
			final Collection<? extends G> inLabels, //
			final Collection<? extends H> outLabels
	) throws ArgumentException {
		final Set<? extends G> gset = new HashSet<>(inLabels);
		for (final Object lbl : inLabels)
			if (!gset.remove(lbl))
				throw new ArgumentException("The input label " + lbl + " is duplicated.");
		final Set<? extends H> hset = new HashSet<>(outLabels);
		for (final Object lbl : outLabels)
			if (!hset.remove(lbl))
				throw new ArgumentException("The input label " + lbl + " is duplicated.");
		gset.retainAll(hset);
		if (gset.size() > 0) {
			throw new ArgumentException("The labels " + gset + " are in both the inputs and the outputs.");
		}
	}

	/***
	 * Unique object labels identifying each of the random variable arguments to the
	 * functions
	 */
	private final List<G> mInputLabels;

	/***
	 * Unique object labels identifying each of the functions (in order)
	 */
	private final List<H> mOutputLabels;

	/***
	 * A mechanism to ensure that input values remain with constrained ranges.
	 */
	private final Map<G, Constraint> mConstraints;

	/***
	 * An alternative way to input quantities into the model
	 */
	private final Map<Object, Double> mAdditionalInputs;

	/**
	 * Constructs an instance of the LabeledMultivariateJacobianFunction class for a
	 * N dimensional function of M input values.
	 *
	 * @param inputLabels  A List containing labels for M input variables
	 * @param outputLabels A list containing labels for N output values
	 * @throws ArgumentException For repeated labels
	 */
	public ExplicitMeasurementModel(
			final List<G> inputLabels, //
			final List<H> outputLabels //
	) throws ArgumentException //
	{
		// assert inputLabels != outputLabels;
		mInputLabels = new FastIndex<>(inputLabels);
		mOutputLabels = new FastIndex<>(outputLabels);
		mConstraints = new HashMap<>();
		mAdditionalInputs = new HashMap<>();
		validateLabels(inputLabels, outputLabels);
	}

	/**
	 * Defines a set of additional input values which are may be used by the
	 * calculation but for which the Jacobian matrix elements are not computed.
	 *
	 * @param inputs A {@link Map} containing labels and the associated value as a {@link Double}.
	 */
	public void addAdditionalInputs(
			final Map<? extends G, Double> inputs
	) {
		for (final Map.Entry<? extends G, Double> me : inputs.entrySet())
			setAdditionalInput(me.getKey(), me.getValue());
	}

	/**
	 * Declares a set of {@link Constraint}s that limits the range of values which
	 * the associated label can take on. For example, the value associated with a
	 * label may be constrained to be positive definite using a
	 * {@link Constraint.Range} instance.
	 *
	 * @param constraints A map from an input label to a {@link Constraint}
	 *                    instance.
	 */
	public void addConstraints(
			final Map<? extends G, Constraint> constraints
	) {
		mConstraints.putAll(constraints);
	}

	private RealVector applyConstraints(
			final RealVector point
	) {
		if (!mConstraints.isEmpty()) {
			final RealVector res = point.copy();
			for (final Entry<G, Constraint> con : mConstraints.entrySet()) {
				final int idx = inputIndex(con.getKey());
				if (idx >= 0)
					res.setEntry(idx, con.getValue().limit(point.getEntry(idx)));
			}
			return res;
		}
		return point;
	}

	public void clearAdditionalInputs() {
		mAdditionalInputs.clear();
	}

	public void clearConstraints() {
		mConstraints.clear();
	}

	/**
	 * Applies {@link Constraint}s to the input variables and then
	 * calls computeValue(...).
	 *
	 * @param inp Evaluation point
	 * @return {@link RealVector} The result of computeValue(...)
	 */
	final public RealVector compute(
			final RealVector inp
	) {
		return computeValue(applyConstraints(inp).toArray());
	}

	@Override
	public boolean equals(
			final Object obj
	) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final ExplicitMeasurementModel<?, ?> other = (ExplicitMeasurementModel<?, ?>) obj;
		return Objects.equals(mInputLabels, other.mInputLabels) && //
				Objects.equals(mOutputLabels, other.mOutputLabels);
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
	 * <p>
	 * Also applies any applicable {@link Constraint}s to the input values.
	 * </p>
	 *
	 * @param point A {@link RealVector} containing the values at which to evaluated the function value(...).
	 * @return Pair&lt;{@link RealVector}, {@link RealMatrix}&gt; As from a call to
	 *         <code>value(x)</code>.
	 */
	final public Pair<RealVector, RealMatrix> evaluate(
			final RealVector point
	) {
		if (point.getDimension() != getInputDimension())
			throw new DimensionMismatchException(point.getDimension(), getInputDimension());
		final Pair<RealVector, RealMatrix> res = value(applyConstraints(point));
		if (res.getFirst().getDimension() != getOutputDimension())
			throw new DimensionMismatchException(res.getFirst().getDimension(), getOutputDimension());
		if (res.getSecond().getRowDimension() != getOutputDimension())
			throw new DimensionMismatchException(res.getSecond().getRowDimension(), getOutputDimension());
		if (res.getSecond().getColumnDimension() != getInputDimension())
			throw new DimensionMismatchException(res.getSecond().getColumnDimension(), getInputDimension());
		return res;
	}

	public Map<Object, Double> getAdditionalInputs() {
		return Collections.unmodifiableMap(mAdditionalInputs);
	}

	public Map<G, Constraint> getConstraints() {
		return Collections.unmodifiableMap(mConstraints);
	}

	/**
	 * Number of input random variables.
	 *
	 * @return int
	 */
	final public int getInputDimension() {
		return mInputLabels.size();
	}

	/**
	 * Returns the <code>idx</code><sup>th</sup> input label.
	 *
	 * @param idx The integer index of the input label
	 * @return G The input label at index idx
	 */
	final public G getInputLabel(
			final int idx
	) {
		return mInputLabels.get(idx);
	}

	/**
	 * Returns a unmodifiable list of the labels associated with the M input random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<G> getInputLabels() {
		return Collections.unmodifiableList(mInputLabels);
	}

	/**
	 * Returns the number of output values.
	 *
	 * @return int The number of output values
	 */
	final public int getOutputDimension() {
		return mOutputLabels.size();
	}

	/**
	 * Returns the <code>idx</code><sup>th</sup> output label.
	 *
	 * @param idx The integer index of the output label
	 * @return G The output label at index idx
	 */
	final public H getOutputLabel(
			final int idx
	) {
		return mOutputLabels.get(idx);
	}

	/**
	 * Returns an array consisting of the labels associated with the N output
	 * function values.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<H> getOutputLabels() {
		return Collections.unmodifiableList(mOutputLabels);
	}

	@Override
	public int hashCode() {
		return Objects.hash(mInputLabels, mOutputLabels);
	}

	/**
	 * A value either with or without uncertainty has been defined for this label.
	 *
	 * @param label A label
	 * @return true if a value has been defined.
	 */
	final public boolean hasValue(
			final G label
	) {
		return (inputIndex(label) != -1);
	}

	/**
	 * Returns the index of the input variable identified by the specified
	 * {@link Object}.
	 *
	 * @param label A label
	 * @return int An integer index or -1 if not found.
	 */
	final public int inputIndex(
			final Object label
	) {
		return mInputLabels.indexOf(label);
	}

	/**
	 * Returns the index of the output value identified by the specified
	 * {@link Object}.
	 *
	 * @param label A label
	 * @return int An integer index or -1 if not found.
	 */
	final public int outputIndex(
			final Object label
	) {
		return mOutputLabels.indexOf(label);
	}

	/**
	 * Specifies a value to associate with an input label. Works in cooperation with
	 * the getArg(...) function. If there is not a value associated with the
	 * specified (inputIndex(lbl)==-1) then the additional inputs will be checked
	 * for a value. Constraints are applied to the input value.
	 *
	 * @param lbl A label
	 * @param value The value to associate with the label
	 */
	public void setAdditionalInput(
			final G lbl, //
			final double value
	) {
		if (mConstraints.containsKey(lbl))
			mAdditionalInputs.put(lbl, mConstraints.get(lbl).limit(value));
		else
			mAdditionalInputs.put(lbl, value);
	}

	/**
	 * Places the specified {@link Constraint} on the specified input label.
	 *
	 * @param lbl A label
	 * @param con A {@link Constraint} derived class
	 */
	public void setConstraint(
			final G lbl, //
			final Constraint con
	) {
		mConstraints.put(lbl, con);
	}

	@Override
	public String toHTML(
			final Mode mode
	) {
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

	/**
	 * A helper to build an zeroed matrix of the correct dimensions to contain the
	 * result Jacobian. Useful for implementing value(...).
	 *
	 * @return {@link RealMatrix}
	 */
	final protected RealMatrix buildJacobian() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
	}

	/**
	 * A helper to build a zeroed vector of the correct dimension to contain the
	 * result vector. Useful for implementing value(...).
	 *
	 * @return {@link RealVector}
	 */
	final protected RealVector buildResult() {
		return new ArrayRealVector(getOutputDimension());
	}

	protected void dumpArguments(
			final RealVector point, //
			final ExplicitMeasurementModel<?, ?> parent
	) {
		if (sDump != null) {
			final StringBuffer sb = new StringBuffer();
			final NumberFormat nf = new HalfUpFormat("0.00E0");
			sb.append("\"Args[" + toString() + "]\"");
			for (int i = 0; i < getInputDimension(); ++i) {
				sb.append(",");
				sb.append(getInputLabel(i));
			}
			sb.append("\n");
			sb.append("\"Args[" + toString() + "]\"");
			for (int i = 0; i < getInputDimension(); ++i) {
				sb.append(",");
				sb.append(nf.format(point.getEntry(i)));
			}
			sb.append("\n");
			sDump.append(sb);
		}
	}

	/**
	 * Checks two places to see if there is a value associated with the specified
	 * label. First it checks the input {@link RealVector} (using inputIndex(label))
	 * and if a value is available here it returns this value. Second it checks the
	 * getAlternativeInputs() to see if a value is available here. The function
	 * throws a NullPointerException if neither place holds a value. Typically used
	 * to implement the value(...) function.
	 *
	 * @param inLabel A label
	 * @param point   The argument to the value(...) function
	 * @return double The value at point.getEntry(inputIndex(inLabel))
	 */
	final protected double getArg(
			final G inLabel, //
			final RealVector point
	) {
		final int idx = inputIndex(inLabel);
		if (idx >= 0)
			return point.getEntry(idx);
		else {
			assert mAdditionalInputs.containsKey(inLabel) : //
			"Missing " + inLabel;
			return mAdditionalInputs.get(inLabel).doubleValue();

		}
	}

	/**
	 * Checks two places to see if there is a value associated with the specified
	 * label. First it checks the input {@link RealVector} (using inputIndex(label))
	 * and if a value is available here it returns this value. Second it checks the
	 * getAlternativeInputs() to see if a value is available here. The function
	 * throws a NullPointerException if neither place holds a value. Typically used
	 * to implement the value(...) function.
	 *
	 * @param inLabel A label
	 * @param point   The argument to the value(...) function
	 * @return double The value at point.getEntry(inputIndex(inLabel))
	 */
	final protected double getArg(
			final G inLabel, //
			final double[] point
	) {
		final int idx = inputIndex(inLabel);
		if (idx >= 0)
			return point[idx];
		else {
			assert mAdditionalInputs.containsKey(inLabel) : //
			"Missing " + inLabel;
			return mAdditionalInputs.get(inLabel).doubleValue();

		}
	}

	/**
	 * <p>
	 * Sets the value associated with the specified row (outLabel) and column
	 * (inLabel) in the Jacobian matrix in jacobian to value. Typically used to
	 * implement the value(...) function.
	 * </p>
	 * <p>
	 * jacobian[outLabel,inLabel] = value = &delta;H/&delta;G
	 * </p>
	 *
	 * @param inLabel  G An input label
	 * @param outLabel H A output label
	 * @param jacobian {@link RealMatrix}
	 * @param value    double The value to be assigned to
	 *                 jacobian[outputIndex(outLabel),inputIndex(inLabel)]
	 */
	final protected void setJacobian(
			final G inlabel, //
			final H outLabel, //
			final RealMatrix jacobian, //
			final double value
	) {
		final int idx = inputIndex(inlabel);
		if (idx >= 0)
			jacobian.setEntry(outputIndex(outLabel), idx, value);
	}

	/**
	 * Sets the value associated with the specified output variable label to the
	 * result {@link RealVector} point which is to be returned by the value(...)
	 * function. Typically used to implement the value(...) function.
	 *
	 * @param outLabel H An output label
	 * @param result   {@link RealVector}
	 * @param value    double The value to be assigned to
	 *                 result[outputIndex(outLabel)]
	 */
	final protected void setResult(
			final H outLabel, final RealVector result, final double value
	) {
		result.setEntry(outputIndex(outLabel), value);
	}

	/**
	 * For internal use. Helps to transfer additional inputs from one step in a
	 * calculation to an embedded step.
	 *
	 * @param inputs A Map from label to value as a Double.
	 */
	protected void applyAdditionalInputs(
			final Map<Object, Double> inputs
	) {
		mAdditionalInputs.putAll(inputs);
	}

	/**
	 * Avoid using this function use value(...) instead.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.ValueAndJacobianFunction#computeJacobian(double[])
	 */
	@Override
	public RealMatrix computeJacobian(
			final double[] params
	) {
		return value(new ArrayRealVector(params)).getSecond();
	}

	/**
	 * Implement this method in derived classes to optimize the calculation of only
	 * the values.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.ValueAndJacobianFunction#computeValue(double[])
	 */
	@Override
	public RealVector computeValue(
			final double[] params
	) {
		return value(new ArrayRealVector(params)).getFirst();
	}

}
