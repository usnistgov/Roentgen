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
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.Constraint;
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
abstract public class ExplicitMeasurementModel<G, H> //
		implements MultivariateJacobianFunction, IToHTML //
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
	) //
			throws ArgumentException {
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
	 * @throws ArgumentException
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

	public void addAdditionalInputs(
			final Map<? extends G, Double> inputs
	) {
		for (Map.Entry<? extends G, Double> me : inputs.entrySet())
			setAdditionalInput(me.getKey(), me.getValue());
	}

	public void addConstraints(
			final Map<? extends G, Constraint> constraints
	) {
		mConstraints.putAll(constraints);
	}

	private void applyConstraints(
			final RealVector point
	) {
		if (!mConstraints.isEmpty()) {
			for (Entry<G, Constraint> con : mConstraints.entrySet()) {
				int idx = inputIndex(con.getKey());
				if (idx >= 0)
					point.setEntry(idx, con.getValue().limit(point.getEntry(idx)));
			}
		}
	}

	public void clearAdditionalInputs() {
		mAdditionalInputs.clear();
	}

	public void clearConstraints() {
		mConstraints.clear();
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
	final public RealVector compute(
			final RealVector inp
	) {
		if (this instanceof ILabeledMultivariateFunction) {
			applyConstraints(inp);
			return ((ILabeledMultivariateFunction<?, ?>) this).optimized(inp);
		} else
			return value(inp).getFirst();
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
	 *
	 * @param x
	 * @return Pair&lt;{@link RealVector}, {@link RealMatrix}&gt; As from a call to
	 *         <code>value(x)</code>.
	 */
	final public Pair<RealVector, RealMatrix> evaluate(
			final RealVector x
	) {
		if (x.getDimension() != getInputDimension())
			throw new DimensionMismatchException(x.getDimension(), getInputDimension());
		for (final Map.Entry<G, Constraint> con : mConstraints.entrySet()) {
			final int idx = inputIndex(con.getKey());
			if (idx >= 0)
				x.setEntry(idx, con.getValue().limit(x.getEntry(idx)));
		}
		final Pair<RealVector, RealMatrix> res = value(x);
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
	 * Number of input random variables expected.
	 *
	 * @return int
	 */
	final public int getInputDimension() {
		return mInputLabels.size();
	}

	final public G getInputLabel(
			final int idx
	) {
		return mInputLabels.get(idx);
	}

	/**
	 * Returns an array consisting of the labels associated with the M input random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	final public List<G> getInputLabels() {
		return Collections.unmodifiableList(mInputLabels);
	}

	/**
	 * Number of output random variable produced.
	 *
	 * @return int
	 */
	final public int getOutputDimension() {
		return mOutputLabels.size();
	}

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
	 * @param label
	 * @return true if a value has been defined.
	 */
	final public boolean hasValue(
			final G label
	) {
		return (inputIndex(label) != -1);
	}

	/**
	 * Returns the index of the input variable identified by the specified instance
	 * of H.
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	final public int inputIndex(
			final Object label
	) {
		return mInputLabels.indexOf(label);
	}

	/**
	 * Returns the index of the output function identified by the specified instance
	 * of H.
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	final public int outputIndex(
			final Object label
	) {
		return mOutputLabels.indexOf(label);
	}

	/**
	 * Specifies a value to associate with an input label. Works with the
	 * getArg(...) function. If there is not a value associated with the specified
	 * (inputIndex(lbl)==-1) then the additional inputs will be checked for a value.
	 * Constraints are applied to the input value.
	 *
	 * @param lbl
	 * @param value
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
	 * @param lbl
	 * @param con
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
	 * A helper to build an zeroed matrix to contain the result Jacobian in
	 * value(...)
	 *
	 * @return {@link RealMatrix}
	 */
	final protected RealMatrix buildJacobian() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
	}

	/**
	 * A helper to build a zeroed vector to contain the result for value(...)
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
			sb.append("\"Args["+toString()+"]\"");
			for (int i = 0; i < getInputDimension(); ++i) {
				sb.append(",");
				sb.append(getInputLabel(i));
			}
			sb.append("\n");
			sb.append("\"Args["+toString()+"]\"");
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
	 * @param inLabel A G class label
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
	 * Sets the value associated with the specified row (outLabel) and column
	 * (inLabel) in the Jacobian matrix in jacobian to value. Typically used to
	 * implement the value(...) function.
	 *
	 * @param inLabel  G
	 * @param outLabel H
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
	 * @param outLabel H A H-class label
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
	 * For internal use.
	 * 
	 * @param inputs
	 */
	void applyAdditionalInputs(
			final Map<Object, Double> inputs
	) {
		mAdditionalInputs.putAll(inputs);
	}

}
