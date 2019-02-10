package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

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
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.utility.FastIndex;

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

	/***
	 * A set of constant values for use evaluating the MulitvariateJacobianFunction.
	 */
	private final Map<Object, Double> mConstants = new HashMap<>();

	/**
	 * Check labels are only used once.
	 *
	 * @param labels
	 */
	public void validateLabels(final List<? extends Object> labels) {
		for (int i = 0; i < labels.size(); ++i)
			for (int j = i + 1; j < labels.size(); ++j)
				if (labels.get(i).equals(labels.get(j)))
					throw new RuntimeException("The label " + labels.get(i).toString() + " is duplicated.");
	}

	public LabeledMultivariateJacobianFunction( //
			final List<? extends Object> inputLabels, //
			final List<? extends Object> outputLabels //
	) {
		validateLabels(inputLabels);
		assert inputLabels != outputLabels;
		mInputLabels = new FastIndex<>(inputLabels);
		validateLabels(outputLabels);
		mOutputLabels = new FastIndex<>(outputLabels);
	}

	/**
	 * Returns an array consisting of the labels associated with the n input random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	public List<? extends Object> getInputLabels() {
		return mInputLabels;
	}

	/**
	 * Returns an array consisting of the labels associated with the m output random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	public List<? extends Object> getOutputLabels() {
		return mOutputLabels;
	}

	public Object getOutputLabel(int idx) {
		return mOutputLabels.get(idx);
	}

	public Object getInputLabel(int idx) {
		return mInputLabels.get(idx);
	}

	/**
	 * Number of input random variables expected.
	 *
	 * @return int
	 */
	public int getInputDimension() {
		return mInputLabels.size();
	}

	public int getConstantDimension() {
		return mConstants.size();
	}

	/**
	 * Number of output random variable produced.
	 *
	 * @return int
	 */
	public int getOutputDimension() {
		return getOutputLabels().size();
	}

	/**
	 * Returns the index of the input variable identified by the specified Object
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	public int inputIndex(final Object label) {
		return mInputLabels.indexOf(label);
	}

	/**
	 * Returns the index of the output variable identified by the specified Object
	 *
	 * @param label
	 * @return int Index or -1 for not found.
	 */
	public int outputIndex(final Object label) {
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
	public Pair<RealVector, RealMatrix> evaluate(final RealVector x) {
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
	public RealVector compute(final RealVector inp) {
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
	public HashMap<? extends Object, UncertainValue> getOutputValues(final UncertainValues uvs)
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
			final UncertainValues uvs, //
			final double tol//
	) throws ArgumentException {
		final UncertainValues ordered = UncertainValues.build(this.getInputLabels(), uvs);
		final Pair<RealVector, RealMatrix> pvm = evaluate(ordered.getValues());
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<Object, UncertainValue> res = new HashMap<>();
		final List<? extends Object> inLabels = getInputLabels();
		final List<? extends Object> outLabels = getOutputLabels();
		for (int i = 0; i < outLabels.size(); ++i)
			res.put(outLabels.get(i), new UncertainValue(vals.getEntry(i)));
		for (int inIdx = 0; inIdx < inLabels.size(); ++inIdx) {
			final Object inLabel = inLabels.get(inIdx);
			final UncertainValues zeroed = UncertainValues.zeroBut(inLabel, ordered);
			final RealMatrix covLabel = jac.multiply(zeroed.getCovariances()).multiply(jac.transpose());
			for (int outIdx = 0; outIdx < outLabels.size(); ++outIdx) {
				final UncertainValue val = res.get(outLabels.get(outIdx));
				final double sigma = Math.sqrt(covLabel.getEntry(outIdx, outIdx));
				if (sigma > Math.abs(tol * val.doubleValue()))
					val.assignComponent(inLabel, sigma);
			}
		}
		return res;
	}

	/**
	 * Implements a {@link LabeledMultivariateJacobianFunction} that normalizes the
	 * values in the inLabels and returns them as the outLabels.
	 *
	 * @param inLabels
	 * @param outLabels
	 * @return {@link LabeledMultivariateJacobianFunction}
	 */
	public static LabeledMultivariateJacobianFunction normalize( //
			final List<? extends Object> inLabels, //
			final List<? extends Object> outLabels //
	) {
		return new LabeledMultivariateJacobianFunction(inLabels, outLabels) {

			private double computePartial(final int f, final int v, final double sum, final RealVector point) {
				double res = (f == v ? 1.0 / sum : 0.0);
				if (f != v)
					res -= point.getEntry(v) / (sum * sum);
				return res;
			}

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				double sum = 0.0;
				for (int i = 0; i < point.getDimension(); ++i)
					sum += point.getEntry(i);
				final double[][] jac = new double[point.getDimension()][point.getDimension()];
				for (int f = 0; f < point.getDimension(); ++f)
					for (int v = 0; v < point.getDimension(); ++v)
						jac[f][v] = computePartial(f, v, sum, point);
				return Pair.create(point.mapMultiply(1.0 / sum), MatrixUtils.createRealMatrix(jac));
			}

		};
	}

	/**
	 * Implements a {@link LabeledMultivariateJacobianFunction} that sum the values
	 * in the inLabels and returns them as the outLabel.
	 *
	 * @param inLabels
	 * @param outLabels
	 * @return {@link LabeledMultivariateJacobianFunction}
	 */
	public static LabeledMultivariateJacobianFunction sum( //
			final List<? extends Object> inLabels, //
			final Object outLabel //
	) {
		return new LabeledMultivariateJacobianFunction(inLabels, Collections.singletonList(outLabel)) {

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				RealVector rv = new ArrayRealVector(1);
				RealMatrix rm = MatrixUtils.createRealMatrix(1, point.getDimension());

				double sum = 0.0;
				for (int i = 0; i < point.getDimension(); ++i) {
					sum += point.getEntry(i);
					rm.setEntry(0, i, 1.0);
				}
				rv.setEntry(0, sum);
				return Pair.create(rv, rm);
			}

		};
	}

	/**
	 * Implements a linear function of the inLabels quantities which is returned as
	 * outLabel.
	 *
	 * @param inLabels Labels for the input values
	 * @param coeffs   A vector of linear parameters of length inLabels.size()
	 * @param outLabel The label for the output value.
	 * @return {@link LabeledMultivariateJacobianFunction}
	 */
	public static LabeledMultivariateJacobianFunction linear(final List<? extends Object> inLabels,
			final RealVector coeffs, final Object outLabel) {
		return new LabeledMultivariateJacobianFunction(inLabels, Collections.singletonList(outLabel)) {

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				RealVector rv = new ArrayRealVector(1);
				RealMatrix rm = MatrixUtils.createRealMatrix(1, point.getDimension());
				for (int v = 0; v < point.getDimension(); ++v)
					rm.setEntry(0, v, coeffs.getEntry(v));
				rv.setEntry(0, point.dotProduct(coeffs));
				return Pair.create(rv, rm);
			}
		};
	}

	/**
	 * Initializes the constant labels with the specified values. The label for
	 * value <code>vals.getEntry(i)</code> is <code>list.get(i)</code>. Checks first
	 * to see is the label is being used as an input variable and won't define it as
	 * a constant if it is an input item.
	 *
	 * @param list
	 * @param vals
	 */
	public void initializeConstants(final List<? extends Object> list, final RealVector vals) {
		for (int i = 0; i < list.size(); ++i)
			if (inputIndex(list.get(i)) == -1)
				mConstants.put(list.get(i), vals.getEntry(i));
	}

	/**
	 * Initializes the constant labels with the associated values.
	 *
	 * @param mod Map&lt;Object,Double&gt; where Object is a label
	 */
	public void initializeConstants(final Map<Object, Double> mod) {
		mConstants.putAll(mod);
	}

	/**
	 * Returns the constant value associated with <code>label</code>.
	 *
	 * @param label
	 * @return double
	 */
	public double getConstant(final Object label) {
		assert inputIndex(label) == -1 : "Label " + label + " is a variable.";
		assert mConstants.containsKey(label) : "Label " + label + " is not a constant";
		return mConstants.get(label).doubleValue();
	}

	/**
	 * Check whether <code>label</code> is defined as a constant value.
	 *
	 * @param label A label
	 * @return true if <code>label</code> is initialized as a constant, false
	 *         otherwise.
	 */
	public boolean isConstant(final Object label) {
		return mConstants.containsKey(label);
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
	public double getValue(final Object label, final RealVector point) {
		final int p = inputIndex(label);
		if (p != -1) {
			// assert !isConstant(label);
			return point.getEntry(p);
		} else {
			assert isConstant(label) : //
			"Can't find the constant " + label + " in " + toString();
			return getConstant(label);
		}
	}

	/**
	 * A value either with or without uncertainty has been defined for this label.
	 * 
	 * @param label
	 * @return true if a value has been defined.
	 */
	public boolean hasValue(final Object label) {
		return (inputIndex(label) != -1) || //
				isConstant(label);
	}

	public RealMatrix computeDelta( //
			final RealVector inp, //
			final RealVector dinp) {
		assert inp.getDimension() == getInputDimension();
		final RealMatrix rm = NullableRealMatrix.build(getOutputDimension(), getInputDimension());
		for (int c = 0; c < getInputDimension(); ++c) {
			final RealVector pt0 = new ArrayRealVector(inp), pt1 = new ArrayRealVector(inp);
			final double deltaX = Math.abs(dinp.getEntry(c));
			pt0.setEntry(c, pt0.getEntry(c) + 0.5 * deltaX);
			pt1.setEntry(c, pt1.getEntry(c) - 0.5 * deltaX);
			final RealVector output0 = compute(pt0), output1 = compute(pt1);
			for (int r = 0; r < getOutputDimension(); ++r)
				rm.setEntry(r, c, (output0.getEntry(r) - output1.getEntry(r)) / deltaX);
		}
		return rm;
	}

	/**
	 * Returns an unmodifiable view of the map of labels and the associated constant
	 * values.
	 *
	 * @return Map<Object, Double>
	 */
	public Map<Object, Double> getConstants() {
		return Collections.unmodifiableMap(mConstants);
	}

	@Override
	public int hashCode() {
		return Objects.hash(mInputLabels, mOutputLabels, mConstants);
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
				Objects.equals(mOutputLabels, other.mOutputLabels) && //
				Objects.equals(mConstants, other.mConstants);
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE: {
			return HTML.escape("V[" + getOutputLabels().size() + " values]=F(" + getInputLabels().size()
					+ " arguments, " + getConstants().size() + " constants)");
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
			{
				final StringBuffer sb = new StringBuffer();
				for (final Object label : getConstants().keySet()) {
					if (sb.length() != 0)
						sb.append("<br/>");
					sb.append(HTML.toHTML(label, Mode.TERSE));
				}
				if (sb.length() > 0)
					t.addRow(Table.td("Constants"), Table.td(sb.toString()));
			}
			return t.toHTML(Mode.NORMAL);
		}
		}
	}
}
