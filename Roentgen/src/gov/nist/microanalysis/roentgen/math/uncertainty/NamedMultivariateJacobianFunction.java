package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Arrays;
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
 * @author nritchie
 * @version $Rev: $
 */
abstract public class NamedMultivariateJacobianFunction implements MultivariateJacobianFunction, IToHTML {

	/***
	 * Unique object tags identifying each of the random variable arguments to the
	 * functions
	 */
	private final List<? extends Object> mInputTags;
	/***
	 * Unique object tags identifying each of the functions (in order)
	 */
	private final List<? extends Object> mOutputTags;

	/***
	 * A set of constant values for use evaluating the MulitvariateJacobianFunction.
	 */
	private final Map<Object, Double> mConstants = new HashMap<>();

	/**
	 * Check tags are only used once.
	 *
	 * @param tags
	 */
	public void validateTags(final List<? extends Object> tags) {
		for (int i = 0; i < tags.size(); ++i)
			for (int j = i + 1; j < tags.size(); ++j)
				if (tags.get(i).equals(tags.get(j)))
					throw new RuntimeException("The tag " + tags.get(i).toString() + " is duplicated.");
	}

	public NamedMultivariateJacobianFunction(final List<? extends Object> inputTags,
			final List<? extends Object> outputTags) {
		validateTags(inputTags);
		assert inputTags != outputTags;
		mInputTags = Collections.unmodifiableList(inputTags);
		validateTags(outputTags);
		mOutputTags = Collections.unmodifiableList(outputTags);
	}

	/**
	 * Validate that the specified set of input tags matches the
	 *
	 * @param tags
	 * @return
	 */
	public boolean validateInputs(final List<? extends Object> tags) {
		for (int i = 0; i < mInputTags.size(); ++i)
			if (tags.get(i).equals(mInputTags.get(i)))
				return false;
		return true;
	}

	/**
	 * Returns an array consisting of the tags associated with the n input random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	public List<? extends Object> getInputTags() {
		return mInputTags;
	}

	/**
	 * Returns an array consisting of the tags associated with the m output random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	public List<? extends Object> getOutputTags() {
		return mOutputTags;
	}

	/**
	 * Number of input random variables expected.
	 *
	 * @return int
	 */
	public int getInputDimension() {
		return mInputTags.size();
	}

	/**
	 * Number of output random variable produced.
	 *
	 * @return int
	 */
	public int getOutputDimension() {
		return getOutputTags().size();
	}

	/**
	 * Returns the index of the input variable identified by the specified Object
	 *
	 * @param tag
	 * @return int Index or -1 for not found.
	 */
	public int inputIndex(final Object tag) {
		for (int i = 0; i < mInputTags.size(); ++i)
			if (mInputTags.get(i).equals(tag))
				return i;
		return -1;
	}

	/**
	 * Returns the index of the output variable identified by the specified Object
	 *
	 * @param tag
	 * @return int Index or -1 for not found.
	 */
	public int outputIndex(final Object tag) {
		for (int i = 0; i < mOutputTags.size(); ++i)
			if (mOutputTags.get(i).equals(tag))
				return i;
		return -1;
	}

	/**
	 * Only writes <code>value</code> to the Jacobian matrix <code>jacob</code> if
	 * <code>tag</code> is an input variable (not a constant.)
	 *
	 *
	 * @param row
	 * @param tag
	 * @param value
	 * @param jacob
	 */
	public void writeJacobian(final int row, final Object tag, final double value, final RealMatrix jacob) {
		final int p = inputIndex(tag);
		if (p != -1)
			jacob.setEntry(row, p, value);
	}

	/**
	 * Extracts from <code>point</code> representing an argument to
	 * <code>nmvj.compute(point)</code> the input RealArray that is suitable as the
	 * argument to <code>this.comptue(...)</code>.
	 *
	 * @param nmvj  The outer {@link NamedMultivariateJacobianFunction}
	 * @param point A RealVector of length nmvj.getInputDimension()
	 * @return RealVector of length this.getInputDimension() containing the
	 *         appropriate input values from point.
	 */
	public RealVector extractArgument(final NamedMultivariateJacobianFunction nmvj, final RealVector point) {
		assert point.getDimension() == nmvj.getInputDimension();
		final int dim = getInputDimension();
		final RealVector res = new ArrayRealVector(dim);
		final List<? extends Object> tags = getInputTags();
		for (int i = 0; i < dim; ++i) {
			final int idx = nmvj.inputIndex(tags.get(i));
			assert idx != -1 : "Can't find " + tags.get(i) + " in the arguments to " + nmvj.toString();
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
	 * then checks the output RealVector and RealMatrix to ensure that they are all
	 * the correct dimensions.
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

	public RealVector compute(final RealVector inp) {
		if (this instanceof INamedMultivariateFunction)
			return ((INamedMultivariateFunction) this).optimized(inp);
		else
			return evaluate(inp).getFirst();
	}

	public static NamedMultivariateJacobianFunction identity(final List<? extends Object> inTags) {
		return new NamedMultivariateJacobianFunction(inTags, inTags) {
			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				return Pair.create(point, MatrixUtils.createRealIdentityMatrix(point.getDimension()));
			}
		};
	}

	public HashMap<? extends Object, UncertainValue> getOutputValues(final UncertainValues uvs) {
		final Pair<RealVector, RealMatrix> pvm = evaluate(uvs.getValues());
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<Object, UncertainValue> res = new HashMap<>();
		final List<? extends Object> inTags = getInputTags();
		final List<? extends Object> outTags = getOutputTags();
		for (int i = 0; i < outTags.size(); ++i)
			res.put(outTags.get(i), new UncertainValue(vals.getEntry(i)));
		for (int inIdx = 0; inIdx < inTags.size(); ++inIdx) {
			final Object inTag = inTags.get(inIdx);
			final UncertainValues uvsTag = UncertainValues.zeroBut(inTag, uvs);
			final RealMatrix covTag = jac.multiply(uvsTag.getCovariances()).multiply(jac.transpose());
			for (int outIdx = 0; outIdx < outTags.size(); ++outIdx) {
				final UncertainValue val = res.get(outTags.get(outIdx));
				final double sigma = Math.sqrt(covTag.getEntry(outIdx, outIdx));
				if (sigma > Math.abs(1.0e-8 * val.doubleValue()))
					val.assignComponent(inTag, sigma);
			}
		}
		return res;
	}

	public static NamedMultivariateJacobianFunction normalize(final List<? extends Object> inTags,
			final List<? extends Object> outTags) {
		return new NamedMultivariateJacobianFunction(inTags, outTags) {

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

	public static NamedMultivariateJacobianFunction sum(final List<? extends Object> inTags, final Object outTag) {
		return new NamedMultivariateJacobianFunction(inTags, Collections.singletonList(outTag)) {

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				double sum = 0.0;
				for (int i = 0; i < point.getDimension(); ++i)
					sum += point.getEntry(i);
				final double[] data = new double[point.getDimension()];
				Arrays.fill(data, 1.0);
				return Pair.create(MatrixUtils.createRealVector(new double[] { sum }),
						MatrixUtils.createRowRealMatrix(data));
			}

		};
	}

	public static NamedMultivariateJacobianFunction linear(final List<? extends Object> inTags, final RealVector coeffs,
			final Object outTag) {
		return new NamedMultivariateJacobianFunction(inTags, Collections.singletonList(outTag)) {

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				final double dot = point.dotProduct(coeffs);
				final double[][] jac = new double[1][point.getDimension()];
				for (int v = 0; v < point.getDimension(); ++v)
					jac[0][v] = coeffs.getEntry(v);
				final RealVector vals = MatrixUtils.createRealVector(new double[] { dot });
				return Pair.create(vals, MatrixUtils.createRealMatrix(jac));
			}
		};
	}
	
	
	/**
	 * Initializes the constant tags with the specified values. The tag for value
	 * <code>vals.getEntry(i)</code> is <code>list.get(i)</code>. Checks first to
	 * see is the tag is being used as an input variable and won't define it as a
	 * constant if it is an input item.
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
	 * Initializes the constant tags with the associated values.
	 * 
	 * @param mod Map&lt;Object,Double&gt; where Object is a tag
	 */
	public void initializeConstants(final Map<Object, Double> mod) {
		mConstants.putAll(mod);
	}

	/**
	 * Returns the constant value associated with <code>tag</code>.
	 * 
	 * @param tag
	 * @return double
	 */
	public double getConstant(final Object tag) {
		assert inputIndex(tag) == -1 : "Tag " + tag + " is not a constant.";
		return mConstants.get(tag).doubleValue();
	}

	/**
	 * Check whether <code>tag</code> is defined as a constant value.
	 * 
	 * @param tag A tag
	 * @return true if <code>tag</code> is initialized as a constant, false
	 *         otherwise.
	 */
	public boolean isConstant(final Object tag) {
		return mConstants.containsKey(tag);
	}

	/**
	 * First checks to see if <code>tag</code> is an input variable in which case it
	 * gets the associated value from <code>point</code>. Otherwise, it returns the
	 * constant value associate with <code>tag</code>.
	 *
	 * @param tag
	 * @param point
	 * @return double
	 */
	public double getValue(final Object tag, final RealVector point) {
		final int p = inputIndex(tag);
		if (p != -1) {
			assert !isConstant(tag);
			return point.getEntry(p);
		} else {
			assert isConstant(tag) : "Can't find the constant " + tag;
			return getConstant(tag);
		}
	}

	/**
	 * Returns an unmodifiable view of the map of tags and the associated constant
	 * values.
	 * 
	 * @return Map<Object, Double>
	 */
	public Map<Object, Double> getConstants() {
		return Collections.unmodifiableMap(mConstants);
	}


	@Override
	public int hashCode() {
		return Objects.hash(mInputTags, mOutputTags, mConstants);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final NamedMultivariateJacobianFunction other = (NamedMultivariateJacobianFunction) obj;
		return Objects.equals(mInputTags, other.mInputTags) && //
				Objects.equals(mOutputTags, other.mOutputTags) && //
				Objects.equals(mConstants, other.mConstants);
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE: {
			return HTML.escape("V[" + getOutputTags().size() + " values]=F(" + getInputTags().size() + " arguments, "
					+ getConstants().size() + " constants)");
		}
		case NORMAL: {
			final StringBuffer sb = new StringBuffer();
			for (final Object tag : getOutputTags())
				if (sb.length() == 0)
					sb.append("<" + tag.toString());
				else
					sb.append("," + tag.toString());
			sb.append(">=F(");
			boolean first = true;
			for (final Object tag : getOutputTags()) {
				if (!first)
					sb.append(",");
				sb.append(tag.toString());
				first = false;
			}
			first = true;
			for (final Object tag : getConstants().keySet()) {
				sb.append(first ? ";" : ",");
				sb.append(tag.toString());
				first = false;
			}
			sb.append(")");
			return HTML.escape(sb.toString());
		}
		default:
		case VERBOSE: {
			final Table res = new Table();
			final Table outs = new Table();
			for (final Object tag : getOutputTags())
				outs.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			final Table args = new Table();
			for (final Object tag : getInputTags())
				args.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			final Table consts = new Table();
			for (final Object tag : getConstants().keySet())
				consts.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			res.addRow(//
					Table.td(outs.toHTML(Mode.NORMAL)), //
					Table.td(" = F("), //
					Table.td(args.toHTML(Mode.NORMAL)), //
					Table.td(";"), //
					Table.td(consts.toHTML(Mode.NORMAL)), //
					Table.td(")")); //
			return res.toHTML(Mode.NORMAL);
		}
		}
	}
}
