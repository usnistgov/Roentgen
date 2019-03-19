package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.EstimateUncertainValues;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;

/**
 * A simple class for facilitating the evaluation of an uncertainty problem
 * using a Monte Carlo approach. The problem is defined by a
 * {@link LabeledMultivariateJacobianFunction}. A {@link UncertainValuesBase}
 * object provides the initial values and associated uncertainties. A
 * {@link MultivariateRealDistribution} object or the default
 * {@link SafeMultivariateNormalDistribution} object provides the random
 * variates.
 *
 * @author Nicholas
 */
public class MCPropagator<H, K> {

	private final LabeledMultivariateJacobianFunction<? extends H, ? extends K> mFunction;
	private final MultivariateRealDistribution mDistribution;
	private final Constraint[] mConstraints;
	private final UncertainValues<? extends H> mValues;
	private final EstimateUncertainValues<K> mOutputs;
	private final EstimateUncertainValues<H> mInputs;

	/**
	 * An interface for limiting the range of values accessible to a random
	 * variable.
	 *
	 * @author Nicholas
	 *
	 */
	public interface Constraint {

		/**
		 * Takes the value <code>val</code> and ensures that it is within the bounds of
		 * the constraint.
		 *
		 * @param val
		 * @return double val constrained
		 */
		public double limit(double val);
	}

	/**
	 * Unconstrained.
	 *
	 * @author Nicholas
	 *
	 */
	public static class None //
			implements Constraint {
		@Override
		public double limit(final double val) {
			return val;
		}
	}

	/**
	 * Constrains a value to within a [min,max] range. Values less than min are set
	 * to min and values larger than max are set to max.
	 *
	 * @author Nicholas
	 *
	 */
	public static class Range //
			implements Constraint {

		private final double mMin;
		private final double mMax;

		public Range(final double min, final double max) {
			mMin = Math.min(min, max);
			mMax = Math.max(min, max);
		}

		@Override
		public double limit(final double val) {
			if (val < mMin)
				return mMin;
			else if (val > mMax)
				return mMax;
			else
				return val;
		}
	}

	private static <J> Constraint[] buildSigma(final double sigma, final UncertainValuesBase<J> uv) {
		final int n = uv.getDimension();
		final Constraint[] res = new Constraint[n];
		for (int i = 0; i < n; ++i) {
			final double va = uv.getEntry(i);
			final double unc = uv.getUncertainty(i);
			res[i] = new Range(va - sigma * unc, va + sigma * unc);
		}
		return res;
	}

	/**
	 * Perform one evaluation of the {@link LabeledMultivariateJacobianFunction}
	 */
	private void evaluate() {
		final RealVector pt = MatrixUtils.createRealVector(mDistribution.sample());
		assert !pt.isNaN();
		for (int i = 0; i < pt.getDimension(); ++i)
			if ((mConstraints != null) && (i < mConstraints.length))
				pt.setEntry(i, mConstraints[i].limit(pt.getEntry(i)));
			else
				pt.setEntry(i, pt.getEntry(i));
		mInputs.add(pt);
		mOutputs.add(mFunction.compute(pt));
	}

	private final static class MultiEvaluate<J, L> extends RecursiveAction {

		private static final int MAX_ITERATIONS = 2000;

		private static final long serialVersionUID = 2121417483581902926L;

		private final MCPropagator<J, L> mMCP;
		private final int mNEvals;

		private MultiEvaluate(final MCPropagator<J, L> mcp, final int nEvals) {
			super();
			mNEvals = nEvals;
			mMCP = mcp;
		}

		@Override
		protected void compute() {
			if (mNEvals <= MAX_ITERATIONS) {
				for (int i = 0; i < mNEvals; ++i)
					mMCP.evaluate();
			} else
				invokeAll(new MultiEvaluate<J, L>(mMCP, mNEvals / 2), new MultiEvaluate<J, L>(mMCP, mNEvals - (mNEvals / 2)));
		}
	}

	/**
	 * Construct a {@link MCPropagator} to evaluate the specified
	 * {@link LabeledMultivariateJacobianFunction} at the specified input values
	 * randomized and constrained.
	 *
	 *
	 * @param nmvjf       {@link LabeledMultivariateJacobianFunction}
	 * @param inputs      The input values and associated covariance matrix
	 * @param constraints The constraints to the input values
	 * @param mrd         A {@link MultivariateRealDistribution} to generate the
	 *                    random variates
	 * @throws ArgumentException
	 */
	private MCPropagator( //
			UncertainValuesCalculator<H, K> uvc, //
			final Constraint[] constraints, //
			final MultivariateRealDistribution mrd//
	) throws ArgumentException {
		mFunction = uvc.getFunction();
		mConstraints = constraints != null ? constraints.clone() : null;
		mValues = uvc.getInputs();
		mOutputs = new EstimateUncertainValues<K>(mFunction.getOutputLabels());
		mInputs = new EstimateUncertainValues<H>(mFunction.getInputLabels());
		mDistribution = mrd;
	}

	/**
	 * Construct a {@link MCPropagator} to evaluate the specified
	 * {@link LabeledMultivariateJacobianFunction} at the specified input values
	 * randomized and unconstrained.
	 *
	 *
	 * @param nmvjf  {@link LabeledMultivariateJacobianFunction}
	 * @param inputs The input values and associated covariance matrix
	 * @param mrd    A {@link MultivariateRealDistribution} to generate the random
	 *               variates
	 * @throws ArgumentException
	 */
	public MCPropagator( //
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> nmvjf, //
			final UncertainValuesBase<H> uv, //
			final MultivariateRealDistribution mrd //
	) throws ArgumentException {
		this(new UncertainValuesCalculator<H,K>(nmvjf, uv), null, mrd);
	}

	/**
	 * Construct a {@link MCPropagator} to evaluate the specified
	 * {@link LabeledMultivariateJacobianFunction} at the specified input values
	 * randomized using a default {@link MultivariateNormalDistribution} and
	 * unconstrained.
	 *
	 *
	 * @param nmvjf  {@link LabeledMultivariateJacobianFunction}
	 * @param inputs The input values and associated covariance matrix in the same
	 *               order as the nmvjf input labels.
	 * @throws ArgumentException
	 */
	public MCPropagator(final UncertainValuesCalculator<H, K> uvc) throws ArgumentException {
		this(uvc, null, new SafeMultivariateNormalDistribution(uvc.getInputs()));
	}

	/**
	 * Construct a {@link MCPropagator} to evaluate the specified
	 * {@link LabeledMultivariateJacobianFunction} at the specified input values
	 * randomized and constrained.
	 *
	 *
	 * @param nmvjf  {@link LabeledMultivariateJacobianFunction}
	 * @param inputs The input values and associated covariance matrix
	 * @param sigma  Constrain the inputs to within sigma of the mean
	 * @param mrd    A {@link MultivariateRealDistribution} to generate the random
	 *               variates
	 * @throws ArgumentException
	 */
	public MCPropagator( //
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> nmvjf, //
			final UncertainValuesBase<H> uv, //
			final double sigma, //
			final MultivariateRealDistribution mrd //
	) throws ArgumentException {
		this(new UncertainValuesCalculator<H, K>(nmvjf, uv), buildSigma(sigma, uv), mrd);
	}

	/**
	 * Construct a {@link MCPropagator} to evaluate the specified
	 * {@link LabeledMultivariateJacobianFunction} at the specified input values
	 * randomized and constrained.
	 *
	 *
	 * @param nmvjf  {@link LabeledMultivariateJacobianFunction}
	 * @param inputs The input values and associated covariance matrix
	 * @param sigma  Constrain the inputs to within sigma of the mean
	 * @throws ArgumentException
	 */

	public MCPropagator(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> nmvjf, //
			final UncertainValuesBase<H> uv, //
			final double sigma //
	) throws ArgumentException {
		this(nmvjf, uv, sigma, new SafeMultivariateNormalDistribution(uv));
	}

	/**
	 * Specify a {@link Constraint} model for a specific input variable.
	 *
	 * @param label
	 * @param con
	 */
	public void setConstraint(final H label, final Constraint con) {
		final int idx = mFunction.inputIndex(label);
		mConstraints[idx] = con;
	}

	/**
	 * Evaluates the {@link LabeledMultivariateJacobianFunction} nEvals times and
	 * computes the average result and the covariance matrix which are returned as
	 * an {@link UncertainValues} object.
	 *
	 * @param nEvals
	 * @return {@link UncertainValues}
	 */
	public EstimateUncertainValues<K> compute(final int nEvals) {
		for (int eval = 0; eval < nEvals; ++eval)
			evaluate();
		return mOutputs;
	}

	/**
	 * Evaluates the {@link LabeledMultivariateJacobianFunction} nEvals times using
	 * a multi-threaded implementation and computes the average result and the
	 * covariance matrix which are returned as an {@link UncertainValues} object.
	 *
	 * @param nEvals
	 * @return {@link UncertainValues}
	 */
	public EstimateUncertainValues<K> computeMT(final int nEvals) {
		final ForkJoinPool fjp = new ForkJoinPool();
		fjp.invoke(new MultiEvaluate<H, K>(this, nEvals));
		return mOutputs;
	}

	/**
	 * Returns the {@link UncertainValues} object used to seed the
	 * {@link MultivariateRealDistribution} to provide the randomized arrays of
	 * input values.
	 *
	 * @return {@link UncertainValues}
	 */
	public UncertainValues<? extends H> getInputDistribution() {
		return mValues;
	}

	/**
	 * Computes the UncertainValues object from the randomized distribution of draws
	 * from the {@link MultivariateRealDistribution} used as the source of
	 * randomness.
	 *
	 * @return UncertainValues
	 */
	public EstimateUncertainValues<? extends H> estimatedInputDistribution() {
		return mInputs;
	}

	/**
	 * Returns the estimated output set of random values that results from
	 * evaluating the {@link LabeledMultivariateJacobianFunction} at the randomized
	 * distribution of values generated by the {@link MultivariateRealDistribution}.
	 *
	 * @return UncertainValues
	 */
	public EstimateUncertainValues<K> estimatedOutputDistribution() {
		return mOutputs;
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * output variable label.
	 *
	 * @param label An output variable label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getOutputStatistics(final K label) {
		return mOutputs.getDescriptiveStatistics(label);
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * input variable label.
	 *
	 * @param label An input variable label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getInputStatistics(final H label) {
		return mInputs.getDescriptiveStatistics(label);
	}
}