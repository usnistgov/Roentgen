package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import gov.nist.microanalysis.roentgen.math.EstimateUncertainValues;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;

/**
 * A simple class for facilitating the evaluation of an uncertainty problem
 * using a Monte Carlo approach. The problem is defined by a
 * {@link LabeledMultivariateJacobianFunction}. A {@link UncertainValues} object
 * provides the initial values and associated uncertainties. A
 * {@link MultivariateRealDistribution} object or the default
 * {@link SafeMultivariateNormalDistribution} object provides the random
 * variates.
 *
 * @author Nicholas
 */
public class MCPropagator {

	private final LabeledMultivariateJacobianFunction mFunction;
	private final MultivariateRealDistribution mDistribution;
	private final Constraint[] mConstraints;
	private final UncertainValues mValues;
	private final EstimateUncertainValues mOutputs;
	private final EstimateUncertainValues mInputs;

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

	private static Constraint[] buildNone(final int n) {
		final Constraint[] res = new Constraint[n];
		for (int i = 0; i < n; ++i)
			res[i] = new None();
		return res;
	}

	private static Constraint[] buildSigma(final double sigma, final UncertainValues uv) {
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
		for (int i = 0; i < mConstraints.length; ++i)
			pt.setEntry(i, mConstraints[i].limit(pt.getEntry(i)));
		mInputs.add(pt);
		mOutputs.add(mFunction.compute(pt));
	}

	private final static class MultiEvaluate extends RecursiveAction {

		private static final int MAX_ITERATIONS = 2000;

		private static final long serialVersionUID = 2121417483581902926L;

		private final MCPropagator mMCP;
		private final int mNEvals;

		private MultiEvaluate(final MCPropagator mcp, final int nEvals) {
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
				invokeAll(new MultiEvaluate(mMCP, mNEvals / 2), new MultiEvaluate(mMCP, mNEvals - (mNEvals / 2)));
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
	 */
	private MCPropagator( //
			final LabeledMultivariateJacobianFunction nmvjf, //
			final UncertainValues inputs, //
			final Constraint[] constraints, //
			final MultivariateRealDistribution mrd//
	) {
		mFunction = nmvjf;
		mConstraints = constraints.clone();
		mValues = inputs;
		mOutputs = new EstimateUncertainValues(nmvjf.getOutputLabels());
		mInputs = new EstimateUncertainValues(nmvjf.getInputLabels());
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
	 */
	public MCPropagator( //
			final LabeledMultivariateJacobianFunction nmvjf, //
			final UncertainValues uv, //
			final MultivariateRealDistribution mrd //
	) {
		this(nmvjf, uv, buildNone(uv.getDimension()), mrd);
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
	 */
	public MCPropagator(final LabeledMultivariateJacobianFunction nmvjf, final UncertainValues uv) {
		this(nmvjf, uv, buildNone(uv.getDimension()),
				new SafeMultivariateNormalDistribution(uv.getValues(), uv.getCovariances()));
		assert nmvjf.getInputLabels().equals(uv.getLabels());
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
	 */
	public MCPropagator( //
			final LabeledMultivariateJacobianFunction nmvjf, //
			final UncertainValues uv, //
			final double sigma, //
			final MultivariateRealDistribution mrd //
	) {
		this(nmvjf, uv, buildSigma(sigma, uv), mrd);
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
	 */

	public MCPropagator(final LabeledMultivariateJacobianFunction nmvjf, final UncertainValues uv, final double sigma) {
		this(nmvjf, uv, sigma, new SafeMultivariateNormalDistribution(uv.getValues(), uv.getCovariances()));
	}

	/**
	 * Specify a {@link Constraint} model for a specific input variable.
	 *
	 * @param label
	 * @param con
	 */
	public void setConstraint(final Object label, final Constraint con) {
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
	public UncertainValues compute(final int nEvals) {
		for (int eval = 0; eval < nEvals; ++eval)
			evaluate();
		return mOutputs.estimateDistribution();
	}

	/**
	 * Evaluates the {@link LabeledMultivariateJacobianFunction} nEvals times using
	 * a multi-threaded implementation and computes the average result and the
	 * covariance matrix which are returned as an {@link UncertainValues} object.
	 *
	 * @param nEvals
	 * @return {@link UncertainValues}
	 */
	public UncertainValues computeMT(final int nEvals) {
		final ForkJoinPool fjp = new ForkJoinPool();
		fjp.invoke(new MultiEvaluate(this, nEvals));
		return mOutputs.estimateDistribution();
	}

	/**
	 * Returns the {@link UncertainValues} object used to seed the
	 * {@link MultivariateRealDistribution} to provide the randomized arrays of
	 * input values.
	 *
	 * @return {@link UncertainValues}
	 */
	public UncertainValues getInputDistribution() {
		return mValues;
	}

	/**
	 * Computes the UncertainValues object from the randomized distribution of draws
	 * from the {@link MultivariateRealDistribution} used as the source of
	 * randomness.
	 *
	 * @return UncertainValues
	 */
	public UncertainValues estimatedInputDistribution() {
		return mInputs.estimateDistribution();
	}

	/**
	 * Returns the estimated output set of random values that results from
	 * evaluating the {@link LabeledMultivariateJacobianFunction} at the randomized
	 * distribution of values generated by the {@link MultivariateRealDistribution}.
	 *
	 * @return UncertainValues
	 */
	public UncertainValues estimatedOutputDistribution() {
		return mOutputs.estimateDistribution();
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * output variable label.
	 *
	 * @param label An output variable label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getOutputStatistics(final Object label) {
		return mOutputs.getDescriptiveStatistics(label);
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * input variable label.
	 *
	 * @param label An input variable label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getInputStatistics(final Object label) {
		return mInputs.getDescriptiveStatistics(label);
	}
}