package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import gov.nist.microanalysis.roentgen.math.EstimateUncertainValues;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;

/**
 * A simple class for facilitating the evaluation of an uncertainty problem
 * using a Monte Carlo approach. The problem is defined by a
 * {@link NamedMultivariateJacobianFunction}. A {@link UncertainValues} object
 * provides the initial values and associated uncertainties. A
 * {@link MultivariateRealDistribution} object or the default
 * {@link SafeMultivariateNormalDistribution} object provides the random
 * variates.
 * 
 * @author Nicholas
 */
public class MCPropagator {

	public interface Constraint {
		public double limit(double val);
	}

	private final NamedMultivariateJacobianFunction mFunction;
	private final MultivariateRealDistribution mDistribution;
	private final Constraint[] mConstraints;
	private final UncertainValues mValues;
	private final EstimateUncertainValues mOutputs;
	private final EstimateUncertainValues mInputs;

	public static class None //
			implements Constraint {
		@Override
		public double limit(final double val) {
			return val;
		}
	}

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

	private MCPropagator( //
			final NamedMultivariateJacobianFunction nmvjf, //
			final UncertainValues uv, //
			Constraint[] constraints, //
			MultivariateRealDistribution mrd//
	) {
		mFunction = nmvjf;
		mConstraints = constraints.clone();
		mValues = uv;
		mOutputs = new EstimateUncertainValues(nmvjf.getOutputTags());
		mInputs = new EstimateUncertainValues(nmvjf.getInputTags());
		mDistribution = mrd;
	}

	private static Constraint[] buildNone(int n) {
		Constraint[] res = new Constraint[n];
		for (int i = 0; i < n; ++i)
			res[i] = new None();
		return res;
	}

	private static Constraint[] buildSigma(double sigma, UncertainValues uv) {
		final int n = uv.getDimension();
		Constraint[] res = new Constraint[n];
		for (int i = 0; i < n; ++i) {
			final double va = uv.getEntry(i);
			final double unc = uv.getUncertainty(i);
			res[i] = new Range(va - sigma * unc, va + sigma * unc);
		}
		return res;
	}

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

		private MultiEvaluate(MCPropagator mcp, int nEvals) {
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

	public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv,
			MultivariateRealDistribution mrd) {
		this(nmvjf, uv, buildNone(uv.getDimension()), mrd);
	}

	public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv) {
		this(nmvjf, uv, buildNone(uv.getDimension()),
				new SafeMultivariateNormalDistribution(uv.getValues(), uv.getCovariances()));
	}

	public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv, final double sigma,
			MultivariateRealDistribution mrd) {
		this(nmvjf, uv, buildSigma(sigma, uv), mrd);
	}

	public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv, final double sigma) {
		this(nmvjf, uv, sigma, new SafeMultivariateNormalDistribution(uv.getValues(), uv.getCovariances()));
	}

	public void setConstraint(final Object tag, final Constraint con) {
		final int idx = mFunction.inputIndex(tag);
		mConstraints[idx] = con;
	}

	public UncertainValues compute(final int nEvals) {
		for (int eval = 0; eval < nEvals; ++eval)
			evaluate();
		return mOutputs.estimateDistribution();
	}

	/**
	 * Evaluates using a multi-threaded algorithm the Monte Carlo model.
	 * 
	 * @param nEvals Number of evaluations to perform
	 * @return UncertainValues
	 */
	public UncertainValues computeMT(final int nEvals) {
		ForkJoinPool fjp = new ForkJoinPool();
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
	 * evaluating the {@link NamedMultivariateJacobianFunction} at the randomized
	 * distribution of values generated by the {@link MultivariateRealDistribution}.
	 * 
	 * @return UncertainValues
	 */
	public UncertainValues estimatedOutputDistribution() {
		return mOutputs.estimateDistribution();
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * output variable tag.
	 * 
	 * @param tag An output variable tag
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getOutputStatistics(final Object tag) {
		return mOutputs.getDescriptiveStatistics(tag);
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object associated with the specified
	 * input variable tag.
	 * 
	 * @param tag An input variable tag
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getInputStatistics(final Object tag) {
		return mInputs.getDescriptiveStatistics(tag);
	}
}