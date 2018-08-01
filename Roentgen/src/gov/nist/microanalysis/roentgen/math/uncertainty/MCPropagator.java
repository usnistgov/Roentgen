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
	private final EstimateUncertainValues mEvals;
	private final EstimateUncertainValues mPts;

	public static class None implements Constraint {
		@Override
		public double limit(final double val) {
			return val;
		}
	}

	public static class Range implements Constraint {

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

	private MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv,
			Constraint[] constraints, MultivariateRealDistribution mrd) {
		mFunction = nmvjf;
		mConstraints = constraints.clone();
		mValues = uv;
		mEvals = new EstimateUncertainValues(nmvjf.getOutputTags());
		mPts = new EstimateUncertainValues(nmvjf.getInputTags());
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
		return mEvals.estimateDistribution();
	}

	private void evaluate() {
		final RealVector pt = MatrixUtils.createRealVector(mDistribution.sample());
		for (int i = 0; i < mConstraints.length; ++i)
			pt.setEntry(i, mConstraints[i].limit(pt.getEntry(i)));
		mPts.add(pt);
		mEvals.add(mFunction.evaluate(pt).getFirst());
	}

	private final static class MultiEvaluate extends RecursiveAction {

		private static final int MAX_ITERATIONS = 2000;

		private static final long serialVersionUID = 2121417483581902926L;

		private final MCPropagator mMCP;
		private final int mNEvals;

		private MultiEvaluate(MCPropagator mcp, int nEvals) {
			assert nEvals >= MAX_ITERATIONS / 2;
			mNEvals = nEvals;
			mMCP = mcp;
		}

		@Override
		protected void compute() {
			if (mNEvals <= MAX_ITERATIONS) {
				for (int i = 0; i < mNEvals; ++i)
					mMCP.evaluate();
				return;
			}
			invokeAll(new MultiEvaluate(mMCP, mNEvals / 2), new MultiEvaluate(mMCP, mNEvals - (mNEvals / 2)));
		}

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
		return mEvals.estimateDistribution();
	}

	public UncertainValues getInputDistribution() {
		return mValues;
	}

	public UncertainValues estimatedInputDistribution() {
		return mPts.estimateDistribution();
	}

	public UncertainValues estimatedOutputDistribution() {
		return mEvals.estimateDistribution();
	}

	public DescriptiveStatistics getOutputStatistics(final Object tag) {
		return mEvals.getDescriptiveStatistics(tag);
	}

	public DescriptiveStatistics getInputStatistics(final Object tag) {
		return mPts.getDescriptiveStatistics(tag);
	}
}