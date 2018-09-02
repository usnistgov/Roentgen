/**
 * 
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

/**
 * @author nicho
 *
 */
public class ParallelNamedMultivariateJacobianFunction extends NamedMultivariateJacobianFunction
		implements INamedMultivariateFunction {

	private final ForkJoinPool mPool;
	final List<NamedMultivariateJacobianFunction> mFuncs;

	private static List<Object> buildInputs(List<NamedMultivariateJacobianFunction> steps) {
		Set<Object> res = new HashSet<>();
		for (NamedMultivariateJacobianFunction func : steps)
			res.addAll(func.getInputTags());
		return new ArrayList<>(res);
	}

	private static List<Object> buildOutputs(List<NamedMultivariateJacobianFunction> steps) {
		Set<Object> res = new HashSet<>();
		for (NamedMultivariateJacobianFunction func : steps)
			res.addAll(func.getOutputTags());
		return new ArrayList<>(res);
	}

	public ParallelNamedMultivariateJacobianFunction(List<NamedMultivariateJacobianFunction> steps, ForkJoinPool fjp) {
		super(buildInputs(steps), buildOutputs(steps));
		mFuncs = new ArrayList<>(steps);
		mPool = fjp;
	}

	private static final class Result {
		private final NamedMultivariateJacobianFunction mFunc;
		private final RealVector mValues;
		private final RealMatrix mJacobian;

		Result(NamedMultivariateJacobianFunction func, RealVector vals, RealMatrix jac) {
			mFunc = func;
			mValues = vals;
			mJacobian = jac;
		}

	}

	private static final class NMJFValue implements Callable<Result> {

		private final NamedMultivariateJacobianFunction mFunc;
		private final RealVector mPoint;

		private NMJFValue(NamedMultivariateJacobianFunction func, RealVector point) {
			mFunc = func;
			mPoint = point;
		}

		@Override
		public Result call() throws Exception {
			Pair<RealVector, RealMatrix> tmp = mFunc.evaluate(mPoint);
			return new Result(mFunc, tmp.getFirst(), tmp.getSecond());
		}
	}

	private static final class NMJFOptimized implements Callable<Result> {

		private final NamedMultivariateJacobianFunction mFunc;
		private final RealVector mPoint;

		private NMJFOptimized(NamedMultivariateJacobianFunction func, RealVector point) {
			mFunc = func;
			mPoint = point;
		}

		@Override
		public Result call() throws Exception {
			final RealVector tmp = mFunc.compute(mPoint);
			return new Result(mFunc, tmp, null);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		RealVector rv = new ArrayRealVector(getOutputDimension());
		RealMatrix rm = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());
		if (mPool!=null) {
			List<Callable<Result>> tasks = new ArrayList<>();
			for (NamedMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				tasks.add(new NMJFValue(func, func.extractArgument(this, point)));
			}
			final List<Future<Result>> res = mPool.invokeAll(tasks);
			for (Future<Result> fp : res) {
				try {
					final Result tmp = fp.get();
					List<? extends Object> outTags = tmp.mFunc.getOutputTags();
					List<? extends Object> inTags = tmp.mFunc.getInputTags();
					for (int r = 0; r < outTags.size(); ++r) {
						final Object outTag = outTags.get(r);
						final int outIdx = outputIndex(outTag);
						rv.setEntry(outIdx, tmp.mValues.getEntry(r));
						for (int c = 0; c < inTags.size(); ++c)
							rm.setEntry(outIdx, inputIndex(inTags.get(c)), tmp.mJacobian.getEntry(r, c));
					}
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
		} else {
			for (NamedMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				Pair<RealVector, RealMatrix> pr = func.evaluate(func.extractArgument(this, point));
				RealVector prv = pr.getFirst();
				RealMatrix prm = pr.getSecond();
				List<? extends Object> outTags = func.getOutputTags();
				List<? extends Object> inTags = func.getInputTags();
				for (int r = 0; r < outTags.size(); ++r) {
					final Object outTag = outTags.get(r);
					final int outIdx = outputIndex(outTag);
					rv.setEntry(outIdx, prv.getEntry(r));
					for (int c = 0; c < inTags.size(); ++c)
						rm.setEntry(outIdx, inputIndex(inTags.get(c)), prm.getEntry(r, c));
				}
			}
		}
		return Pair.create(rv, rm);
	}

	@Override
	public RealVector optimized(RealVector point) {
		RealVector rv = new ArrayRealVector(getOutputDimension());
		if (mPool!=null) {
			List<Callable<Result>> tasks = new ArrayList<>();
			for (NamedMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				tasks.add(new NMJFOptimized(func, func.extractArgument(this, point)));
			}
			List<Future<Result>> res = mPool.invokeAll(tasks);
			for (Future<Result> fp : res) {
				try {
					final Result tmp = fp.get();
					List<? extends Object> outTags = tmp.mFunc.getOutputTags();
					for (int r = 0; r < outTags.size(); ++r)
						rv.setEntry(outputIndex(outTags.get(r)), tmp.mValues.getEntry(r));
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
		} else {
			for (NamedMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				RealVector prv = func.compute(func.extractArgument(this, point));
				List<? extends Object> outTags = func.getOutputTags();
				for (int r = 0; r < outTags.size(); ++r) {
					final int outIdx = outputIndex(outTags.get(r));
					rv.setEntry(outIdx, prv.getEntry(r));
				}
			}
		}
		return rv;
	}
	
	public String toString() {
		StringBuffer sb=new StringBuffer();
		for(NamedMultivariateJacobianFunction func: mFuncs) {
			if(sb.length()>0)
				sb.append(",");
			sb.append(func.toString());
		}
		return "Parallelized["+sb.toString()+"]";
	}

}
