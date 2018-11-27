package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;

/**
 * <p>
 * Takes multiple independent {@link LabeledMultivariateJacobianFunction}
 * objects and combines them into a single
 * {@link LabeledMultivariateJacobianFunction}. The input
 * {@link LabeledMultivariateJacobianFunction} objects usually take a similar
 * list of arguments. None of the input
 * {@link LabeledMultivariateJacobianFunction}s may return the same output
 * label.
 * </p>
 * 
 * <p>
 * This class wraps a set of {@link LabeledMultivariateJacobianFunction} that
 * could be calculated in parallel. The output from one is not used as input to
 * another. This class is different from
 * {@link SerialLabeledMultivariateJacobianFunction} which implements a sequential
 * set of steps in which the output from earlier steps can become input to
 * subsequent steps.
 * </p>
 * 
 *
 * @author Nicholas
 */
public class LabeledMultivariateJacobianFunctionBuilder implements IToHTML {

	private final String mName;
	private final List<LabeledMultivariateJacobianFunction> mFuncs = new ArrayList<>();

	private final SimplyLazy<List<? extends Object>> mInputLabels = new SimplyLazy<List<? extends Object>>() {
		@Override
		protected List<? extends Object> initialize() {
			final Set<Object> res = new HashSet<>();
			final Set<Object> out = new HashSet<>();
			for (final LabeledMultivariateJacobianFunction func : mFuncs) {
				final List<? extends Object> inp = func.getInputLabels();
				for (final Object tag : inp)
					if (!out.contains(tag))
						res.add(tag);
				out.addAll(func.getOutputLabels());
			}
			return new ArrayList<>(res);
		}

	};

	private final SimplyLazy<List<? extends Object>> mOutputLabels = new SimplyLazy<List<? extends Object>>() {
		@Override
		protected List<? extends Object> initialize() {
			final Set<Object> res = new HashSet<>();
			for (final LabeledMultivariateJacobianFunction func : mFuncs) {
				final List<? extends Object> out = func.getOutputLabels();
				for (final Object tag : out) {
					assert !res.contains(tag) : "Duplicate output tag";
					res.add(tag);
				}
			}
			return new ArrayList<>(res);
		}
	};

	/**
	 * Ensure that output tags are not repeated.
	 *
	 * @throws ArgumentException
	 */
	private void validate() throws ArgumentException {
		final List<Object> all = new ArrayList<>();
		for (final LabeledMultivariateJacobianFunction func : mFuncs) {
			for (final Object tag : func.getOutputLabels()) {
				if (all.contains(tag))
					throw new ArgumentException("The output tag " + tag.toString() + " is repeated.");
				all.add(tag);
			}
		}
	}

	/**
	 * Constructs the builder.
	 * 
	 * @param name The name of the resulting function.
	 * @throws ArgumentException
	 */
	public LabeledMultivariateJacobianFunctionBuilder( //
			final String name //
	) throws ArgumentException {
		mName = name;
	}

	/**
	 * Constructs the builder.
	 * 
	 * @param name  The name of the resulting function.
	 * @param funcs A list of {@link LabeledMultivariateJacobianFunction} instances
	 *              that are combined.
	 * @throws ArgumentException
	 */
	public LabeledMultivariateJacobianFunctionBuilder( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> funcs //
	) throws ArgumentException {
		this(name);
		mFuncs.addAll(funcs);
		validate();
	}

	/**
	 * Add an additional {@link LabeledMultivariateJacobianFunction} to the builder.
	 * 
	 * @param func
	 * @return {@link LabeledMultivariateJacobianFunctionBuilder}
	 * @throws ArgumentException
	 */
	public LabeledMultivariateJacobianFunctionBuilder add( //
			final LabeledMultivariateJacobianFunction func //
	) throws ArgumentException {
		mFuncs.add(func);
		validate();
		mInputLabels.reset();
		mOutputLabels.reset();
		return this;
	}

	/**
	 * Static means to build the combination
	 * {@link LabeledMultivariateJacobianFunction} from a list of constituent
	 * {@link LabeledMultivariateJacobianFunction} instances.
	 * 
	 * @param name  A name for the resulting
	 *              {@link LabeledMultivariateJacobianFunction}
	 * @param funcs A list of constituent
	 *              {@link LabeledMultivariateJacobianFunction} instances.
	 * @return {@link LabeledMultivariateJacobianFunction}
	 * @throws ArgumentException
	 */
	public static LabeledMultivariateJacobianFunction join( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> funcs //
	) throws ArgumentException {
		final LabeledMultivariateJacobianFunctionBuilder j = new LabeledMultivariateJacobianFunctionBuilder(name,
				funcs);
		return j.build();
	}

	@Override
	public String toHTML(final Mode mode) {
		final Report rep = new Report("Combined NMJF - " + mName);
		rep.addHeader("Combined NMJF - " + mName);
		for (int i = 0; i < mFuncs.size(); ++i) {
			rep.addSubHeader("Function " + Integer.toString(i + 1));
			rep.add(mFuncs.get(i));
		}
		if (mode == Mode.VERBOSE) {
			rep.addHeader("Result");
			rep.add(build());
		}
		return rep.toHTML(mode);
	}

	/**
	 * Returns a list of input labelsthat is constructed from specified
	 * constituent {@link LabeledMultivariateJacobianFunction} instances.
	 * 
	 * @return List&lt;? extends Object&gt;
	 */
	public List<? extends Object> getInputLabels() {
		return mInputLabels.get();
	}

	/**
	 * Returns a list of output labels that is constructed from specified
	 * constituent {@link LabeledMultivariateJacobianFunction} instances.
	 * 
	 * @return List&lt;? extends Object&gt;
	 */
	public List<? extends Object> getOutputLabels() {
		return mOutputLabels.get();
	}

	@Override
	public String toString() {
		return mName;
	}

	private static class Implementation //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction, IToHTML {

		private final List<LabeledMultivariateJacobianFunction> mFuncs;
		private final String mName;

		private Implementation(final String name, final List<? extends Object> inpLabels,
				final List<? extends Object> outLabels, final List<LabeledMultivariateJacobianFunction> funcs) {
			super(inpLabels, outLabels);
			mFuncs = funcs;
			mName = name;
		}

		/*
		 * This function evaluates each of the individual
		 * NamedMultivariateJacobianFunction instances and then combines the outputs
		 * into a single RealVector and RealMatrix corresponding to the values and
		 * Jacobian elements.
		 *
		 * @see org.apache.commons.math3.fitting.leastsquares.
		 * MultivariateJacobianFunction#
		 * value(org.apache.commons.math3.linear.RealVector)
		 */
		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final List<? extends Object> output = getOutputLabels();
			final int oDim = output.size();
			assert oDim > 0 : "Output dimensions is zero in " + toString();
			final int iDim = getInputDimension();
			assert iDim > 0 : "Input dimensions is zero in " + toString();
			final RealVector vals = new ArrayRealVector(oDim);
			final RealMatrix cov = new NullableRealMatrix(oDim, iDim);
			for (final LabeledMultivariateJacobianFunction func : mFuncs) {
				final RealVector funcPoint = func.extractArgument(this, point);
				Map<Object, Double> consts = new HashMap<>(getConstants());
				if (func instanceof ImplicitMeasurementModel) {
					final List<? extends Object> fOut = func.getOutputLabels();
					for (Object tag : fOut)
						consts.put(tag, getValue(tag, point));
				}
				func.initializeConstants(consts);
				final Pair<RealVector, RealMatrix> v = func.value(funcPoint);
				final RealVector fVals = v.getFirst();
				final RealMatrix fJac = v.getSecond();
				final List<? extends Object> fout = func.getOutputLabels();
				final List<? extends Object> fin = func.getInputLabels();
				final int[] findx = new int[fin.size()];
				for (int c = 0; c < findx.length; ++c)
					findx[c] = inputIndex(fin.get(c));
				for (int r = 0; r < fout.size(); ++r) {
					final int idxr = outputIndex(fout.get(r));
					vals.setEntry(idxr, fVals.getEntry(r));
					for (int c = 0; c < findx.length; ++c)
						cov.setEntry(idxr, findx[c], fJac.getEntry(r, c));
				}
			}
			return Pair.create(vals, cov);
		}

		/*
		 * This function evaluates each of the individual
		 * NamedMultivariateJacobianFunction instances and then combines the outputs
		 * into a single RealVector and RealMatrix corresponding to the values and
		 * Jacobian elements.
		 */
		@Override
		public RealVector optimized(final RealVector point) {
			final List<? extends Object> output = getOutputLabels();
			final int oDim = output.size();
			assert oDim > 0 : "Output dimensions is zero in " + toString();
			final int iDim = getInputDimension();
			assert iDim > 0 : "Input dimensions is zero in " + toString();
			final RealVector vals = new ArrayRealVector(oDim);
			for (final LabeledMultivariateJacobianFunction func : mFuncs) {
				final RealVector funcPoint = func.extractArgument(this, point);
				func.initializeConstants(getConstants());
				final RealVector fVals = func.compute(funcPoint);
				final List<? extends Object> fout = func.getOutputLabels();
				final int fOutDim = fout.size();
				for (int r = 0; r < fOutDim; ++r) {
					final int idxr = outputIndex(fout.get(r));
					vals.setEntry(idxr, fVals.getEntry(r));
				}
			}
			return vals;
		}

		@Override
		public String toHTML(final Mode mode) {
			return LabeledMultivariateJacobianFunctionBuilder.toHTML(mName, mFuncs, getInputLabels(), getOutputLabels(),
					mode);
		}

	}

	private static class ParallelImplementation //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final ForkJoinPool mPool;
		private final List<LabeledMultivariateJacobianFunction> mFuncs;
		private final String mName;

		private static List<Object> buildInputs(List<LabeledMultivariateJacobianFunction> steps) {
			Set<Object> res = new HashSet<>();
			for (LabeledMultivariateJacobianFunction func : steps)
				res.addAll(func.getInputLabels());
			return new ArrayList<>(res);
		}

		private static List<Object> buildOutputs(List<LabeledMultivariateJacobianFunction> steps) {
			Set<Object> res = new HashSet<>();
			for (LabeledMultivariateJacobianFunction func : steps)
				res.addAll(func.getOutputLabels());
			return new ArrayList<>(res);
		}

		public ParallelImplementation(String name, List<LabeledMultivariateJacobianFunction> steps, ForkJoinPool fjp) {
			super(buildInputs(steps), buildOutputs(steps));
			mFuncs = new ArrayList<>(steps);
			mPool = fjp;
			mName = name;
		}

		private static final class Result {
			private final LabeledMultivariateJacobianFunction mFunc;
			private final RealVector mValues;
			private final RealMatrix mJacobian;

			Result(LabeledMultivariateJacobianFunction func, RealVector vals, RealMatrix jac) {
				mFunc = func;
				mValues = vals;
				mJacobian = jac;
			}

		}

		private static final class NMJFValue implements Callable<Result> {

			private final LabeledMultivariateJacobianFunction mFunc;
			private final RealVector mPoint;

			private NMJFValue(LabeledMultivariateJacobianFunction func, RealVector point) {
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

			private final LabeledMultivariateJacobianFunction mFunc;
			private final RealVector mPoint;

			private NMJFOptimized(LabeledMultivariateJacobianFunction func, RealVector point) {
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
			List<Callable<Result>> tasks = new ArrayList<>();
			for (LabeledMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				tasks.add(new NMJFValue(func, func.extractArgument(this, point)));
			}
			final List<Future<Result>> res = mPool.invokeAll(tasks);
			for (Future<Result> fp : res) {
				try {
					final Result tmp = fp.get();
					List<? extends Object> outLabels = tmp.mFunc.getOutputLabels();
					List<? extends Object> inLabels = tmp.mFunc.getInputLabels();
					final int[] inIdx = new int[inLabels.size()];
					for (int c = 0; c < inIdx.length; ++c)
						inIdx[c] = inputIndex(inLabels.get(c));
					for (int r = 0; r < outLabels.size(); ++r) {
						final Object outLabel = outLabels.get(r);
						final int outIdx = outputIndex(outLabel);
						rv.setEntry(outIdx, tmp.mValues.getEntry(r));
						for (int c = 0; c < inLabels.size(); ++c)
							rm.setEntry(outIdx, inIdx[c], tmp.mJacobian.getEntry(r, c));
					}
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(RealVector point) {
			RealVector rv = new ArrayRealVector(getOutputDimension());
			List<Callable<Result>> tasks = new ArrayList<>();
			for (LabeledMultivariateJacobianFunction func : mFuncs) {
				func.initializeConstants(getConstants());
				tasks.add(new NMJFOptimized(func, func.extractArgument(this, point)));
			}
			List<Future<Result>> res = mPool.invokeAll(tasks);
			for (Future<Result> fp : res) {
				try {
					final Result tmp = fp.get();
					List<? extends Object> outLabels = tmp.mFunc.getOutputLabels();
					for (int r = 0; r < outLabels.size(); ++r)
						rv.setEntry(outputIndex(outLabels.get(r)), tmp.mValues.getEntry(r));
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
			return rv;
		}

		public String toString() {
			return "Parallelized[" + mName + "]";
		}

		@Override
		public String toHTML(final Mode mode) {
			return LabeledMultivariateJacobianFunctionBuilder.toHTML(mName, mFuncs, getInputLabels(), getOutputLabels(),
					mode);
		}
	}

	protected static String toHTML(//
			String name, //
			List<LabeledMultivariateJacobianFunction> funcs, //
			List<? extends Object> inpLabels, //
			List<? extends Object> outLabels, //
			final Mode mode) {
		switch (mode) {
		case TERSE:
			return HTML.escape(name) + ": " + Integer.toString(funcs.size()) + "Consolidated Calculations";
		case NORMAL: {
			final Table tbl = new Table();
			for (int si = 0; si < funcs.size(); ++si) {
				final String html = HTML.toHTML(funcs.get(si), Mode.NORMAL);
				tbl.addRow(Table.td("Function " + si), Table.td(html));
			}
			return HTML.subHeader(HTML.escape(name)) + tbl.toHTML(mode);
		}
		default:
		case VERBOSE: {
			final Report report = new Report(name + " Consolidated Calculation Report");
			final Map<Object, Integer> valLabels = new HashMap<>();
			final Map<Object, Integer> usage = new HashMap<>();
			for (final Object inp : inpLabels) {
				valLabels.put(inp, 0);
				usage.put(inp, 0);
			}
			for (int i = 0; i < funcs.size(); ++i) {
				report.addSubHeader("Step " + Integer.toString(i + 1));
				final LabeledMultivariateJacobianFunction func = funcs.get(i);
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Input Variable"), Table.th("Source"));
					for (final Object tag : func.getInputLabels())
						tbl.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), Table.td("Input"));
					report.add(tbl);
				}
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Output Variable"), Table.th("Source"));
					for (final Object tag : func.getOutputLabels())
						tbl.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), Table.td("Output"));
					report.add(tbl);
				}
			}
			return report.toHTML(Mode.VERBOSE);
		}
		}
	}

	/**
	 * @return NamedMultivariateJacobianFunction
	 */
	public LabeledMultivariateJacobianFunction build() {
		final LabeledMultivariateJacobianFunction res = new Implementation( //
				mName, //
				mInputLabels.get(), //
				mOutputLabels.get(), //
				mFuncs);
		return res;
	}

	public LabeledMultivariateJacobianFunction buildParallel(ForkJoinPool pool) {
		final LabeledMultivariateJacobianFunction res = new ParallelImplementation( //
				mName, //
				mFuncs, //
				pool);
		return res;
	}

}
