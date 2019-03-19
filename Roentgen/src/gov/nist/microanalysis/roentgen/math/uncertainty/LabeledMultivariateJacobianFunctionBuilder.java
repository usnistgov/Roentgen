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
import gov.nist.microanalysis.roentgen.utility.FastIndex;

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
 * {@link CompositeLabeledMultivariateJacobianFunction} which implements a
 * sequential set of steps in which the output from earlier steps can become
 * input to subsequent steps.
 * </p>
 *
 *
 * @author Nicholas
 */
public class LabeledMultivariateJacobianFunctionBuilder<H, K> implements IToHTML {

	private final String mName;
	private final List<LabeledMultivariateJacobianFunction<? extends H, ? extends K>> mFuncs = new ArrayList<>();

	private final SimplyLazy<List<H>> mInputLabels = new SimplyLazy<List<H>>() {
		@Override
		protected FastIndex<H> initialize() {
			final FastIndex<H> res = new FastIndex<>();
			for (final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func : mFuncs)
				res.addMissing(func.getInputLabels());
			return res;
		}
	};

	private final SimplyLazy<List<K>> mOutputLabels = new SimplyLazy<List<K>>() {
		@Override
		protected FastIndex<K> initialize() {
			final FastIndex<K> res = new FastIndex<>();
			for (final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func : mFuncs) {
				final List<? extends K> out = func.getOutputLabels();
				for (final K tag : out) {
					assert (func.inputIndex(tag) != -1) || (!res.contains(tag)) : "Duplicate output tag";
					res.addIfMissing(tag);
				}
			}
			return res;
		}
	};

	/**
	 * Ensure that output tags are not repeated.
	 *
	 * @throws ArgumentException
	 */
	private void validate() throws ArgumentException {
		final List<K> all = new ArrayList<>();
		for (final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func : mFuncs) {
			for (final K tag : func.getOutputLabels()) {
				if (all.contains(tag) && (func.inputIndex(tag) == -1))
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
			final List<LabeledMultivariateJacobianFunction<? extends H, ? extends K>> funcs //
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
	public LabeledMultivariateJacobianFunctionBuilder<? extends H, ? extends K> add( //
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func //
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
	public static <H, K> LabeledMultivariateJacobianFunction<H, K> join( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction<? extends H, ? extends K>> funcs //
	) throws ArgumentException {
		final LabeledMultivariateJacobianFunctionBuilder<H, K> j = new LabeledMultivariateJacobianFunctionBuilder<H, K>(
				name, funcs);
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
	 * Returns a list of input labelsthat is constructed from specified constituent
	 * {@link LabeledMultivariateJacobianFunction} instances.
	 *
	 * @return List&lt;H&gt;
	 */
	public List<? extends H> getInputLabels() {
		return mInputLabels.get();
	}

	/**
	 * Returns a list of output labels that is constructed from specified
	 * constituent {@link LabeledMultivariateJacobianFunction} instances.
	 *
	 * @return List&lt;H&gt;
	 */
	public List<K> getOutputLabels() {
		return mOutputLabels.get();
	}

	@Override
	public String toString() {
		return mName;
	}

	private static class Implementation<J, L> //
			extends LabeledMultivariateJacobianFunction<J, L> //
			implements ILabeledMultivariateFunction<J, L>, IToHTML {

		private final List<LabeledMultivariateJacobianFunction<? extends J, ? extends L>> mFuncs;
		private final String mName;

		private Implementation( //
				final String name, //
				final List<J> inpLabels, //
				final List<L> outLabels, //
				final List<LabeledMultivariateJacobianFunction<? extends J, ? extends L>> funcs //
		) {
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
			final int iDim = getInputDimension(), oDim = getOutputDimension();
			final RealVector vals = new ArrayRealVector(oDim);
			final RealMatrix cov = NullableRealMatrix.build(oDim, iDim);
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func : mFuncs) {
				final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
				for (int i = 0; i < func.getInputDimension(); ++i) {
					final J fLbl = func.getInputLabel(i);
					funcPoint.setEntry(i, point.getEntry(inputIndex(fLbl)));
				}
				if (func instanceof ImplicitMeasurementModel) {
					@SuppressWarnings("unchecked")
					final ImplicitMeasurementModel<L> imm = (ImplicitMeasurementModel<L>) func;
					final Map<L, Double> consts = new HashMap<>();
					for (final L tag : func.getOutputLabels())
						consts.put(tag, point.getEntry(inputIndex(tag)));
					imm.initializeConstants(consts);
				}
				func.dumpArguments(funcPoint, this);
				final Pair<RealVector, RealMatrix> v = func.value(funcPoint);
				final RealVector fVals = v.getFirst();
				final RealMatrix fJac = v.getSecond();
				final List<? extends L> fout = func.getOutputLabels();
				final List<? extends J> fin = func.getInputLabels();
				final int[] idxc = new int[fin.size()];
				for (int c = 0; c < idxc.length; ++c)
					idxc[c] = inputIndex(fin.get(c));
				for (int r = 0; r < fout.size(); ++r) {
					final int idxr = outputIndex(fout.get(r));
					vals.setEntry(idxr, fVals.getEntry(r));
					for (int c = 0; c < idxc.length; ++c)
						cov.setEntry(idxr, idxc[c], fJac.getEntry(r, c));
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
			final int oDim = getOutputDimension();
			final RealVector vals = new ArrayRealVector(oDim);
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func : mFuncs) {
				final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
				for (int i = 0; i < func.getInputDimension(); ++i) {
					final J fLbl = func.getInputLabel(i);
					funcPoint.setEntry(i, point.getEntry(inputIndex(fLbl)));
				}
				@SuppressWarnings("unchecked")
				ILabeledMultivariateFunction<L, L> altFunc = func instanceof ILabeledMultivariateFunction ? //
						(ILabeledMultivariateFunction<L, L>) func : null;
				if (func instanceof ImplicitMeasurementModel) {
					@SuppressWarnings("unchecked")
					final ImplicitMeasurementModel<L> imm = (ImplicitMeasurementModel<L>) func;
					if (imm.hasAlternativeModel()) {
						altFunc = imm.getAlternativeModel();
					} else {
						// Add in output values as constants...
						final Map<L, Double> consts = new HashMap<>();
						for (final L label : func.getOutputLabels())
							consts.put(label, point.getEntry(inputIndex(label)));
						imm.initializeConstants(consts);
					}
				}
				func.dumpArguments(funcPoint, this);
				final RealVector fVals = altFunc != null ? //
						altFunc.optimized(funcPoint) : func.evaluate(funcPoint).getFirst();
				final List<? extends L> fout = func.getOutputLabels();
				for (int r = 0; r < fout.size(); ++r)
					vals.setEntry(outputIndex(fout.get(r)), fVals.getEntry(r));
			}
			return vals;
		}

		@Override
		public String toHTML(final Mode mode) {
			return LabeledMultivariateJacobianFunctionBuilder.toHTML(mName, mFuncs, getInputLabels(), getOutputLabels(),
					mode);
		}

		@Override
		public String toString() {
			return mName;
		}

	}

	private static class ParallelImplementation<J, L> //
			extends LabeledMultivariateJacobianFunction<J, L> //
			implements ILabeledMultivariateFunction<J, L> {

		private final ForkJoinPool mPool;
		private final List<LabeledMultivariateJacobianFunction<? extends J, ? extends L>> mFuncs;
		private final String mName;

		private static <J, L> List<J> buildInputs(
				final List<LabeledMultivariateJacobianFunction<? extends J, ? extends L>> steps) {
			final Set<J> res = new HashSet<>();
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func : steps)
				res.addAll(func.getInputLabels());
			return new ArrayList<>(res);
		}

		private static <J, K> List<K> buildOutputs(
				final List<LabeledMultivariateJacobianFunction<? extends J, ? extends K>> steps) {
			final Set<K> res = new HashSet<>();
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends K> func : steps)
				res.addAll(func.getOutputLabels());
			return new ArrayList<>(res);
		}

		public ParallelImplementation(final String name,
				final List<LabeledMultivariateJacobianFunction<? extends J, ? extends L>> steps,
				final ForkJoinPool fjp) {
			super(buildInputs(steps), buildOutputs(steps));
			mFuncs = new ArrayList<>(steps);
			mPool = fjp;
			mName = name;
		}

		private static final class Result<J, L> {
			private final LabeledMultivariateJacobianFunction<? extends J, ? extends L> mFunc;
			private final RealVector mValues;
			private final RealMatrix mJacobian;

			Result(final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func, final RealVector vals,
					final RealMatrix jac) {
				mFunc = func;
				mValues = vals;
				mJacobian = jac;
			}

		}

		private static final class NMJFValue<J, L> implements Callable<Result<J, L>> {

			private final LabeledMultivariateJacobianFunction<? extends J, ? extends L> mFunc;
			private final RealVector mPoint;

			private NMJFValue(final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func,
					final RealVector point) {
				mFunc = func;
				mPoint = point;
			}

			@Override
			public Result<J, L> call() throws Exception {
				final Pair<RealVector, RealMatrix> tmp = mFunc.evaluate(mPoint);
				return new Result<J, L>(mFunc, tmp.getFirst(), tmp.getSecond());
			}
		}

		private static final class NMJFOptimized<J, L> implements Callable<Result<J, L>> {

			private final LabeledMultivariateJacobianFunction<? extends J, ? extends L> mFunc;
			private final RealVector mPoint;

			private NMJFOptimized(final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func,
					final RealVector point) {
				mFunc = func;
				mPoint = point;
			}

			@Override
			public Result<J, L> call() throws Exception {
				final RealVector tmp = mFunc.compute(mPoint);
				return new Result<J, L>(mFunc, tmp, null);
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
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());
			final List<Callable<Result<J, L>>> tasks = new ArrayList<>();
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func : mFuncs) {
				final RealVector fPoint = func.extractArgument(this, point);
				func.dumpArguments(fPoint, this);
				tasks.add(new NMJFValue<J, L>(func, fPoint));
			}
			final List<Future<Result<J, L>>> res = mPool.invokeAll(tasks);
			for (final Future<Result<J, L>> fp : res) {
				try {
					final Result<J, L> tmp = fp.get();
					final List<? extends L> outLabels = tmp.mFunc.getOutputLabels();
					final List<? extends J> inLabels = tmp.mFunc.getInputLabels();
					final int[] inIdx = new int[inLabels.size()];
					for (int c = 0; c < inIdx.length; ++c)
						inIdx[c] = inputIndex(inLabels.get(c));
					for (int r = 0; r < outLabels.size(); ++r) {
						final L outLabel = outLabels.get(r);
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
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final List<Callable<Result<J, L>>> tasks = new ArrayList<>();
			for (final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func : mFuncs) {
				final RealVector fPoint = func.extractArgument(this, point);
				func.dumpArguments(fPoint, this);
				tasks.add(new NMJFOptimized<J, L>(func, fPoint));
			}
			final List<Future<Result<J, L>>> res = mPool.invokeAll(tasks);
			for (final Future<Result<J, L>> fp : res) {
				try {
					final Result<J, L> tmp = fp.get();
					final List<? extends L> outLabels = tmp.mFunc.getOutputLabels();
					for (int r = 0; r < outLabels.size(); ++r)
						rv.setEntry(outputIndex(outLabels.get(r)), tmp.mValues.getEntry(r));
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
			return rv;
		}

		@Override
		public String toString() {
			return "Parallelized[" + mName + "]";
		}

		@Override
		public String toHTML(final Mode mode) {
			return LabeledMultivariateJacobianFunctionBuilder.toHTML(mName, mFuncs, getInputLabels(), getOutputLabels(),
					mode);
		}
	}

	protected static <H, K> String toHTML(//
			final String name, //
			final List<LabeledMultivariateJacobianFunction<? extends H, ? extends K>> funcs, //
			final List<H> inpLabels, //
			final List<K> outLabels, //
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
			final Map<H, Integer> valLabels = new HashMap<>();
			final Map<H, Integer> usage = new HashMap<>();
			for (final H inp : inpLabels) {
				valLabels.put(inp, 0);
				usage.put(inp, 0);
			}
			for (int i = 0; i < funcs.size(); ++i) {
				report.addSubHeader("Step " + Integer.toString(i + 1));
				final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func = funcs.get(i);
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Input Variable"), Table.th("Source"));
					for (final H tag : func.getInputLabels())
						tbl.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), Table.td("Input"));
					report.add(tbl);
				}
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Output Variable"), Table.th("Source"));
					for (final K tag : func.getOutputLabels())
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
	public LabeledMultivariateJacobianFunction<H, K> build() {
		final LabeledMultivariateJacobianFunction<H, K> res = new Implementation<H, K>( //
				mName, //
				mInputLabels.get(), //
				mOutputLabels.get(), //
				mFuncs);
		return res;
	}

	public LabeledMultivariateJacobianFunction<H, K> buildParallel(final ForkJoinPool pool) {
		final LabeledMultivariateJacobianFunction<H, K> res = new ParallelImplementation<H, K>( //
				mName, //
				mFuncs, //
				pool);
		return res;
	}

}
