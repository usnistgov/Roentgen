package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
 * Takes multiple {@link NamedMultivariateJacobianFunction} objects and combines
 * them into a single {@link NamedMultivariateJacobianFunction}. The input
 * {@link NamedMultivariateJacobianFunction} objects usually take a similar list
 * of arguments. None of the input {@link NamedMultivariateJacobianFunction} may
 * return the same output tag.
 *
 * @author Nicholas
 */
public class NamedMultivariateJacobianFunctionBuilder implements IToHTML {

	private final String mName;
	private final List<NamedMultivariateJacobianFunction> mFuncs = new ArrayList<>();

	private final SimplyLazy<List<? extends Object>> mInputTags = new SimplyLazy<List<? extends Object>>() {
		@Override
		protected List<? extends Object> initialize() {
			final List<Object> res = new ArrayList<>();
			for (final NamedMultivariateJacobianFunction func : mFuncs) {
				final List<? extends Object> inp = func.getInputTags();
				for (final Object tag : inp)
					if (!res.contains(tag))
						res.add(tag);
			}
			return res;
		}

	};

	private final SimplyLazy<List<? extends Object>> mOutputTags = new SimplyLazy<List<? extends Object>>() {
		@Override
		protected List<? extends Object> initialize() {
			final List<Object> res = new ArrayList<>();
			for (final NamedMultivariateJacobianFunction func : mFuncs) {
				final List<? extends Object> out = func.getOutputTags();
				for (final Object tag : out) {
					assert !res.contains(tag) : "Duplicate output tag";
					res.add(tag);
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
		final List<Object> all = new ArrayList<>();
		for (final NamedMultivariateJacobianFunction func : mFuncs) {
			for (final Object tag : func.getOutputTags()) {
				if (all.contains(tag))
					throw new ArgumentException("The output tag " + tag.toString() + " is repeated.");
				all.add(tag);
			}
		}
	}

	public NamedMultivariateJacobianFunctionBuilder(final String name) throws ArgumentException {
		mName = name;
	}

	public NamedMultivariateJacobianFunctionBuilder(final String name,
			final List<NamedMultivariateJacobianFunction> funcs) throws ArgumentException {
		this(name);
		mFuncs.addAll(funcs);
		validate();
	}

	public NamedMultivariateJacobianFunctionBuilder add(final NamedMultivariateJacobianFunction func)
			throws ArgumentException {
		mFuncs.add(func);
		validate();
		mInputTags.reset();
		return this;
	}

	public static NamedMultivariateJacobianFunction join(final String name,
			final List<NamedMultivariateJacobianFunction> funcs) throws ArgumentException {
		final NamedMultivariateJacobianFunctionBuilder j = new NamedMultivariateJacobianFunctionBuilder(name, funcs);
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

	public List<? extends Object> getInputTags() {
		return mInputTags.get();
	}

	public List<? extends Object> getOutputTags() {
		return mOutputTags.get();
	}

	@Override
	public String toString() {
		return "Combined " + mName;
	}

	private static class Implementation //
			extends NamedMultivariateJacobianFunction //
			implements INamedMultivariateFunction, IToHTML {

		private final List<NamedMultivariateJacobianFunction> mFuncs;
		private final String mName;

		private Implementation(String name, List<? extends Object> inpTags, List<? extends Object> outTags,
				List<NamedMultivariateJacobianFunction> funcs) {
			super(inpTags, outTags);
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
			final List<? extends Object> output = getOutputTags();
			final int oDim = output.size();
			assert oDim > 0 : "Output dimensions is zero in " + toString();
			final int iDim = getInputDimension();
			assert iDim > 0 : "Input dimensions is zero in " + toString();
			final RealVector vals = new ArrayRealVector(oDim);
			final RealMatrix cov = new NullableRealMatrix(oDim, iDim);
			for (final NamedMultivariateJacobianFunction func : mFuncs) {
				final RealVector funcPoint = func.extractArgument(this, point);
				func.initializeConstants(getConstants());
				final Pair<RealVector, RealMatrix> v = func.value(funcPoint);
				final RealVector fVals = v.getFirst();
				final RealMatrix fJac = v.getSecond();
				final List<? extends Object> fout = func.getOutputTags();
				final List<? extends Object> fin = func.getInputTags();
				final int fOutDim = fout.size();
				final int fInDim = fin.size();
				for (int r = 0; r < fOutDim; ++r) {
					final int idxr = outputIndex(fout.get(r));
					vals.setEntry(idxr, fVals.getEntry(r));
					for (int c = 0; c < fInDim; ++c) {
						final int idxc = inputIndex(fin.get(c));
						cov.setEntry(idxr, idxc, fJac.getEntry(r, c));
					}
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
			final List<? extends Object> output = getOutputTags();
			final int oDim = output.size();
			assert oDim > 0 : "Output dimensions is zero in " + toString();
			final int iDim = getInputDimension();
			assert iDim > 0 : "Input dimensions is zero in " + toString();
			final RealVector vals = new ArrayRealVector(oDim);
			for (final NamedMultivariateJacobianFunction func : mFuncs) {
				final RealVector funcPoint = func.extractArgument(this, point);
				func.initializeConstants(getConstants());
				final RealVector fVals = func.compute(funcPoint);
				final List<? extends Object> fout = func.getOutputTags();
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
			switch (mode) {
			case TERSE:
				return HTML.escape(mName) + ": A " + Integer.toString(mFuncs.size()) + "-Step Calculation";
			case NORMAL: {
				final Table tbl = new Table();
				for (int si = 0; si < mFuncs.size(); ++si) {
					String html = HTML.toHTML(mFuncs.get(si), Mode.NORMAL);
					tbl.addRow(Table.td("Step " + si), Table.td(html));
				}
				return HTML.subHeader(HTML.escape(mName)) + tbl.toHTML(mode);
			}
			default:
			case VERBOSE: {
				final Report report = new Report(mName + " Multi-Step Calculation Report");
				final Map<Object, Integer> valTags = new HashMap<>();
				final Map<Object, Integer> usage = new HashMap<>();
				for (final Object inp : getInputTags()) {
					valTags.put(inp, 0);
					usage.put(inp, 0);
				}
				for (int i = 0; i < mFuncs.size(); ++i) {
					report.addSubHeader("Step " + Integer.toString(i + 1));
					final NamedMultivariateJacobianFunction func = mFuncs.get(i);
					final Table tbl = new Table();
					tbl.addRow(Table.th("Input Variable"), Table.th("Source"));
					for (final Object intag : func.getInputTags()) {
						final Integer ss = valTags.get(intag);
						String item = null;
						if (ss == null)
							item = HTML.error("Not defined.");
						else {
							if (ss.intValue() == 0)
								item = "Input variable";
							else
								item = "Calculated in Step " + Integer.toString(ss.intValue());
							usage.put(intag, usage.get(intag) + 1);
						}
						tbl.addRow(Table.td(HTML.toHTML(intag, Mode.TERSE)), Table.td(item));
					}
					report.add(tbl);
					for (final Object outtag : func.getOutputTags()) {
						if (!valTags.containsKey(outtag)) {
							valTags.put(outtag, i + 1);
							usage.put(outtag, 0);
						} else {
							report.addHTML(HTML.error("Output error in Step " + Integer.toString(i + 1)
									+ ": Attempting to redefine " + HTML.toHTML(outtag, Mode.TERSE) + "."));
						}
					}
				}
				report.addSubHeader("Input Utilization Table");
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Input"), Table.th("Use Count"));
					for (final Object tag : getInputTags())
						tbl.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)),
								usage.containsKey(tag) ? Table.td(usage.get(tag)) : Table.td("Missing"));
					report.add(tbl);
				}
				report.addSubHeader("Output Utilization Table");
				{
					final Table tbl = new Table();
					tbl.addRow(Table.th("Output"), Table.th("Defined"), Table.th("Use Count"));
					for (final Object out : getOutputTags())

						tbl.addRow(Table.td(HTML.toHTML(out, Mode.TERSE)), //
								Table.td(valTags.getOrDefault(out, 0).toString()), //
								Table.td(usage.getOrDefault(out, 0).toString()));
					report.add(tbl);
				}
				return report.toHTML(Mode.VERBOSE);
			}
			}

		}

	}

	/**
	 * @return
	 */
	public NamedMultivariateJacobianFunction build() {
		final NamedMultivariateJacobianFunction res = new Implementation( //
				mName, //
				mInputTags.get(), //
				mOutputTags.get(), mFuncs);
		return res;
	}

}
