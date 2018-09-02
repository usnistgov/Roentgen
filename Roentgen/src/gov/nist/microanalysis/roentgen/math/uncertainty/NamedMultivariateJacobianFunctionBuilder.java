package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
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

	private static class Implementation extends NamedMultivariateJacobianFunction implements INamedMultivariateFunction {

		private final List<NamedMultivariateJacobianFunction> mFuncs;

		private Implementation(List<? extends Object> inpTags, List<? extends Object> outTags,
				List<NamedMultivariateJacobianFunction> funcs) {
			super(inpTags, outTags);
			mFuncs = funcs;
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

	}

	/**
	 * @return
	 */
	public NamedMultivariateJacobianFunction build() {
		final NamedMultivariateJacobianFunction res = new Implementation(//
				mInputTags.get(), //
				mOutputTags.get(), mFuncs);
		return res;
	}

}
