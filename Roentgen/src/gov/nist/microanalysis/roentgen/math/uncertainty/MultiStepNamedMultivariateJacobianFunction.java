package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.google.common.collect.Sets;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>
 * Turns a calculation based on a sequential set of
 * {@link NamedMultivariateJacobianFunction} steps into a single
 * {@link MultiStepNamedMultivariateJacobianFunction}.
 * </p>
 * <p>
 * Starting with mStep.get(0), the {@link NamedMultivariateJacobianFunction} is
 * evaluated against the input variables. The input variables plus the output of
 * step 0 become the input to step 1 and then the input variables plus the
 * output of step 0 and 1 become the input to step 2. In this way, complex
 * multistage calculations can be build from a series of simple steps.
 * </p>
 * <p>
 * The utility of this depends upon the chain rule which allows us to calculate
 * df(g(x))/dx = (df/dy)(dy/dx) = (df(y)/dy)dg(x)/dx. At each step, all we need
 * to compute are the partial derivatives for each function with respect to the
 * input variables - ie. the Jacobian.
 * </p>
 *
 * @author Nicholas
 */
public class MultiStepNamedMultivariateJacobianFunction extends NamedMultivariateJacobianFunction
		implements INamedMultivariateFunction, IToHTML {

	/**
	 * Steps in order from inner-most to outer-most. The output of the inner steps
	 * become the input of the subsequent outer.
	 */
	private final List<NamedMultivariateJacobianFunction> mSteps = new ArrayList<>();

	private final String mName;

	private Set<? extends Object> mFinalOutputs;

	private final static List<? extends Object> computeOutputs(final List<NamedMultivariateJacobianFunction> steps)
			throws ArgumentException {
		final Set<Object> sob = new HashSet<>();
		final List<Object> res = new ArrayList<>();
		for (final NamedMultivariateJacobianFunction nmjf : steps) {
			final List<? extends Object> tags = nmjf.getOutputTags();
			for (final Object tag : tags) {
				if (sob.contains(tag))
					throw new ArgumentException("The output " + tag + " has already been defined in an earlier step");
				sob.add(tag);
				res.add(tag);
			}
		}
		return res;
	}

	/**
	 * A list containing the tags associated with all required input values that are
	 * not computed as output values. (Values for these tags must be provided to the
	 * evaluate function.)
	 *
	 * @throws ArgumentException
	 */
	private final static List<? extends Object> computeInputs(final List<NamedMultivariateJacobianFunction> steps)
			throws ArgumentException {
		final Set<Object> sob = new HashSet<>();
		for (final NamedMultivariateJacobianFunction nmjf : steps)
			sob.addAll(nmjf.getInputTags());
		sob.removeAll(computeOutputs(steps));
		return new ArrayList<>(sob);
	}

	/**
	 * Create a MultiStepNamedMultivariateJacobianFunction object to compute the
	 * steps in order, where the output values for earlier steps become available as
	 * input values to later steps.
	 *
	 * @param steps
	 * @throws ArgumentException
	 */
	public MultiStepNamedMultivariateJacobianFunction(final String name,
			final List<NamedMultivariateJacobianFunction> steps) throws ArgumentException {
		super(computeInputs(steps), computeOutputs(steps));
		mName = name;
		mSteps.addAll(steps);
	}

	/**
	 * Returns the number of steps in the calculation
	 *
	 * @return int In range [0, ...)
	 */
	public int getStepCount() {
		return mSteps.size();
	}

	/**
	 * Returns the NamedMultivariateJacobianFunction associated with the specific
	 * step by index.
	 *
	 * @param idx On range [0, getStepCount())
	 * @return NamedMultivariateJacobianFunction
	 */
	public NamedMultivariateJacobianFunction getStep(final int idx) {
		return mSteps.get(idx);
	}

	/**
	 * Returns the index of the step in which the value associated with the
	 * specified tag is requested.
	 *
	 * @param tag A tag associated with an NamedMultivariateJacobianFunction input
	 * @returns int The index of the NamedMultivariateJacobianFunction in which this
	 *          tag is requested as input.
	 */
	final public int findOutput(final Object tag) {
		for (int i = 0; i < mSteps.size(); ++i)
			if (mSteps.get(i).getInputTags().contains(tag))
				return i;
		return -1;
	}

	/**
	 * Returns an array consisting of the tags associated with the m output random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	@Override
	public List<? extends Object> getOutputTags() {
		return mFinalOutputs != null ? new ArrayList<>(mFinalOutputs) : super.getOutputTags();
	}

	/**
	 * A mechanism for optimizing the calculation by keeping a subset of the outputs
	 * at each step. <code>finalOutput</code> specifies which outputs to keep
	 * regardless of whether they are required by a subsequent step. At each step in
	 * the calculation determine which outputs are required by subsequent steps and
	 * and keep only these outputs.
	 *
	 * @param finalOutputs
	 */
	public void trimOutputs(final Set<? extends Object> finalOutputs) {
		mFinalOutputs = finalOutputs;
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE:
			return HTML.escape(mName) + ": A " + Integer.toString(mSteps.size()) + "-Step Calculation";
		case NORMAL: {
			final Table tbl = new Table();
			tbl.addRow(Table.td("Steps"), Table.td(Integer.toString(mSteps.size())));
			final List<? extends Object> inp = getInputTags();
			tbl.addRow(Table.th("Inputs", 2));
			for (int i = 0; i < inp.size(); ++i)
				tbl.addRow(Table.td("Input " + Integer.toString(i + 1)), Table.td(HTML.toHTML(inp.get(i), Mode.TERSE)));
			final List<? extends Object> out = getOutputTags();
			tbl.addRow(Table.th("Outputs", 2));
			for (int i = 0; i < out.size(); ++i)
				tbl.addRow(Table.td("Output " + Integer.toString(i + 1)),
						Table.td(HTML.toHTML(out.get(i), Mode.TERSE)));
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
			// report.addHeader(HTML.escape(mName) + ": A Multi-Step
			// Calculation");
			for (int i = 0; i < mSteps.size(); ++i) {
				report.addSubHeader("Step " + Integer.toString(i + 1));
				final NamedMultivariateJacobianFunction func = mSteps.get(i);
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.fitting.leastsquares.
	 * MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final List<? extends Object> inpTags = getInputTags();
		List<Double> vals = new ArrayList<>(); // give & calculated
		List<Object> valTags = new ArrayList<>(); // assoc. w. vals
		Map<Object, Integer> valMap = new HashMap<>();
		// Add the input tags and values
		for (int i = 0; i < point.getDimension(); ++i) {
			vals.add(point.getEntry(i));
			valTags.add(inpTags.get(i));
			valMap.put(inpTags.get(i), i);
		}
		// Start with a identity matrix for the input values
		RealMatrix expJac = MatrixUtils.createRealIdentityMatrix(point.getDimension());
		for (int step = 0; step < mSteps.size(); ++step) {
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			final List<? extends Object> fin = func.getInputTags();
			final List<? extends Object> fout = func.getOutputTags();
			// Build the vector argument to func and call
			final RealVector pt = new ArrayRealVector(fin.size());
			for (int i = 0; i < pt.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				final Integer idx = valMap.get(tag);
				assert idx != null : "The value associated with " + fin.get(i) + " has not been assigned in Step "
						+ Integer.toString(mSteps.indexOf(func) + 1) + ".";
				pt.setEntry(i, vals.get(idx));
			}
			final int fullInpSize = vals.size();
			final Pair<RealVector, RealMatrix> fres = func.evaluate(pt);
			final RealVector ovals = fres.getFirst();
			final RealMatrix ojac = fres.getSecond();
			if (mFinalOutputs != null) {
				// Build up a list of computed inputs required by subsequent steps
				// and final output
				final Set<Object> reqOutputs = new HashSet<>();
				for (int subStep = step + 1; subStep < mSteps.size(); ++subStep)
					reqOutputs.addAll(mSteps.get(subStep).getInputTags());
				reqOutputs.addAll(mFinalOutputs);
				reqOutputs.removeAll(inpTags);
				// List of computed items available (prev+this-input)
				final Set<Object> avail = new HashSet<>(valTags);
				avail.addAll(fout);
				avail.removeAll(inpTags);
				// Create lists of the tags and values to carry on to the next step
				final List<Object> outTagsx = new ArrayList<>(inpTags);
				outTagsx.addAll(Sets.intersection(avail, reqOutputs));
				final List<Double> outValsx = new ArrayList<>();
				final Map<Object, Integer> outValsxMap = new HashMap<>();
				for (int j = 0; j < outTagsx.size(); ++j) {
					final Object outTag = outTagsx.get(j);
					final Integer idx = valMap.get(outTag);
					final Double val = idx != null ? vals.get(idx) : ovals.getEntry(fout.indexOf(outTag));
					outValsx.add(val);
					outValsxMap.put(outTag, Integer.valueOf(j));
				}
				// Fill jacobian with combination of old and new
				final RealMatrix ojacx = MatrixUtils.createRealMatrix(outTagsx.size(), valTags.size());
				for (int inIdx = 0; inIdx < valTags.size(); ++inIdx) {
					final int finIdx = func.inputIndex(valTags.get(inIdx)); // May be
																			// -1
					for (int outIdx = 0; outIdx < outTagsx.size(); ++outIdx) {
						final int foutIdx = func.outputIndex(outTagsx.get(outIdx)); // May
																					// be
																					// -1
						if ((finIdx != -1) && (foutIdx != -1))
							ojacx.setEntry(outIdx, inIdx, ojac.getEntry(foutIdx, finIdx));
						else // the output value was carried along from the input
							ojacx.setEntry(outIdx, inIdx, valTags.get(inIdx).equals(outTagsx.get(outIdx)) ? 1.0 : 0.0);
					}
				}
				// Outputs of step become inputs of step+1
				vals = outValsx;
				valTags = outTagsx;
				valMap = outValsxMap;
				// Multiply this Jacobian times the inner combined Jacobian
				expJac = ojacx.multiply(expJac);
			} else {
				// Add output values to the end of vals and valTags
				for (int i = 0; i < ovals.getDimension(); ++i) {
					assert !valTags.contains(fout.get(i)) : fout.get(i) + " computed in Step " + i
							+ " has previously been defined.";
					valTags.add(fout.get(i));
					vals.add(ovals.getEntry(i));
					valMap.put(fout.get(i), valTags.size() - 1);
				}
				assert valTags.size() == vals.size();
				assert valMap.size() == vals.size();
				// Build the expanded jacobian matrix
				assert fullInpSize == expJac.getRowDimension();
				final int outSize = ojac.getRowDimension();
				final int rowDim = outSize + fullInpSize;
				final int colDim = fullInpSize;
				// Build the expanded output Jacobian
				final RealMatrix ojacx = MatrixUtils.createRealMatrix(rowDim, colDim);
				// Prepend an identity matrix to "carry" previous inputs
				for (int i = 0; i < fullInpSize; ++i)
					ojacx.setEntry(i, i, 1.0);
				// Append the output Jacobian
				for (int i = 0; i < fin.size(); ++i) {
					final int col = valMap.get(fin.get(i));
					for (int r = 0; r < outSize; ++r)
						ojacx.setEntry(r + fullInpSize, col, ojac.getEntry(r, i));
				}
				// Multiply this Jacobian times the inner combined Jacobian
				expJac = ojacx.multiply(expJac);
			}
		}
		final List<? extends Object> outs = getOutputTags();
		final int outValSize = outs.size();
		final int inpValSize = inpTags.size();
		assert expJac.getRowDimension() == outValSize + inpValSize : expJac.getRowDimension() + "!=" + outValSize + "+"
				+ inpValSize;
		assert vals.size() == expJac.getRowDimension();
		assert expJac.getColumnDimension() == inpValSize;
		if (mFinalOutputs == null) {
			final RealVector resVals = new ArrayRealVector(outValSize);
			for (int r = 0; r < outValSize; ++r) {
				assert outs.get(r) == valTags.get(r + inpValSize) : toString();
				resVals.setEntry(r, vals.get(r + inpValSize));
			}
			final RealMatrix resJac = expJac.getSubMatrix(inpValSize, inpValSize + outValSize - 1, 0, inpValSize - 1);
			return Pair.create(resVals, resJac);
		} else {
			final RealVector resVals = new ArrayRealVector(outValSize);
			final RealMatrix resJac = MatrixUtils.createRealMatrix(outValSize, inpValSize);
			for (int row = 0; row < outValSize; ++row) {
				final int valIdx = valMap.get(outs.get(row));
				resVals.setEntry(row, vals.get(valIdx));
				for (int col = 0; col < inpValSize; ++col)
					resJac.setEntry(row, col, expJac.getEntry(valIdx, col));
			}
			return Pair.create(resVals, resJac);
		}
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final List<? extends Object> inpTags = getInputTags();
		List<Double> vals = new ArrayList<>(); // give & calculated
		List<Object> valTags = new ArrayList<>(); // assoc. w. vals
		Map<Object, Integer> valMap = new HashMap<>();
		// Add the input tags and values
		for (int i = 0; i < point.getDimension(); ++i) {
			vals.add(point.getEntry(i));
			valTags.add(inpTags.get(i));
			valMap.put(inpTags.get(i), i);
		}
		for (int step = 0; step < mSteps.size(); ++step) {
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			final List<? extends Object> fin = func.getInputTags();
			final List<? extends Object> fout = func.getOutputTags();
			// Build the vector argument to func and call
			final RealVector pt = new ArrayRealVector(fin.size());
			for (int i = 0; i < pt.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				final Integer idx = valMap.get(tag);
				assert idx != null : "The value associated with " + fin.get(i) + " has not been assigned in Step "
						+ Integer.toString(mSteps.indexOf(func) + 1) + ".";
				pt.setEntry(i, vals.get(idx));
			}
			final RealVector ovals = func.compute(pt);
			if (mFinalOutputs != null) {
				// Build up a list of computed inputs required by subsequent steps
				// and final output
				final Set<Object> reqOutputs = new HashSet<>();
				for (int subStep = step + 1; subStep < mSteps.size(); ++subStep)
					reqOutputs.addAll(mSteps.get(subStep).getInputTags());
				reqOutputs.addAll(mFinalOutputs);
				reqOutputs.removeAll(inpTags);
				// List of computed items available (prev+this-input)
				final Set<Object> avail = new HashSet<>(valTags);
				avail.addAll(fout);
				avail.removeAll(inpTags);
				// Create lists of the tags and values to carry on to the next step
				final List<Object> outTagsx = new ArrayList<>(inpTags);
				outTagsx.addAll(Sets.intersection(avail, reqOutputs));
				final List<Double> outValsx = new ArrayList<>();
				final Map<Object, Integer> outValsxMap = new HashMap<>();
				for (int j = 0; j < outTagsx.size(); ++j) {
					final Object outTag = outTagsx.get(j);
					final Integer idx = valMap.get(outTag);
					final Double val = idx != null ? vals.get(idx) : ovals.getEntry(fout.indexOf(outTag));
					outValsx.add(val);
					outValsxMap.put(outTag, Integer.valueOf(j));
				}
				// Outputs of step become inputs of step+1
				vals = outValsx;
				valTags = outTagsx;
				valMap = outValsxMap;
			} else {
				// Add output values to the end of vals and valTags
				for (int i = 0; i < ovals.getDimension(); ++i) {
					assert !valTags.contains(fout.get(i)) : fout.get(i) + " computed in Step " + i
							+ " has previously been defined.";
					valTags.add(fout.get(i));
					vals.add(ovals.getEntry(i));
					valMap.put(fout.get(i), valTags.size() - 1);
				}
				assert valTags.size() == vals.size();
				assert valMap.size() == vals.size();
			}
		}
		final List<? extends Object> outs = getOutputTags();
		final int outValSize = outs.size();
		final int inpValSize = inpTags.size();
		final RealVector resVals = new ArrayRealVector(outValSize);
		if (mFinalOutputs == null) {
			for (int r = 0; r < outValSize; ++r) {
				assert outs.get(r) == valTags.get(r + inpValSize) : toString();
				resVals.setEntry(r, vals.get(r + inpValSize));
			}
		} else {
			for (int row = 0; row < outValSize; ++row) {
				final int valIdx = valMap.get(outs.get(row));
				resVals.setEntry(row, vals.get(valIdx));
			}
		}
		return resVals;
	}
}
