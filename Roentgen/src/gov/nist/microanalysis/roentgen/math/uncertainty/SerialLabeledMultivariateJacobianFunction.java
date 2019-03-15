package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.text.NumberFormat;
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

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.utility.FastIndex;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * Turns a calculation based on a sequential set of
 * {@link LabeledMultivariateJacobianFunction} steps into a single
 * {@link SerialLabeledMultivariateJacobianFunction}.
 * </p>
 * <p>
 * Starting with mStep.get(0), the {@link LabeledMultivariateJacobianFunction}
 * is evaluated against the input variables. The input variables plus the output
 * of step 0 become the input to step 1 and then the input variables plus the
 * output of step 0 and 1 become the input to step 2. In this way, complex
 * multistage calculations can be build from a series of simple steps.
 * </p>
 * <p>
 * The utility of this depends upon the chain rule which allows us to calculate
 * &delta;f(g(x))/&delta;x = (&delta;f/&delta;y)(&delta;y/&delta;x) =
 * (&delta;f(y)/&delta;y)&delta;g(x)/&delta;x. At each step, all we need to
 * compute are the partial derivatives for each function with respect to the
 * input variables - ie. the Jacobian.
 * </p>
 *
 * @author Nicholas
 */
public class SerialLabeledMultivariateJacobianFunction //
		extends LabeledMultivariateJacobianFunction //
		implements ILabeledMultivariateFunction, IToHTML {

	private final String mName;

	/**
	 * Should we keep the function inputs as outputs?
	 */
	private final boolean mKeepInputs;

	/**
	 * Steps in order from inner-most to outer-most. The output of the inner steps
	 * become the input of the subsequent outer.
	 */
	private final List<LabeledMultivariateJacobianFunction> mSteps = new ArrayList<>();

	/**
	 * A list of the subset of outputs required to be retained from each step.
	 */
	private final List<FastIndex<? extends Object>> mOutputs = new ArrayList<>();

	public static List<? extends Object> allOutputs( //
			final List<LabeledMultivariateJacobianFunction> steps, //
			final boolean keepInputs //
	) {
		final FastIndex<Object> res = new FastIndex<>();
		for (int i = 0; i < steps.size(); ++i) {
			final LabeledMultivariateJacobianFunction step = steps.get(i);
			if (keepInputs)
				res.addMissing(step.getInputLabels());
			res.addMissing(step.getOutputLabels());
		}
		return new FastIndex<>(res);
	}

	/**
	 * Determines the smallest list of labels required as input to step number
	 * 'step' and all subsequent steps.
	 *
	 * @param step  First step to include
	 * @param steps The list of steps
	 * @return Set&lt;Object&gt; The set of labels
	 */
	private static FastIndex<Object> minimumInputRequirements(//
			final int step, //
			final List<LabeledMultivariateJacobianFunction> steps, //
			final List<? extends Object> outputs //
	) {
		final Set<Object> inputs = new HashSet<>();
		inputs.addAll(outputs);
		for (int st = steps.size() - 1; st >= step; --st) {
			final LabeledMultivariateJacobianFunction nmjf = steps.get(st);
			inputs.removeAll(nmjf.getOutputLabels());
			inputs.addAll(nmjf.getInputLabels());
			if (nmjf instanceof ImplicitMeasurementModel)
				inputs.addAll(nmjf.getOutputLabels());
		}
		return new FastIndex<>(new ArrayList<>(inputs));
	}

	public SerialLabeledMultivariateJacobianFunction( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> steps, //
			final List<? extends Object> outputs, //
			final boolean keepInputs //
	) throws ArgumentException {
		super(minimumInputRequirements(0, steps, outputs), outputs);
		mName = name;
		mKeepInputs = keepInputs;
		for (int step = 0; step < steps.size(); ++step) {
			if (step != 0)
				mOutputs.add(minimumInputRequirements(step, steps, outputs));
			mSteps.add(steps.get(step));
		}
		mOutputs.add(new FastIndex<>(outputs));
	}

	public SerialLabeledMultivariateJacobianFunction( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> steps, //
			final List<? extends Object> outputs //
	) throws ArgumentException {
		this(name, steps, outputs, false);
	}

	public SerialLabeledMultivariateJacobianFunction( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> steps, //
			final boolean keepInputs //
	) throws ArgumentException {
		this(name, steps, allOutputs(steps, keepInputs), keepInputs);
	}

	public SerialLabeledMultivariateJacobianFunction( //
			final String name, //
			final List<LabeledMultivariateJacobianFunction> steps //
	) throws ArgumentException {
		this(name, steps, false);
	}

	public SerialLabeledMultivariateJacobianFunction(//
			final SerialLabeledMultivariateJacobianFunction func, //
			final List<? extends Object> outputs //
	) throws ArgumentException {
		this(func.mName, func.mSteps, outputs, func.mKeepInputs);
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
	public LabeledMultivariateJacobianFunction getStep(final int idx) {
		return mSteps.get(idx);
	}

	/**
	 * Returns the index of the step in which the value associated with the
	 * specified label is required as input.
	 *
	 * @param label A label associated with an NamedMultivariateJacobianFunction
	 *              input
	 * @returns int The index of the NamedMultivariateJacobianFunction in which this
	 *          label is requested as input.
	 */
	final public int findAsInput(final Object label) {
		for (int i = 0; i < mSteps.size(); ++i)
			if (mSteps.get(i).getInputLabels().contains(label))
				return i;
		return -1;
	}

	@Override
	public String toString() {
		return mName.toString();
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE:
			return HTML.escape(mName) + ": A " + Integer.toString(mSteps.size()) + "-Step Calculation";
		case NORMAL: {
			final Table tbl = new Table();
			for (int si = 0; si < mSteps.size(); ++si) {
				final String html = HTML.toHTML(mSteps.get(si), Mode.NORMAL);
				tbl.addRow(Table.td("Step " + si), Table.td(html));
			}
			return HTML.subHeader(HTML.escape(mName)) + tbl.toHTML(mode);
		}
		default:
		case VERBOSE: {
			final Report report = new Report(mName + " Multi-Step Calculation Report");
			final Map<Object, Integer> valLabels = new HashMap<>();
			final Map<Object, Integer> usage = new HashMap<>();
			for (final Object inp : getInputLabels()) {
				valLabels.put(inp, 0);
				usage.put(inp, 0);
			}
			// report.addHeader(HTML.escape(mName) + ": A Multi-Step
			// Calculation");
			for (int i = 0; i < mSteps.size(); ++i) {
				report.addSubHeader("Step " + Integer.toString(i + 1));
				final LabeledMultivariateJacobianFunction func = mSteps.get(i);
				final Table tbl = new Table();
				tbl.addRow(Table.th("Input Variable"), Table.th("Source"));
				for (final Object inlabel : func.getInputLabels()) {
					final Integer ss = valLabels.get(inlabel);
					String item = null;
					if (ss == null)
						item = HTML.error("Not defined.");
					else {
						if (ss.intValue() == 0)
							item = "Input variable";
						else
							item = "Calculated in Step " + Integer.toString(ss.intValue());
						usage.put(inlabel, usage.get(inlabel) + 1);
					}
					tbl.addRow(Table.td(HTML.toHTML(inlabel, Mode.TERSE)), Table.td(item));
				}
				report.add(tbl);
				for (final Object outlabel : func.getOutputLabels()) {
					if (!valLabels.containsKey(outlabel)) {
						valLabels.put(outlabel, i + 1);
						usage.put(outlabel, 0);
					} else {
						report.addHTML(
								HTML.error("Output error in Step " + Integer.toString(i + 1) + " - " + func.toString()
										+ ": Attempting to redefine " + HTML.toHTML(outlabel, Mode.TERSE) + "."));
					}
				}
			}
			report.addSubHeader("Input Utilization Table");
			{
				final Table tbl = new Table();
				tbl.addRow(Table.th("Input"), Table.th("Use Count"));
				for (final Object label : getInputLabels())
					tbl.addRow(Table.td(HTML.toHTML(label, Mode.TERSE)),
							usage.containsKey(label) ? Table.td(usage.get(label)) : Table.td("Missing"));
				report.add(tbl);
			}
			report.addSubHeader("Output Utilization Table");
			{
				final Table tbl = new Table();
				tbl.addRow(Table.th("Output"), Table.th("Defined"), Table.th("Use Count"));
				for (final Object out : getOutputLabels())

					tbl.addRow(Table.td(HTML.toHTML(out, Mode.TERSE)), //
							Table.td(valLabels.getOrDefault(out, 0).toString()), //
							Table.td(usage.getOrDefault(out, 0).toString()));
				report.add(tbl);
			}
			return report.toHTML(Mode.VERBOSE);
		}
		}
	}

	protected void dumpCurrentValues(final int step, final List<? extends Object> index, final RealVector vals) {
		if (sDump != null) {
			final StringBuffer sb = new StringBuffer();
			final NumberFormat nf = new HalfUpFormat("0.00E0");
			sb.append("Step[" + step + "] -> ");
			sb.append(toString() + " : ");
			for (int i = 0; i < index.size(); ++i) {
				if (i != 0)
					sb.append(",");
				sb.append(index.get(i));
				sb.append("=");
				sb.append(nf.format(vals.getEntry(i)));
			}
			sb.append("\n");
			sDump.append(sb.toString());
		}
	}

	/**
	 * This function builds the output vector and Jacobian by calling in sequence
	 * the steps in the calculation and building cumulatively the result Jacobian.
	 * After each step, the outputs from that step become available as inputs to
	 * subsequent steps.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 *      value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		// The current set of output values
		RealVector currVals = point;
		// Cumulative Jacobian
		RealMatrix cumJac = null;
		List<? extends Object> currInputs = getInputLabels();
		// The rows and cols associated with the cumJac
		for (int step = 0; step < mSteps.size(); ++step) {
			dumpCurrentValues(step, currInputs, currVals);
			// Initialize and call the 'func' associated with this step
			final LabeledMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			if (func instanceof ImplicitMeasurementModel) {
				final ImplicitMeasurementModel imm = (ImplicitMeasurementModel) func;
				// Add in output values as constants...
				final Map<Object, Double> consts = new HashMap<>();
				for (final Object label : func.getOutputLabels())
					consts.put(label, currVals.getEntry(currInputs.indexOf(label)));
				imm.initializeConstants(consts);
			}
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputLabels();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object label = fin.get(i);
				assert currInputs.indexOf(label) != -1 : //
				label + " not assigned for " + func + " (step " + step + ")";
				funcPoint.setEntry(i, currVals.getEntry(currInputs.indexOf(label)));
			}
			func.dumpArguments(funcPoint, this);
			final Pair<RealVector, RealMatrix> fres = func.evaluate(funcPoint);
			final RealVector vres = fres.getFirst();
			final RealMatrix mres = fres.getSecond();

			final FastIndex<? extends Object> nextInps = mOutputs.get(step);
			final RealVector newVals = new ArrayRealVector(nextInps.size());
			final RealMatrix newJac = MatrixUtils.createRealMatrix(nextInps.size(), currInputs.size());

			final int[] cidx = new int[func.getInputDimension()];
			for (int cI = 0; cI < cidx.length; ++cI)
				cidx[cI] = currInputs.indexOf(func.getInputLabel(cI));

			for (int r = 0; r < nextInps.size(); ++r) {
				final Object nextVal = nextInps.get(r);
				final int rI = currInputs.indexOf(nextVal);
				if (rI != -1) { // One of the current inputs...
					assert currInputs.get(rI).equals(nextVal) : "1: " + nextVal + " for " + func + " in " + this;
					newVals.setEntry(r, currVals.getEntry(rI));
					newJac.setEntry(r, rI, 1.0);
				} else { // function inputs
					final int frI = func.outputIndex(nextVal);
					assert func.getOutputLabel(frI).equals(nextVal) : //
					"2: " + nextVal + " for " + func + " in " + this;
					assert frI != -1 : //
					nextVal + " is missing in SerialLMJF";
					newVals.setEntry(r, vres.getEntry(frI));
					for (int cI = 0; cI < cidx.length; ++cI)
						newJac.setEntry(r, cidx[cI], mres.getEntry(frI, cI));
				}
			}
			currVals = newVals;
			cumJac = cumJac == null ? newJac : newJac.multiply(cumJac);
			currInputs = nextInps;
			assert cumJac.getRowDimension() == nextInps.size();
		}
		return Pair.create(currVals, cumJac);
	}

	/**
	 * Implements the optimized version of the evaluate function. (non-Javadoc)
	 *
	 * @see gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction#optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(final RealVector point) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		// The current set of input values
		RealVector currVals = point;
		List<? extends Object> currInputs = getInputLabels();
		for (int step = 0; step < mSteps.size(); ++step) {
			dumpCurrentValues(step, currInputs, currVals);
			// Initialize and call the 'func' associated with this step
			final LabeledMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			if (func instanceof ImplicitMeasurementModel) {
				final ImplicitMeasurementModel imm = (ImplicitMeasurementModel) func;
				// Add in output values as constants...
				final Map<Object, Double> consts = new HashMap<>();
				for (final Object label : func.getOutputLabels())
					consts.put(label, currVals.getEntry(currInputs.indexOf(label)));
				imm.initializeConstants(consts);
			}
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputLabels();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object label = fin.get(i);
				assert currInputs.indexOf(label) != -1 : //
				label + " not assigned for " + func + " (step " + step + ")";
				funcPoint.setEntry(i, currVals.getEntry(currInputs.indexOf(label)));
			}
			func.dumpArguments(funcPoint, this);
			final RealVector vres = func.compute(funcPoint);

			final FastIndex<? extends Object> nextInps = mOutputs.get(step);
			final RealVector newVals = new ArrayRealVector(nextInps.size());

			for (int r = 0; r < nextInps.size(); ++r) {
				final Object nextVal = nextInps.get(r);
				final int rI = currInputs.indexOf(nextVal);
				if (rI != -1) { // One of the current inputs...
					assert currInputs.get(rI).equals(nextVal): //
					"1: " + nextVal + " for " + func + " in " + this;
					newVals.setEntry(r, currVals.getEntry(rI));
				} else { // function inputs
					final int frI = func.outputIndex(nextVal);
					assert func.getOutputLabel(frI).equals(nextVal): //
					"2: " + nextVal + " for " + func + " in " + this;
					assert frI != -1 : nextVal + " is missing in SerialLMJF";
					newVals.setEntry(r, vres.getEntry(frI));
				}
			}
			currVals = newVals;
			currInputs = nextInps;
		}
		return currVals;
	}

}
