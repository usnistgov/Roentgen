package gov.nist.juncertainty;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * Turns a calculation based on a sequential set of
 * {@link ExplicitMeasurementModel} steps into a single
 * {@link CompositeMeasurementModel}.
 * </p>
 * <p>
 * Starting with mStep.get(0), the {@link ExplicitMeasurementModel} is evaluated
 * against the input variables. The input variables plus the output of step 0
 * become the input to step 1 and then the input variables plus the output of
 * step 0 and 1 become the input to step 2. In this way, complex multistage
 * calculations can be build from a series of simple steps.
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
public class CompositeMeasurementModel<G> //
		extends ExplicitMeasurementModel<G, G> //
		implements ILabeledMultivariateFunction<G, G>, IToHTML {

	public static <J, K> List<K> allOutputs(
			final List<? extends ExplicitMeasurementModel<? extends J, ? extends K>> steps //
	) throws ArgumentException {
		final List<K> res = new ArrayList<>();
		for (int i = 0; i < steps.size(); ++i) {
			final ExplicitMeasurementModel<? extends J, ? extends K> step = steps.get(i);
			for (final K label : step.getOutputLabels()) {
				if (res.contains(label))
					for (int j = 0; j < i; ++j) {
						final ExplicitMeasurementModel<? extends J, ? extends K> stepj = steps.get(j);
						if (stepj.outputIndex(label) != 1)
							throw new ArgumentException(
									label + "defined in " + stepj + " is also defined in step " + step.toString());
					}
				res.add(label);
			}
		}
		return res;
	}

	/**
	 * If outputs is empty, then builds a list containing all the outputs in steps.
	 * If outputs is not empty, it returns the intersection of the desired list of
	 * outputs and the outputs in steps.
	 *
	 * @param steps   The calculation steps
	 * @param outputs A list of desired output labels or the empty list for all
	 *                outputs
	 * @return List&lt;H&gt; A list of available output labels
	 * @throws ArgumentException
	 */
	public static <H> List<H> buildOutputs(
			final List<? extends ExplicitMeasurementModel<? extends H, ? extends H>> steps, //
			final List<? extends H> outputs //
	) throws ArgumentException {
		final List<H> res = allOutputs(steps);
		if (outputs.size() > 0)
			res.retainAll(outputs);
		res.sort((final H o1, final H o2) -> o1.toString().compareTo(o2.toString()));
		return res;
	}

	/**
	 * Determines the smallest list of labels required as input to step number
	 * 'step' and all subsequent steps plus those labels computed in this and
	 * previous steps that are required as the final output.
	 *
	 * @param step  First step to include
	 * @param steps The list of steps
	 * @return Set&lt;Object&gt; The set of labels
	 */
	private static <G> List<G> labelsToRetain(
			//
			final int step, //
			final List<? extends ExplicitMeasurementModel<? extends G, ? extends G>> steps, //
			final List<? extends G> outputs //
	) {
		final Set<G> inputs = new HashSet<>();
		// Build a list of all the inputs required by this and subsequent steps
		for (int st = step; st < steps.size(); ++st) {
			final ExplicitMeasurementModel<? extends G, ? extends G> func = steps.get(st);
			inputs.addAll(func.getInputLabels());
		}
		// Remove all the outputs from this and subsequent steps
		for (int st = step; st < steps.size() - 1; ++st) {
			final ExplicitMeasurementModel<? extends G, ? extends G> func = steps.get(st);
			inputs.removeAll(func.getOutputLabels());
		}
		// Add all the outputs from previous steps that are in the final output list
		final Set<G> retain = new HashSet<>();
		for (int st = step - 1; st >= 0; --st) {
			final ExplicitMeasurementModel<? extends G, ? extends G> func = steps.get(st);
			retain.addAll(func.getOutputLabels());
		}
		retain.retainAll(outputs);
		inputs.addAll(retain);
		final ArrayList<G> result = new ArrayList<G>(inputs);
		result.sort((final G o1, final G o2) -> o1.toString().compareTo(o2.toString()));
		return Collections.unmodifiableList(result);
	}

	private final String mName;

	/**
	 * Steps in order from inner-most to outer-most. The output of the inner steps
	 * become the input of the subsequent outer.
	 */
	private final List<ExplicitMeasurementModel<? extends G, ? extends G>> mSteps = new ArrayList<>();

	/**
	 * A list of the subset of outputs required to be retained from each step.
	 */
	private final List<List<G>> mOutputs;

	public CompositeMeasurementModel(
			//
			final CompositeMeasurementModel<G> func, //
			final List<G> outputs //
	) throws ArgumentException {
		this(func.mName, func.mSteps, outputs);
	}

	public CompositeMeasurementModel(
			//
			final String name, //
			final List<? extends ExplicitMeasurementModel<? extends G, ? extends G>> steps //
	) throws ArgumentException {
		this(name, steps, Collections.emptyList());
	}

	public CompositeMeasurementModel(
			final String name, //
			final List<? extends ExplicitMeasurementModel<? extends G, ? extends G>> steps, //
			final List<? extends G> outputLabels //
	) throws ArgumentException {
		super(labelsToRetain(0, steps, buildOutputs(steps, outputLabels)), buildOutputs(steps, outputLabels));
		mName = name;
		mSteps.addAll(steps);
		final List<List<G>> outList = new ArrayList<>();
		for (int step = 1; step < steps.size(); ++step)
			outList.add(labelsToRetain(step, steps, getOutputLabels()));
		outList.add(getOutputLabels());
		mOutputs = Collections.unmodifiableList(outList);
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
	final public int findAsInput(
			final Object label
	) {
		for (int i = 0; i < mSteps.size(); ++i)
			if (mSteps.get(i).getInputLabels().contains(label))
				return i;
		return -1;
	}

	/**
	 * Returns the NamedMultivariateJacobianFunction associated with the specific
	 * step by index.
	 *
	 * @param idx On range [0, getStepCount())
	 * @return NamedMultivariateJacobianFunction
	 */
	public ExplicitMeasurementModel<? extends G, ? extends G> getStep(
			final int idx
	) {
		return mSteps.get(idx);
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
	 * Implements the optimized version of the evaluate function. (non-Javadoc)
	 *
	 * @see gov.nist.juncertainty.ILabeledMultivariateFunction#optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(
			final RealVector point
	) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		// The current set of input values. Starts as point and evolves.
		RealVector currVals = point;
		List<? extends G> currInputs = getInputLabels();
		for (int step = 0; step < mSteps.size(); ++step) {
			dumpCurrentValues(step, currInputs, currVals);
			// Initialize and call the 'func' associated with this step
			final ExplicitMeasurementModel<? extends G, ? extends G> func = mSteps.get(step);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends G> fin = func.getInputLabels();
			for (int i = 0; i < fin.size(); ++i) {
				final int cidx = currInputs.indexOf(fin.get(i));
				assert cidx != -1 : fin.get(i) + " not assigned for " + func + " (step " + step + ")";
				funcPoint.setEntry(i, currVals.getEntry(cidx));
			}
			// Initialize constants in ImplicitMeasurementModels
			func.dumpArguments(funcPoint, this);
			func.applyAdditionalInputs(getAdditionalInputs());
			RealVector vres = null;
			if (func instanceof ILabeledMultivariateFunction)
				vres = ((ILabeledMultivariateFunction<?, ?>) func).optimized(funcPoint);
			else
				vres = func.evaluate(funcPoint).getFirst();

			final List<G> nextInputs = mOutputs.get(step);
			if ((nextInputs.size() == 0) || (currInputs.size() == 0)) {
				System.err.println("Warning: ");
				System.err.println("Step " + step + " in " + toString() + " at " + func);
				System.err.println("Current inputs = " + currInputs);
				System.err.println("Next inputs = " + nextInputs);
			}
			final RealVector newVals = new ArrayRealVector(nextInputs.size());

			for (int r = 0; r < nextInputs.size(); ++r) {
				final G nextVal = nextInputs.get(r);
				final int frI = func.outputIndex(nextVal);
				if (frI != -1) {
					newVals.setEntry(r, vres.getEntry(frI));
				} else {
					final int rI = currInputs.indexOf(nextVal);
					assert rI != -1 : //
					nextVal + " is missing in " + toString() + " at " + func;
					newVals.setEntry(r, currVals.getEntry(rI));
				}
			}
			currVals = newVals;
			currInputs = nextInputs;
		}
		return currVals;
	}

	@Override
	public String toHTML(
			final Mode mode
	) {
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
				final ExplicitMeasurementModel<? extends G, ? extends G> func = mSteps.get(i);
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

	@Override
	public String toString() {
		return mName.toString();
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
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
			) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		// The current set of output values
		RealVector currVals = point;
		// Cumulative Jacobian
		RealMatrix cumJac = null;
		List<G> currInputs = getInputLabels();
		// The rows and cols associated with the cumJac
		for (int step = 0; step < mSteps.size(); ++step) {
			dumpCurrentValues(step, currInputs, currVals);
			// Initialize and call the 'func' associated with this step
			final ExplicitMeasurementModel<? extends G, ? extends G> func = mSteps.get(step);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends G> fin = func.getInputLabels();

			final int[] cidx = new int[func.getInputDimension()];
			for (int cI = 0; cI < cidx.length; ++cI) {
				cidx[cI] = currInputs.indexOf(fin.get(cI));
				assert cidx[cI] != -1 : fin.get(cI) + " not assigned for " + func + " (step " + step + ")";
				funcPoint.setEntry(cI, currVals.getEntry(cidx[cI]));
			}
			func.dumpArguments(funcPoint, this);
			func.applyAdditionalInputs(getAdditionalInputs());
			final Pair<RealVector, RealMatrix> fres = func.evaluate(funcPoint);
			final RealVector vres = fres.getFirst();
			final RealMatrix mres = fres.getSecond();

			final List<G> nextInputs = mOutputs.get(step);
			final RealVector newVals = new ArrayRealVector(nextInputs.size());
			if ((nextInputs.size() == 0) || (currInputs.size() == 0)) {
				System.err.println("Warning: ");
				System.err.println("Step " + step + " in " + toString() + " at " + func);
				System.err.println("Current inputs = " + currInputs);
				System.err.println("Next inputs = " + nextInputs);
			}
			final RealMatrix newJac = NullableRealMatrix.build(nextInputs.size(), currInputs.size());

			for (int r = 0; r < nextInputs.size(); ++r) {
				final G nextVal = nextInputs.get(r);
				final int frI = func.outputIndex(nextVal);
				if (frI != -1) {
					newVals.setEntry(r, vres.getEntry(frI));
					for (int cI = 0; cI < cidx.length; ++cI)
						newJac.setEntry(r, cidx[cI], mres.getEntry(frI, cI));
				} else {
					final int rI = currInputs.indexOf(nextVal);
					assert rI != -1 : nextVal + " is missing in " + toString() + " at " + func;
					newVals.setEntry(r, currVals.getEntry(rI));
					newJac.setEntry(r, rI, 1.0);
				}
			}
			currVals = newVals;
			cumJac = cumJac == null ? newJac : newJac.multiply(cumJac);
			currInputs = nextInputs;
			assert cumJac.getRowDimension() == nextInputs.size();
		}
		return Pair.create(currVals, cumJac);
	}

	protected void dumpCurrentValues(
			final int step, final List<? extends G> index, final RealVector vals
			) {
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

}
