package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collections;
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

/**
 * <p>
 * Turns a calculation based on a sequential set of
 * {@link LabeledMultivariateJacobianFunction} steps into a single
 * {@link SerialNamedMultivariateJacobianFunction}.
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
public class SerialNamedMultivariateJacobianFunction extends LabeledMultivariateJacobianFunction //
		implements ILabeledMultivariateFunction, IToHTML {

	/**
	 * Steps in order from inner-most to outer-most. The output of the inner steps
	 * become the input of the subsequent outer.
	 */
	private final List<LabeledMultivariateJacobianFunction> mSteps = new ArrayList<>();

	private final String mName;

	private Set<? extends Object> mFinalOutputs;

	private final static List<? extends Object> computeOutputs( //
			final List<LabeledMultivariateJacobianFunction> steps //
	) throws ArgumentException {
		final Set<Object> sob = new HashSet<>();
		final List<Object> res = new ArrayList<>();
		for (final LabeledMultivariateJacobianFunction nmjf : steps) {
			final List<? extends Object> labels = nmjf.getOutputLabels();
			for (final Object label : labels) {
				if (sob.contains(label))
					throw new ArgumentException("The output " + label + " has already been defined in an earlier step");
				sob.add(label);
				res.add(label);
			}
		}
		return res;
	}

	/**
	 * A list containing the labels associated with all required input values that
	 * are not computed as output values. (Values for these labels must be provided
	 * to the evaluate function.)
	 *
	 * @throws ArgumentException
	 */
	private final static List<? extends Object> computeInputs( //
			final List<LabeledMultivariateJacobianFunction> steps //
	) throws ArgumentException {
		return new ArrayList<>(minimumInputRequirements(0, steps));
	}

	/**
	 * Determines the smallest list of labels required as input to step number
	 * 'step' and all subsequent steps.
	 *
	 * @param step  First step to include
	 * @param steps The list of steps
	 * @return Set&lt;Object&gt; The set of labels
	 */
	private static Set<Object> minimumInputRequirements(//
			final int step, //
			final List<LabeledMultivariateJacobianFunction> steps //
	) {
		final Set<Object> inputs = new HashSet<>();
		for (int st = steps.size() - 1; st >= step; --st) {
			final LabeledMultivariateJacobianFunction nmjf = steps.get(st);
			inputs.removeAll(nmjf.getOutputLabels());
			inputs.addAll(nmjf.getInputLabels());
			if (nmjf instanceof ImplicitMeasurementModel)
				inputs.addAll(nmjf.getOutputLabels());
		}
		return inputs;
	}

	/**
	 * Create a MultiStepNamedMultivariateJacobianFunction object to compute the
	 * steps in order, where the output values for earlier steps become available as
	 * input values to later steps.
	 *
	 * @param steps
	 * @throws ArgumentException
	 */
	public SerialNamedMultivariateJacobianFunction(final String name,
			final List<LabeledMultivariateJacobianFunction> steps) throws ArgumentException {
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

	/**
	 * Returns an array consisting of the labels associated with the m output random
	 * variables.
	 *
	 * @return List&lt;Object&gt;
	 */
	@Override
	public List<? extends Object> getOutputLabels() {
		return mFinalOutputs != null ? new ArrayList<>(mFinalOutputs) : super.getOutputLabels();
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
		mFinalOutputs = new HashSet<>(getOutputLabels());
		mFinalOutputs.retainAll(finalOutputs);
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
						report.addHTML(HTML.error("Output error in Step " + Integer.toString(i + 1)
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
		final List<? extends Object> inpLabels = getInputLabels();
		// inputs contains labels and the associated values (inputs and calculated
		// values)
		final Map<Object, Double> inputs = new HashMap<>(); // label -> index in valLabels
		// Add the input labels and values
		for (int i = 0; i < point.getDimension(); ++i)
			inputs.put(inpLabels.get(i), point.getEntry(i));
		// Cumulative Jacobian: Start with an identity matrix representing the partials
		// the input values relative to themselves and accumulate.
		RealMatrix cumJac = MatrixUtils.createRealIdentityMatrix(point.getDimension());
		// prevRows keeps track of the rows in the previous step's Jacobian
		List<? extends Object> prevRows = Collections.unmodifiableList(inpLabels);

		for (int step = 0; step < mSteps.size(); ++step) {
			// Initialize and call the 'func' associated with this step
			final LabeledMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			final List<? extends Object> fout = func.getOutputLabels();
			if (func instanceof ImplicitMeasurementModel) {
				// Add in output values as constants...
				for (final Object label : fout)
					consts.put(label, inputs.get(label));
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputLabels();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object label = fin.get(i);
				assert label != null : "Missmatch in length between input labels and point dimension.";
				assert inputs.containsKey(label) : "The value for label " + label + " has not been assigned by Step "
						+ step + ".";
				funcPoint.setEntry(i, inputs.get(label));
			}
			final Pair<RealVector, RealMatrix> fres = func.evaluate(funcPoint);
			final RealVector ovals = fres.getFirst();
			final RealMatrix ojac = fres.getSecond();
			// (ovals, ojac) = func.evaluate(funcPoint)
			// ovals[fout], ojac[fout][fin]
			// expJac[expRows][expCols]
			for (int i = 0; i < ovals.getDimension(); ++i)
				inputs.put(fout.get(i), ovals.getEntry(i));
			{ // Build the resulting Jacobian
				final List<? extends Object> currCols = prevRows; // output becomes input
				// Build currRows from prevRows + newItems
				final List<Object> currRows = new ArrayList<>(prevRows); // old + new output
				// foutIdx is where to place the out
				final int[] foutRow = new int[fout.size()];
				if (func instanceof ImplicitMeasurementModel) {
					// Output label will be one of the previous steps output labels also
					for (int i = 0; i < fout.size(); ++i) {
						final int prevIdx = prevRows.indexOf(fout.get(i));
						assert prevIdx != -1 : "Jacobian elements should alread exist for " + fout.get(i);
						foutRow[i] = prevIdx;
					}
				} else {
					// Add the output label to the end off currRows
					for (int i = 0; i < fout.size(); ++i) {
						final Object label = fout.get(i);
						assert prevRows.indexOf(label) == -1 : label + " was already calculated.";
						currRows.add(label);
						foutRow[i] = currRows.size() - 1;
					}
				}
				// Build the extended Jacobian previous inputs and current outputs
				final RealMatrix newJac = MatrixUtils.createRealMatrix(currRows.size(), currCols.size());
				for (int rc = 0; rc < currCols.size(); ++rc)
					newJac.setEntry(rc, rc, fout.contains(currCols.get(rc)) ? 0.0 : 1.0);
				for (int c = 0; c < fin.size(); ++c) {
					final int col = currCols.indexOf(fin.get(c));
					assert col != -1 : fin.get(c) + " was not a row in the previous Jacobian.";
					for (int r = 0; r < foutRow.length; ++r) {
						assert foutRow[r] < newJac.getRowDimension() : fout.get(foutRow[r]);
						assert foutRow[r] != -1;
						newJac.setEntry(foutRow[r], col, ojac.getEntry(r, c));
					}
				}
				/*
				 * Keep the Jacobian row count as small as possible each step by eliminating
				 * rows associated with values that are not required as input to subsequent
				 * steps.
				 */
				if (mFinalOutputs != null) {
					// Determine the minimum set if outputs to retain
					final Set<Object> minOutputs = minimumInputRequirements(step + 1, mSteps);
					minOutputs.addAll(mFinalOutputs);
					// Build a list of rows to keep and their index
					final List<Object> trimmedRows = new ArrayList<>();
					final List<Integer> trimmedIdx = new ArrayList<>();
					for (int r = 0; r < currRows.size(); ++r) {
						final Object label = currRows.get(r);
						if (minOutputs.contains(label)) {
							trimmedRows.add(label);
							trimmedIdx.add(r);
						} else
							inputs.remove(label);
					}
					// Build the trimmed Jacobian
					final RealMatrix trimmedJac = MatrixUtils.createRealMatrix(trimmedIdx.size(), currCols.size());
					for (int r = trimmedIdx.size() - 1; r >= 0; --r)
						trimmedJac.setRow(r, newJac.getRow(trimmedIdx.get(r)));
					// Calculate the cumulative Jacobian
					cumJac = trimmedJac.multiply(cumJac);
					prevRows = Collections.unmodifiableList(trimmedRows);
				} else {
					// Calculate the cumulative Jacobian
					cumJac = newJac.multiply(cumJac);
					prevRows = Collections.unmodifiableList(currRows);
				}
			}
		}
		// Build the output structures
		final List<? extends Object> outs = getOutputLabels();
		final RealVector resVals = new ArrayRealVector(outs.size());
		final RealMatrix resJac = MatrixUtils.createRealMatrix(outs.size(), inpLabels.size());
		{ // Ensure the output rows are ordered correctly.
			assert inputs.size() == cumJac.getRowDimension();
			assert cumJac.getColumnDimension() == inpLabels.size();
			for (int row = outs.size() - 1; row >= 0; --row) {
				final Object rowLabel = outs.get(row);
				resVals.setEntry(row, inputs.get(rowLabel));
				resJac.setRow(row, cumJac.getRow(prevRows.indexOf(rowLabel)));
			}
		}
		return Pair.create(resVals, resJac);
	}

	/**
	 * Implements the optimized version of the evaluate function. (non-Javadoc)
	 * 
	 * @see gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction#optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(final RealVector point) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		final List<? extends Object> inpLabels = getInputLabels();
		final Map<Object, Double> inputs = new HashMap<>(); // label -> index in valLabels
		// Add the input labels and values
		for (int i = 0; i < point.getDimension(); ++i)
			inputs.put(inpLabels.get(i), point.getEntry(i));
		// Start with a identity matrix for the input values
		for (int step = 0; step < mSteps.size(); ++step) {
			// Initialize and call the 'func' associated with this step
			final LabeledMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			final List<? extends Object> fout = func.getOutputLabels();
			if (func instanceof ImplicitMeasurementModel) {
				// Add in output values as constants...
				for (final Object label : fout)
					consts.put(label, inputs.get(label));
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputLabels();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object label = fin.get(i);
				assert label != null : "Missmatch in length between input labels and point dimension.";
				assert inputs.containsKey(label) : "The value for label " + label + " has not been assigned by Step "
						+ step + ".";
				funcPoint.setEntry(i, inputs.get(label));
			}
			final RealVector ovals = func.compute(funcPoint);
			for (int i = 0; i < ovals.getDimension(); ++i)
				inputs.put(fout.get(i), ovals.getEntry(i));
		}
		// Build the outputs in the correct order.
		final List<? extends Object> outs = getOutputLabels();
		final RealVector resVals = new ArrayRealVector(outs.size());
		for (int r = outs.size() - 1; r >= 0; --r)
			resVals.setEntry(r, inputs.get(outs.get(r)));
		return resVals;
	}
}
