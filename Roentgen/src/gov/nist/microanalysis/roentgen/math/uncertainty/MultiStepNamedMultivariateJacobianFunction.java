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
public class MultiStepNamedMultivariateJacobianFunction //
		extends NamedMultivariateJacobianFunction //
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
		return new ArrayList<>(minimumInputRequirements(0, steps));
	}

	private static Set<Object> minimumInputRequirements(final int step,
			final List<NamedMultivariateJacobianFunction> steps) {
		final Set<Object> sob = new HashSet<>();
		final Set<Object> outs = new HashSet<>();
		for (int st = step; st < steps.size(); ++st) {
			final NamedMultivariateJacobianFunction nmjf = steps.get(st);
			for (final Object tag : nmjf.getInputTags())
				if (!outs.contains(tag))
					sob.add(tag);
			outs.addAll(nmjf.getOutputTags());
		}
		return sob;
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
		mFinalOutputs = new HashSet<>(getOutputTags());
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
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		final List<? extends Object> inpTags = getInputTags();
		final Map<Object, Double> inputs = new HashMap<>(); // tag -> index in valTags
		// Add the input tags and values
		for (int i = 0; i < point.getDimension(); ++i)
			inputs.put(inpTags.get(i), point.getEntry(i));
		// Start with a identity matrix for the input values
		RealMatrix expJac = MatrixUtils.createRealIdentityMatrix(point.getDimension());
		List<? extends Object> prevCols = Collections.unmodifiableList(inpTags);
		List<? extends Object> prevRows = Collections.unmodifiableList(inpTags);

		for (int step = 0; step < mSteps.size(); ++step) {
			// Initialize and call the 'func' associated with this step
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			final List<? extends Object> fout = func.getOutputTags();
			if (func instanceof ImplicitMeasurementModel2) {
				// Add in output values as constants...
				for (final Object tag : fout)
					consts.put(tag, inputs.get(tag));
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputTags();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				assert inputs.containsKey(tag) : "The value for tag " + tag + " has not been assigned by Step " + step
						+ ".";
				funcPoint.setEntry(i, inputs.get(tag));
			}
			final Pair<RealVector, RealMatrix> fres = func.evaluate(funcPoint);
			final RealVector ovals = fres.getFirst();
			final RealMatrix ojac = fres.getSecond();
			// (ovals, ojac) = func.evaluate(funcPoint)
			// ovals[fout], ojac[fout][fin]
			// expJac[expRows][expCols]
			for (int i = 0; i < ovals.getDimension(); ++i)
				inputs.put(fout.get(i), ovals.getEntry(i));
			{
				final List<? extends Object> currCols = prevRows; // output becomes input
				// Build currRows from prevRows + newItems
				final List<Object> currRows = new ArrayList<>(prevRows); // old + new output
				final int[] foutIdx = new int[fout.size()];
				if (func instanceof ImplicitMeasurementModel2) {
					for (int i = 0; i < fout.size(); ++i) {
						final int prevIdx = prevCols.indexOf(fout.get(i));
						assert prevIdx != -1;
						foutIdx[i] = prevIdx;
					}
				} else {
					for (int i = 0; i < fout.size(); ++i) {
						final Object tag = fout.get(i);
						assert prevCols.indexOf(tag) == -1;
						currRows.add(tag);
						foutIdx[i] = currRows.size() - 1;
					}
				}
				final RealVector newVals = new ArrayRealVector(currRows.size());
				// Build the extended output
				for (int row = 0; row < currRows.size(); ++row)
					newVals.setEntry(row, inputs.get(currRows.get(row)));
				// Build the extended Jacobian
				final RealMatrix newJac = MatrixUtils.createRealMatrix(currRows.size(), currCols.size());
				for (int rc = 0; rc < currCols.size(); ++rc)
					newJac.setEntry(rc, rc, fout.contains(currCols.get(rc)) ? 0.0 : 1.0);
				for (int c = 0; c < fin.size(); ++c) {
					final int col = currCols.indexOf(fin.get(c));
					if (col == -1)
						continue;
					for (int r = 0; r < foutIdx.length; ++r) {
						final int row = foutIdx[r];
						assert row != -1;
						newJac.setEntry(row, col, ojac.getEntry(r, c));
					}
				}
				expJac = newJac.multiply(expJac);
				prevRows = Collections.unmodifiableList(currRows);
				prevCols = Collections.unmodifiableList(currCols);
			}
		}
		final List<? extends Object> outs = getOutputTags();
		final RealVector resVals = new ArrayRealVector(outs.size());
		for (int r = 0; r < outs.size(); ++r)
			resVals.setEntry(r, inputs.get(outs.get(r)));
		final int outValSize = outs.size();
		final int inpValSize = inpTags.size();
		assert inputs.size() == expJac.getRowDimension();
		assert expJac.getColumnDimension() == inpValSize;
		final RealMatrix resJac = MatrixUtils.createRealMatrix(outValSize, inpValSize);
		for(int row=0;row<outValSize;++row) {
			final int outIdx=prevRows.indexOf(outs.get(row));
			for(int c=0;c<inpValSize;++c)
				resJac.setEntry(row, c, expJac.getEntry(outIdx, c));
		}
		return Pair.create(resVals, resJac);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		assert point.getDimension() > 0 : "Zero dimension input in " + this.toString();
		final List<? extends Object> inpTags = getInputTags();
		final Map<Object, Double> inputs = new HashMap<>(); // tag -> index in valTags
		// Add the input tags and values
		for (int i = 0; i < point.getDimension(); ++i)
			inputs.put(inpTags.get(i), point.getEntry(i));
		// Start with a identity matrix for the input values
		for (int step = 0; step < mSteps.size(); ++step) {
			// Initialize and call the 'func' associated with this step
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			final List<? extends Object> fout = func.getOutputTags();
			if (func instanceof ImplicitMeasurementModel2) {
				// Add in output values as constants...
				for (final Object tag : fout)
					consts.put(tag, inputs.get(tag));
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputTags();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				assert inputs.containsKey(tag) : "The value for tag " + tag + " has not been assigned by Step " + step
						+ ".";
				funcPoint.setEntry(i, inputs.get(tag));
			}
			final RealVector ovals = func.compute(funcPoint);
			for (int i = 0; i < ovals.getDimension(); ++i)
				inputs.put(fout.get(i), ovals.getEntry(i));
		}

		final List<? extends Object> outs = getOutputTags();
		final RealVector resVals = new ArrayRealVector(outs.size());
		for (int r = 0; r < outs.size(); ++r)
			resVals.setEntry(r, inputs.get(outs.get(r)));
		return resVals;
	}
}
