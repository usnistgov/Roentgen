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

	private static Set<Object> minimumInputRequirements(int step, final List<NamedMultivariateJacobianFunction> steps) {
		final Set<Object> sob = new HashSet<>();
		final Set<Object> outs = new HashSet<>();
		for (int st = step; st < steps.size(); ++st) {
			final NamedMultivariateJacobianFunction nmjf = steps.get(st);
			for (Object tag : nmjf.getInputTags())
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
				String html = HTML.toHTML(mSteps.get(si), Mode.NORMAL);
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
	
	
	private static RealMatrix expand( //
			List<? extends Object> colTags, // Fixed column tags
			List<? extends Object> inTags, // Input tags to step i+1
			List<? extends Object> outTags, // Output tags to step i+1
			RealMatrix ip1Jac) {
		final List<Object> rOutTags = new ArrayList<>(colTags);
		for (Object outTag : outTags)
			if (!rOutTags.contains(outTag))
				rOutTags.add(outTag);
		final RealMatrix res = MatrixUtils.createRealMatrix(rOutTags.size(), colTags.size());
		for (int rc = 0; rc < colTags.size(); ++rc)
			res.setEntry(rc, rc, outTags.contains(colTags.get(rc)) ? 0.0 : 1.0);
		int[] colIdx = new int[colTags.size()];
		for (int i = 0; i < colIdx.length; ++i) {
			Object colTag = colTags.get(i);
			colIdx[i] = inTags.indexOf(colTag);
		}
		for (int r = 0; r < outTags.size(); ++r) {
			final Object outTag = outTags.get(r);
			final int row = rOutTags.indexOf(outTag);
			for (int col = 0; col < colIdx.length; ++col)
				if (colIdx[col] != -1)
					res.setEntry(row, col, ip1Jac.getEntry(r, colIdx[col]));
		}
		return res;
	}

	private static RealMatrix trim(//
			List<? extends Object> curRowTags, //
			List<? extends Object> resRowTags, //
			RealMatrix input) {
		assert input.getRowDimension() == curRowTags.size();
		RealMatrix res = MatrixUtils.createRealMatrix(resRowTags.size(), input.getColumnDimension());
		for (int r=0;r<res.getRowDimension();++r) {
			final Object rowTag = resRowTags.get(r);
			final int currRow=curRowTags.indexOf(rowTag);
			for(int col=0;col<input.getColumnDimension();++col)
				res.setEntry(r, col, input.getEntry(currRow, col));
		}
		return res;
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
		List<Double> inpVals = new ArrayList<>(); // give & calculated
		List<Object> inpValTags = new ArrayList<>(); // assoc. w. vals
		Map<Object, Integer> inpValMap = new HashMap<>(); // tag -> index in valTags
		// Add the input tags and values
		for (int i = 0; i < point.getDimension(); ++i) {
			inpVals.add(point.getEntry(i));
			inpValTags.add(inpTags.get(i));
			inpValMap.put(inpTags.get(i), i);
		}
		// Start with a identity matrix for the input values
		RealMatrix expJac = MatrixUtils.createRealIdentityMatrix(point.getDimension());
		for (int step = 0; step < mSteps.size(); ++step) {
			// Initialize and call the 'func' associated with this step
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants based on this->getConstants()
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			if (func instanceof ImplicitMeasurementModel2) {
				// Add in output values as constants...
				List<? extends Object> fOut = func.getOutputTags();
				for (Object tag : fOut) {
					consts.put(tag, inpVals.get(inpValMap.get(tag)));
					// It's a constant so remove it from the inputs
					inpVals.remove(inpValMap.get(tag).intValue());
					inpValTags.remove(tag);
					inpValMap.remove(tag);
				}
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call evaluate
			final RealVector funcPoint = new ArrayRealVector(func.getInputDimension());
			final List<? extends Object> fin = func.getInputTags();
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				final Integer idx = inpValMap.get(tag);
				assert idx != null : "The value for tag " + tag + " has not been assigned by Step " + step + ".";
				funcPoint.setEntry(i, inpVals.get(idx));
			}
			final Pair<RealVector, RealMatrix> fres = func.evaluate(funcPoint);
			final RealVector ovals = fres.getFirst();
			final RealMatrix ojac = fres.getSecond();
			final int valsSize = inpVals.size();
			final List<? extends Object> fout = func.getOutputTags();
			if (mFinalOutputs != null) {
				// Build up a list of computed inputs required by subsequent steps
				// and final output
				final Set<Object> reqOutputs = minimumInputRequirements(step + 1, mSteps);
				reqOutputs.addAll(mFinalOutputs);
				reqOutputs.removeAll(inpTags);
				// List of computed items available (prev+this-input)
				final Set<Object> avail = new HashSet<>(inpValTags);
				avail.addAll(fout);
				avail.removeAll(inpTags);
				// Create lists of the tags and values to carry on to the next step
				final List<Object> outTagsExt = new ArrayList<>(inpTags);
				outTagsExt.addAll(Sets.intersection(avail, reqOutputs));
				final List<Double> outValsExt = new ArrayList<>();
				final Map<Object, Integer> outValsExtMap = new HashMap<>();
				for (int j = 0; j < outTagsExt.size(); ++j) {
					final Object outTag = outTagsExt.get(j);
					final Integer idx = inpValMap.get(outTag);
					final Double val = idx != null ? inpVals.get(idx) : ovals.getEntry(func.outputIndex(outTag));
					outValsExt.add(val);
					outValsExtMap.put(outTag, j);
				}
				// Fill jacobian with combination of old and new
				final RealMatrix ojacx = MatrixUtils.createRealMatrix(outTagsExt.size(), inpValTags.size());
				final int[] foutIdxs = new int[outTagsExt.size()];
				for (int outIdx = 0; outIdx < outTagsExt.size(); ++outIdx)
					foutIdxs[outIdx] = func.outputIndex(outTagsExt.get(outIdx)); // May be -1
				for (int ojacxCol = 0; ojacxCol < inpValTags.size(); ++ojacxCol) {
					final int ojacCol = func.inputIndex(inpValTags.get(ojacxCol)); // May be -1??
					for (int ojacxRow = 0; ojacxRow < outTagsExt.size(); ++ojacxRow) {
						final int ojacRow = foutIdxs[ojacxRow]; // May be -1
						if ((ojacCol != -1) && (ojacRow != -1))
							ojacx.setEntry(ojacxRow, ojacxCol, ojac.getEntry(ojacRow, ojacCol));
						else {
							// the output value was carried over from one of the the inputs
							if (inpValTags.get(ojacxCol).equals(outTagsExt.get(ojacxRow)))
								ojacx.setEntry(ojacxRow, ojacxCol, 1.0); // dp/dp=1
							else
								ojacx.setEntry(ojacxRow, ojacxCol, 0.0); // dp/dq=0
						}
					}
				}
				// Outputs of step become inputs of step+1
				inpVals = outValsExt;
				inpValTags = outTagsExt;
				inpValMap = outValsExtMap;
				// Multiply this Jacobian times the inner combined Jacobian
				expJac = ojacx.multiply(expJac);
			} else {
				if (func instanceof ImplicitMeasurementModel2) {
					// assert false: "Does not work correctly...";
					// IMPLICIT MEASUREMENT MODEL
					List<Object> newTags = new ArrayList<>();
					for (int i = 0; i < ovals.getDimension(); ++i) {
						final Object fOutTag = fout.get(i);
						final int fOutIdx = inpValTags.indexOf(fOutTag);
						if (fOutIdx != -1) {
							assert func instanceof ImplicitMeasurementModel2;
							assert inpVals.get(fOutIdx) == ovals.getEntry(i) : //
							inpVals.get(fOutIdx) + " != " + ovals.getEntry(i) + " for " + fOutTag;
							inpVals.set(fOutIdx, ovals.getEntry(i));
						} else {
							newTags.add(fOutTag);
							inpValTags.add(fOutTag);
							inpVals.add(ovals.getEntry(i));
							inpValMap.put(fOutTag, inpValTags.size() - 1);
						}
					}
					assert inpValTags.size() == inpVals.size();
					assert inpValMap.size() == inpVals.size();
					final int rowDim = inpVals.size();
					final int colDim = expJac.getRowDimension();
					assert colDim <= rowDim;
					final RealMatrix ojacx = MatrixUtils.createRealMatrix(rowDim, colDim);
					// Fill the diagonal with 1.0 unless part of the implicit model
					for (int rc = 0; rc < colDim; ++rc) {
						final Object inpTag = inpTags.get(rc);
						ojacx.setEntry(rc, rc, func.outputIndex(inpTag) == -1 ? 1.0 : 0.0);
					}
					for (int rOut = 0; rOut < fout.size(); ++rOut) {
						final int rowX = inpTags.indexOf(fout.get(rOut));
						if (rowX == -1)
							continue;
						for (int cOut = 0; cOut < fin.size(); ++cOut) {
							final int colX = inpTags.indexOf(fin.get(cOut));
							if (colX != -1)
								ojacx.setEntry(rowX, colX, ojac.getEntry(rOut, cOut));
						}
					}
					expJac = ojacx.multiply(expJac);
				} else {
					// EXPLICIT MEASUREMENT MODEL
					// Add output values to the end of vals and valTags
					for (int i = 0; i < ovals.getDimension(); ++i) {
						final Object fOutTag = fout.get(i);
						assert !inpValTags.contains(fOutTag) : "Duplicate value defined for " + fOutTag;
						inpValTags.add(fOutTag);
						inpVals.add(ovals.getEntry(i));
						inpValMap.put(fOutTag, inpValTags.size() - 1);
					}
					assert inpValTags.size() == inpVals.size();
					assert inpValMap.size() == inpVals.size();
					// Build the expanded jacobian matrix
					assert valsSize == expJac.getRowDimension();
					final int outSize = ojac.getRowDimension();
					final int rowDim = outSize + valsSize;
					final int colDim = valsSize;
					// Build the expanded output Jacobian
					final RealMatrix ojacx = MatrixUtils.createRealMatrix(rowDim, colDim);
					// Prepend an identity matrix to "carry" previous inputs
					for (int i = 0; i < valsSize; ++i)
						ojacx.setEntry(i, i, 1.0);
					// Append the output Jacobian
					for (int i = 0; i < fin.size(); ++i) {
						if (!inpValMap.containsKey(fin.get(i)))
							assert inpValMap.containsKey(fin.get(i)) : fin.get(i) + " is missing...";
						final int col = inpValMap.get(fin.get(i));
						for (int r = 0; r < outSize; ++r)
							ojacx.setEntry(r + valsSize, col, ojac.getEntry(r, i));
					}
					// Multiply this Jacobian times the inner combined Jacobian
					expJac = ojacx.multiply(expJac);
				}
			}
		}
		final List<? extends Object> outs = getOutputTags();
		final int outValSize = outs.size();
		final int inpValSize = inpTags.size();
		assert expJac.getRowDimension() == outValSize + inpValSize : expJac.getRowDimension() + "!=" + outValSize + "+"
				+ inpValSize;
		assert inpVals.size() == expJac.getRowDimension();
		assert expJac.getColumnDimension() == inpValSize;
		if (mFinalOutputs == null) {
			final RealVector resVals = new ArrayRealVector(outValSize);
			for (int r = 0; r < outValSize; ++r) {
				assert outs.get(r) == inpValTags.get(r + inpValSize) : toString();
				resVals.setEntry(r, inpVals.get(r + inpValSize));
			}
			final RealMatrix resJac = expJac.getSubMatrix(inpValSize, inpValSize + outValSize - 1, 0, inpValSize - 1);
			return Pair.create(resVals, resJac);
		} else {
			final RealVector resVals = new ArrayRealVector(outValSize);
			final RealMatrix resJac = MatrixUtils.createRealMatrix(outValSize, inpValSize);
			for (int row = 0; row < outValSize; ++row) {
				final int valIdx = inpValMap.get(outs.get(row));
				resVals.setEntry(row, inpVals.get(valIdx));
				for (int col = 0; col < inpValSize; ++col)
					resJac.setEntry(row, col, expJac.getEntry(valIdx, col));
			}
			return Pair.create(resVals, resJac);
		}
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final List<? extends Object> inpTags = getInputTags();
		// Keep track of the inputs available to the next step
		List<Double> inpVals = new ArrayList<>(); // give & calculated
		List<Object> inpValTags = new ArrayList<>(); // assoc. w. vals
		Map<Object, Integer> inpValMap = new HashMap<>();
		// Add the initial input tags and values
		for (int i = 0; i < point.getDimension(); ++i) {
			inpVals.add(point.getEntry(i));
			inpValTags.add(inpTags.get(i));
			inpValMap.put(inpTags.get(i), i);
		}
		for (int step = 0; step < mSteps.size(); ++step) {
			final NamedMultivariateJacobianFunction func = mSteps.get(step);
			// Initialize constants
			final Map<Object, Double> consts = new HashMap<>(getConstants());
			if (func instanceof ImplicitMeasurementModel2) {
				// Specify the output values as constants since they must be provided.
				final List<? extends Object> fOut = func.getOutputTags();
				for (Object tag : fOut) {
					consts.put(tag, inpVals.get(inpValMap.get(tag)));
					// It's a constant so remove it from the inputs
					inpVals.remove(inpValMap.get(tag).intValue());
					inpValTags.remove(tag);
					inpValMap.remove(tag);
				}
			}
			func.initializeConstants(consts);
			// Build the vector argument to func and call
			final List<? extends Object> fin = func.getInputTags();
			final RealVector funcPoint = new ArrayRealVector(fin.size());
			for (int i = 0; i < funcPoint.getDimension(); ++i) {
				final Object tag = fin.get(i);
				assert tag != null : "Missmatch in length between input tags and point dimension.";
				final Integer idx = inpValMap.get(tag);
				assert idx != null : "The value for tag " + tag + " has not been assigned by Step " + step + ".";
				funcPoint.setEntry(i, inpVals.get(idx));
			}
			final RealVector ovals = func.compute(funcPoint);
			final List<? extends Object> fout = func.getOutputTags();
			if (mFinalOutputs != null) {
				// Build up a list of computed inputs required by subsequent steps
				// and final output
				final Set<Object> reqOutputs = new HashSet<>();
				for (int subStep = step + 1; subStep < mSteps.size(); ++subStep)
					reqOutputs.addAll(mSteps.get(subStep).getInputTags());
				reqOutputs.addAll(mFinalOutputs);
				reqOutputs.removeAll(inpTags);
				// List of computed items available (prev+this-input)
				final Set<Object> avail = new HashSet<>(inpValTags);
				avail.addAll(fout);
				avail.removeAll(inpTags);
				// Create lists of the tags and values to carry on to the next step
				final List<Object> outTagsExt = new ArrayList<>(inpTags);
				outTagsExt.addAll(Sets.intersection(avail, reqOutputs));
				final List<Double> outValsExt = new ArrayList<>();
				final Map<Object, Integer> outValsExtMap = new HashMap<>();
				for (int j = 0; j < outTagsExt.size(); ++j) {
					final Object outTag = outTagsExt.get(j);
					final Integer idx = inpValMap.get(outTag);
					final Double val = idx != null ? inpVals.get(idx) : ovals.getEntry(fout.indexOf(outTag));
					outValsExt.add(val);
					outValsExtMap.put(outTag, Integer.valueOf(j));
				}
				// Outputs of step become inputs of step+1
				inpVals = outValsExt;
				inpValTags = outTagsExt;
				inpValMap = outValsExtMap;
			} else {
				// Add output values to the end of vals and valTags
				for (int i = 0; i < ovals.getDimension(); ++i) {
					final Object fOutTag = fout.get(i);
					final int fOutIdx = inpValTags.indexOf(fOutTag);
					if (fOutIdx != -1) {
						assert func instanceof ImplicitMeasurementModel2;
						assert inpVals.get(fOutIdx) == ovals.getEntry(i) : //
						inpVals.get(fOutIdx) + " != " + ovals.getEntry(i) + " for " + fOutTag;
						inpVals.set(fOutIdx, ovals.getEntry(i));
					} else {
						inpValTags.add(fOutTag);
						inpVals.add(ovals.getEntry(i));
						inpValMap.put(fOutTag, inpValTags.size() - 1);
					}
				}
				assert inpValTags.size() == inpVals.size();
				assert inpValMap.size() == inpVals.size();
			}
		}
		// The final inpVals / inpValTags represent the outputs
		final List<? extends Object> outs = getOutputTags();
		final int outValSize = outs.size();
		final int inpValSize = inpTags.size();
		final RealVector resVals = new ArrayRealVector(outValSize);
		if (mFinalOutputs == null) {
			// Skip the values associated with the initial inputs.
			for (int r = 0; r < outValSize; ++r) {
				assert outs.get(r) == inpValTags.get(r + inpValSize) : toString();
				resVals.setEntry(r, inpVals.get(r + inpValSize));
			}
		} else {
			// Output only the requested mFinalOutput values
			for (int row = 0; row < outValSize; ++row) {
				final int valIdx = inpValMap.get(outs.get(row));
				resVals.setEntry(row, inpVals.get(valIdx));
			}
		}
		return resVals;
	}
}
