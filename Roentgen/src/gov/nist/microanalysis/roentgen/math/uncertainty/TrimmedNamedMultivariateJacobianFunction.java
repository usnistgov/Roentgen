package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * This class serves to reduce the number of input and output variables when the
 * Jacobian associated with a fraction of the input variables and output
 * variables is required.
 *
 * @author Nicholas
 *
 */
public class TrimmedNamedMultivariateJacobianFunction
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	private final LabeledMultivariateJacobianFunction mBase;
	private RealVector mBaseInputs;
	

	/**
	 * @param base         The base {@link LabeledMultivariateJacobianFunction}
	 * @param inputLabels  The input labels for the
	 *                     {@link TrimmedNamedMultivariateJacobianFunction}
	 * @param outputLabels The output labels for the
	 *                     {@link TrimmedNamedMultivariateJacobianFunction}
	 * @throws ArgumentException
	 */
	public TrimmedNamedMultivariateJacobianFunction(//
			final LabeledMultivariateJacobianFunction base, //
			final List<? extends Object> inputLabels, //
			final List<? extends Object> outputLabels //
	) throws ArgumentException {
		super(inputLabels, outputLabels);
		assert base.getInputLabels().containsAll(inputLabels);
		assert base.getOutputLabels().containsAll(outputLabels);
		if (!base.getInputLabels().containsAll(inputLabels))
			throw new ArgumentException("Some of the requested input labels are not in the base function's inputs.");
		if (!base.getOutputLabels().containsAll(outputLabels))
			throw new ArgumentException("Some of the requested output labels are not in the base function's outputs.");
		mBase = base;
	}
	
	/**
	 * The set of inputs from which the full set of inputs required
	 * by the base {@link LabeledMultivariateJacobianFunction} is
	 * constructed.  
	 * 
	 * @param point
	 */
	public void setBaseInputs(RealVector point) {
		assert point.getDimension()==mBase.getInputDimension();
		mBaseInputs = point.copy();
	}
	
	public void setBaseInputs(UncertainValuesBase uvb) {
		RealVector basePt = new ArrayRealVector(mBase.getInputDimension());
		for(int i=0;i<mBase.getInputDimension();++i) {
			Object lbk = mBase.getInputLabel(i);
			basePt.setEntry(i, uvb.getEntry(lbk));
		}
		setBaseInputs(basePt);
	}
	
	private final RealVector buildBaseInput(final RealVector point) {
		RealVector basePt = mBaseInputs.copy();
		assert point.getDimension()==getInputDimension();
		for(int i=0;i<getInputDimension();++i) {
			Object lbl = getInputLabel(i);
			basePt.setEntry(mBase.inputIndex(lbl), point.getEntry(i));
		}
		return basePt;	
	}
	

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector basePoint = buildBaseInput(point);
		final Pair<RealVector, RealMatrix> tmp = mBase.evaluate(basePoint);
		final RealVector rv1 = tmp.getFirst();
		final RealMatrix rm1 = tmp.getSecond();
		final List<? extends Object> outLabels = getOutputLabels();
		final int[] rIdx = new int[outLabels.size()];
		for (int i = 0; i < rIdx.length; ++i)
			rIdx[i] = outputIndex(outLabels.get(i));
		final List<? extends Object> inLabels = getInputLabels();
		final int[] cIdx = new int[inLabels.size()];
		for (int i = 0; i < cIdx.length; ++i)
			cIdx[i] = inputIndex(inLabels.get(i));
		final RealVector rv2 = new ArrayRealVector(rIdx.length);
		final RealMatrix rm2 = MatrixUtils.createRealMatrix(rIdx.length, cIdx.length);
		for (int r = 0; r < rIdx.length; ++r) {
			rv2.setEntry(r, rv1.getEntry(rIdx[r]));
			for (int c = 0; c < cIdx.length; ++c)
				rm2.setEntry(r, c, rm1.getEntry(rIdx[r], cIdx[c]));
		}
		return Pair.create(rv2, rm2);
	}

	@Override
	public RealVector optimized(RealVector point) {
		final RealVector basePoint = buildBaseInput(point);
		final RealVector rv1 = mBase.compute(basePoint);
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (int r = 0; r < rv.getDimension(); ++r)
			rv.setEntry(r, rv1.getEntry(outputIndex(getOutputLabel(r))));
		return rv;
	}

}
