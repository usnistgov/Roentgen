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
public class TrimmedNamedMultivariateJacobianFunction //
		extends NamedMultivariateJacobianFunction {

	private final NamedMultivariateJacobianFunction mBase;

	public TrimmedNamedMultivariateJacobianFunction(//
			NamedMultivariateJacobianFunction base, //
			List<? extends Object> inputs, //
			List<? extends Object> outputs //
	) throws ArgumentException {
		super(inputs, outputs);
		assert base.getInputTags().containsAll(inputs);
		assert base.getOutputTags().containsAll(outputs);
		if (!base.getInputTags().containsAll(inputs))
			throw new ArgumentException("Some of the requested input tags are not in the base function's inputs.");
		if (base.getOutputTags().containsAll(outputs))
			throw new ArgumentException("Some of the requested output tags are not in the base function's outputs.");
		mBase = base;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		List<? extends Object> baseInputs = mBase.getInputTags();
		RealVector basePoint = new ArrayRealVector(baseInputs.size());
		for (int i = 0; i < baseInputs.size(); ++i) {
			int idx = inputIndex(baseInputs.get(i));
			assert !((idx == -1) && (!getConstants().containsKey(baseInputs.get(i))));
			final double value = idx != -1 ? point.getEntry(idx) : getConstant(baseInputs.get(i));
			basePoint.setEntry(i, value);
		}
		Pair<RealVector, RealMatrix> tmp = mBase.value(point);
		List<? extends Object> outTags = getOutputTags();
		int[] outIdx = new int[outTags.size()];
		for (int i = 0; i < outIdx.length; ++i)
			outIdx[i] = outputIndex(outTags.get(i));
		List<? extends Object> inTags = getInputTags();
		int[] inIdx = new int[inTags.size()];
		for (int i = 0; i < inIdx.length; ++i)
			inIdx[i] = inputIndex(inTags.get(i));
		RealVector rv = new ArrayRealVector(outIdx.length);
		RealMatrix rm = MatrixUtils.createRealMatrix(outIdx.length, inIdx.length);
		for (int r = 0; r < outIdx.length; ++r) {
			rv.setEntry(r, tmp.getFirst().getEntry(outIdx[r]));
			for (int c = 0; c < inIdx.length; ++c)
				rm.setEntry(r, c, tmp.getSecond().getEntry(outIdx[r], inIdx[c]));
		}
		return Pair.create(rv, rm);
	}

}
