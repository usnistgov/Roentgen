package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;

public class Trimmed<G,H> //
		extends LabeledMultivariateJacobianFunction<G,H> implements ILabeledMultivariateFunction<G,H> {

	final private LabeledMultivariateJacobianFunction<G,H> mBase;
	final private boolean mUntrimmed;

	/**
	 * Ensure the output labels are in the same order as in base
	 *
	 * @param labels
	 * @param func
	 * @return Same items as in labels but in same order as base.getOutputLabels()
	 */
	private static <G, H> List<H> sorted(final Collection<? extends H> labels,
			final LabeledMultivariateJacobianFunction<G, H> base) {
		LabeledMultivariateJacobianFunction<G, H> func;
		if (base instanceof Trimmed)
			func = ((Trimmed<G, H>) base).mBase;
		else
			func = base;
		final List<H> res = new ArrayList<>();
		for (int i = 0; i < func.getOutputDimension(); ++i) {
			final H lbl = func.getOutputLabel(i);
			if (labels.contains(lbl))
				res.add(lbl);
		}
		return res;
	}

	public Trimmed(final LabeledMultivariateJacobianFunction<G,H> base, final List<? extends H> outputLabels)
			throws ArgumentException {
		super(base.getInputLabels(), sorted(outputLabels, base));
		// Apply recursively to SerialLabeledMultivariateJacobianFunction bases
		mBase = base;
		boolean trim = (base.getOutputDimension() == outputLabels.size());
		if (trim) {
			int i = 0;
			for (Object outLbl : outputLabels) {
				if (!base.getOutputLabel(i).equals(outLbl)) {
					trim = false;
					break;
				}
				++i;
			}
		}
		mUntrimmed = !trim;
		for (Object lbl : outputLabels)
			if (mBase.outputIndex(lbl) == -1)
				throw new ArgumentException("The label " + lbl + " is not present in the base LMJF.");
	}

	public LabeledMultivariateJacobianFunction<G,H> getBase() {
		return mBase;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final Pair<RealVector, RealMatrix> baseRM = mBase.value(point);
		if (mUntrimmed)
			return baseRM;
		final RealVector rv1 = baseRM.getFirst();
		final RealMatrix rm1 = baseRM.getSecond();
		final RealVector rv2 = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm2 = NullableRealMatrix.build(getOutputDimension(), getInputDimension());
		for (int r = 0; r < rv2.getDimension(); ++r) {
			final Object rLbl = getOutputLabel(r);
			final int rIdx = mBase.outputIndex(rLbl);
			rv2.setEntry(r, rv1.getEntry(rIdx));
			rm2.setRow(r, rm1.getRow(rIdx));
		}
		return Pair.create(rv2, rm2);
	}

	@SuppressWarnings("unchecked")
	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector rv1;
		if (mBase instanceof ILabeledMultivariateFunction)
			rv1 = ((ILabeledMultivariateFunction<G,H>) mBase).optimized(point);
		else
			rv1 = mBase.value(point).getFirst();
		if (mUntrimmed)
			return rv1;
		final RealVector rv2 = new ArrayRealVector(getOutputDimension());
		for (int r = 0; r < rv2.getDimension(); ++r) {
			final Object rLbl = getOutputLabel(r);
			final int rIdx = mBase.outputIndex(rLbl);
			rv2.setEntry(r, rv1.getEntry(rIdx));
		}
		return rv2;
	}

}