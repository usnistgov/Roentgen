package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.IntermediateMassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

public class ElementByDifference //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	private final Element mElementByDifference;

	public ElementByDifference(final Material mat, final Element elmByDiff) {
		super(buildInputs(mat, elmByDiff), MaterialLabel.buildMassFractionTags(mat));
		mElementByDifference = elmByDiff;
	}

	private static List<? extends Object> buildInputs(final Material mat, final Element elmByDif) {
		final List<Object> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			if (!elm.equals(elmByDif))
				res.add(MaterialLabel.buildIntermediateMFTag(mat, elm));
		return res;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		double sum = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			sum += point.getEntry(i);
		for (int oIdx = 0; oIdx < getOutputDimension(); ++oIdx) {
			final MassFraction mfl = (MassFraction) getOutputLabel(oIdx);
			if (mfl.getElement().equals(mElementByDifference)) {
				for (int iIdx = 0; iIdx < getInputDimension(); ++iIdx)
					writeJacobian(oIdx, getInputLabel(iIdx), -1.0, rm);
				rv.setEntry(oIdx, 1.0 - sum);

			} else {
				final IntermediateMassFraction ebdl = MaterialLabel.buildIntermediateMFTag(mfl.getMaterial(),
						mfl.getElement());
				final double cx = getValue(ebdl, point);
				for (int iIdx = 0; iIdx < getInputDimension(); ++iIdx) {
					final IntermediateMassFraction inputLabel = (IntermediateMassFraction) getInputLabel(iIdx);
					if (inputLabel.getElement().equals(mfl.getElement()))
						writeJacobian(oIdx, getInputLabel(iIdx), 1.0, rm);
					else
						writeJacobian(oIdx, getInputLabel(iIdx), 0.0, rm);
				}
				rv.setEntry(oIdx, cx);
			}
		}
		return Pair.create(rv, rm);

	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		double sum = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			sum += Math.max(0, point.getEntry(i));
		for (int oIdx = 0; oIdx < getOutputDimension(); ++oIdx) {
			final MassFraction mfl = (MassFraction) getOutputLabel(oIdx);
			if (mfl.getElement().equals(mElementByDifference))
				rv.setEntry(oIdx, 1.0 - sum);
			else {
				final IntermediateMassFraction label = MaterialLabel.buildIntermediateMFTag(mfl.getMaterial(),
						mfl.getElement());
				rv.setEntry(oIdx, getValue(label, point));
			}
		}
		return rv;
	}

}
