package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.List;
import java.util.Map;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;

/**
 * A matrix correction model is a model that computes the values associated with
 * {@link MatrixCorrectionLabel} for a set of {@link ElementXRaySet} objects
 * associated with {@link MatrixCorrectionDatum}s associated with Standards
 * relative to a {@link MatrixCorrectionDatum} associated with an unknown.
 * 
 * A matrix correction model may also calculate values for {@link KRatioLabel}
 * and other related tags.
 * 
 * @author nicholas
 *
 */
abstract public class MatrixCorrectionModel //
		extends SerialLabeledMultivariateJacobianFunction {

	protected final MatrixCorrectionDatum mUnknown;
	protected final Map<ElementXRaySet, MatrixCorrectionDatum> mStandards;

	public MatrixCorrectionModel(//
			String name, //
			MatrixCorrectionDatum unk, //
			Map<ElementXRaySet, MatrixCorrectionDatum> stds, //
			List<LabeledMultivariateJacobianFunction> steps //
	) throws ArgumentException {
		super(name, steps);
		mUnknown = unk;
		mStandards = stds;
		// Validate the inputs...
		final List<? extends Object> outputTags = getOutputLabels();
		final List<? extends Object> inputTags = getInputLabels();
		for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			MatrixCorrectionLabel mct = new MatrixCorrectionLabel(unk, me.getValue(), me.getKey());
			if (!outputTags.contains(mct))
				throw new ArgumentException(toString() + " does not calculate the required output " + mct.toString());
		}
		for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum std = me.getValue();
			MatrixCorrectionLabel mct = new MatrixCorrectionLabel(unk, std, me.getKey());
			if (!outputTags.contains(mct))
				throw new ArgumentException(toString() + " does not calculate the required output " + mct.toString());
			for (Element elm : std.getComposition().getElementSet()) {
				MassFractionTag mft = Composition.buildMassFractionTag(std.getComposition(), elm);
				if (!inputTags.contains(mft))
					throw new ArgumentException(toString() + " must take " + mft.toString() + " as an argument.");
			}
		}
		for (Element elm : unk.getComposition().getElementSet()) {
			MassFractionTag mft = Composition.buildMassFractionTag(unk.getComposition(), elm);
			if (!inputTags.contains(mft))
				throw new ArgumentException(toString() + " must take " + mft.toString() + " as an argument.");
		}

	}

	public MatrixCorrectionDatum getUnknown() {
		return mUnknown;
	}

	public Map<ElementXRaySet, MatrixCorrectionDatum> getStandards() {
		return mStandards;
	}

}
