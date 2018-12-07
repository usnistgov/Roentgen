package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class KRatioCorrectionModel extends ImplicitMeasurementModel {

	static List<? extends Object> buildOutputs(//
			Composition unk //
	) {
		List<Object> res = new ArrayList<>();
		for (Element elm : unk.getElementSet())
			res.add(Composition.buildMassFractionTag(unk, elm));
		return res;
	}

	public KRatioCorrectionModel(//
			final UnknownMatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds //
	) {
		super(new KRatioHModel(unk, stds), buildOutputs(unk.getComposition()));
	}
	
	
	static public LabeledMultivariateJacobianFunction buildXPPModel( //
			final UnknownMatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds //
	) throws ArgumentException {
		List<LabeledMultivariateJacobianFunction> steps = new ArrayList<>();
		steps.add(new XPPMatrixCorrection(unk, stds));
		steps.add(new KRatioCorrectionModel(unk, stds));
		return new SerialLabeledMultivariateJacobianFunction("Full K-ratio correction", steps);
	}
		

	public UncertainValues getInputs(MatrixCorrectionModel mcm, UncertainValues inputs, UncertainValues krs) //
			throws ArgumentException {
		return UncertainValues.extract(getInputLabels(), UncertainValues.propagate(mcm, inputs), krs, inputs);
	}

	public String toString() {
		return "K-implicit model";
	}

}
