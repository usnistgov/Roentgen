package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * @author nicholas
 *
 */
public class KRatioCorrectionModel
		extends ImplicitMeasurementModel {

	static List<? extends Object> buildOutputs(//
			Composition unk //
	) {
		List<Object> res = new ArrayList<>();
		for (Element elm : unk.getElementSet())
			res.add(Composition.buildMassFractionTag(unk, elm));
		return res;
	}

	public KRatioCorrectionModel(//
			final MatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, MatrixCorrectionDatum> stds //
	) {
		super(new KRatioHModel(unk, stds), buildOutputs(unk.getComposition()));
	}
	
	public String toString() {
		return "K-implicit model";
	}
	
}
