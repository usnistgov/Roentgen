package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class KRatioCorrectionModel2 //
		extends ImplicitMeasurementModel {

	static List<? extends Object> buildOutputs(//
			Set<KRatioLabel> krs //
	) {
		assert KRatioLabel.areAllSameUnknownElements(krs);
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		Set<Element> elms = unk.getElementSet();
		final List<Object> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(Composition.buildMassFractionTag(unk.getComposition(), elm));
		return res;
	}

	public KRatioCorrectionModel2(//
			Set<KRatioLabel> krs //
	) throws ArgumentException {
		super(new KRatioHModel2(krs), buildOutputs(krs));
	}

	static public Pair<LabeledMultivariateJacobianFunction,UncertainValues> buildXPPModel( //
			Set<KRatioLabel> krs,
			Set<MatrixCorrectionModel2.Variates> variates //
	) throws ArgumentException {
		final List<LabeledMultivariateJacobianFunction> steps = new ArrayList<>();
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(krs, variates);
		steps.add(xpp);
		
		steps.add(new KRatioCorrectionModel2(krs));
		return Pair.create(new SerialLabeledMultivariateJacobianFunction("Full K-ratio correction", steps),
				xpp.buildInput());
	}

	public UncertainValues getInputs(final MatrixCorrectionModel mcm, final UncertainValues inputs,
			final UncertainValues krs) //
			throws ArgumentException {
		return UncertainValues.extract(getInputLabels(), UncertainValues.propagate(mcm, inputs), krs, inputs);
	}

	@Override
	public String toString() {
		return "K-implicit model";
	}

}
