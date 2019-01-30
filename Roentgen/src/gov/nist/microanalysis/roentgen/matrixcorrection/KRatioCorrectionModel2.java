package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class KRatioCorrectionModel2 //
		extends ImplicitMeasurementModel {

	static List<? extends Object> buildOutputs(//
			final Set<KRatioLabel> krs //
	) {
		assert KRatioLabel.areAllSameUnknownElements(krs);
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		return unk.getComposition().massFractionTags();
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs //
	) throws ArgumentException {
		super(new KRatioHModel2(krs), buildOutputs(krs));
	}

	static public Pair<LabeledMultivariateJacobianFunction, UncertainValues> buildXPPModel( //
			final UncertainValues krs, final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		Set<KRatioLabel> keySet = new HashSet<>(krs.extractTypeOfLabel(KRatioLabel.class));
		final List<LabeledMultivariateJacobianFunction> steps = new ArrayList<>();
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(keySet, variates);
		steps.add(xpp);
		steps.add(new KRatioCorrectionModel2(keySet));
		final UncertainValues input = UncertainValues.combine(xpp.buildInput(), krs);
		return Pair.create( //
				new SerialLabeledMultivariateJacobianFunction("Full K-ratio correction", steps), //
				input //
		);
	}

	public UncertainValues getInputs(final MatrixCorrectionModel2 mcm, final UncertainValues inputs,
			final UncertainValues krs) //
			throws ArgumentException {
		return UncertainValues.extract(getInputLabels(), UncertainValues.propagate(mcm, inputs), krs, inputs);
	}

	@Override
	public String toString() {
		return "K-implicit model";
	}

}
