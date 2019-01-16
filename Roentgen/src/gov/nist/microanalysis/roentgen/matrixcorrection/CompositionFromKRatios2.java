package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;

public class CompositionFromKRatios2 //
		extends SerialLabeledMultivariateJacobianFunction {

	private final XPPMatrixCorrection2 mXPP;

	private static List<LabeledMultivariateJacobianFunction> buildSteps(final Set<KRatioLabel> kratios,
			final Set<MatrixCorrectionModel2.Variate> variates) //
			throws ArgumentException {
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(kratios, variates);
		final KRatioHModel2 hModel = new KRatioHModel2(kratios);

		// Build the full model with xpp and hModel
		final List<LabeledMultivariateJacobianFunction> steps = new ArrayList<>();
		steps.add(xpp);
		steps.add(hModel);
		return steps;
	}

	public CompositionFromKRatios2(final Set<KRatioLabel> kratios, final Set<MatrixCorrectionModel2.Variate> variates) //
			throws ArgumentException {
		super("Composition From K-Ratios", buildSteps(kratios, variates));
		mXPP = (XPPMatrixCorrection2) getStep(0);
	}

	public UncertainValues buildInputs(final UncertainValues kratios) //
			throws ArgumentException {
		return UncertainValues.build(getInputLabels(), mXPP.buildInput(), kratios);
	}
}
