package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
abstract public class MatrixCorrectionModel2 //
		extends SerialLabeledMultivariateJacobianFunction {

	protected final Set<KRatioLabel> mKRatios;

	public MatrixCorrectionModel2(//
			String name, //
			Set<KRatioLabel> kratios, //
			List<LabeledMultivariateJacobianFunction> steps //
	) throws ArgumentException {
		super(name, steps);
		mKRatios = kratios;
		// Validate the inputs...
		final List<? extends Object> outputTags = getOutputLabels();
		final List<? extends Object> inputTags = getInputLabels();
		for(KRatioLabel krl : mKRatios){
			MatrixCorrectionLabel mct = new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet());
			if (!outputTags.contains(mct))
				throw new ArgumentException(toString() + " does not calculate the required output " + mct.toString());
			final Composition stdComp = krl.getStandard().getComposition();
			for (Element elm : stdComp.getElementSet()) {
				MassFractionTag mft = Composition.buildMassFractionTag(stdComp, elm);
				if (!inputTags.contains(mft))
					throw new ArgumentException(toString() + " must take " + mft.toString() + " as an argument.");
			}
			final Composition unkComp = krl.getUnknown().getComposition();
			for (Element elm : unkComp.getElementSet()) {
				MassFractionTag mft = Composition.buildMassFractionTag(unkComp, elm);
				if (!inputTags.contains(mft))
					throw new ArgumentException(toString() + " must take " + mft.toString() + " as an argument.");
			}
		}
	}

	public Set<KRatioLabel> getKRatios() {
		return Collections.unmodifiableSet(mKRatios);
	}

	public Set<KRatioLabel> getKRatios(Element elm) {
		Set<KRatioLabel> res = new HashSet<>();
		for(KRatioLabel krl : mKRatios)
			if(krl.getXRaySet().getElement().equals(elm))
				res.add(krl);
		return Collections.unmodifiableSet(res);
	}
	
	public Set<Element> getElementSet() {
		Set<Element> res = new HashSet<>();
		for(KRatioLabel krl : mKRatios)
			res.add(krl.getXRaySet().getElement());
		return Collections.unmodifiableSet(res);
	}
	

	
	
}
