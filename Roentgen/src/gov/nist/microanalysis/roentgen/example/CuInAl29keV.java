package gov.nist.microanalysis.roentgen.example;

import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.duckandcover.html.Report;
import com.duckandcover.html.IToHTML.Mode;
import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class CuInAl29keV {

	public static void main(String[] args) throws ArgumentException, IOException {
		Map<Element,Double> men = new HashMap<>();
		men.put(Element.Aluminum, 1.0-1e-8);
		men.put(Element.Copper,   1e-8); // Must be a little bit of Cu 
		Composition unkComp = Composition.massFraction("Al", men);
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(), UncertainValue.valueOf(29.0), UncertainValue.toRadians(40.0, 0.0));
		StandardMatrixCorrectionDatum std = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Copper), UncertainValue.valueOf(29.0), UncertainValue.toRadians(40.0, 0.0));
		KRatioLabel krl = new KRatioLabel(unk, std, CharacteristicXRay.create(Element.Copper, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		Set<KRatioLabel> kratios = Collections.singleton(krl);
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, true);
		Report r=new Report("Cu K-L3 in Al at 29.0 keV");
		r.addHeader("Cu K-L3 in Al at 29.0 keV");
		// r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.inBrowser(Mode.VERBOSE);
	}

}
