package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioCorrectionModel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioTag;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionTag;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class KRatioCorrectionModelTest {

	@Test
	public void test1() throws ParseException, ArgumentException, IOException {

		final Composition std0 = Composition.massFraction("K411", XPPMatrixCorrectionTest.buildK411());

		final Composition std1 = Composition.parse("Al");

		final Composition unk = Composition.massFraction("K412", XPPMatrixCorrectionTest.buildK412());

		MatrixCorrectionDatum stdK411Mcd = new MatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), //
						Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum stdAlMcd = new MatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), //
						Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		final CharacteristicXRay feTr = CharacteristicXRay.create(Element.Iron, XRayTransition.KA1);
		final CharacteristicXRay mgTr = CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1);
		final CharacteristicXRay caTr = CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1);
		final CharacteristicXRay oTr = CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1);
		final CharacteristicXRay alTr = CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1);
		final CharacteristicXRay siTr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
		{
			final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
			scxr.add(siTr);
			scxr.add(feTr);
			scxr.add(mgTr);
			scxr.add(caTr);
			scxr.add(oTr);

			stds.put(stdK411Mcd, scxr);
			stds.put(stdAlMcd, CharacteristicXRaySet.build(alTr));
		}

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<MatrixCorrectionDatum, CharacteristicXRaySet> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getKey();
			for (final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zaTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, meStd, cxr));
			}
		}

		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		xpp.trimOutputs(outputs);

		KRatioCorrectionModel krcm = new KRatioCorrectionModel(unkMcd, stds);
		Map<Object, Double> constants = new HashMap<>();
		
		constants.put(Composition.buildMassFractionTag(unk, Element.Silicon), 0.211982);
		constants.put(Composition.buildMassFractionTag(unk, Element.Iron), 0.0774196);
		constants.put(Composition.buildMassFractionTag(unk, Element.Magnesium), 0.116567);
		constants.put(Composition.buildMassFractionTag(unk, Element.Calcium), 0.10899);
		constants.put(Composition.buildMassFractionTag(unk, Element.Oxygen), 0.42758);
		constants.put(Composition.buildMassFractionTag(unk, Element.Aluminum), 0.0490615);

		krcm.initializeConstants(constants);

		Report r = new Report("KRationCorrectionModel - test 1");
		try {
			r.add(krcm, Mode.VERBOSE);
			
			Map<Object, UncertainValue> mouv = new HashMap<>();
			
			mouv.put(new KRatioTag(unkMcd, stdK411Mcd, siTr), new UncertainValue(0.794983));
			mouv.put(new KRatioTag(unkMcd, stdK411Mcd, feTr), new UncertainValue(0.687101));
			mouv.put(new KRatioTag(unkMcd, stdK411Mcd, mgTr), new UncertainValue(1.364719));
			mouv.put(new KRatioTag(unkMcd, stdK411Mcd, caTr), new UncertainValue(0.980070));
			mouv.put(new KRatioTag(unkMcd, stdK411Mcd, oTr), new UncertainValue(1.011527));
			mouv.put(new KRatioTag(unkMcd, stdAlMcd, alTr), new UncertainValue(0.032001));

			mouv.put(Composition.buildMassFractionTag(std0, Element.Silicon), new UncertainValue(0.25190067871134, 0.00448737523776));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Iron), new UncertainValue(0.11255374113608, 0.00209872307367));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Magnesium), new UncertainValue(0.09117902759616, 0.0012060717936));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Calcium), new UncertainValue(0.11070559976183, 0.00107203615005));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Oxygen), new UncertainValue(0.42346095279459, 0.00693579374492));
			mouv.put(Composition.buildMassFractionTag(std1, Element.Aluminum), new UncertainValue(1.0,0.001));

			mouv.put(new MatrixCorrectionTag(unkMcd, stdK411Mcd, siTr), new UncertainValue(0.9519));
			mouv.put(new MatrixCorrectionTag(unkMcd, stdK411Mcd, feTr), new UncertainValue(0.9948));
			mouv.put(new MatrixCorrectionTag(unkMcd, stdK411Mcd, mgTr), new UncertainValue(1.0357));
			mouv.put(new MatrixCorrectionTag(unkMcd, stdK411Mcd, caTr), new UncertainValue(0.9942));
			mouv.put(new MatrixCorrectionTag(unkMcd, stdK411Mcd, oTr), new UncertainValue(1.0023));
			mouv.put(new MatrixCorrectionTag(unkMcd, stdAlMcd, alTr), new UncertainValue(0.6523));

			UncertainValues uvs = UncertainValues.extract(krcm.getInputTags(), new UncertainValues(mouv));
			
			NamedMultivariateJacobian nmvj = new NamedMultivariateJacobian(krcm, uvs.getValues());
			r.addHeader("Jacobian");

			r.addHTML(nmvj.toHTML(Mode.NORMAL, new BasicNumberFormat("0.000")));

			r.addHeader("Input");
			r.add(uvs);
			
			r.addHeader("Result");
			r.add(UncertainValues.propagate(nmvj, uvs));
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

}
