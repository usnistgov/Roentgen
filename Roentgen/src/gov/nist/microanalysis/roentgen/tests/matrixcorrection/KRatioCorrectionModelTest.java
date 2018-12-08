package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.CompositionFromKRatios2;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class KRatioCorrectionModelTest {

	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);

	@Test
	public void test1() throws ParseException, ArgumentException, IOException {

		final Composition std0 = Composition.massFraction("K411", XPPMatrixCorrectionTest.buildK411());

		final Composition std1 = Composition.parse("Al");

		final Composition unk = Composition.massFraction("K412", XPPMatrixCorrectionTest.buildK412());

		final StandardMatrixCorrectionDatum stdK411Mcd = new StandardMatrixCorrectionDatum( //
				std0, new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), //
						Math.toRadians(0.9)) //
		);

		final StandardMatrixCorrectionDatum stdAlMcd = new StandardMatrixCorrectionDatum( //
				std1, new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), //
						Math.toRadians(0.7)) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk, new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final ElementXRaySet feTr = ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1);
		final ElementXRaySet mgTr = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet caTr = ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1);
		final ElementXRaySet oTr = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet alTr = ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1);
		final ElementXRaySet siTr = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);

		final Set<KRatioLabel> kratios = new HashSet<>();
		kratios.add(new KRatioLabel(unkMcd, stdK411Mcd, siTr, Method.Measured));
		kratios.add(new KRatioLabel(unkMcd, stdK411Mcd, feTr, Method.Measured));
		kratios.add(new KRatioLabel(unkMcd, stdK411Mcd, mgTr, Method.Measured));
		kratios.add(new KRatioLabel(unkMcd, stdK411Mcd, caTr, Method.Measured));
		kratios.add(new KRatioLabel(unkMcd, stdK411Mcd, oTr, Method.Measured));
		kratios.add(new KRatioLabel(unkMcd, stdAlMcd, alTr, Method.Measured));

		final Set<Object> outputs = new HashSet<>();
		for (final KRatioLabel krl : kratios) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection2.zafTag(krl.getUnknown(), meStd, cxr));
				outputs.add(XPPMatrixCorrection2.tagFxF(krl.getUnknown(), cxr));
				outputs.add(XPPMatrixCorrection2.tagFxF(meStd, cxr));
			}
		}

		final CompositionFromKRatios2 krcm = new CompositionFromKRatios2(kratios,
				XPPMatrixCorrection2.defaultVariates());
		final Map<Object, Double> constants = new HashMap<>();

		constants.put(Composition.buildMassFractionTag(unk, Element.Silicon), 0.211982);
		constants.put(Composition.buildMassFractionTag(unk, Element.Iron), 0.0774196);
		constants.put(Composition.buildMassFractionTag(unk, Element.Magnesium), 0.116567);
		constants.put(Composition.buildMassFractionTag(unk, Element.Calcium), 0.10899);
		constants.put(Composition.buildMassFractionTag(unk, Element.Oxygen), 0.42758);
		constants.put(Composition.buildMassFractionTag(unk, Element.Aluminum), 0.0490615);

		krcm.initializeConstants(constants);

		final Report r = new Report("KRationCorrectionModel - test 1");
		try {
			r.add(krcm, Mode.VERBOSE);

			final Map<Object, Number> mouv = new HashMap<>();

			mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, siTr, Method.Measured), new UncertainValue(0.794983));
			mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, feTr, Method.Measured), new UncertainValue(0.687101));
			mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, mgTr, Method.Measured), new UncertainValue(1.364719));
			mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, caTr, Method.Measured), new UncertainValue(0.980070));
			mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, oTr, Method.Measured), new UncertainValue(1.011527));
			mouv.put(new KRatioLabel(unkMcd, stdAlMcd, alTr, Method.Measured), new UncertainValue(0.032001));

			mouv.put(Composition.buildMassFractionTag(std0, Element.Silicon),
					new UncertainValue(0.25190067871134, 0.00448737523776));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Iron),
					new UncertainValue(0.11255374113608, 0.00209872307367));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Magnesium),
					new UncertainValue(0.09117902759616, 0.0012060717936));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Calcium),
					new UncertainValue(0.11070559976183, 0.00107203615005));
			mouv.put(Composition.buildMassFractionTag(std0, Element.Oxygen),
					new UncertainValue(0.42346095279459, 0.00693579374492));
			mouv.put(Composition.buildMassFractionTag(std1, Element.Aluminum), new UncertainValue(1.0, 0.001));

			mouv.put(new MatrixCorrectionLabel(unkMcd, stdK411Mcd, siTr), new UncertainValue(0.9519));
			mouv.put(new MatrixCorrectionLabel(unkMcd, stdK411Mcd, feTr), new UncertainValue(0.9948));
			mouv.put(new MatrixCorrectionLabel(unkMcd, stdK411Mcd, mgTr), new UncertainValue(1.0357));
			mouv.put(new MatrixCorrectionLabel(unkMcd, stdK411Mcd, caTr), new UncertainValue(0.9942));
			mouv.put(new MatrixCorrectionLabel(unkMcd, stdK411Mcd, oTr), new UncertainValue(1.0023));
			mouv.put(new MatrixCorrectionLabel(unkMcd, stdAlMcd, alTr), new UncertainValue(0.6523));

			final UncertainValues uvs = UncertainValues.extract(krcm.getInputLabels(), new UncertainValues(mouv));

			final LabeledMultivariateJacobian nmvj = new LabeledMultivariateJacobian(krcm, uvs.getValues());
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

	public static Map<Element, Number> buildK240() {
		final Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Magnesium, new UncertainValue(0.030154, 0.00030154));
		res.put(Element.Silicon, new UncertainValue(0.186986, 0.00186986));
		res.put(Element.Titanium, new UncertainValue(0.059950, 0.00059950));
		res.put(Element.Zinc, new UncertainValue(0.040168, 0.00040168));
		res.put(Element.Zirconium, new UncertainValue(0.074030, 0.00074030));
		res.put(Element.Barium, new UncertainValue(0.268689, 0.000268689));
		res.put(Element.Oxygen, new UncertainValue(0.340023, 0.00340023));
		return res;
	}

	@Test
	public void test2() throws ParseException, ArgumentException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds
		final Composition unk = Composition.massFraction("K240", buildK240());
		final Composition std0 = Composition.parse("Mg2SiO4");
		final Composition std1 = Composition.parse("BaTiSi3O9");
		final Composition std2 = Composition.parse("Zn");
		final Composition std3 = Composition.parse("Zr");

		// Mg
		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.2));
		// Ba, Ti, Si, O
		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)), //
				MatrixCorrectionDatum.roughness(10.0, 3.6));
		// Zn
		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		// Zr
		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 5.62));

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.5));

		final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds = new HashMap<>();
		// Mg
		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		stds.put(mgTrs, std0Mcd);
		// Ba, Ti, Si, O
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		stds.put(baTrs, std1Mcd);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		stds.put(tiTrs, std1Mcd);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		stds.put(siTrs, std1Mcd);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		stds.put(oTrs, std1Mcd);
		// Zn
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		stds.put(znTrs, std2Mcd);
		// Zr
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);
		stds.put(zrTrs, std3Mcd);

		final Map<KRatioLabel, Number> lkr = new HashMap<>();
		lkr.put(new KRatioLabel(unkMcd, std0Mcd, mgTrs, Method.Measured), new UncertainValue(0.066685, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, baTrs, Method.Measured), new UncertainValue(0.800181, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, tiTrs, Method.Measured), new UncertainValue(0.511920, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, siTrs, Method.Measured), new UncertainValue(0.914851, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, oTrs, Method.Measured), new UncertainValue(0.991157, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std2Mcd, znTrs, Method.Measured), new UncertainValue(0.035440, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std3Mcd, zrTrs, Method.Measured), new UncertainValue(0.051319, 0.001));
		final UncertainValues kratios = new UncertainValues(lkr);

		final CompositionFromKRatios2 cfk = new CompositionFromKRatios2(lkr.keySet(),
				XPPMatrixCorrection2.defaultVariates());
		final UncertainValues msInp = cfk.buildInputs(kratios);

		final Report report = new Report("K-Ratio (2)");
		try {
			report.addHeader("K-ratio Test2");
			report.addSubHeader("Inputs");
			report.add(msInp);
			report.addSubHeader("Function");
			report.add(cfk);
			// UncertainValues res = UncertainValues.propagate(ms, msInp);
			final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp.getValues());
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<? extends Object, UncertainValue> outVals = eval.getOutputValues(msInp);
			t.addRow(Table.th("Element"), Table.td("Mass Fraction"));
			for (final Element elm : unk.getElementSet()) {
				final MassFractionTag mft = Composition.buildMassFractionTag(unk, elm);
				final UncertainValue uv = outVals.get(mft);
				t.addRow(Table.td(elm), Table.td(HTML.toHTML(uv, Mode.VERBOSE)));
			}
			report.add(t);
		} catch (final Throwable e) {
			e.printStackTrace();
			report.addNote(e.getMessage());
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	@Test
	public void test3() throws ParseException, ArgumentException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds
		final Composition unk = Composition.massFraction("K240", buildK240());

		// Mg
		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("MgO"), //
				new UncertainValue(14.8, 0.2), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.6));
		// Ba
		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("BaSi2O5"), //
				new UncertainValue(14.9, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)), //
				MatrixCorrectionDatum.roughness(10.0, 3.74));
		// Zn
		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Zn"), //
				new UncertainValue(15.0, 0.05), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)), //
				MatrixCorrectionDatum.roughness(10.0, 7.14));

		// Zr
		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Zr"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 6.52));
		// Si, O
		final StandardMatrixCorrectionDatum std4Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("SiO2"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 2.65));

		// Ti
		final StandardMatrixCorrectionDatum std5Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Ti"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 4.5));

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.5));

		final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds = new HashMap<>();
		// Mg
		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		stds.put(mgTrs, std0Mcd);
		// Ba
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		stds.put(baTrs, std1Mcd);
		// Ti
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		stds.put(tiTrs, std5Mcd);
		// Si
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		stds.put(siTrs, std4Mcd);
		// O
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		stds.put(oTrs, std4Mcd);
		// Zn
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		stds.put(znTrs, std2Mcd);
		// Zr
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);
		stds.put(zrTrs, std3Mcd);

		final Map<KRatioLabel, Number> lkr = new HashMap<>();
		lkr.put(new KRatioLabel(unkMcd, std0Mcd, mgTrs, Method.Measured), new UncertainValue(0.036570, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, baTrs, Method.Measured), new UncertainValue(0.508579, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std2Mcd, znTrs, Method.Measured), new UncertainValue(0.035440, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std3Mcd, zrTrs, Method.Measured), new UncertainValue(0.051319, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std4Mcd, siTrs, Method.Measured), new UncertainValue(0.364811, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std4Mcd, oTrs, Method.Measured), new UncertainValue(0.571425, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std5Mcd, tiTrs, Method.Measured), new UncertainValue(0.055222, 0.001));

		final UncertainValues kuv = new UncertainValues(lkr);

		final CompositionFromKRatios2 cfk = new CompositionFromKRatios2(lkr.keySet(),
				XPPMatrixCorrection2.allVariates());

		final UncertainValues msInp = cfk.buildInputs(kuv);

		final Report report = new Report("K-Ratio (3)");
		try {
			report.addHeader("K-ratio Test3");
			report.addSubHeader("Inputs");
			report.add(msInp);
			report.addSubHeader("Function");
			report.addHTML(HTML.toHTML(cfk, Mode.NORMAL));
			// UncertainValues res = UncertainValues.propagate(ms, msInp);
			final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp.getValues());
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<? extends Object, UncertainValue> outVals = eval.getOutputValues(msInp, 1.0e-3);
			t.addRow(Table.th("Element"), Table.th("Mass Fraction"));
			for (final Element elm : unk.getElementSet()) {
				final MassFractionTag mft = Composition.buildMassFractionTag(unk, elm);
				final UncertainValue uv = outVals.get(mft);
				t.addRow(Table.td(elm), Table.td(HTML.toHTML(uv, Mode.VERBOSE)));
			}
			report.add(t);

			final UncertainValues results = UncertainValues.propagate(eval, msInp);
			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.00E0");
				for (final Object tag : msInp.sort().getLabels())
					t2.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
							Table.td(bnf.formatHTML(msInp.getEntry(tag))), //
							Table.td(bnf.formatHTML(msInp.getUncertainty(tag))));
				report.addSubHeader("Inputs");
				report.add(t2);
			}

			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("#,##0.00000");
				for (final Object tag : results.getLabels())
					t2.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
							Table.td(bnf.formatHTML(results.getEntry(tag))), //
							Table.td(bnf.formatHTML(results.getUncertainty(tag))));
				report.addSubHeader("Outputs");
				report.add(t2);
			}

			report.addImage(results.asCovarianceBitmap(4, V2L3, L2C), "Correlation matrix");

			final Set<Object> tags = new HashSet<>();
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Oxygen));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Magnesium));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Silicon));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Titanium));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Zinc));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Zirconium));
			tags.add(Composition.buildMassFractionTag(unkMcd.getComposition(), Element.Barium));
			cfk.trimOutputs(tags);

			final UncertainValues resultsTr = UncertainValues.propagate(cfk, msInp);
			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("#,##0.00000");
				for (final Object tag : resultsTr.getLabels())
					t2.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), Table.td(bnf.formatHTML(resultsTr.getEntry(tag))),
							Table.td(bnf.formatHTML(resultsTr.getUncertainty(tag))));
				report.add(t2);
			}
			report.addImage(resultsTr.asCovarianceBitmap(16, V2L3, L2C), "Correlation matrix (Trimmed)");

		} catch (final Throwable e) {
			e.printStackTrace();
			report.addNote(e.getMessage());
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

}
