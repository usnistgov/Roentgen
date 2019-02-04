package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.CompositionFromKRatios2;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
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

		final Composition std0 = XPPMatrixCorrection2Test.buildK411(false);

		final Composition std1 = Composition.parse("Al");

		final Composition unk = XPPMatrixCorrection2Test.buildK412(false);

		final StandardMatrixCorrectionDatum stdK411Mcd = new StandardMatrixCorrectionDatum( //
				std0, new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum stdAlMcd = new StandardMatrixCorrectionDatum( //
				std1, new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk, new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final ElementXRaySet feTr = ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1);
		final ElementXRaySet mgTr = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet caTr = ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1);
		final ElementXRaySet oTr = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet alTr = ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1);
		final ElementXRaySet siTr = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);

		final Map<KRatioLabel, Number> mouv = new HashMap<>();
		mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, siTr, Method.Measured), new UncertainValue(0.794983));
		mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, feTr, Method.Measured), new UncertainValue(0.687101));
		mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, mgTr, Method.Measured), new UncertainValue(1.364719));
		mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, caTr, Method.Measured), new UncertainValue(0.980070));
		mouv.put(new KRatioLabel(unkMcd, stdK411Mcd, oTr, Method.Measured), new UncertainValue(1.011527));
		mouv.put(new KRatioLabel(unkMcd, stdAlMcd, alTr, Method.Measured), new UncertainValue(0.032001));
		final UncertainValues kuv = new UncertainValues(mouv);

		final Set<Object> outputs = new HashSet<>();
		for (final KRatioLabel krl : mouv.keySet()) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(krl.getUnknown(), meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(krl.getUnknown(), cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}

		final Pair<SerialLabeledMultivariateJacobianFunction, UncertainValues> pr = KRatioCorrectionModel2
				.buildXPPModel(kuv, MatrixCorrectionModel2.defaultVariates());

		final LabeledMultivariateJacobianFunction krcm = pr.getFirst();
		final UncertainValues uvs = pr.getSecond();

		final Report r = new Report("KRationCorrectionModel - test 1");
		try {
			r.add(krcm, Mode.VERBOSE);

			final LabeledMultivariateJacobian nmvj = new LabeledMultivariateJacobian(krcm, uvs);
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
		// K240 glass using Mg2SiO4, BaTiSi3O9
		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Mg2SiO4"), //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.2));
		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		// Ba, Ti, Si, O
		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("BaTiSi3O9"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)), //
				MatrixCorrectionDatum.roughness(10.0, 3.6));
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		// Zn
		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Zn"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		// Zr
		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				Composition.parse("Zr"), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 5.62));
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);
		// Unknown
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				Composition.massFraction("K240", buildK240()), //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.5));

		final Map<KRatioLabel, Number> lkr = new HashMap<>();
		lkr.put(new KRatioLabel(unkMcd, std0Mcd, mgTrs, Method.Measured), new UncertainValue(0.066685, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, baTrs, Method.Measured), new UncertainValue(0.800181, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, tiTrs, Method.Measured), new UncertainValue(0.511920, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, siTrs, Method.Measured), new UncertainValue(0.914851, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, oTrs, Method.Measured), new UncertainValue(0.991157, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std2Mcd, znTrs, Method.Measured), new UncertainValue(0.035440, 0.001));
		lkr.put(new KRatioLabel(unkMcd, std3Mcd, zrTrs, Method.Measured), new UncertainValue(0.051319, 0.001));
		final UncertainValues kratios = new UncertainValues(lkr);

		final Pair<SerialLabeledMultivariateJacobianFunction, UncertainValues> pr = //
				KRatioCorrectionModel2.buildXPPModel(kratios, MatrixCorrectionModel2.defaultVariates());

		final LabeledMultivariateJacobianFunction krcm = pr.getFirst();
		final UncertainValues uvs = pr.getSecond();

		final Report report = new Report("K-Ratio (2)");
		try {
			report.addHeader("K-ratio Test2");
			report.addSubHeader("Inputs");
			report.add(uvs);
			report.addSubHeader("Function");
			report.add(krcm);
			// UncertainValues res = UncertainValues.propagate(ms, msInp);
			final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(krcm, uvs);
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<? extends Object, UncertainValue> outVals = eval.getOutputValues(uvs);
			t.addRow(Table.th("Element"), Table.td("Mass Fraction"));
			for (final Element elm : unkMcd.getElementSet()) {
				final MassFractionTag mft = Composition.buildMassFractionTag(unkMcd.getComposition(), elm);
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
		// K240 using MgO, BaSi2O5, Zn, Zr, SiO2 and Ti as standards
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

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

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
				MatrixCorrectionModel2.allVariates());

		final UncertainValues msInp = cfk.buildInputs(kuv);

		final Report report = new Report("K-Ratio (3)");
		try {
			report.addHeader("K-ratio Test3");
			report.addSubHeader("Inputs");
			report.add(msInp);
			report.addSubHeader("Function");
			report.addHTML(HTML.toHTML(cfk, Mode.NORMAL));
			// UncertainValues res = UncertainValues.propagate(ms, msInp);
			final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<? extends Object, UncertainValue> outVals = eval.getOutputValues(msInp, 1.0e-6);
			t.addRow(Table.th("Element"), Table.th("Mass Fraction"));
			for (final Map.Entry<? extends Object, UncertainValue> me : outVals.entrySet()) {
				if (me.getKey() instanceof KRatioLabel) {
					final UncertainValue uv = me.getValue();
					t.addRow(Table.td(me.getKey()), Table.td(HTML.toHTML(uv, Mode.VERBOSE)));
				}
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
				for (final Object tag : results.sort().getLabels())
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

	@Test
	public void test4() //
			throws ParseException, ArgumentException, IOException {
		// K240 using MgO, BaSi2O5, Zn, Zr, SiO2 and Ti as standards
		final Composition unk = Composition.massFraction("K240", buildK240());

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mgo = Composition.parse("MgO");
		final Composition basi2o5 = Composition.parse("BaSi2O5");
		final Composition zn = Composition.parse("Zn");
		final Composition zr = Composition.parse("Zr");
		final Composition ti = Composition.parse("Ti");
		final Composition sio2 = Composition.parse("SiO2");

		// KRatios computed using DTSA-II. 10 keV to 30 keV by 1 keV
		final double[] kOk = { 0.610759, 0.599076, 0.588755, 0.579709, 0.571847, 0.565075, 0.559301, 0.554432, 0.550382,
				0.547066, 0.544407, 0.542332, 0.540773, 0.539669, 0.538963, 0.538604, 0.538544, 0.538741, 0.539158,
				0.539761, 0.540519 };
		final double[] kMgK = { 0.043289, 0.041490, 0.039785, 0.038176, 0.036665, 0.035249, 0.033928, 0.032699,
				0.031557, 0.030499, 0.029521, 0.028618, 0.027786, 0.027019, 0.026314, 0.025667, 0.025073, 0.024529,
				0.024030, 0.023575, 0.023159 };
		final double[] kSiK = { 0.403768, 0.393957, 0.384338, 0.374920, 0.365717, 0.356747, 0.348022, 0.339558,
				0.331364, 0.323450, 0.315821, 0.308482, 0.301434, 0.294677, 0.288208, 0.282024, 0.276119, 0.270489,
				0.265126, 0.260023, 0.255171 };
		final double[] kTiK = { 0.055908, 0.055749, 0.055580, 0.055398, 0.055205, 0.055001, 0.054787, 0.054563,
				0.054329, 0.054086, 0.053836, 0.053577, 0.053312, 0.053040, 0.052761, 0.052477, 0.052188, 0.051893,
				0.051595, 0.051292, 0.050985 };
		final double[] kZnK = { 0.035080, 0.035080, 0.035194, 0.035282, 0.035346, 0.035391, 0.035418, 0.035431,
				0.035431, 0.035420, 0.035399, 0.035368, 0.035330, 0.035284, 0.035231, 0.035172, 0.035108, 0.035038,
				0.034963, 0.034884, 0.034800 };
		final double[] kZrL = { 0.054396, 0.053683, 0.052905, 0.052079, 0.051219, 0.050334, 0.049435, 0.048530,
				0.047623, 0.046722, 0.045831, 0.044953, 0.044092, 0.043249, 0.042427, 0.041628, 0.040853, 0.040101,
				0.039375, 0.038674, 0.037998 };
		final double[] kBaL = { 0.503054, 0.504528, 0.505786, 0.506870, 0.507813, 0.508641, 0.509374, 0.510028,
				0.510614, 0.511143, 0.511623, 0.512058, 0.512456, 0.512819, 0.513153, 0.513458, 0.513739, 0.513998,
				0.514236, 0.514455, 0.514656 };

		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final double roughness = MatrixCorrectionDatum.roughness(10.0, 3.6);

		final Map<Integer, Map<? extends Object, UncertainValue>> outVals = new TreeMap<>();
		final int MIN_E = 10, MAX_E = 31;
		List<? extends Object> outLabels = null, inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgoMcd = new StandardMatrixCorrectionDatum(mgo, e0, toa, roughness);
				final StandardMatrixCorrectionDatum basi2o5Mcd = new StandardMatrixCorrectionDatum(basi2o5, e0, toa,
						roughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, toa, roughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, toa, roughness);
				final StandardMatrixCorrectionDatum sio2Mcd = new StandardMatrixCorrectionDatum(sio2, e0, toa,
						roughness);
				final StandardMatrixCorrectionDatum tiMcd = new StandardMatrixCorrectionDatum(ti, e0, toa, roughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk, e0, toa, roughness);

				final Map<KRatioLabel, Number> lkr = new HashMap<>();
				lkr.put(new KRatioLabel(unkMcd, mgoMcd, mgTrs, Method.Measured),
						new UncertainValue(kMgK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, basi2o5Mcd, baTrs, Method.Measured),
						new UncertainValue(kBaL[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured),
						new UncertainValue(kZnK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured),
						new UncertainValue(kZrL[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, sio2Mcd, siTrs, Method.Measured),
						new UncertainValue(kSiK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, sio2Mcd, oTrs, Method.Measured),
						new UncertainValue(kOk[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured),
						new UncertainValue(kTiK[ie0 - MIN_E], 0.001));

				final UncertainValues kuv = new UncertainValues(lkr);

				final Pair<SerialLabeledMultivariateJacobianFunction, UncertainValues> pair = //
						KRatioCorrectionModel2.buildXPPModel(kuv, MatrixCorrectionModel2.allVariates());
				final SerialLabeledMultivariateJacobianFunction cfk = pair.getFirst();
				final Set<Object> finalOutputs = new HashSet<>();
				for (final Object output : cfk.getOutputLabels())
					if (output instanceof Composition.MassFractionTag)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final UncertainValues msInp = pair.getSecond();
				final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Jacobian");
					initReport.addHTML(eval.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
					initReport.addSubHeader("Output");
					initReport.addHTML(
							UncertainValues.propagate(cfk, msInp).toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, eval.getOutputValues(msInp, 0.0));
				if (ie0 == MIN_E)
					for (final Object label : cfk.getOutputLabels())
						System.out.println(label);
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = cfk.getInputLabels();
				}
			}
			initReport.inBrowser(Mode.NORMAL);
		}
		final Report report = new Report("K-Ratio (4)");
		try {
			report.addHeader("K240 against simple standards");
			final Table valTable = new Table();
			{
				final List<Item> row = new ArrayList<>();
				row.add(Table.td("Output"));
				row.add(Table.td("Abbrev."));
				row.add(Table.td("Input"));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				valTable.addRow(row);
			}
			final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
			for (int oi = 0; oi < outLabels.size(); ++oi) {
				final Object outTag = outLabels.get(oi);
				if (outTag instanceof MassFractionTag) {
					final MassFractionTag mft = (MassFractionTag) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final Object inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						boolean addRow = false;
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table
										.td(bnf.format(100.0 * getComponentByName(tmp, inTag) / tmp.doubleValue()))); //
							else
								row.add(Table.td("---"));
						}
						if (addRow)
							valTable.addRow(row);
					}
				}
			}
			report.addHTML(valTable.toHTML(Mode.NORMAL));
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	@Test
	public void test5() //
			throws ParseException, ArgumentException, IOException {
		// K240 using MgO, BaSi2O5, Zn, Zr, SiO2 and Ti as standards
		final Composition unk = Composition.massFraction("K240", buildK240());

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mgo = Composition.parse("MgO");
		final Composition benitoite = Composition.parse("BaTiSi3O9");
		final Composition zn = Composition.parse("Zn");
		final Composition zr = Composition.parse("Zr");

		final double[] kOk = { 0.986842, 0.98807, 0.989192, 0.990207, 0.991117, 0.991925, 0.992635, 0.993254, 0.993788,
				0.994244, 0.994629, 0.99495, 0.995213, 0.995425, 0.995591, 0.995717, 0.995808, 0.995869, 0.995903,
				0.995914, 0.995906 };
		final double[] kMgK = { 0.043289, 0.04149, 0.039785, 0.038176, 0.036665, 0.035249, 0.033928, 0.032699, 0.031557,
				0.030499, 0.029521, 0.028618, 0.027786, 0.027019, 0.026314, 0.025667, 0.025073, 0.024529, 0.02403,
				0.023575, 0.023159 };
		final double[] kSiK = { 0.915356, 0.91524, 0.915111, 0.914972, 0.914825, 0.914671, 0.914512, 0.914349, 0.914183,
				0.914016, 0.913847, 0.913679, 0.913511, 0.913344, 0.913179, 0.913015, 0.912854, 0.912696, 0.912541,
				0.912388, 0.912239 };
		final double[] kTiK = { 0.513686, 0.513278, 0.512855, 0.512407, 0.51193, 0.511425, 0.510893, 0.510335, 0.509753,
				0.509148, 0.508522, 0.507876, 0.507211, 0.506528, 0.505829, 0.505115, 0.504386, 0.503645, 0.502891,
				0.502126, 0.50135 };
		final double[] kZnK = { 0.023714, 0.03508, 0.035194, 0.035282, 0.035346, 0.035391, 0.035418, 0.035431, 0.035431,
				0.03542, 0.035399, 0.035368, 0.03533, 0.035284, 0.035231, 0.035172, 0.035108, 0.035038, 0.034963,
				0.034884, 0.0348 };
		final double[] kZrL = { 0.054396, 0.053683, 0.052905, 0.052079, 0.051219, 0.050334, 0.049435, 0.04853, 0.047623,
				0.046722, 0.045831, 0.044953, 0.044092, 0.043249, 0.042427, 0.041628, 0.040853, 0.040101, 0.039375,
				0.038674, 0.037998 };
		final double[] kBaL = { 0.802829, 0.802215, 0.801573, 0.800888, 0.800156, 0.799378, 0.798557, 0.797694,
				0.796793, 0.795857, 0.794887, 0.793885, 0.792854, 0.791796, 0.790712, 0.789604, 0.788474, 0.787324,
				0.786154, 0.784967, 0.783764 };

		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final double roughness = MatrixCorrectionDatum.roughness(10.0, 3.6);

		final Map<Integer, Map<? extends Object, UncertainValue>> outVals = new TreeMap<>();
		final int MIN_E = 10, MAX_E = 31;
		List<? extends Object> outLabels = null, inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgoMcd = new StandardMatrixCorrectionDatum(mgo, e0, toa, roughness);
				final StandardMatrixCorrectionDatum benitoiteMcd = new StandardMatrixCorrectionDatum(benitoite, e0, toa,
						roughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, toa, roughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, toa, roughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk, e0, toa, roughness);

				final Map<KRatioLabel, Number> lkr = new HashMap<>();
				lkr.put(new KRatioLabel(unkMcd, mgoMcd, mgTrs, Method.Measured),
						new UncertainValue(kMgK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, benitoiteMcd, baTrs, Method.Measured),
						new UncertainValue(kBaL[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured),
						new UncertainValue(kZnK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured),
						new UncertainValue(kZrL[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, benitoiteMcd, siTrs, Method.Measured),
						new UncertainValue(kSiK[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, benitoiteMcd, oTrs, Method.Measured),
						new UncertainValue(kOk[ie0 - MIN_E], 0.001));
				lkr.put(new KRatioLabel(unkMcd, benitoiteMcd, tiTrs, Method.Measured),
						new UncertainValue(kTiK[ie0 - MIN_E], 0.001));

				final UncertainValues kuv = new UncertainValues(lkr);

				final Pair<SerialLabeledMultivariateJacobianFunction, UncertainValues> pair = //
						KRatioCorrectionModel2.buildXPPModel(kuv, MatrixCorrectionModel2.allVariates());
				final SerialLabeledMultivariateJacobianFunction cfk = pair.getFirst();
				final Set<Object> finalOutputs = new HashSet<>();
				for (final Object output : cfk.getOutputLabels())
					if (output instanceof Composition.MassFractionTag)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final UncertainValues msInp = pair.getSecond();
				final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Jacobian");
					initReport.addHTML(eval.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
					initReport.addSubHeader("Output");
					initReport.addHTML(
							UncertainValues.propagate(cfk, msInp).toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, eval.getOutputValues(msInp, 0.0));
				if (ie0 == MIN_E)
					for (final Object label : cfk.getOutputLabels())
						System.out.println(label);
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = cfk.getInputLabels();
				}
			}
			initReport.inBrowser(Mode.NORMAL);
		}
		final Report report = new Report("K-Ratio (5)");
		try {
			report.addHeader("K240 against simple standards");
			final Table valTable = new Table();
			{
				final List<Item> row = new ArrayList<>();
				row.add(Table.td("Output"));
				row.add(Table.td("Abbrev."));
				row.add(Table.td("Input"));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				for (int i = MIN_E; i < MAX_E; ++i)
					row.add(Table.td(i));
				valTable.addRow(row);
			}
			final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
			for (int oi = 0; oi < outLabels.size(); ++oi) {
				final Object outTag = outLabels.get(oi);
				if (outTag instanceof MassFractionTag) {
					final MassFractionTag mft = (MassFractionTag) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final Object inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						boolean addRow = false;
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<? extends Object, UncertainValue> oVals = outVals.get(ie0);
							final UncertainValue tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table
										.td(bnf.format(100.0 * getComponentByName(tmp, inTag) / tmp.doubleValue()))); //
							else
								row.add(Table.td("---"));
						}
						if (addRow)
							valTable.addRow(row);
						else 
							System.out.println("Skipping "+outTag +" - "+ inTag);
					}
				}
			}
			report.addHTML(valTable.toHTML(Mode.NORMAL));
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	public double getComponentByName(final UncertainValue uv, final Object tag) {
		final String name = tag.toString();
		for (final Map.Entry<Object, Double> me : uv.getComponents().entrySet())
			if (me.getKey().toString().equals(name))
				return me.getValue().doubleValue();
		return 0.0;
	}

	public UncertainValue getByMFT(final Map<? extends Object, UncertainValue> oVals, final MassFractionTag mft) {
		for (final Map.Entry<? extends Object, UncertainValue> me : oVals.entrySet()) {
			if (me.getKey() instanceof MassFractionTag) {
				final MassFractionTag mft2 = (MassFractionTag) me.getKey();
				if (mft2.toString().equals(mft.toString()))
					return me.getValue();
			}
		}
		return null;
	}

}
