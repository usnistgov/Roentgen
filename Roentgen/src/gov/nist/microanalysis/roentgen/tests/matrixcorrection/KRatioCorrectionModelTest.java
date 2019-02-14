package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioCorrectionModel3;
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
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
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
				unk.getMaterial(), new UncertainValue(15.0, 0.12), //
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

		KRatioCorrectionModel2 krcm = new KRatioCorrectionModel2(mouv.keySet(),
				MatrixCorrectionModel2.defaultVariates());
		final UncertainValues uvs = UncertainValues.combine(krcm.buildInput(unk), kuv);

		final Report r = new Report("KRatioCorrectionModel - test 1");
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

	@Test
	public void test2vs3() throws ParseException, ArgumentException, IOException {

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
				unk.getMaterial(), new UncertainValue(15.0, 0.12), //
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

		KRatioCorrectionModel2 krcm2 = new KRatioCorrectionModel2(mouv.keySet(),
				MatrixCorrectionModel2.defaultVariates());
		final UncertainValues uvs2 = UncertainValues.combine(krcm2.buildInput(unk), kuv);

		KRatioCorrectionModel3 krcm3 = new KRatioCorrectionModel3(mouv.keySet(),
				MatrixCorrectionModel2.defaultVariates());
		final UncertainValues uvs3 = UncertainValues.combine(krcm2.buildInput(unk), kuv);

		final Report r = new Report("KRatioCorrectionModel - test 2 to 3 comparison");
		try {
			r.add(krcm2, Mode.VERBOSE);

			final LabeledMultivariateJacobian nmvj2 = new LabeledMultivariateJacobian(krcm2, uvs2);
			r.addHeader("Jacobian - 2");
			r.addHTML(nmvj2.toHTML(Mode.NORMAL, new BasicNumberFormat("0.000")));

			final LabeledMultivariateJacobian nmvj3 = new LabeledMultivariateJacobian(krcm3, uvs3);
			r.addHeader("Jacobian - 3");
			r.addHTML(nmvj3.toHTML(Mode.NORMAL, new BasicNumberFormat("0.000")));

			r.addHeader("Results");

			UncertainValues res2 = UncertainValues.propagate(nmvj2, uvs2);
			UncertainValues res3 = UncertainValues.propagate(nmvj3, uvs3);

			for (Object label : res2.getLabels()) {
				if (label instanceof MassFraction) {
					r.addSubHeader(HTML.toHTML(label, Mode.NORMAL));
					r.add(res2.getUncertainValue(label));
					r.add(res3.getUncertainValue(label));
				}
			}
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

		final Composition unk = Composition.massFraction("K240", buildK240());
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
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

		KRatioCorrectionModel2 krcm = new KRatioCorrectionModel2(lkr.keySet(),
				MatrixCorrectionModel2.defaultVariates());
		final UncertainValues uvs = UncertainValues.combine(krcm.buildInput(unk), kratios);

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
				final MaterialLabel.MassFraction mft = MaterialLabel
						.buildMassFractionTag(unkMcd.getMaterial(), elm);
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
				unk.getMaterial(), //
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
		lkr.put(new KRatioLabel(unkMcd, std0Mcd, mgTrs, Method.Measured),
				new UncertainValue(0.036570, 0.036570 * 0.014));
		lkr.put(new KRatioLabel(unkMcd, std1Mcd, baTrs, Method.Measured),
				new UncertainValue(0.508579, 0.508579 * 0.014));
		lkr.put(new KRatioLabel(unkMcd, std2Mcd, znTrs, Method.Measured),
				new UncertainValue(0.035440, 0.035440 * 0.111));
		lkr.put(new KRatioLabel(unkMcd, std3Mcd, zrTrs, Method.Measured),
				new UncertainValue(0.051319, 0.051319 * 0.0175));
		lkr.put(new KRatioLabel(unkMcd, std4Mcd, siTrs, Method.Measured),
				new UncertainValue(0.364811, 0.364811 * 0.0045));
		lkr.put(new KRatioLabel(unkMcd, std4Mcd, oTrs, Method.Measured),
				new UncertainValue(0.571425, 0.571425 * 0.005));
		lkr.put(new KRatioLabel(unkMcd, std5Mcd, tiTrs, Method.Measured),
				new UncertainValue(0.055222, 0.055222 * 0.0295));

		final UncertainValues kuv = new UncertainValues(lkr);

		KRatioCorrectionModel2 krcm = new KRatioCorrectionModel2(lkr.keySet(),
				MatrixCorrectionModel2.defaultVariates());
		final UncertainValues uvs = UncertainValues.combine(krcm.buildInput(unk), kuv);

		final Report report = new Report("K-Ratio (3)");
		try {
			report.addHeader("K-ratio Test3");
			report.addSubHeader("Inputs");
			report.add(uvs);
			report.addSubHeader("Function");
			report.addHTML(HTML.toHTML(krcm, Mode.NORMAL));
			// UncertainValues res = UncertainValues.propagate(ms, msInp);
			final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(krcm, uvs);
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<? extends Object, UncertainValue> outVals = eval.getOutputValues(uvs, 1.0e-6);
			t.addRow(Table.th("Quantity"), Table.th("Verbose"), Table.th("Normal"));
			for (final Map.Entry<? extends Object, UncertainValue> me : outVals.entrySet()) {
				if ((me.getKey() instanceof MaterialLabel.MassFraction) //
						|| (me.getKey() instanceof KRatioLabel)) {
					final UncertainValue uv = me.getValue();
					t.addRow(Table.td(me.getKey()), Table.td(uv.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.0000"))), Table.td(uv.toHTML(Mode.TERSE, new BasicNumberFormat("0.0000"))));
				}
			}
			report.add(t);

			final UncertainValues results = UncertainValues.propagate(eval, uvs);
			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.00E0");
				for (final Object tag : uvs.sort().getLabels())
					t2.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
							Table.td(bnf.formatHTML(uvs.getEntry(tag))), //
							Table.td(bnf.formatHTML(uvs.getUncertainty(tag))));
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
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Oxygen));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Magnesium));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Silicon));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Titanium));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Zinc));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Zirconium));
			tags.add(MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), Element.Barium));
			krcm.trimOutputs(tags);

			final UncertainValues resultsTr = UncertainValues.propagate(krcm, uvs);
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

		final Composition mgo = Composition.parse("Mg");
		final Composition basi2o5 = Composition.parse("BaSi2O5");
		final Composition zn = Composition.parse("Zn");
		final Composition zr = Composition.parse("Zr");
		final Composition ti = Composition.parse("Ti");
		final Composition sio2 = Composition.parse("SiO2");

		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final double roughness = MatrixCorrectionDatum.roughness(10.0, 3.6);

		final Map<Integer, Map<? extends Object, UncertainValue>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValues> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
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
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0, toa,
						roughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgoMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, basi2o5Mcd, baTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, sio2Mcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, sio2Mcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured));

				KRatioCorrectionModel2 cfk = new KRatioCorrectionModel2(lkr, MatrixCorrectionModel2.allVariates());
				final UncertainValues input = cfk.buildInput(unk);
				MatrixCorrectionModel2 mcm = cfk.getModel();
				// Calculate the optimal k-ratios
				RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				UncertainValues krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(), Method.Measured);
				final UncertainValues msInp = UncertainValues.combine(cfk.buildInput(unk), krs);

				final Set<Object> finalOutputs = new HashSet<>();
				for (final Object output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
				final UncertainValues res = UncertainValues.propagate(cfk, msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Jacobian");
					initReport.addHTML(eval.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, eval.getOutputValues(msInp, 0.0));
				resVals.put(ie0, res);
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
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
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
					final List<Item> row = new ArrayList<>();
					row.add(Table.td(mft));
					row.add(Table.td(mft.getElement()));
					row.add(Table.td("Overall"));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getEntry(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getUncertainty(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
						final UncertainValues uvs = resVals.get(ie0);
						row.add(Table.td(bnf.format(100.0 * uvs.getUncertainty(outTag) / uvs.getEntry(outTag))));
					}
					valTable.addRow(row);
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
		// K240 using Benitoite, Mg, Zn and Zr as standards

		final double dkO[] = { 0.0014, 0.0015, 0.0015, 0.0015, 0.0016, 0.0014, 0.0016, 0.0017, 0.0017, 0.0017, 0.0017,
				0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0018, 0.0020, 0.0020 };
		final double dkMg[] = { 0.0052, 0.0055, 0.0059, 0.0063, 0.0067, 0.0071, 0.0076, 0.0081, 0.0085, 0.0090, 0.0096,
				0.0101, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100 };
		final double dkSi[] = { 0.0015, 0.0014, 0.0014, 0.0014, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013,
				0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013, 0.0013 };
		final double dkTi[] = { 0.0138, 0.0126, 0.0113, 0.0102, 0.0096, 0.0091, 0.0087, 0.0081, 0.0077, 0.0073, 0.0068,
				0.0067, 0.0064, 0.0063, 0.0062, 0.0060, 0.0058, 0.0056, 0.0055 };
		final double dkZn[] = { 0.0971, 0.0590, 0.0405, 0.0328, 0.0290, 0.0224, 0.0196, 0.0193, 0.0170, 0.0141, 0.0141,
				0.0142, 0.0111, 0.0113, 0.0113, 0.0114, 0.0086, 0.0087, 0.0087 };
		final double dkZr[] = { 0.0056, 0.0057, 0.0058, 0.0059, 0.0061, 0.0062, 0.0063, 0.0043, 0.0043, 0.0045, 0.0045,
				0.0046, 0.0047, 0.0048, 0.0049, 0.0050, 0.0052, 0.0052, 0.0052 };
		final double dkBa[] = { 0.0083, 0.0072, 0.0066, 0.0061, 0.0057, 0.0052, 0.0050, 0.0047, 0.0045, 0.0043, 0.0042,
				0.0040, 0.0039, 0.0037, 0.0036, 0.0035, 0.0034, 0.0033, 0.0033 };

		final Composition unk = Composition.massFraction("K240", buildK240());

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mg = Composition.parse("Mg");
		final Composition zn = Composition.parse("Zn");
		final Composition zr = Composition.parse("Zr");
		final Composition benitoite = Composition.parse("BaTiSi3O9");

		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final double roughness = MatrixCorrectionDatum.roughness(10.0, 3.6);

		final Map<Integer, Map<? extends Object, UncertainValue>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValues> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<? extends Object> outLabels = null, inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, toa, roughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, toa, roughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, toa, roughness);
				final StandardMatrixCorrectionDatum benitoiteMcd = new StandardMatrixCorrectionDatum(//
						benitoite, e0, toa, roughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(//
						unk.getMaterial(), e0, toa, roughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, tiTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, baTrs, Method.Measured));

				KRatioCorrectionModel2 cfk = new KRatioCorrectionModel2(lkr, MatrixCorrectionModel2.allVariates());
				final UncertainValues input = cfk.buildInput(unk);
				MatrixCorrectionModel2 mcm = cfk.getModel();
				// Calculate the optimal k-ratios
				RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				UncertainValues krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(), Method.Measured);
				for (Object label : krs.getLabels()) {
					KRatioLabel krl = (KRatioLabel) label;
					final double v = krs.getEntry(krl);
					double dk = 0.0;
					switch (krl.getXRaySet().getElement()) {
					case Oxygen:
						dk = v * dkO[ie0 - MIN_E];
						break;
					case Magnesium:
						dk = v * dkMg[ie0 - MIN_E];
						break;
					case Silicon:
						dk = v * dkSi[ie0 - MIN_E];
						break;
					case Titanium:
						dk = v * dkTi[ie0 - MIN_E];
						break;
					case Zinc:
						dk = v * dkZn[ie0 - MIN_E];
						break;
					case Zirconium:
						dk = v * dkZr[ie0 - MIN_E];
						break;
					case Barium:
						dk = v * dkBa[ie0 - MIN_E];
						break;
					default:
						assert false;
						break;

					}
					krs.set(krl, new UncertainValue(v, dk));
				}

				final UncertainValues msInp = UncertainValues.combine(cfk.buildInput(unk), krs);

				final Set<Object> finalOutputs = new HashSet<>();
				for (final Object output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
				final UncertainValues res = UncertainValues.propagate(cfk, msInp);
				resVals.put(ie0, res);
				if ((ie0 - MIN_E) % 5 == 0) {
					initReport.addHeader("E<sub>0</sub> = " + ie0);
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Jacobian");
					initReport.addHTML(eval.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, eval.getOutputValues(msInp, 0.0));
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = cfk.getInputLabels();
				}
			}
			initReport.inBrowser(Mode.NORMAL);
		}
		final Report report = new Report("K-Ratio (5)");
		try {
			report.addHeader("K240 against benitoite");
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
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
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
					}
					final List<Item> row = new ArrayList<>();
					row.add(Table.td(mft));
					row.add(Table.td(mft.getElement()));
					row.add(Table.td("Overall"));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getEntry(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getUncertainty(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
						final UncertainValues uvs = resVals.get(ie0);
						row.add(Table.td(bnf.format(100.0 * uvs.getUncertainty(outTag) / uvs.getEntry(outTag))));
					}
					valTable.addRow(row);
				}
			}
			report.addHTML(valTable.toHTML(Mode.NORMAL));
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	@Test
	public void test6() //
			throws ParseException, ArgumentException, IOException {
		// K240 using simple elements
		final Composition unk = Composition.massFraction("K240", buildK240());

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mg = Composition.parse("Mg");
		final Composition ba = Composition.parse("Ba");
		final Composition zn = Composition.parse("Zn");
		final Composition zr = Composition.parse("Zr");
		final Composition ti = Composition.parse("Ti");
		final Composition si = Composition.parse("Si");
		final Composition o = Composition.parse("O");

		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final double roughness = MatrixCorrectionDatum.roughness(10.0, 3.6);

		final Map<Integer, Map<? extends Object, UncertainValue>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValues> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<? extends Object> outLabels = null, inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, toa, roughness);
				final StandardMatrixCorrectionDatum baMcd = new StandardMatrixCorrectionDatum(ba, e0, toa, roughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, toa, roughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, toa, roughness);
				final StandardMatrixCorrectionDatum tiMcd = new StandardMatrixCorrectionDatum(ti, e0, toa, roughness);
				final StandardMatrixCorrectionDatum siMcd = new StandardMatrixCorrectionDatum(si, e0, toa, roughness);
				final StandardMatrixCorrectionDatum oMcd = new StandardMatrixCorrectionDatum(o, e0, toa, roughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0, toa,
						roughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, baMcd, baTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, siMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, oMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured));

				final KRatioCorrectionModel2 cfk = new KRatioCorrectionModel2(lkr,
						MatrixCorrectionModel2.allVariates());
				final UncertainValues input = cfk.buildInput(unk);
				MatrixCorrectionModel2 mcm = cfk.getModel();
				// Calculate the optimal k-ratios
				RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				UncertainValues krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(), Method.Measured);
				final UncertainValues msInp = UncertainValues.combine(cfk.buildInput(unk), krs);

				final Set<Object> finalOutputs = new HashSet<>();
				for (final Object output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final LabeledMultivariateJacobian eval = new LabeledMultivariateJacobian(cfk, msInp);
				final UncertainValues res = UncertainValues.propagate(cfk, msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Jacobian");
					initReport.addHTML(eval.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, eval.getOutputValues(msInp, 0.0));
				resVals.put(ie0, res);
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
		final Report report = new Report("K-Ratio (6)");
		try {
			report.addHeader("K240 using simple elements");
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
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
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
							System.out.println("Skipping " + outTag + " - " + inTag);
					}
					final List<Item> row = new ArrayList<>();
					row.add(Table.td(mft));
					row.add(Table.td(mft.getElement()));
					row.add(Table.td("Overall"));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getEntry(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0)
						row.add(Table.td(bnf.format(resVals.get(ie0).getUncertainty(outTag))));
					for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
						final UncertainValues uvs = resVals.get(ie0);
						row.add(Table.td(bnf.format(100.0 * uvs.getUncertainty(outTag) / uvs.getEntry(outTag))));
					}
					valTable.addRow(row);
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

	public UncertainValue getByMFT(final Map<? extends Object, UncertainValue> oVals,
			final MaterialLabel.MassFraction mft) {
		for (final Map.Entry<? extends Object, UncertainValue> me : oVals.entrySet()) {
			if (me.getKey() instanceof MaterialLabel.MassFraction) {
				final MaterialLabel.MassFraction mft2 = (MaterialLabel.MassFraction) me.getKey();
				if (mft2.toString().equals(mft.toString()))
					return me.getValue();
			}
		}
		return null;
	}

	@Test
	public void iterationTest1() throws ArgumentException, ParseException {
		// Material K240 = [O(0.3400 mass frac),Mg(0.0302 mass frac),Si(0.1870 mass
		// frac),Ti(0.0600 mass frac),Zn(0.0402 mass frac),Zr(0.0740 mass
		// frac),Ba(0.2687 mass frac),Σ=1.0001]
		// Detector Si(Li) - FWHM[Mn Kα]=154.1 eV - 2014-05-29 00:00
		// Algorithm XPP - Pouchou & Pichoir Simplified
		// MAC NIST-Chantler 2005
		// E0 15 keV
		// Take-off 40°
		// IUPAC Seigbahn Standard Energy ZAF Z A F k-ratio
		// O K-L3 O Kα1 SiO2 0.5249 0.8951 1.1356 0.7885 0.9996 0.571425
		// Mg K-L3 Mg Kα1 Pure Mg 1.2536 0.5677 1.0961 0.5160 1.0038 0.017146
		// Si K-L3 Si Kα1 SiO2 1.7397 0.9119 1.1390 0.7985 1.0027 0.364811
		// Ti K-L3 Ti Kα1 Pure Ti 4.5109 0.9204 0.9441 0.9684 1.0067 0.055222
		// Zn K-L3 Zn Kα1 Pure Zn 8.6389 0.8816 0.8920 0.9884 1.0000 0.035440
		// Zr L3-M5 Zr Lα1 Pure Zr 2.0423 0.6935 0.8393 0.8248 1.0017 0.051319
		// Ba L3-M5 Ba Lα1 BaF2 4.4663 0.8279 0.8194 1.0107 0.9996 0.283990

		final UncertainValue beamEnergy = new UncertainValue(15.0, 0.1);
		final UncertainValue takeOffAngle = new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9));
		final StandardMatrixCorrectionDatum mcd_sio2 = new StandardMatrixCorrectionDatum( //
				Composition.parse("SiO2"), beamEnergy, takeOffAngle);
		final StandardMatrixCorrectionDatum mcd_mg = new StandardMatrixCorrectionDatum( //
				Composition.pureElement(Element.Magnesium), beamEnergy, takeOffAngle);
		final StandardMatrixCorrectionDatum mcd_ti = new StandardMatrixCorrectionDatum( //
				Composition.pureElement(Element.Titanium), beamEnergy, takeOffAngle);
		final StandardMatrixCorrectionDatum mcd_zn = new StandardMatrixCorrectionDatum( //
				Composition.pureElement(Element.Zinc), beamEnergy, takeOffAngle);
		final StandardMatrixCorrectionDatum mcd_zr = new StandardMatrixCorrectionDatum( //
				Composition.pureElement(Element.Zirconium), beamEnergy, takeOffAngle);
		final StandardMatrixCorrectionDatum mcd_baf2 = new StandardMatrixCorrectionDatum( //
				Composition.parse("BaF2"), beamEnergy, takeOffAngle);

		final Set<KRatioLabel> skrl = new HashSet<>();

		final Set<Element> elms = new HashSet<>(Arrays.asList(Element.Oxygen, Element.Silicon, Element.Titanium,
				Element.Zinc, Element.Zirconium, Element.Barium));

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(new Material("Unknown", elms),
				beamEnergy, takeOffAngle);

		skrl.add(new KRatioLabel(unkMcd, mcd_sio2, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_sio2, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_mg, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_ti, ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_zn, ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_zr, ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, mcd_baf2, ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1),
				Method.Measured));

		final KRatioCorrectionModel2 kri = new KRatioCorrectionModel2(skrl, MatrixCorrectionModel2.defaultVariates());

		final Map<Object, Number> vals = new HashMap<>();
		vals.put(new KRatioLabel(unkMcd, mcd_sio2, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
				KRatioLabel.Method.Measured), 0.571425);
		vals.put(new KRatioLabel(unkMcd, mcd_mg, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1),
				KRatioLabel.Method.Measured), 0.017146);
		vals.put(new KRatioLabel(unkMcd, mcd_sio2, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
				KRatioLabel.Method.Measured), 0.364811);
		vals.put(new KRatioLabel(unkMcd, mcd_ti, ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1),
				KRatioLabel.Method.Measured), 0.05522);
		vals.put(new KRatioLabel(unkMcd, mcd_zn, ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1),
				KRatioLabel.Method.Measured), 0.035440);
		vals.put(new KRatioLabel(unkMcd, mcd_zr, ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1),
				KRatioLabel.Method.Measured), 0.051319);
		vals.put(new KRatioLabel(unkMcd, mcd_baf2, ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1),
				KRatioLabel.Method.Measured), 0.283990);

		final UncertainValues kratios = new UncertainValues(vals);

		final Composition res = kri.optimize(kratios);

		final Report rep = new Report("K240 Iteration");
		rep.addSubHeader("Unknown");
		rep.add(unkMcd);
		rep.addSubHeader("K-ratios");
		rep.add(vals, Mode.NORMAL, Mode.VERBOSE);
		rep.addSubHeader("Results");
		rep.add(res);

	}

}
