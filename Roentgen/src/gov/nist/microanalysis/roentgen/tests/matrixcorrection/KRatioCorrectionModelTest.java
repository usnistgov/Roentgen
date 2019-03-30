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
import java.util.Map.Entry;
import java.util.Set;

import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValueEx;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.DataStore.CompositionFactory;
import gov.nist.microanalysis.roentgen.math.Utility;
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
import gov.nist.microanalysis.roentgen.physics.composition.ElementByDifference;
import gov.nist.microanalysis.roentgen.physics.composition.ElementByStoichiometry;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class KRatioCorrectionModelTest {

	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);

	@Test
	public void test1() throws Exception {

		final Composition std0 = CompositionFactory.instance().findComposition("K411").getObject();

		final Composition std1 = Composition.parse("Al");

		final Composition unk = CompositionFactory.instance().findComposition("K412").getObject();

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
		final UncertainValues<KRatioLabel> kuv = new UncertainValues<KRatioLabel>(mouv);

		final Set<Object> outputs = new HashSet<>();
		for (final KRatioLabel krl : mouv.keySet()) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(krl.getUnknown(), meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(krl.getUnknown(), cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}

		final KRatioCorrectionModel2 krcm = KRatioCorrectionModel2.buildXPPModel(mouv.keySet(), null);
		final UncertainValuesBase<EPMALabel> uvs = krcm.buildInput(kuv);
		krcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final Report r = new Report("KRatioCorrectionModel - test 1");
		try {
			r.addHeader("KRatioCorrectionModel - Test 1");

			r.add(krcm, Mode.VERBOSE);

			final UncertainValuesCalculator<EPMALabel> nmvj = new UncertainValuesCalculator<>(krcm, uvs);
			r.addHeader("Jacobian");

			r.addHTML(nmvj.toHTML(Mode.NORMAL, new BasicNumberFormat("0.000")));

			r.addHeader("Input");
			r.add(uvs);

			r.addHeader("Result");
			r.add(nmvj);
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
	}


	@Test
	public void test2() throws Exception {
		// K240 glass using Mg2SiO4, BaTiSi3O9
		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				CompositionFactory.instance().findComposition("Forsterite").getObject(), //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.2));
		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		// Ba, Ti, Si, O
		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				CompositionFactory.instance().findComposition("Benitoite").getObject(), //
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
		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();
		
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
		final UncertainValues<KRatioLabel> kratios = new UncertainValues<KRatioLabel>(lkr);

		final KRatioCorrectionModel2 krcm = KRatioCorrectionModel2.buildXPPModel(lkr.keySet(), null);
		final UncertainValuesBase<EPMALabel> uvs = krcm.buildInput(kratios);
		krcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final Report report = new Report("K-Ratio (2)");
		try {
			report.addHeader("KRatioCorrectionModel - Test 2");

			report.addSubHeader("Inputs");
			report.add(uvs);
			report.addSubHeader("Function");
			report.add(krcm);
			final UncertainValuesCalculator<EPMALabel> eval = new UncertainValuesCalculator<EPMALabel>(krcm, uvs);
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = eval.getOutputValues();
			t.addRow(Table.th("Element"), Table.td("Mass Fraction"));
			for (final Element elm : unkMcd.getElementSet()) {
				final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(unkMcd.getMaterial(), elm);
				final UncertainValue uv = outVals.get(mft);
				t.addRow(Table.td(elm), Table.td(HTML.toHTML(uv, Mode.VERBOSE)));
			}
			report.add(t);
		} catch (final Throwable e) {
			report.addThrowable(e);
			report.inBrowser(Mode.NORMAL);
			throw e;
		}
		report.inBrowser(Mode.NORMAL);
	}

	@Test
	public void test3() throws Exception {
		// K240 using MgO, BaSi2O5, Zn, Zr, SiO2 and Ti as standards
		
		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

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

		final UncertainValuesBase<KRatioLabel> kuv = new UncertainValues<KRatioLabel>(lkr);

		final List<EPMALabel> tags = new ArrayList<>();
		tags.addAll(MaterialLabel.buildMassFractionTags(unkMcd.getMaterial()));
		for (final KRatioLabel krl : lkr.keySet()) {
			tags.add(MatrixCorrectionModel2.zafLabel(krl));
			tags.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
		}
		final KRatioCorrectionModel2 krcm = new KRatioCorrectionModel2(lkr.keySet(), null, tags);
		final UncertainValuesBase<EPMALabel> uvs = krcm.buildInput(kuv);
		krcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final Report report = new Report("K-Ratio (3)");
		try {
			report.addHeader("KRatioCorrectionModel - Test 3");

			report.addSubHeader("Inputs");
			report.add(uvs);
			report.addSubHeader("Function");
			report.addHTML(HTML.toHTML(krcm, Mode.NORMAL));
			// UncertainValuesBase res = UncertainValues.propagate(ms, msInp);
			final UncertainValuesCalculator<EPMALabel> eval = new UncertainValuesCalculator<>(krcm, uvs);
			report.addSubHeader("Jacobian");
			report.add(eval);
			report.addSubHeader("Result");
			final Table t = new Table();
			final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = eval.getOutputValues(1.0e-6);
			t.addRow(Table.th("Quantity"), Table.th("Verbose"), Table.th("Normal"));
			for (final Entry<EPMALabel, UncertainValueEx<EPMALabel>> me : outVals.entrySet()) {
				if ((me.getKey() instanceof MaterialLabel.MassFraction) //
						|| (me.getKey() instanceof KRatioLabel)) {
					final UncertainValue uv = me.getValue();
					t.addRow(Table.td(me.getKey()), Table.td(uv.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.0000"))),
							Table.td(uv.toHTML(Mode.TERSE, new BasicNumberFormat("0.0000"))));
				}
			}
			report.add(t);

			final UncertainValuesBase<EPMALabel> results = eval;
			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.00E0");
				for (final EPMALabel tag : uvs.sort().getLabels())
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
				for (final EPMALabel tag : results.sort().getLabels())
					t2.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
							Table.td(bnf.formatHTML(results.getEntry(tag))), //
							Table.td(bnf.formatHTML(results.getUncertainty(tag))));
				report.addSubHeader("Outputs");
				report.add(t2);
			}

			report.addImage(results.asCovarianceBitmap(4, V2L3, L2C), "Correlation matrix");

			final UncertainValuesBase<EPMALabel> resultsTr = UncertainValuesBase.propagateAnalytical(krcm, uvs);
			{
				final Table t2 = new Table();
				t2.addRow(Table.th("Tag"), Table.th("Value"), Table.th("Uncertainty"));
				final BasicNumberFormat bnf = new BasicNumberFormat("#,##0.00000");
				for (final EPMALabel tag : resultsTr.getLabels())
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
	public void iterationTest1() throws ArgumentException, ParseException, IOException {
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

		final Set<Element> elms = new HashSet<>(Arrays.asList(Element.Oxygen, Element.Magnesium, Element.Silicon,
				Element.Titanium, Element.Zinc, Element.Zirconium, Element.Barium));

		final Material unkMat = new Material("Unknown", elms);
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unkMat, beamEnergy, takeOffAngle);

		final Map<KRatioLabel, Number> vals = new HashMap<>();
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

		final UncertainValuesBase<KRatioLabel> kratios = new UncertainValues<KRatioLabel>(vals);

		final KRatioCorrectionModel2 iter = KRatioCorrectionModel2.buildXPPModel(vals.keySet(), null);

		// Takes 71 steps using basic iteration
		final UncertainValues<EPMALabel> uvs = UncertainValues.asUncertainValues(iter.iterate(kratios));

		final Composition comp = Composition.massFraction(unkMat, uvs);

		final Report rep = new Report("K240 Iteration");

		rep.addHeader("Iteration - O-by-K-ratio");
		rep.addSubHeader("Unknown");
		rep.add(unkMcd);
		rep.addSubHeader("K-ratios");
		rep.add(vals, Mode.NORMAL, Mode.VERBOSE);
		rep.addSubHeader("All Uncertainties");
		rep.add(uvs);
		rep.addSubHeader("Results");
		rep.add(comp);
		rep.inBrowser(Mode.NORMAL);
	}

	@Test
	public void iterationTest2() throws ArgumentException, ParseException, IOException {
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

		final Set<Element> elms = new HashSet<>(Arrays.asList(Element.Oxygen, Element.Magnesium, Element.Silicon,
				Element.Titanium, Element.Zinc, Element.Zirconium, Element.Barium));

		final Material unkMat = new Material("Unknown", elms);
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unkMat, beamEnergy, takeOffAngle);

		final Map<KRatioLabel, Number> vals = new HashMap<>();
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

		final UncertainValuesBase<KRatioLabel> kratios = new UncertainValues<KRatioLabel>(vals);

		final ExplicitMeasurementModel<MassFraction, MassFraction> elmByDiff = new ElementByDifference(unkMat,
				Element.Oxygen);

		final KRatioCorrectionModel2 iter = KRatioCorrectionModel2.buildXPPModel(vals.keySet(), elmByDiff);
		final UncertainValues<EPMALabel> uvs = UncertainValues.asUncertainValues(iter.iterate(kratios));
		final Composition comp = Composition.massFraction(iter.getUnknownMaterial(), uvs);

		final Report rep = new Report("K240 Iteration - O by differences");
		rep.addHeader("K240 Iteration - O-by-Difference");

		rep.addSubHeader("Unknown");
		rep.add(unkMcd);
		rep.addSubHeader("K-ratios");
		rep.add(vals, Mode.NORMAL, Mode.VERBOSE);
		rep.addSubHeader("All Uncertainties");
		rep.add(uvs);
		rep.addSubHeader("Results");
		rep.add(comp);
		rep.inBrowser(Mode.NORMAL);
	}

	@Test
	public void iterationTest3() throws ArgumentException, ParseException, IOException {
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

		final Set<Element> elms = new HashSet<>(Arrays.asList(Element.Oxygen, Element.Magnesium, Element.Silicon,
				Element.Titanium, Element.Zinc, Element.Zirconium, Element.Barium));

		final Material unkMat = new Material("Unknown", elms);
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unkMat, beamEnergy, takeOffAngle);

		final Map<KRatioLabel, Number> vals = new HashMap<>();
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

		// Takes 4 steps using basic iteration
		final UncertainValuesBase<KRatioLabel> kratios = new UncertainValues<KRatioLabel>(vals);

		final KRatioCorrectionModel2 iter = KRatioCorrectionModel2.buildXPPModel(//
				vals.keySet(), ElementByStoichiometry.buildDefaultOxygen(unkMat));
		// Takes 77 steps using basic iteration
		final UncertainValues<EPMALabel> uvs = UncertainValues.asUncertainValues(iter.iterate(kratios));
		final Composition comp = Composition.massFraction(iter.getUnknownMaterial(), uvs);

		final Report rep = new Report("K240 Iteration - O by stoichiometry");
		rep.addSubHeader("Unknown");
		rep.add(unkMcd);
		rep.addSubHeader("K-ratios");
		rep.add(vals, Mode.NORMAL, Mode.VERBOSE);
		rep.addSubHeader("All Uncertainties");
		rep.add(uvs);
		rep.addSubHeader("Results");
		rep.add(comp);
		rep.inBrowser(Mode.NORMAL);
	}

	/**
	 * Computes the calculated k-ratio for the unknown relative to the standards in
	 * kratios. Then converts the calculated k-ratios to 'measured' k-ratios and
	 * hands them as input to the iteration algorithm. The resulting 'measured'
	 * composition is compared to the input composition.
	 *
	 * @param estUnk
	 * @param kratios
	 * @param r
	 * @return
	 * @throws ArgumentException
	 */
	public boolean testComposition(
			final Composition estUnk, //
			final Set<KRatioLabel> kratios, //
			final Report r
		) throws ArgumentException {
		final Composition result = KRatioCorrectionModel2.roundTripXPP(estUnk, kratios);
		final Table t = new Table();
		t.addRow(Table.th("Input"), Table.th("Output"));
		t.addRow(Table.td(estUnk), Table.td(result));
		final Map<MassFraction, Double> diff = estUnk.differences(MassFraction.class, result);
		final MassFraction maxMF = Utility.largest(diff);
		t.addRow(Table.td("Max:<br/>&nbsp;&nbsp" + maxMF.toHTML(Mode.TERSE) + //
				"<br/>Value:<br/>&nbsp;&nbsp" + diff.get(maxMF)), //
				Table.td("RMS:<br/>&nbsp;&nbsp;" + Utility.rms(diff.values())));
		r.add(t, Mode.NORMAL);
		return (diff.get(maxMF) < 0.001) && (Utility.rms(diff.values()) < 0.002);
	}

}
