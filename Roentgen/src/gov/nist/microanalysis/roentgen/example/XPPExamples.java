package gov.nist.microanalysis.roentgen.example;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.duckandcover.html.Report;
import com.d3x.morpheus.frame.DataFrame;
import com.d3x.morpheus.util.text.Formats;
import com.d3x.morpheus.viz.chart.Chart;
import com.duckandcover.html.IToHTML.Mode;
import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.EPMALabel.MaterialMAC;
import gov.nist.microanalysis.roentgen.DataStore.CompositionFactory;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.ZAFLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.ZAFMultiLineLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class XPPExamples {

	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);

	public static void main(String[] args) throws Exception {
		// CuInAl29keV();
		// K240at15keV();
		// K240at20keV();
		// K240at20keVAlt();
		// K240at20keVSimple();
		// K240at30keV();
		K240atMultikeV();
		K240atMultikeVSimilar();
	}

	public static void CuInAl29keV() throws ArgumentException, IOException {
		Map<Element, Double> men = new HashMap<>();
		men.put(Element.Aluminum, 1.0 - 1e-8);
		men.put(Element.Copper, 1e-8); // Must be a little bit of Cu
		Composition unkComp = Composition.massFraction("Al", men);
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(29.0), UncertainValue.toRadians(40.0, 0.0));
		StandardMatrixCorrectionDatum std = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Copper),
				UncertainValue.valueOf(29.0), UncertainValue.toRadians(40.0, 0.0));
		KRatioLabel krl = new KRatioLabel(unk, std, CharacteristicXRay.create(Element.Copper, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		Set<KRatioLabel> kratios = Collections.singleton(krl);
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, true);
		Report r = new Report("Cu K-L3 in Al at 29.0 keV");
		r.addHeader("Cu K-L3 in Al at 29.0 keV");
		// r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
				.computePhiRhoZCurve(xpp.getValueMap(), 0.003, 0.003 / 100, 0.1);
		final String userHome = System.getProperty("user.home");
		prz.write().csv(new File(userHome + "\\Desktop\\prz_Cu29keV.csv")).apply(options -> {
			options.setSeparator(",");
			options.setIncludeRowHeader(true);
			options.setIncludeColumnHeader(true);
		});
		Chart.create().asHtml().withLinePlot(prz, chart -> {
			chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
			chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
			chart.legend().on().bottom();
			chart.show();
		});
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240at20keV() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(20.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 5.0));
		StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(20.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 1.738));
		StandardMatrixCorrectionDatum stdSi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Silicon), UncertainValue.valueOf(20.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 2.3290));
		StandardMatrixCorrectionDatum stdTi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Titanium), UncertainValue.valueOf(20.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 4.506));
		StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Zinc),
				UncertainValue.valueOf(20.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(20.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 6.52));
		StandardMatrixCorrectionDatum stdBa = new StandardMatrixCorrectionDatum(Composition.parse("BaF2"),
				UncertainValue.valueOf(20.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 4.893));
		StandardMatrixCorrectionDatum stdO = new StandardMatrixCorrectionDatum(Composition.parse("SiO2"),
				UncertainValue.valueOf(20.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 2.196));

		KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlSi = new KRatioLabel(unk, stdSi, CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlTi = new KRatioLabel(unk, stdTi, CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
				CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
		KRatioLabel krlBa = new KRatioLabel(unk, stdBa, CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlO = new KRatioLabel(unk, stdO, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);

		Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
				Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
		Report r = new Report("K240 at 20.0 keV");
		r.addHeader("K240 at 20.0 keV");
		r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.addImage(xpp.asCovarianceBitmap(32, L2C, L2C, true), "Covariances");

		// r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));

		// List<ZAFMultiLineLabel> zafs =
		// kratios.stream().flatMap(krl->MatrixCorrectionModel2.zafLabel(krl)).collect(Collectors.toList());
		ArrayList<EPMALabel> labels = new ArrayList<>();
		labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
		// labels.addAll(xpp.getLabels(MaterialMAC.class));
		labels.addAll(xpp.getLabels(ZAFLabel.class));
		UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
		r.addHTML(zafMac.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addImage(zafMac.asCovarianceBitmap(64, L2C, L2C, true), "MACS -> ZAF");

		File file = new File(System.getProperty("user.home") + "\\Desktop", "K240_ZAF_Elm.csv");
		try (PrintWriter pw = new PrintWriter(file, Charset.forName("UTF-8"))) {
			pw.write(zafMac.toCSV());
		}
		if (false) {
			final double range = 0.0015;
			final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
					.computePhiRhoZCurve(xpp.getValueMap(), range, range / 100, 0.1);
			prz.write().csv(new File(System.getProperty("user.home") + "\\Desktop", "prz_K240.csv")).apply(options -> {
				options.setSeparator(",");
				options.setIncludeRowHeader(true);
				options.setIncludeColumnHeader(true);
			});
			Chart.create().asHtml().withLinePlot(prz, chart -> {
				chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
				chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
				chart.legend().on().bottom();
				chart.show();
			});
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240at20keVAlt() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		final Composition benitoite = Composition.parse("BaTi(Si3O9)");

		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdMulti = new StandardMatrixCorrectionDatum(benitoite,
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Zinc),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));

		KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlSi = new KRatioLabel(unk, stdMulti, CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlTi = new KRatioLabel(unk, stdMulti, CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
				CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
		KRatioLabel krlBa = new KRatioLabel(unk, stdMulti, CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlO = new KRatioLabel(unk, stdMulti, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);

		Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
				Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
		Report r = new Report("K240-similar");
		r.addHeader("K240 at 20.0 keV (Ba,Si,Ti,O => BaTi(Si3O9))");
		r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.addImage(xpp.asCovarianceBitmap(32, L2C, L2C, true), "Covariances");

		// r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));

		// List<ZAFMultiLineLabel> zafs =
		// kratios.stream().flatMap(krl->MatrixCorrectionModel2.zafLabel(krl)).collect(Collectors.toList());
		ArrayList<EPMALabel> labels = new ArrayList<>();
		labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
		// labels.addAll(xpp.getLabels(MaterialMAC.class));
		labels.addAll(xpp.getLabels(ZAFMultiLineLabel.class));
		UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
		r.addHTML(zafMac.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addImage(zafMac.asCovarianceBitmap(64, L2C, L2C, true), "MACS -> ZAF");

		File file = new File(System.getProperty("user.home") + "\\Desktop", "K240_ZAF_Elm.csv");
		try (PrintWriter pw = new PrintWriter(file, Charset.forName("UTF-8"))) {
			pw.write(zafMac.toCSV());
		}
		if (false) {
			final double range = 0.0015;
			final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
					.computePhiRhoZCurve(xpp.getValueMap(), range, range / 100, 0.1);
			prz.write().csv(new File(System.getProperty("user.home") + "\\Desktop", "prz_K240.csv")).apply(options -> {
				options.setSeparator(",");
				options.setIncludeRowHeader(true);
				options.setIncludeColumnHeader(true);
			});
			Chart.create().asHtml().withLinePlot(prz, chart -> {
				chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
				chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
				chart.legend().on().bottom();
				chart.show();
			});
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240at20keVSimple() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdSi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Silicon), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdTi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Titanium), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Zinc),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(20.0, 0.0),
				UncertainValue.toRadians(40.0, 0.0), MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdBa = new StandardMatrixCorrectionDatum(Composition.parse("BaF2"),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));
		StandardMatrixCorrectionDatum stdO = new StandardMatrixCorrectionDatum(Composition.parse("SiO2"),
				UncertainValue.valueOf(20.0, 0.0), UncertainValue.toRadians(40.0, 0.0),
				MatrixCorrectionDatum.roughness(0.1, 5.0));

		KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlSi = new KRatioLabel(unk, stdSi, CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlTi = new KRatioLabel(unk, stdTi, CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
				CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
		KRatioLabel krlBa = new KRatioLabel(unk, stdBa, CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlO = new KRatioLabel(unk, stdO, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);

		Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
				Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
		Report r = new Report("K240-simple");
		r.addHeader("\"K240 at 20.0 keV (O=SiO2, Ba=BaF2, pure otherwise)\"");
		r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.addImage(xpp.asCovarianceBitmap(32, L2C, L2C, true), "Covariances");

		// r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));

		// List<ZAFMultiLineLabel> zafs =
		// kratios.stream().flatMap(krl->MatrixCorrectionModel2.zafLabel(krl)).collect(Collectors.toList());
		ArrayList<EPMALabel> labels = new ArrayList<>();
		labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
		// labels.addAll(xpp.getLabels(MaterialMAC.class));
		labels.addAll(xpp.getLabels(ZAFMultiLineLabel.class));
		UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
		r.addHTML(zafMac.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addImage(zafMac.asCovarianceBitmap(64, L2C, L2C, true), "MACS -> ZAF");

		File file = new File(System.getProperty("user.home") + "\\Desktop", "K240_ZAF_Elm.csv");
		try (PrintWriter pw = new PrintWriter(file, Charset.forName("UTF-8"))) {
			pw.write(zafMac.toCSV());
		}
		if (false) {
			final double range = 0.0015;
			final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
					.computePhiRhoZCurve(xpp.getValueMap(), range, range / 100, 0.1);
			prz.write().csv(new File(System.getProperty("user.home") + "\\Desktop", "prz_K240.csv")).apply(options -> {
				options.setSeparator(",");
				options.setIncludeRowHeader(true);
				options.setIncludeColumnHeader(true);
			});
			Chart.create().asHtml().withLinePlot(prz, chart -> {
				chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
				chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
				chart.legend().on().bottom();
				chart.show();
			});
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240at30keV() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(30.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 5.0));
		StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(30.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 1.738));
		StandardMatrixCorrectionDatum stdSi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Silicon), UncertainValue.valueOf(30.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 2.3290));
		StandardMatrixCorrectionDatum stdTi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Titanium), UncertainValue.valueOf(30.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 4.506));
		StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Zinc),
				UncertainValue.valueOf(30.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(30.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 6.52));
		StandardMatrixCorrectionDatum stdBa = new StandardMatrixCorrectionDatum(Composition.parse("BaF2"),
				UncertainValue.valueOf(30.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 4.893));
		StandardMatrixCorrectionDatum stdO = new StandardMatrixCorrectionDatum(Composition.parse("SiO2"),
				UncertainValue.valueOf(30.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 2.196));

		KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlSi = new KRatioLabel(unk, stdSi, CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlTi = new KRatioLabel(unk, stdTi, CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
				CharacteristicXRay.create(Element.Zirconium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlBa = new KRatioLabel(unk, stdBa, CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlO = new KRatioLabel(unk, stdO, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);

		Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
				Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
		Report r = new Report("K240 at 30.0 keV");
		r.addHeader("K240 at 30.0 keV");
		r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.addImage(xpp.asCovarianceBitmap(16, L2C, L2C, true), "Covariances");

		// r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));

		// List<ZAFMultiLineLabel> zafs =
		// kratios.stream().flatMap(krl->MatrixCorrectionModel2.zafLabel(krl)).collect(Collectors.toList());
		ArrayList<EPMALabel> labels = new ArrayList<>();
		labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
		// labels.addAll(xpp.getLabels(MaterialMAC.class));
		labels.addAll(xpp.getLabels(ZAFMultiLineLabel.class));
		UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
		r.addHTML(zafMac.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addImage(zafMac.asCovarianceBitmap(64, L2C, L2C, true), "MACS -> ZAF");

		File file = new File(System.getProperty("user.home") + "\\Desktop", "K240_ZAF_Elm.csv");
		try (PrintWriter pw = new PrintWriter(file, Charset.forName("UTF-8"))) {
			pw.write(zafMac.toCSV());
		}
		if (false) {
			final double range = 0.0015;
			final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
					.computePhiRhoZCurve(xpp.getValueMap(), range, range / 100, 0.1);
			prz.write().csv(new File(System.getProperty("user.home") + "\\Desktop", "prz_K240.csv")).apply(options -> {
				options.setSeparator(",");
				options.setIncludeRowHeader(true);
				options.setIncludeColumnHeader(true);
			});
			Chart.create().asHtml().withLinePlot(prz, chart -> {
				chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
				chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
				chart.legend().on().bottom();
				chart.show();
			});
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240at15keV() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
				UncertainValue.valueOf(15.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 5.0));
		StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(15.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 1.738));
		StandardMatrixCorrectionDatum stdSi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Silicon), UncertainValue.valueOf(15.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 2.3290));
		StandardMatrixCorrectionDatum stdTi = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Titanium), UncertainValue.valueOf(15.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 4.506));
		StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(Composition.pureElement(Element.Zinc),
				UncertainValue.valueOf(15.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
				Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(15.0, 0.1),
				UncertainValue.toRadians(40.0, 0.5), MatrixCorrectionDatum.roughness(10.0, 6.52));
		StandardMatrixCorrectionDatum stdBa = new StandardMatrixCorrectionDatum(Composition.parse("BaF2"),
				UncertainValue.valueOf(15.0, 0.1), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 4.893));
		StandardMatrixCorrectionDatum stdO = new StandardMatrixCorrectionDatum(Composition.parse("SiO2"),
				UncertainValue.valueOf(15.0, 0.25), UncertainValue.toRadians(40.0, 0.5),
				MatrixCorrectionDatum.roughness(10.0, 2.196));

		KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
		KRatioLabel krlSi = new KRatioLabel(unk, stdSi, CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlTi = new KRatioLabel(unk, stdTi, CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
				CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
		KRatioLabel krlBa = new KRatioLabel(unk, stdBa, CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5),
				KRatioLabel.Method.Calculated);
		KRatioLabel krlO = new KRatioLabel(unk, stdO, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
				KRatioLabel.Method.Calculated);

		Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
				Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
		Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
		UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
		Report r = new Report("K240 at 15.0 keV");
		r.addHeader("K240 at 15.0 keV");
		r.addHTML(xpp.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
		r.addImage(xpp.asCovarianceBitmap(16, L2C, L2C, true), "Covariances");

		// r.addHTML(xpp.toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));

		// List<ZAFMultiLineLabel> zafs =
		// kratios.stream().flatMap(krl->MatrixCorrectionModel2.zafLabel(krl)).collect(Collectors.toList());
		ArrayList<EPMALabel> labels = new ArrayList<>();
		labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
		// labels.addAll(xpp.getLabels(MaterialMAC.class));
		labels.addAll(xpp.getLabels(ZAFMultiLineLabel.class));
		UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
		r.addHTML(zafMac.toHTML(Mode.VERBOSE, new BasicNumberFormat("0.00E0")));
		r.addImage(zafMac.asCovarianceBitmap(64, L2C, L2C, true), "MACS -> ZAF");

		File file = new File(System.getProperty("user.home") + "\\Desktop", "K240_ZAF_Elm.csv");
		try (PrintWriter pw = new PrintWriter(file, Charset.forName("UTF-8"))) {
			pw.write(zafMac.toCSV());
		}
		if (false) {
			final double range = 0.0015;
			final DataFrame<Double, String> prz = ((XPPMatrixCorrection2) xpp.getFunction())
					.computePhiRhoZCurve(xpp.getValueMap(), range, range / 100, 0.1);
			prz.write().csv(new File(System.getProperty("user.home") + "\\Desktop", "prz_K240.csv")).apply(options -> {
				options.setSeparator(",");
				options.setIncludeRowHeader(true);
				options.setIncludeColumnHeader(true);
			});
			Chart.create().asHtml().withLinePlot(prz, chart -> {
				chart.plot().axes().domain().label().withText("\\u03C1\\u00B7z  (g/cm\\u00B2)");
				chart.plot().axes().range(0).label().withText("	\\u03D5(\\u03C1\\u00B7z)");
				chart.legend().on().bottom();
				chart.show();
			});
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static void K240atMultikeV() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();

		Report r = new Report("K240 at Multi keV");
		for (double e0 : new double[] { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0 }) {

			UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
					UncertainValue.valueOf(e0, 0.1), UncertainValue.toRadians(40.0, 0.5),
					MatrixCorrectionDatum.roughness(1.0, 5.0));
			StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 1.738));
			StandardMatrixCorrectionDatum stdSi = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Silicon), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 2.3290));
			StandardMatrixCorrectionDatum stdTi = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Titanium), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 4.506));
			StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Zinc), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 7.14));
			StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 6.52));
			StandardMatrixCorrectionDatum stdBa = new StandardMatrixCorrectionDatum(Composition.parse("BaF2"),
					UncertainValue.valueOf(e0, 0.1), UncertainValue.toRadians(40.0, 0.1),
					MatrixCorrectionDatum.roughness(1.0, 4.893));
			StandardMatrixCorrectionDatum stdO = new StandardMatrixCorrectionDatum(Composition.parse("SiO2"),
					UncertainValue.valueOf(e0, 0.1), UncertainValue.toRadians(40.0, 0.1),
					MatrixCorrectionDatum.roughness(1.0, 2.196));

			KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
					CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlSi = new KRatioLabel(unk, stdSi,
					CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlTi = new KRatioLabel(unk, stdTi,
					CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
					KRatioLabel.Method.Calculated);
			KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
					CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
			KRatioLabel krlBa = new KRatioLabel(unk, stdBa,
					CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
			KRatioLabel krlO = new KRatioLabel(unk, stdO, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
					KRatioLabel.Method.Calculated);

			Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
					Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
			Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
			UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
			r.addHeader("K240 at " + e0 + " keV");
			ArrayList<EPMALabel> labels = new ArrayList<>();
			for(KRatioLabel krl : kratios)
				labels.add(MatrixCorrectionModel2.zaLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet().getBrightest()));
			r.addHTML(xpp.extract(labels).toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
			labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
			labels.addAll(xpp.getLabels(MaterialMAC.class));
			UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
			r.addImage(zafMac.asCovarianceBitmap(8, L2C, L2C, true), "MACS -> ZAF");
		}
		r.inBrowser(Mode.NORMAL);
	}
	
	public static void K240atMultikeVSimilar() throws Exception {
		final Composition unkComp = CompositionFactory.instance().findComposition("K240").getObject();
		final Composition benitoite = Composition.parse("BaTi(Si3O9)");
		Report r = new Report("K240 at Multi keV");
		for (double e0 : new double[] { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0 }) {

			UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unkComp.getMaterial(),
					UncertainValue.valueOf(e0, 0.1), UncertainValue.toRadians(40.0, 0.1),
					MatrixCorrectionDatum.roughness(1.0, 5.0));
			StandardMatrixCorrectionDatum stdMg = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Magnesium), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 1.738));
			StandardMatrixCorrectionDatum stdMulti = new StandardMatrixCorrectionDatum(
					benitoite, UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 3.65));
			StandardMatrixCorrectionDatum stdZn = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Zinc), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 7.14));
			StandardMatrixCorrectionDatum stdZr = new StandardMatrixCorrectionDatum(
					Composition.pureElement(Element.Zirconium), UncertainValue.valueOf(e0, 0.1),
					UncertainValue.toRadians(40.0, 0.1), MatrixCorrectionDatum.roughness(1.0, 6.52));

			KRatioLabel krlMg = new KRatioLabel(unk, stdMg,
					CharacteristicXRay.create(Element.Magnesium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlSi = new KRatioLabel(unk, stdMulti,
					CharacteristicXRay.create(Element.Silicon, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlTi = new KRatioLabel(unk, stdMulti,
					CharacteristicXRay.create(Element.Titanium, XRayTransition.KL3), KRatioLabel.Method.Calculated);
			KRatioLabel krlZn = new KRatioLabel(unk, stdZn, CharacteristicXRay.create(Element.Zinc, XRayTransition.KL3),
					KRatioLabel.Method.Calculated);
			KRatioLabel krlZr = new KRatioLabel(unk, stdZr,
					CharacteristicXRay.create(Element.Zirconium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
			KRatioLabel krlBa = new KRatioLabel(unk, stdMulti,
					CharacteristicXRay.create(Element.Barium, XRayTransition.L3M5), KRatioLabel.Method.Calculated);
			KRatioLabel krlO = new KRatioLabel(unk, stdMulti, CharacteristicXRay.create(Element.Oxygen, XRayTransition.KL3),
					KRatioLabel.Method.Calculated);

			Set<KRatioLabel> kratios = new HashSet<KRatioLabel>(
					Arrays.asList(krlMg, krlSi, krlTi, krlZn, krlZr, krlBa, krlO));
			Map<MassFraction, Double> estUnknown = unkComp.getValueMap(MassFraction.class);
			UncertainValuesCalculator<EPMALabel> xpp = XPPMatrixCorrection2.buildAnalytical(kratios, estUnknown, false);
			r.addHeader("K240 at " + e0 + " keV (Benitoite)");
			ArrayList<EPMALabel> labels = new ArrayList<>();
			for(KRatioLabel krl : kratios)
				labels.add(MatrixCorrectionModel2.zaLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet().getBrightest()));
			r.addHTML(xpp.extract(labels).toHTML(Mode.NORMAL, new BasicNumberFormat("0.00E0")));
			labels.addAll(xpp.getLabels(ElementalMAC.ElementMAC.class));
			labels.addAll(xpp.getLabels(MaterialMAC.class));
			UncertainValuesBase<EPMALabel> zafMac = xpp.extract(labels);
			r.addImage(zafMac.asCovarianceBitmap(8, L2C, L2C, true), "MACS -> ZAF");
		}
		r.inBrowser(Mode.NORMAL);
	}

	
}
