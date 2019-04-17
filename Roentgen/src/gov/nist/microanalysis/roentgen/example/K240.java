package gov.nist.microanalysis.roentgen.example;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.linear.RealVector;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValueEx;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.DataStore.CompositionFactory;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioCorrectionModel3;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.Layer;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * The K240 class serves as an example of how to use
 * {@link KRatioCorrectionModel3} to calculate the matrix correction with
 * associated uncertainties for K240 glass (O, Mg, Si, Ti, Zn, Zr, Ba). The
 * results are output as HTML pages.
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class K240 {

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

	public static void main(
			final String[] args
	) throws Exception {
		final K240 k240 = new K240();
		try {
			// k240.simple();
			// k240.benitoite();
			// k240.elements();
			// k240.k240();
			k240.surfaceRoughness();
			k240.surfaceCoating();
		} catch (IOException | ArgumentException | ParseException e) {
			e.printStackTrace();
		}
	}

	private final UncertainValue mTOA;

	private final double mRoughness;

	public K240() {
		mTOA = UncertainValue.toRadians(40.0, 0.5);
		mRoughness = MatrixCorrectionDatum.roughness(10.0, 3.6);
	}

	public void benitoite() throws Exception {
		// K240 using Benitoite, Mg, Zn and Zr as standards

		final double[] dkO = { 0.002175978085581, 0.002197549502352, 0.002222917439207, 0.002259155722945,
				0.002296599040709, 0.002348395478965, 0.002387912718087, 0.002436731785879, 0.002482245465167,
				0.002533966940594, 0.002580622126306, 0.00263742877099, 0.002689424504791, 0.002749521302058,
				0.002802880345387, 0.002851611085874, 0.002911168928636, 0.002954109640922, 0.002999994065294 };
		final double[] dkMg = { 0.004962906114096, 0.005001404889014, 0.004924646968079, 0.00491029272899,
				0.004888829359764, 0.004858183896549, 0.004937708903069, 0.004937272359369, 0.00497256515775,
				0.004936277149524, 0.005032282567414, 0.005146316851665, 0.005192328070361, 0.005230358335188,
				0.005201109570042, 0.00534889375152, 0.005327245053272, 0.005479085928104, 0.005534800055348 };
		final double[] dkSi = { 0.001540612407138, 0.001481555994616, 0.001425181252983, 0.001389009082235,
				0.001361726594707, 0.001334054483399, 0.001313935477242, 0.001298536909002, 0.001284616309707,
				0.001280450577798, 0.001268442997156, 0.00126374058475, 0.001265849189384, 0.001263734641962,
				0.001261893758378, 0.001263952257026, 0.001269622687434, 0.001275024133209, 0.001275896158657 };
		final double[] dkTi = { 0.01404057687607, 0.01229778662795, 0.011155590526459, 0.010357134615828,
				0.009809937458651, 0.009039481437831, 0.008514949836098, 0.008046240928013, 0.007626096177276,
				0.00729970338153, 0.006968410778483, 0.006682418598052, 0.006411703616743, 0.00626132149355,
				0.006114523370149, 0.005900629072074, 0.005751532974722, 0.005623044191438, 0.005445931668182 };
		final double[] dkZn = { 0.078306155151287, 0.054715770563326, 0.042693216334796, 0.034162025607166,
				0.028721862336203, 0.024233475201946, 0.021164636450853, 0.018596859827944, 0.016726233023588,
				0.015237771546719, 0.013915939819911, 0.012983363309353, 0.012256379344987, 0.01127724633209,
				0.011014041461235, 0.010472874963737, 0.009827448291624, 0.009332191780822, 0.009008492081708 };
		final double[] dkZr = { 0.006400090620752, 0.006149588746253, 0.005870045333023, 0.005768652979469,
				0.005500637755102, 0.005427769698708, 0.005344231373782, 0.005287733737021, 0.005202011444425,
				0.00519393669193, 0.005174378846592, 0.005165240966671, 0.005112859902926, 0.00507493269162,
				0.005097188010633, 0.005160348386447, 0.005107329776744, 0.005119366124717, 0.00516414607156 };
		final double[] dkBa = { 0.008182058447819, 0.007360938324068, 0.00668716148458, 0.006090709449536,
				0.005639271864227, 0.005289641813193, 0.004971912248055, 0.004714049877284, 0.004507028922697,
				0.004308293402505, 0.004128654614498, 0.0039938302209, 0.003861512410919, 0.003728404071926,
				0.003611433547977, 0.003519063068446, 0.003433831396981, 0.003344645059098, 0.003344645059098 };

		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

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
		final Composition benitoite = CompositionFactory.instance().findComposition("Benitoite").getObject();

		final Map<Integer, Map<EPMALabel, UncertainValueEx<String>>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValuesBase<EPMALabel>> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<EPMALabel> outLabels = null;
		List<String> inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum benitoiteMcd = new StandardMatrixCorrectionDatum(//
						benitoite, e0, mTOA, mRoughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(//
						unk.getMaterial(), e0, mTOA, mRoughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, tiTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, benitoiteMcd, baTrs, Method.Measured));

				final KRatioCorrectionModel3 cfk = KRatioCorrectionModel3.buildXPPModel(lkr, null);
				final XPPMatrixCorrection2 mcm = (XPPMatrixCorrection2) cfk.getModel();
				final UncertainValues<EPMALabel> input = mcm.buildInput(unk.getMaterial());
				mcm.addConstraints(mcm.buildConstraints(input));
				mcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));
				// Calculate the optimal k-ratios
				final RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				final UncertainValues<KRatioLabel> krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(),
						Method.Measured);
				for (final KRatioLabel label : krs.getLabels()) {
					final KRatioLabel krl = label;
					final double v = krs.getEntry(krl);
					double dk = 0.0;
					switch (krl.getElement()) {
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

				final UncertainValuesBase<EPMALabel> msInp = cfk.buildInput(krs);

				final Set<EPMALabel> finalOutputs = new HashSet<>();
				for (final EPMALabel output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final UncertainValuesCalculator<EPMALabel> res = UncertainValuesBase.propagateAnalytical(cfk, msInp);
				final Map<String, Collection<? extends EPMALabel>> labels = extractLabelBlocks(msInp);
				resVals.put(ie0, res);
				if ((ie0 - MIN_E) % 5 == 0) {
					initReport.addHeader("E<sub>0</sub> = " + ie0);
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, res.getOutputValues(labels, 0.0));
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = new ArrayList<>(labels.keySet());
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
				final EPMALabel outTag = outLabels.get(oi);
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final String inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						boolean addRow = false;
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
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
						final UncertainValuesBase<EPMALabel> uvs = resVals.get(ie0);
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

	public void surfaceRoughness() throws Exception {
		// K240 using Benitoite, Zircon, Mg, Zn as standards

		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mg = Composition.parse("Mg");
		final Composition zn = Composition.parse("Zn");
		final Composition benitoite = CompositionFactory.instance().findComposition("Benitoite").getObject();
		final Composition zircon = CompositionFactory.instance().findComposition("Zircon").getObject();

		final UncertainValue e0 = new UncertainValue(15.0, 0.1);

		MatrixCorrectionDatum.roughness(10.0, 3.6);

		Report report = new Report("K240 roughness");
		report.addHeader("K240 roughness");
		try {
			for (double roughNm : Arrays.asList(1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0)) {
				final double roughness = MatrixCorrectionDatum.roughness(roughNm, 3.6);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, mTOA, 0.0);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, mTOA, 0.0);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zircon, e0, mTOA, 0.0);
				final StandardMatrixCorrectionDatum benitoiteMcd = new StandardMatrixCorrectionDatum(//
						benitoite, e0, mTOA, 0.0);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(//
						unk.getMaterial(), e0, mTOA, roughness);
				final Set<KRatioLabel> skr = new HashSet<>();
				skr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, benitoiteMcd, siTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, benitoiteMcd, oTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, benitoiteMcd, tiTrs, Method.Measured));
				skr.add(new KRatioLabel(unkMcd, benitoiteMcd, baTrs, Method.Measured));

				Composition res = KRatioCorrectionModel3.roundTripXPP(unk, skr);

				report.addSubHeader("Roughness = " + roughNm + " nm");
				report.add(res, Mode.NORMAL);
			}
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	public void surfaceCoating() throws Exception {
		// K240 using Benitoite, Zircon, Mg, Zn as standards

		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		final Composition mg = Composition.parse("Mg");
		final Composition zn = Composition.parse("Zn");
		final Composition benitoite = CompositionFactory.instance().findComposition("Benitoite").getObject();
		final Composition zircon = CompositionFactory.instance().findComposition("Zircon").getObject();

		final UncertainValue e0 = new UncertainValue(15.0, 0.1);

		MatrixCorrectionDatum.roughness(10.0, 3.6);

		Report report = new Report("K240 coating");
		report.addHeader("K240 coating");
		try {
			for (boolean same : Arrays.asList(true, false)) {
				for (double coatingDelta : Arrays.asList(0.0, 1.0, 2.0, 4.0, 8.0)) {

					final Layer coating = Layer.carbonCoating(UncertainValue.valueOf(10.0, coatingDelta));
					final Layer unkCoating = same ? coating : Layer.carbonCoating(UncertainValue.valueOf(10.0, coatingDelta));
					final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, mTOA, 0.0,
							coating);
					final StandardMatrixCorrectionDatum znMcd = //
							new StandardMatrixCorrectionDatum(zn, e0, mTOA, 0.0, coating);
					final StandardMatrixCorrectionDatum zrMcd = //
							new StandardMatrixCorrectionDatum(zircon, e0, mTOA, 0.0, coating);
					final StandardMatrixCorrectionDatum benitoiteMcd = //
							new StandardMatrixCorrectionDatum(benitoite, e0, mTOA, 0.0, coating);
					final UnknownMatrixCorrectionDatum unkMcd = //
							new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0, mTOA, 0.0, unkCoating);
					final Set<KRatioLabel> skr = new HashSet<>();
					skr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, benitoiteMcd, siTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, benitoiteMcd, oTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, benitoiteMcd, tiTrs, Method.Measured));
					skr.add(new KRatioLabel(unkMcd, benitoiteMcd, baTrs, Method.Measured));

					Composition res = KRatioCorrectionModel3.roundTripXPP(unk, skr);
					
					report.addSubHeader("Coating: same="+Boolean.valueOf(same)+"  delta= "+coatingDelta+" nm");
					report.add(res, Mode.NORMAL);
				}
			}
		} finally {
			report.inBrowser(Mode.NORMAL);
		}
	}

	public void elements() throws Exception {

		final double[] dkO = { 0.000676038081616, 0.000681498814916, 0.000687566080938, 0.000702437667351,
				0.000712176738527, 0.000725866898102, 0.000737679495025, 0.000753980573912, 0.000766820966007,
				0.00078475226352, 0.000799975385373, 0.00081455044318, 0.000834733831206, 0.000844741196755,
				0.000861408882083, 0.000874772066433, 0.000899668009833, 0.000911921418429, 0.000933457794373 };
		final double[] dkMg = { 0.002244376011274, 0.002226923505178, 0.002176470588235, 0.002210154079313,
				0.002134614101794, 0.002207819955844, 0.002190994258084, 0.002172164119067, 0.002215216835648,
				0.002256725040621, 0.002288111354753, 0.002311325494925, 0.002336448598131, 0.002330743618202,
				0.002329373398556, 0.002445884798826, 0.002435273006921, 0.002404488378306, 0.00250626566416 };
		final double[] dkSi = { 0.000691058714631, 0.000662890596798, 0.000639925931977, 0.000620146298149,
				0.000606321818089, 0.000598530191681, 0.000590272312293, 0.000579653350968, 0.000574936165176,
				0.000570165435694, 0.00056422623652, 0.000567188164674, 0.000568521549907, 0.00056847020607,
				0.000567978627175, 0.000566393272991, 0.000563710568446, 0.000570663250451, 0.000566094549834 };
		final double[] dkTi = { 0.005341630546213, 0.004820421056476, 0.00432137993645, 0.0039335270024,
				0.003669155364086, 0.003414776669954, 0.003227817921725, 0.003052929709052, 0.002911667068489,
				0.002777984935491, 0.002655545789464, 0.002564102564103, 0.00248077869939, 0.002401829965688,
				0.00233318671231, 0.002260016854363, 0.002200559791526, 0.002148789128289, 0.002088253088468 };
		final double[] dkZn = { 0.041532110864621, 0.027013681122375, 0.019346382881728, 0.014706693891066,
				0.012368339310423, 0.010599807790152, 0.009270809707417, 0.008209207853758, 0.007409911813597,
				0.006772773450728, 0.006207516431661, 0.005807036428531, 0.005411651963133, 0.005077860528098,
				0.004855894198635, 0.004595353586929, 0.004405286343612, 0.004198320671731, 0.004030632809351 };
		final double[] dkZr = { 0.002869548801208, 0.002735323935041, 0.00261592419633, 0.00254935673208,
				0.002473306388369, 0.002429533241368, 0.002382196217741, 0.002353040934433, 0.002320839840361,
				0.002293173399189, 0.002297556031085, 0.00226819712695, 0.002273952692407, 0.002273107936736,
				0.00225181598063, 0.002272278205888, 0.002260681721132, 0.002269829125223, 0.00227855311877 };
		final double[] dkBa = { 0.002693749587821, 0.002392916965781, 0.002186211430051, 0.002005810952318,
				0.001856345192504, 0.001742025394859, 0.001638975017612, 0.001556868551731, 0.001481264830339,
				0.001414937428843, 0.001363339247362, 0.001310773819725, 0.001266289082286, 0.001224543880251,
				0.001192882769673, 0.00115847588017, 0.001128023980621, 0.00110585296185, 0.00110585296185 };

		// K240 using elements
		final Composition unk = CompositionFactory.instance().findComposition("K412").getObject();

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

		final Map<Integer, Map<EPMALabel, UncertainValueEx<String>>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValuesBase<EPMALabel>> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<? extends EPMALabel> outLabels = null;
		List<String> inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(mg, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum baMcd = new StandardMatrixCorrectionDatum(ba, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum tiMcd = new StandardMatrixCorrectionDatum(ti, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum siMcd = new StandardMatrixCorrectionDatum(si, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum oMcd = new StandardMatrixCorrectionDatum(o, e0, mTOA, mRoughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0,
						mTOA, mRoughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, baMcd, baTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, siMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, oMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured));

				final KRatioCorrectionModel3 cfk = KRatioCorrectionModel3.buildXPPModel(lkr, null);
				final XPPMatrixCorrection2 mcm = (XPPMatrixCorrection2) cfk.getModel();
				final UncertainValues<EPMALabel> input = mcm.buildInput(unk.getMaterial());
				mcm.addConstraints(mcm.buildConstraints(input));
				mcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));
				// Calculate the optimal k-ratios
				final RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				final UncertainValues<KRatioLabel> krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(),
						Method.Measured);
				for (final KRatioLabel label : krs.getLabels()) {
					final KRatioLabel krl = label;
					final double v = krs.getEntry(krl);
					double dk = 0.0;
					switch (krl.getElement()) {
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

				final UncertainValuesBase<EPMALabel> msInp = cfk.buildInput(krs);

				final Set<EPMALabel> finalOutputs = new HashSet<>();
				for (final EPMALabel output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final UncertainValuesCalculator<EPMALabel> res = new UncertainValuesCalculator<>(cfk, msInp);
				final Map<String, Collection<? extends EPMALabel>> labels = extractLabelBlocks(msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, res.getOutputValues(labels, 0.0));
				resVals.put(ie0, res);
				if (ie0 == MIN_E)
					for (final EPMALabel label : cfk.getOutputLabels())
						System.out.println(label);
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = new ArrayList<String>(labels.keySet());
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
				final EPMALabel outTag = outLabels.get(oi);
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final String inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						boolean addRow = false;
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
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
						final UncertainValuesBase<EPMALabel> uvs = resVals.get(ie0);
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

	public List<? extends EPMALabel> filter(
			final List<? extends EPMALabel> labels, final String name
			) {
		final List<EPMALabel> res = new ArrayList<>();
		for (final EPMALabel label : labels)
			if (label.toString().startsWith(name))
				res.add(label);
		return res;
	}

	public <H> UncertainValueEx<String> getByMFT(
			final Map<H, UncertainValueEx<String>> oVals, final MaterialLabel.MassFraction mft
			) {
		for (final Map.Entry<H, UncertainValueEx<String>> me : oVals.entrySet()) {
			if (me.getKey() instanceof MaterialLabel.MassFraction) {
				final MaterialLabel.MassFraction mft2 = (MaterialLabel.MassFraction) me.getKey();
				if (mft2.toString().equals(mft.toString()))
					return me.getValue();
			}
		}
		return null;
	}

	public double getComponentByName(
			final UncertainValueEx<String> uv, final String name
			) {
		for (final Map.Entry<String, Double> me : uv.getComponents().entrySet())
			if (me.getKey().toString().equals(name))
				return me.getValue().doubleValue();
		return 0.0;
	}

	public void k240() throws Exception {

		final double[] dkO = { 0.00213730595774, 0.002153401156992, 0.002183789128298, 0.002214099174781,
				0.002250268857424, 0.002293120836805, 0.002328292172395, 0.002377052394699, 0.002428879224355,
				0.002476552320898, 0.002531283217462, 0.002573524601466, 0.002630804180641, 0.002677525966007,
				0.0027308818302, 0.00278074566441, 0.002831060723404, 0.002884710931689, 0.002929106637445 };
		final double[] dkMg = { 0.007139753389642, 0.007014356496, 0.006946733436851, 0.006903426435649,
				0.006873264680701, 0.006869265257026, 0.006891680540841, 0.006907549970205, 0.006973728024607,
				0.007030108361156, 0.007117328683991, 0.007157737087594, 0.007296049071621, 0.007381388667523,
				0.007473734887447, 0.007535762942779, 0.007651372243953, 0.007796973859683, 0.007833621005858 };
		final double[] dkSi = { 0.002177443609023, 0.002089698766817, 0.002021163957383, 0.001967589129893,
				0.001925435159412, 0.001888781516296, 0.001856676609086, 0.001838303135489, 0.001816653561454,
				0.00180758528516, 0.001796799737775, 0.001794689882398, 0.001785956163806, 0.001784417615563,
				0.001785857131429, 0.001791426488176, 0.001793583888538, 0.001800684260019, 0.001807426059161 };
		final double[] dkTi = { 0.016961989316451, 0.01495454603787, 0.013410677710168, 0.012471673946599,
				0.011525439427514, 0.010850440764694, 0.010307448054625, 0.009710203159639, 0.00916241225594,
				0.008845861616106, 0.008482746682983, 0.008073816600677, 0.007818605142337, 0.007617207699947,
				0.007384478193792, 0.007160468876314, 0.006922178237133, 0.006806180207432, 0.006593433317357 };
		final double[] dkZn = { 0.096655960199049, 0.058637910729642, 0.04396613701324, 0.033422593126996,
				0.028349822589152, 0.024607737034644, 0.021091390843811, 0.018260884560413, 0.017092353256542,
				0.015352806136363, 0.013921969582544, 0.013135552177497, 0.012197495295152, 0.011752344708962,
				0.010745531746816, 0.010335373998332, 0.00975601530013, 0.009284794272665, 0.009033216679099 };
		final double[] dkZr = { 0.00904643210491, 0.008630414259884, 0.008346199769945, 0.008125177155089,
				0.007912995433412, 0.007734786681564, 0.007500526299961, 0.007436826229482, 0.007393904083465,
				0.007271531003301, 0.007272100775691, 0.007190120512302, 0.007161698428888, 0.007138116165722,
				0.00713058742647, 0.007094644626686, 0.007149815409671, 0.007211895463139, 0.007238790132689 };
		final double[] dkBa = { 0.008484743161996, 0.007590751011964, 0.006900126919433, 0.006361762365905,
				0.005895470137267, 0.005487648516967, 0.005147871821179, 0.004913872616228, 0.004709654842562,
				0.00447519440929, 0.004315111303428, 0.004159453334672, 0.003995562248976, 0.003883656254936,
				0.003759892408307, 0.003665472596181, 0.003574832224132, 0.003492056863254, 0.003407377435567 };
		// K240 using simple elements
		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

		final ElementXRaySet mgTrs = ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1);
		final ElementXRaySet baTrs = ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1);
		final ElementXRaySet tiTrs = ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1);
		final ElementXRaySet siTrs = ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1);
		final ElementXRaySet oTrs = ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1);
		final ElementXRaySet znTrs = ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1);
		final ElementXRaySet zrTrs = ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1);

		// Must be a different instance of K240
		final Composition std = Composition.massFraction("K240s", buildK240());

		final Map<Integer, Map<EPMALabel, UncertainValueEx<String>>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValuesBase<EPMALabel>> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<? extends EPMALabel> outLabels = null;
		List<String> inLabels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum baMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum tiMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum siMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum oMcd = new StandardMatrixCorrectionDatum(std, e0, mTOA, mRoughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0,
						mTOA, mRoughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, baMcd, baTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, siMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, oMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured));

				final KRatioCorrectionModel3 cfk = KRatioCorrectionModel3.buildXPPModel(lkr, null);
				final XPPMatrixCorrection2 mcm = (XPPMatrixCorrection2) cfk.getModel();
				final UncertainValues<EPMALabel> input = mcm.buildInput(unk.getMaterial());
				mcm.addConstraints(mcm.buildConstraints(input));

				// Calculate the optimal k-ratios
				final RealVector calculated = mcm.optimized(input.extractValues(mcm.getInputLabels()));
				final UncertainValues<KRatioLabel> krs = KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(),
						Method.Measured);
				for (final KRatioLabel label : krs.getLabels()) {
					final KRatioLabel krl = label;
					final double v = krs.getEntry(krl);
					double dk = 0.0;
					switch (krl.getElement()) {
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

				final UncertainValuesBase<EPMALabel> msInp = cfk.buildInput(krs);

				final Set<EPMALabel> finalOutputs = new HashSet<>();
				for (final EPMALabel output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				// cfk.trimOutputs(finalOutputs);
				final UncertainValuesCalculator<EPMALabel> res = UncertainValuesBase.propagateAnalytical(cfk, msInp);
				final Map<String, Collection<? extends EPMALabel>> labels = extractLabelBlocks(msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, res.getOutputValues(labels, 0.0));
				resVals.put(ie0, res);
				if (ie0 == MIN_E)
					for (final EPMALabel label : cfk.getOutputLabels())
						System.out.println(label);
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = new ArrayList<>(labels.keySet());
				}
			}
			initReport.inBrowser(Mode.NORMAL);
		}
		final Report report = new Report("K-Ratio (7)");
		try {
			report.addHeader("K240 using K240");
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
				final EPMALabel outTag = outLabels.get(oi);
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final String inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						boolean addRow = false;
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
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
						final UncertainValuesBase<EPMALabel> uvs = resVals.get(ie0);
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

	public void simple() throws Exception {

		final double[] dkO = { 0.000675032606014, 0.000681693829774, 0.000690067208914, 0.000700899426437,
				0.00071072153042, 0.000726035676794, 0.00073750326059, 0.000754022397678, 0.000765121084199,
				0.000783567577461, 0.000797175937492, 0.000815312565192, 0.000832393534182, 0.000846830898216,
				0.000862673828776, 0.000879601860387, 0.000896661232938, 0.000909711166705, 0.000927583463386 };
		final double[] dkMg = { 0.002241100745296, 0.002236886254334, 0.002191423833215, 0.002199736031676,
				0.00214563497385, 0.002210969260395, 0.002202141392665, 0.002175489485134, 0.002208818282219,
				0.002236536052961, 0.002284408909195, 0.002317146886964, 0.002335704427222, 0.00234794275492,
				0.002331002331002, 0.002442002442002, 0.002431844361961, 0.00240706071142, 0.002504522053708 };
		final double[] dkSi = { 0.000688465475405, 0.000661280430899, 0.000639096306897, 0.000620824466077,
				0.000609424812214, 0.000597065847834, 0.000587809730512, 0.000581183264395, 0.000576244783718,
				0.000569625374221, 0.000570322032999, 0.000566658749279, 0.00056578362768, 0.000563290466575,
				0.000563816022495, 0.000564109370852, 0.000568579755548, 0.000568044430286, 0.000571345832502 };
		final double[] dkTi = { 0.005281752970986, 0.004749824580342, 0.004244239960054, 0.003952999202147,
				0.00371767737186, 0.003446253116293, 0.00324639392039, 0.00307191235797, 0.002925883821912,
				0.002776173352462, 0.002691562047093, 0.002577319587629, 0.002474172285493, 0.002394526795895,
				0.002342504944222, 0.00226809672087, 0.002194586686174, 0.002151454654701, 0.002115863095821 };
		final double[] dkZn = { 0.039833243679398, 0.026112941965034, 0.019295181930091, 0.015049930691109,
				0.012277350115188, 0.010704193189279, 0.009309078486634, 0.00810368273527, 0.007467560123797,
				0.006749535969402, 0.006271994551027, 0.005840788648926, 0.005482894511394, 0.005169591315087,
				0.004826802952868, 0.004604868003318, 0.004428978171465, 0.004242056964765, 0.004036326942482 };
		final double[] dkZr = { 0.002860422665086, 0.002754290336871, 0.002613746369797, 0.002547795860325,
				0.002473853580048, 0.002440174708307, 0.002379411826094, 0.002364922447588, 0.002325733040624,
				0.002292314135202, 0.002295591114712, 0.002277642295127, 0.002263340877803, 0.0022676278226,
				0.002261673151751, 0.002282481951026, 0.002254226675016, 0.002271915045694, 0.002290831467694 };
		final double[] dkBa = { 0.002703601371197, 0.00240966713384, 0.00218648065737, 0.002003608469243,
				0.0018559493269, 0.001736046477452, 0.00163714575891, 0.001552033183017, 0.00147972078822,
				0.001416344696043, 0.001356934108761, 0.001310682551123, 0.001266480864021, 0.001224589284844,
				0.001191258386089, 0.001156587094513, 0.001129237081176, 0.001103578440835, 0.001103578440835 };

		// K240 using MgO, BaSi2O5, Zn, Zr, SiO2 and Ti as standards
		final Composition unk = CompositionFactory.instance().findComposition("K240").getObject();

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
		final Composition si = Composition.parse("Si");

		final Map<Integer, Map<EPMALabel, UncertainValueEx<String>>> outVals = new TreeMap<>();
		final Map<Integer, UncertainValuesBase<EPMALabel>> resVals = new TreeMap<>();
		final int MIN_E = 12, MAX_E = 31;
		List<? extends EPMALabel> outLabels = null;
		List<String> inLabels = null;
		Map<String, Collection<? extends EPMALabel>> labels = null;
		{
			final Report initReport = new Report("K-Ratio - Initialization");
			for (int ie0 = MIN_E; ie0 < MAX_E; ie0++) {
				initReport.addHeader("E<sub>0</sub> = " + ie0);
				final UncertainValue e0 = new UncertainValue(ie0, 0.1);
				final StandardMatrixCorrectionDatum mgoMcd = new StandardMatrixCorrectionDatum(mgo, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum basi2o5Mcd = new StandardMatrixCorrectionDatum(basi2o5, e0, mTOA,
						mRoughness);
				final StandardMatrixCorrectionDatum znMcd = new StandardMatrixCorrectionDatum(zn, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum zrMcd = new StandardMatrixCorrectionDatum(zr, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum siMcd = new StandardMatrixCorrectionDatum(si, e0, mTOA, mRoughness);
				final StandardMatrixCorrectionDatum tiMcd = new StandardMatrixCorrectionDatum(ti, e0, mTOA, mRoughness);
				final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0,
						mTOA, mRoughness);

				final Set<KRatioLabel> lkr = new HashSet<>();
				lkr.add(new KRatioLabel(unkMcd, mgoMcd, mgTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, basi2o5Mcd, baTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, siMcd, siTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, mgoMcd, oTrs, Method.Measured));
				lkr.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Measured));

				final List<EPMALabel> defOut = KRatioCorrectionModel3.buildDefaultOutputs(lkr);
				defOut.add(new KRatioLabel(unkMcd, mgoMcd, mgTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, basi2o5Mcd, baTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, znMcd, znTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, zrMcd, zrTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, siMcd, siTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, mgoMcd, oTrs, Method.Calculated));
				defOut.add(new KRatioLabel(unkMcd, tiMcd, tiTrs, Method.Calculated));

				final KRatioCorrectionModel3 cfk = KRatioCorrectionModel3.buildXPPModel(lkr, null, defOut);
				final XPPMatrixCorrection2 mcm = (XPPMatrixCorrection2) cfk.getModel();
				final UncertainValues<EPMALabel> input = mcm.buildInput(unk.getMaterial());
				mcm.addConstraints(mcm.buildConstraints(input));
				mcm.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				// Calculate the optimal k-ratios
				final RealVector inputs = input.extractValues(mcm.getInputLabels());
				final RealVector calculated = mcm.optimized(inputs);
				final UncertainValues<KRatioLabel> krs = UncertainValues.asUncertainValues(//
						KRatioLabel.extractKRatios(calculated, mcm.getOutputLabels(), Method.Calculated));
				final Map<KRatioLabel, UncertainValue> measKrs = new HashMap<>();
				for (final KRatioLabel label : krs.getLabels()) {
					final KRatioLabel krl = label;
					final double v = krs.getEntry(krl);
					double dk = 0.0;
					switch (krl.getElement()) {
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
					measKrs.put(krl.as(Method.Measured), new UncertainValue(v, dk));
				}

				final UncertainValues<KRatioLabel> uncKrs = new UncertainValues<KRatioLabel>(measKrs);
				final UncertainValuesBase<EPMALabel> msInp = cfk.buildInput(uncKrs);

				final Set<EPMALabel> finalOutputs = new HashSet<>();
				for (final EPMALabel output : cfk.getOutputLabels())
					if (output instanceof MaterialLabel.MassFraction)
						finalOutputs.add(output);
				final UncertainValuesCalculator<EPMALabel> res = UncertainValuesBase.propagateAnalytical(cfk, msInp);
				labels = extractLabelBlocks(msInp);
				if (ie0 % 5 == 0) {
					initReport.addSubHeader("Inputs");
					initReport.add(msInp.sort(), Mode.VERBOSE);
					initReport.addSubHeader("Output");
					initReport.addHTML(res.toHTML(Mode.NORMAL, new BasicNumberFormat("0.0E0")));
				}
				outVals.put(ie0, res.getOutputValues(labels, 0.0));
				resVals.put(ie0, res);
				if (ie0 == MIN_E)
					for (final EPMALabel label : cfk.getOutputLabels())
						System.out.println(label);
				if (ie0 == MIN_E) {
					outLabels = cfk.getOutputLabels();
					inLabels = new ArrayList<String>(labels.keySet());
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
				final EPMALabel outTag = outLabels.get(oi);
				if (outTag instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) outTag;
					for (int ii = 0; ii < labels.size(); ++ii) {
						final String inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						row.add(Table.td(mft));
						row.add(Table.td(mft.getElement()));
						row.add(Table.td(inTag));
						boolean addRow = false;
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int ie0 = MIN_E; ie0 < MAX_E; ++ie0) {
							final Map<EPMALabel, UncertainValueEx<String>> oVals = outVals.get(ie0);
							final UncertainValueEx<String> tmp = getByMFT(oVals, mft);
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
						final UncertainValuesBase<EPMALabel> uvs = resVals.get(ie0);
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

	private Map<String, Collection<? extends EPMALabel>> extractLabelBlocks(
			final UncertainValuesBase<EPMALabel> res
			) {
		final Map<String, Collection<? extends EPMALabel>> labels = new TreeMap<>();
		labels.put("[µ/ρ]", filter(res.getLabels(), "[μ/ρ]"));
		labels.put("dz", filter(res.getLabels(), "dz"));
		labels.put("TOA", filter(res.getLabels(), "TOA"));
		labels.put("E0", filter(res.getLabels(), "E0"));
		labels.put("k", filter(res.getLabels(), "k"));
		labels.put("Fs", filter(res.getLabels(), "Fs"));
		labels.put("E0", filter(res.getLabels(), "E0"));
		labels.put("C", filter(res.getLabels(), "C["));
		labels.put("J", filter(res.getLabels(), "J["));
		// labels.put("w", filter(res.getLabels(), "w["));
		labels.put("m", filter(res.getLabels(), "m["));
		return labels;
	}
}
