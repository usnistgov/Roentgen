package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import static org.junit.Assert.assertEquals;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.MathUtilities;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;
import gov.nist.microanalysis.roentgen.math.uncertainty.MCPropagator;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioTag;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionTag;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioTag.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection.Variates;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.Shell.Principle;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat.OutputMode;
import joinery.DataFrame;

public class XPPMatrixCorrectionTest {

	private static final double SIGMA = 3.0;
	private static final double DELTA_JAC = 0.0001;
	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);
	public static boolean DUMP = false;
	public static int MC_ITERATIONS = 16 * 8000;

	/**
	 * Computes Si and O in Al2SiO5 using SiO2
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP1() throws ArgumentException, ParseException, IOException {
		// final Composition unk = Composition.parse("Al2SiO5").asMassFraction();
		final List<Element> elmsU = Arrays.asList(Element.Aluminum, Element.Silicon, Element.Oxygen);
		final RealVector valsU = new ArrayRealVector(new double[] { 0.3330, 0.1733, 0.4937 });
		final RealVector varsU = new ArrayRealVector(new double[] { 1.0e-6, 0.4e-6, 4.0e-6 });
		final Composition unk = Composition.massFraction("Al<sub>2</sub>SiO<sub>5</sub>", elmsU, valsU, varsU);
		final List<Element> elmsS = Arrays.asList(Element.Silicon, Element.Oxygen);
		final RealVector valsS = new ArrayRealVector(new double[] { 0.4674, 0.5326 });
		final RealVector varsS = new ArrayRealVector(new double[] { 2.0e-6, 0.9e-6 });
		final Composition std = Composition.massFraction("SiO<sub>2</sub>", elmsS, valsS, varsS);

		MatrixCorrectionDatum stdMcd = new MatrixCorrectionDatum( //
				std, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), stdMcd);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		final Report r = new Report("XPP Report - test1");
		UncertainValues resultsD = null;
		try {
			{
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
				r.addHeader("test1()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				Object tagAu = XPPMatrixCorrection.tagShell("A", unkMcd, cxr.getInner());
				assertEquals(results.getEntry(tagAu), 401.654, 0.001);
				Object tagau = XPPMatrixCorrection.tagShell("a", unkMcd, cxr.getInner());
				assertEquals(results.getEntry(tagau), 11255.385, 0.001);
				Object tagBu = XPPMatrixCorrection.tagShell("B", unkMcd, cxr.getInner());
				assertEquals(results.getEntry(tagBu), -529730.331, 0.001);
				Object tagbu = XPPMatrixCorrection.tagShell("b", unkMcd, cxr.getInner());
				assertEquals(results.getEntry(tagbu), 12643.340, 0.001);
				Object tagPhi0u = XPPMatrixCorrection.tagPhi0(unkMcd, cxr.getInner());
				assertEquals(results.getEntry(tagPhi0u), 1.252, 0.001);

				Object tagAs = XPPMatrixCorrection.tagShell("A", stdMcd, cxr.getInner());
				assertEquals(results.getEntry(tagAs), 396.744, 0.001);
				Object tagas = XPPMatrixCorrection.tagShell("a", stdMcd, cxr.getInner());
				assertEquals(results.getEntry(tagas), 11382.116, 0.001);
				Object tagBs = XPPMatrixCorrection.tagShell("B", stdMcd, cxr.getInner());
				assertEquals(results.getEntry(tagBs), -532506.458, 0.001);
				Object tagbs = XPPMatrixCorrection.tagShell("b", stdMcd, cxr.getInner());
				assertEquals(results.getEntry(tagbs), 12795.314, 0.001);
				Object tagPhi0s = XPPMatrixCorrection.tagPhi0(stdMcd, cxr.getInner());
				assertEquals(results.getEntry(tagPhi0s), 1.254, 0.001);

				Object tagChiu = XPPMatrixCorrection.tagChi(unkMcd, cxr);
				assertEquals(results.getEntry(tagChiu), 2542.429, 0.001);
				Object tagChis = XPPMatrixCorrection.tagChi(stdMcd, cxr);
				assertEquals(results.getEntry(tagChis), 1038.418, 0.001);
				Object tagFChiFu = XPPMatrixCorrection.tagFxF(unkMcd, cxr);
				assertEquals(results.getEntry(tagFChiFu), 0.635, 0.001);
				Object tagFChiFs = XPPMatrixCorrection.tagFxF(stdMcd, cxr);
				assertEquals(results.getEntry(tagFChiFs), 0.822, 0.001);
				Object tagZA = XPPMatrixCorrection.zafTag(unkMcd, stdMcd, cxr);
				assertEquals(results.getEntry(tagZA), 0.781, 0.001);

				// Check that INamedMultivariateFunction works...
				RealVector quick = xpp.optimized(inputs.getValues());

				assertEquals(results.getEntry(tagAu), quick.getEntry(xpp.outputIndex(tagAu)), 0.001);
				assertEquals(results.getEntry(tagau), quick.getEntry(xpp.outputIndex(tagau)), 0.001);
				assertEquals(results.getEntry(tagBu), quick.getEntry(xpp.outputIndex(tagBu)), 0.001);
				assertEquals(results.getEntry(tagbu), quick.getEntry(xpp.outputIndex(tagbu)), 0.001);
				assertEquals(results.getEntry(tagPhi0u), quick.getEntry(xpp.outputIndex(tagPhi0u)), 0.001);

				assertEquals(results.getEntry(tagAs), quick.getEntry(xpp.outputIndex(tagAs)), 0.001);
				assertEquals(results.getEntry(tagas), quick.getEntry(xpp.outputIndex(tagas)), 0.001);
				assertEquals(results.getEntry(tagBs), quick.getEntry(xpp.outputIndex(tagBs)), 0.001);
				assertEquals(results.getEntry(tagbs), quick.getEntry(xpp.outputIndex(tagbs)), 0.001);
				assertEquals(results.getEntry(tagPhi0s), quick.getEntry(xpp.outputIndex(tagPhi0s)), 0.001);

				assertEquals(results.getEntry(tagChiu), quick.getEntry(xpp.outputIndex(tagChiu)), 0.001);
				assertEquals(results.getEntry(tagChis), quick.getEntry(xpp.outputIndex(tagChis)), 0.001);
				assertEquals(results.getEntry(tagFChiFu), quick.getEntry(xpp.outputIndex(tagFChiFu)), 0.001);
				assertEquals(results.getEntry(tagFChiFs), quick.getEntry(xpp.outputIndex(tagFChiFs)), 0.001);
				assertEquals(results.getEntry(tagZA), quick.getEntry(xpp.outputIndex(tagZA)), 0.001);

				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					final UncertainValue uv = outVals.get(outTag);
					valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
							Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
							Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final NamedMultivariateJacobian jac = NamedMultivariateJacobian.compute(xpp, inputs.getValues());
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > //
							0.01 * Math.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx)))) {
								System.out.print(jac.getOutputTags().get(oIdx));
								System.out.print(jac.getInputTags().get(iIdx));
								assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
										.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
							}
						}
				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}

				final Object unkCompTag = new XPPMatrixCorrection.CompositionTag("J", unk);
				assertEquals(jac.getEntry(unkCompTag, Composition.buildMassFractionTag(unk, Element.Oxygen)), -0.027565,
						0.00001);
				assertEquals(jac.getEntry(unkCompTag, XPPMatrixCorrection.meanIonizationTag(Element.Oxygen)), 0.609601,
						0.00001);

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1);
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();

				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("A", unkMcd, cxr.getInner())), 2366.373,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("a", unkMcd, cxr.getInner())), 11402.291,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("B", unkMcd, cxr.getInner())), -1506725.664,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("b", unkMcd, cxr.getInner())), 12050.502,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagPhi0(unkMcd, cxr.getInner())), 1.258, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("A", stdMcd, cxr.getInner())), 2307.215,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("a", stdMcd, cxr.getInner())), 11531.967,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("B", stdMcd, cxr.getInner())), -1505332.755,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagShell("b", stdMcd, cxr.getInner())), 12196.382,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagPhi0(stdMcd, cxr.getInner())), 1.26, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tagChi(unkMcd, cxr)), 5836.018, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagChi(stdMcd, cxr)), 6414.025, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagFxF(unkMcd, cxr)), 0.376, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tagFxF(stdMcd, cxr)), 0.353, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.zafTag(unkMcd, stdMcd, cxr)), 1.078, 0.001);

				r.addHeader("Analyic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("Monte Carlo Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : resultsMc.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				r.addSubHeader("Phi0");
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object tag : xpp.getOutputTags()) {
					r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
					r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
				}

				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					BasicNumberFormat bnf2 = new BasicNumberFormat("0.0000");
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(results.getUncertainty(tag), bnf2),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
						}
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof KRatioTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(results.getUncertainty(tag), bnf2),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (Exception e) {
			e.printStackTrace();
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	public static Map<Element, Number> buildK411() {
		Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Silicon, new UncertainValue(0.25190067871134, 0.00448737523776));
		res.put(Element.Iron, new UncertainValue(0.11255374113608, 0.00209872307367));
		res.put(Element.Magnesium, new UncertainValue(0.09117902759616, 0.0012060717936));
		res.put(Element.Calcium, new UncertainValue(0.11070559976183, 0.00107203615005));
		res.put(Element.Oxygen, new UncertainValue(0.42346095279459, 0.00693579374492));
		return res;
	}

	public static Map<Element, Number> buildK412() {
		Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Silicon, new UncertainValue(0.21226219744446, 0.00359924888862));
		res.put(Element.Iron, new UncertainValue(0.07726410130474, 0.00139914871578));
		res.put(Element.Magnesium, new UncertainValue(0.11855685731088, 0.001507589742));
		res.put(Element.Calcium, new UncertainValue(0.11034825437848, 0.00107203615005));
		res.put(Element.Aluminum, new UncertainValue(0.0494320434, 0.0015348279));
		res.put(Element.Oxygen, new UncertainValue(0.43003654616144, 0.00728714860355));
		return res;
	}

	/**
	 * Computes K412 vs K411 as in SP 260-74
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP2() throws ArgumentException, ParseException, IOException {
		// K411 and K412 as in SP 260-74
		final boolean combined = false;
		final Composition std = combined ? Composition.combine("K411", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015))) : //
				Composition.massFraction("K411", buildK411());

		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum stdMcd = new MatrixCorrectionDatum( //
				std, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), stdMcd);

		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		final Report r = new Report("XPP Report - test2");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test2()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);

				final long start = System.currentTimeMillis();
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				System.out.println("Full Timing(2) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final NamedMultivariateJacobian jac = NamedMultivariateJacobian.compute(xpp, inputs.getValues());
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8)
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();
				r.addHeader("Analytic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}

				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	/**
	 * Calculates the correction for Mg in K412 against pure Mg
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP3() throws ArgumentException, ParseException, IOException {
		// K411 and K412 as in SP 260-74

		final Composition std = Composition.parse("Mg").asMassFraction();

		final boolean combined = false;
		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum stdMcd = new MatrixCorrectionDatum( //
				std, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);
		// final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		// stds.put(ElementXRaySet.singleton((Element.Silicon,
		// XRayTransition.KA1));
		// stds.put(ElementXRaySet.singleton((Element.Iron, XRayTransition.KA1));
		// stds.put(ElementXRaySet.singleton((Element.Iron, XRayTransition.LA1));
		// stds.put(ElementXRaySet.singleton((Element.Magnesium, XRayTransition.KA1));
		// stds.put(ElementXRaySet.singleton((Element.Magnesium, XRayTransition.KA2));
		// stds.put(ElementXRaySet.singleton((Element.Calcium,
		// XRayTransition.KA1));
		// stds.put(ElementXRaySet.singleton((Element.Calcium,
		// XRayTransition.L3M1));
		// stds.put(ElementXRaySet.singleton((Element.Oxygen,
		// XRayTransition.KA1));

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), stdMcd);
		stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2)), stdMcd);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(meStd, cxr));
			}
		}

		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final Report r = new Report("XPP Report - test3");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test3()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);
				final long start = System.currentTimeMillis();
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				System.out.println("Timing(3) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final NamedMultivariateJacobian jac = NamedMultivariateJacobian.compute(xpp, inputs.getValues());
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInput();
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();
				r.addHeader("Analytic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(16 * MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsMc, L2C, 8), "Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}

				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	/**
	 * Compute ZAF for K412 as in SP 260-74 using K411 and Al
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP4() throws Exception {
		final boolean combined = false;
		final Composition std0 = combined ? Composition.combine("K411", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015))) : //
				Composition.massFraction("K411", buildK411());

		final Composition std1 = Composition.parse("Al");

		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		{
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), std0Mcd);
			stds.put(new ElementXRaySet(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)), std1Mcd);
		}

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(meStd, cxr));
			}
		}

		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInput();
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds);
		final UncertainValues inputs2 = xpp2.buildInput();
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs2).sort();
		System.out.println("Full Timing (4) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		for (final Object outTag : outputs) {
			for (final Object outTag2 : outputs)
				if (outTag == outTag2) {
					final int oi1 = results.indexOf(outTag);
					final double v1 = results.getEntry(oi1);
					final int oi2 = results2.indexOf(outTag);
					final double v2 = results2.getEntry(oi2);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (results.getVariance(oi1) > 1.0e-8)
						assertEquals(results.getVariance(oi1), results2.getVariance(oi2),
								0.01 * Math.abs(results.getVariance(oi1)));
				} else {
					final double cov1 = results.getCovariance(outTag, outTag2);
					final double cov2 = results2.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}
		final Report r = new Report("XPP Report - test4()");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test4()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				System.out.println(
						"Trimmed Delta Timing (4) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("Monte Carlo Results");
				r.add(resultsMc);
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsMc, L2C, 8), "Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}

				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	/**
	 * Compute ZAF for K412 as in SP 260-74 using elements and simple compounds
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP5() throws ArgumentException, ParseException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds

		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al");
		final Composition std2 = Composition.parse("Mg");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final boolean combined = false;
		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std4Mcd = new MatrixCorrectionDatum( //
				std4, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1), std0Mcd);
		stds.put(ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1), std0Mcd);
		stds.put(ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1), std1Mcd);

		stds.put(ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1), std2Mcd);
		stds.put(ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA2), std2Mcd);

		stds.put(ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1), std3Mcd);
		stds.put(ElementXRaySet.singleton(Element.Calcium, XRayTransition.L3M1), std3Mcd);

		stds.put(ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1), std4Mcd);
		stds.put(ElementXRaySet.singleton(Element.Iron, XRayTransition.LA1), std4Mcd);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(meStd, cxr));
			}
		}
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds, XPPMatrixCorrection.defaultVariates());
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInput();
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (5) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds);
		final UncertainValues inputs2 = xpp2.buildInput();
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs2).sort();
		System.out.println("Full Timing (5) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		for (final Object outTag : outputs) {
			for (final Object outTag2 : outputs)
				if (outTag == outTag2) {
					final int oi1 = results.indexOf(outTag);
					final double v1 = results.getEntry(oi1);
					final int oi2 = results2.indexOf(outTag);
					final double v2 = results2.getEntry(oi2);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (results.getVariance(oi1) > 1.0e-8)
						assertEquals(results.getVariance(oi1), results2.getVariance(oi2),
								0.01 * Math.abs(results.getVariance(oi1)));
				} else {
					final double cov1 = results.getCovariance(outTag, outTag2);
					final double cov2 = results2.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test5()");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test5()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				System.out.println(
						"Trimmed Delta Timing (5) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsMc, L2C, 8), "Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	/**
	 * Compute ZAF for K412 as in SP 260-74 using elements and simple compounds Only
	 * calculate the uncertainties for a subset of the input parameters.
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP6() throws ArgumentException, ParseException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds

		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al");
		final Composition std2 = Composition.parse("Mg");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final boolean combined = false;
		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std4Mcd = new MatrixCorrectionDatum( //
				std4, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1), std0Mcd);
		stds.put(ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1), std0Mcd);

		stds.put(ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1), std1Mcd);

		stds.put(ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1), std2Mcd);
		stds.put(ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA2), std2Mcd);

		stds.put(ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1), std3Mcd);
		stds.put(ElementXRaySet.singleton(Element.Calcium, XRayTransition.L3M1), std3Mcd);

		stds.put(ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1), std4Mcd);
		stds.put(ElementXRaySet.singleton(Element.Iron, XRayTransition.LA1), std4Mcd);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(meStd, cxr));
			}
		}

		Set<XPPMatrixCorrection.Variates> variates = new HashSet<>();
		variates.add(Variates.UnknownComposition);
		variates.add(Variates.StandardComposition);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds, variates);
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInput();
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (6) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds, variates);
		final UncertainValues inputs2 = xpp2.buildInput();
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs2).sort();
		System.out.println("Full Timing (6) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		for (final Object outTag : outputs) {
			for (final Object outTag2 : outputs)
				if (outTag == outTag2) {
					final int oi1 = results.indexOf(outTag);
					final double v1 = results.getEntry(oi1);
					final int oi2 = results2.indexOf(outTag);
					final double v2 = results2.getEntry(oi2);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (results.getVariance(oi1) > 1.0e-8)
						assertEquals(results.getVariance(oi1), results2.getVariance(oi2),
								0.01 * Math.abs(results.getVariance(oi1)));
				} else {
					final double cov1 = results.getCovariance(outTag, outTag2);
					final double cov2 = results2.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test6()");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test6()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				System.out.println(
						"Trimmed Delta Timing (6) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsMc, L2C, 8), "Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	public static Map<Element, Number> buildK240() {
		Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Magnesium, new UncertainValue(0.030154, 0.00030154));
		res.put(Element.Silicon, new UncertainValue(0.186986, 0.00186986));
		res.put(Element.Titanium, new UncertainValue(0.059950, 0.00059950));
		res.put(Element.Zinc, new UncertainValue(0.040168, 0.00040168));
		res.put(Element.Zirconium, new UncertainValue(0.074030, 0.00074030));
		res.put(Element.Barium, new UncertainValue(0.268689, 0.000268689));
		res.put(Element.Oxygen, new UncertainValue(0.340023, 0.00340023));
		return res;
	}

	/**
	 * Compute ZAF for K240 using elements and simple compounds Only calculate the
	 * uncertainties for a subset of the input parameters.
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP7() throws ArgumentException, ParseException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds
		final Composition unk = Composition.massFraction("K240", buildK240());
		final Composition std0 = Composition.parse("Mg2SiO4");
		final Composition std1 = Composition.parse("BaTiSi3O9");
		final Composition std2 = Composition.parse("Zn");
		final Composition std3 = Composition.parse("Zr");

		// Mg
		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.2));
		// Ba, Ti, Si, O
		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)), //
				MatrixCorrectionDatum.roughness(10.0, 3.6));
		// Zn
		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		// Zr
		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 5.62));

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, false, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)), //
				MatrixCorrectionDatum.roughness(10.0, 3.5));

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		// Mg
		stds.put(ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1), std0Mcd);
		// Ba, Ti, Si, O
		stds.put(ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1), std1Mcd);
		stds.put(ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1), std1Mcd);
		stds.put(ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1), std1Mcd);
		stds.put(ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1), std1Mcd);
		// Zn
		stds.put(ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1), std2Mcd);
		// Zr
		stds.put(ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1), std3Mcd);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.zTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.aTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagFxF(meStd, cxr));
				outputs.add(new KRatioTag(unkMcd, meStd, cxr, Method.Calculated));
			}
		}
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds, XPPMatrixCorrection.allVariates());
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInput();
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (7) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds, XPPMatrixCorrection.allVariates());
		final UncertainValues inputs2 = xpp2.buildInput();
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs2).sort();
		System.out.println("Full Timing (7) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		for (final Object outTag : outputs) {
			for (final Object outTag2 : outputs)
				if (outTag == outTag2) {
					final int oi1 = results.indexOf(outTag);
					final double v1 = results.getEntry(oi1);
					final int oi2 = results2.indexOf(outTag);
					final double v2 = results2.getEntry(oi2);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (results.getVariance(oi1) > 1.0e-8)
						assertEquals(results.getVariance(oi1), results2.getVariance(oi2),
								0.01 * Math.abs(results.getVariance(oi1)));
				} else {
					final double cov1 = results.getCovariance(outTag, outTag2);
					final double cov2 = results2.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test7()");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test7()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				System.out.println(
						"Trimmed Delta Timing (7) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

				Table res = new Table();
				res.addRow(Table.th("Line"), Table.th("ZA"), Table.th("Z"), Table.th("A"));
				BasicNumberFormat bnf2 = new BasicNumberFormat("0.000");

				for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
					final MatrixCorrectionDatum meStd = me.getValue();
					for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
						Object zaTag = XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr);
						Object zTag = XPPMatrixCorrection.zTag(unkMcd, meStd, cxr);
						Object aTag = XPPMatrixCorrection.aTag(unkMcd, meStd, cxr);
						UncertainValue za = results.getUncertainValue(zaTag);
						UncertainValue a = results.getUncertainValue(aTag);
						UncertainValue z = results.getUncertainValue(zTag);
						res.addRow(Table.td(cxr), //
								Table.td(bnf2.formatHTML(za, OutputMode.ValuePlusUncertainty)), //
								Table.td(bnf2.formatHTML(z, OutputMode.ValuePlusUncertainty)), //
								Table.td(bnf2.formatHTML(a, OutputMode.ValuePlusUncertainty)));
					}
				}
				r.addHeader("ZAF results");
				r.add(res, Mode.NORMAL);

				Table tk = new Table();
				for (final Object tag : xpp.getOutputTags())
					if (tag instanceof KRatioTag) {
						tk.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
								MathUtilities.td(results.getValue(tag).doubleValue(), bnf2),
								MathUtilities.td(results.getUncertainty(tag), bnf2),
								MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
								MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
					}
				r.addHeader("K-ratio");
				r.add(tk, Mode.NORMAL);

				for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
					final MatrixCorrectionDatum meStd = me.getValue();
					for (final CharacteristicXRay cxr : me.getKey().getSetOfCharacteristicXRay()) {
						Object zaTag = XPPMatrixCorrection.zafTag(unkMcd, meStd, cxr);
						Object fUnkTag = XPPMatrixCorrection.tagFxF(unkMcd, cxr);
						Object fStdTag = XPPMatrixCorrection.tagFxF(meStd, cxr);
						double za = results.getEntry(zaTag);
						double a = results.getEntry(fUnkTag) / results.getEntry(fStdTag);
						double z = za / a;
						// IUPAC Seigbahn Standard Energy ZAF Z A F k-ratio
						// O K-L3 O Kα1 BaTiSi3O9 0.5249 1.0163 0.9979 1.0184 1.0001 0.992287
						// Mg K-L3 Mg Kα1 Mg2SiO4 1.2536 0.7388 1.1343 0.6515 0.9996 0.064471
						// Si K-L3 Si Kα1 BaTiSi3O9 1.7397 0.9969 0.9977 0.9980 1.0012 0.914605
						// Ti K-L3 Ti Kα1 BaTiSi3O9 4.5109 0.9868 0.9973 0.9906 0.9989 0.510866
						// Zn K-L3 Zn Kα1 Pure Zn 8.6389 0.8804 0.8920 0.9870 1.0000 0.035367
						// Zr L3-M5 Zr Lα1 Pure Zr 2.0423 0.6802 0.8393 0.8091 1.0016 0.050359
						// Ba L3-M5 Ba Lα1 BaTiSi3O9 4.4663 0.9882 0.9972 0.9906 1.0003 0.799365
						if (cxr.getElement() == Element.Oxygen) {
							assertEquals(0.9979 * 1.0184, za, 0.002);
							assertEquals(0.9979, z, 0.002);
							assertEquals(1.0184, a, 0.002);
						} else if (cxr.getElement() == Element.Magnesium) {
							assertEquals(1.1343 * 0.6515, za, 0.025); // ??????
							assertEquals(1.1343, z, 0.005); // ??????
							assertEquals(0.6515, a, 0.025); // ??????
						} else if (cxr.getElement() == Element.Silicon) {
							assertEquals(0.9977 * 0.9980, za, 0.002);
							assertEquals(0.9977, z, 0.002);
							assertEquals(0.9980, a, 0.002);
						} else if (cxr.getElement() == Element.Titanium) {
							assertEquals(0.9973 * 0.9906, za, 0.003);
							assertEquals(0.9973, z, 0.002);
							assertEquals(0.9906, a, 0.002);
						} else if (cxr.getElement() == Element.Zinc) {
							assertEquals(0.8920 * 0.9870, za, 0.005);
							assertEquals(0.8920, z, 0.005);
							assertEquals(0.9870, a, 0.005);
						} else if (cxr.getElement() == Element.Zirconium) {
							assertEquals(0.8393 * 0.8091, za, 0.02); // ??????
							assertEquals(0.8393, z, 0.008); // ??????
							assertEquals(0.8091, a, 0.02);// ??????
						} else if (cxr.getElement() == Element.Barium) {
							assertEquals(0.9972 * 0.9906, za, 0.002);
							assertEquals(0.9972, z, 0.002);
							assertEquals(0.9906, a, 0.002);
						}
					}
				}
				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}

	/**
	 * Compute ZAF for K412 as in SP 260-74 using elements and simple compounds
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP8() throws ArgumentException, ParseException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds

		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al");
		final Composition std2 = Composition.parse("Mg");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final boolean combined = false;
		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std4Mcd = new MatrixCorrectionDatum( //
				std4, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(ElementXRaySet.build(Element.Silicon, Principle.K, 0.001), std0Mcd);
		stds.put(ElementXRaySet.build(Element.Oxygen, Principle.K, 0.001), std0Mcd);
		stds.put(ElementXRaySet.build(Element.Aluminum, Principle.K, 0.001), std1Mcd);

		stds.put(ElementXRaySet.build(Element.Magnesium, Principle.K, 0.001), std2Mcd);

		stds.put(ElementXRaySet.build(Element.Calcium, Principle.K, 0.001), std3Mcd);
		stds.put(ElementXRaySet.build(Element.Calcium, Principle.L, 0.001), std3Mcd);

		stds.put(ElementXRaySet.build(Element.Iron, Principle.K, 0.001), std4Mcd);
		stds.put(ElementXRaySet.build(Element.Iron, Principle.L, 0.001), std4Mcd);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			final ElementXRaySet exrs = me.getKey();
			outputs.add(new MatrixCorrectionTag(unkMcd, meStd, exrs));
		}
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds, XPPMatrixCorrection.defaultVariates());
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInput();
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (8) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds);
		final UncertainValues inputs2 = xpp2.buildInput();
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs2).sort();
		System.out.println("Full Timing (8) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		for (final Object outTag : outputs) {
			for (final Object outTag2 : outputs)
				if (outTag == outTag2) {
					final int oi1 = results.indexOf(outTag);
					final double v1 = results.getEntry(oi1);
					final int oi2 = results2.indexOf(outTag);
					final double v2 = results2.getEntry(oi2);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (results.getVariance(oi1) > 1.0e-8)
						assertEquals(results.getVariance(oi1), results2.getVariance(oi2),
								0.01 * Math.abs(results.getVariance(oi1)));
				} else {
					final double cov1 = results.getCovariance(outTag, outTag2);
					final double cov2 = results2.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test8()");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test8()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs, trimmed)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof MatrixCorrectionTag) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(results.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final NamedMultivariateJacobian djac = NamedMultivariateJacobian.computeDelta(xpp, inputs, DELTA_JAC);
				System.out.println(
						"Trimmed Delta Timing (8) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
				for (int oIdx = 0; oIdx < jac.getOutputDimension(); ++oIdx)
					for (int iIdx = 0; iIdx < jac.getInputDimension(); ++iIdx)
						if (Math.abs(jac.getEntry(oIdx, iIdx)) > 1.0e-8) {
							if (Math.abs(jac.getEntry(oIdx, iIdx) - djac.getEntry(oIdx, iIdx)) > 0.01
									* Math.abs(jac.getEntry(oIdx, iIdx)))
								System.out.println(xpp.getInputTags().get(iIdx) + ", " + xpp.getOutputTags().get(oIdx)
										+ "=[ " + jac.getEntry(oIdx, iIdx) + " ?=? " + djac.getEntry(oIdx, iIdx) + "]");
							assertEquals(jac.getEntry(oIdx, iIdx), djac.getEntry(oIdx, iIdx), 0.01 * Math
									.max(Math.abs(jac.getEntry(oIdx, iIdx)), Math.abs(djac.getEntry(oIdx, iIdx))));
						}

				resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(results.toCSV());

					System.out.println("Jacobian");
					System.out.println(jac.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(djac.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances()));
				final UncertainValues resultsMc = mcp.computeMT(MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final Object tag : results.getTags()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(results.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsMc, L2C, 8), "Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final Object tag : xpp.getOutputTags()) {
					if (tag instanceof MatrixCorrectionTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof MatrixCorrectionTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (Exception e) {
			r.addHTML(HTML.error(HTML.escape(e.getMessage())));
			throw e;
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}
	}
	
	@Test
	public void testXPP9() throws ArgumentException, ParseException, IOException {
		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al");
		final Composition std2 = Composition.parse("Mg");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final boolean combined = false;
		final Composition unk = combined ? Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)))
				: Composition.massFraction("K412", buildK412());

		MatrixCorrectionDatum std0Mcd = new MatrixCorrectionDatum( //
				std0, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std4Mcd = new MatrixCorrectionDatum( //
				std4, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, true, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<ElementXRaySet, MatrixCorrectionDatum> stds = new HashMap<>();
		stds.put(ElementXRaySet.build(Element.Silicon, Principle.K, 0.001), std0Mcd);
		stds.put(ElementXRaySet.build(Element.Oxygen, Principle.K, 0.001), std0Mcd);
		stds.put(ElementXRaySet.build(Element.Aluminum, Principle.K, 0.001), std1Mcd);

		stds.put(ElementXRaySet.build(Element.Magnesium, Principle.K, 0.001), std2Mcd);

		stds.put(ElementXRaySet.build(Element.Calcium, Principle.K, 0.001), std3Mcd);
		stds.put(ElementXRaySet.build(Element.Calcium, Principle.L, 0.001), std3Mcd);

		stds.put(ElementXRaySet.build(Element.Iron, Principle.K, 0.001), std4Mcd);
		stds.put(ElementXRaySet.build(Element.Iron, Principle.L, 0.001), std4Mcd);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getValue();
			final ElementXRaySet exrs = me.getKey();
			outputs.add(new MatrixCorrectionTag(unkMcd, meStd, exrs));
		}
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds, XPPMatrixCorrection.minimalVariates());
		final UncertainValues results = UncertainValues.propagate(xpp, xpp.buildInput());
		
		DataFrame<Double> df=xpp.computePhiRhoZCurve(results.getValueMap(), 1.201e-3, 2.0e-5, 0.9);
		df.writeCsv("C:\\Users\\nicho\\OneDrive\\Desktop\\prz412.csv");
	}
	
	
	
}

// Unmodifiable FastIndex
// Full Timing(2) = 142 ms
// Timing(3) = 12 ms
// Trimmed Timing (4) = 21 ms
// Full Timing (4) = 117 ms
// Trimmed Delta Timing (4) = 205 ms
// Trimmed Timing (5) = 15 ms
// Full Timing (5) = 111 ms
// Trimmed Delta Timing (5) = 226 ms
// Trimmed Timing (6) = 9 ms
// Full Timing (6) = 48 ms
// Trimmed Delta Timing (6) = 66 ms
// Trimmed Timing (7) = 26 ms
// Full Timing (7) = 103 ms
// Trimmed Delta Timing (7) = 210 ms
// Trimmed Timing (8) = 271 ms
// Full Timing (8) = 1744 ms
// Trimmed Delta Timing (8) = 2682 ms

// ArrayList
// Full Timing(2) = 109 ms
// Timing(3) = 14 ms
// Trimmed Timing (4) = 27 ms
// Full Timing (4) = 98 ms
// Trimmed Delta Timing (4) = 209 ms
// Trimmed Timing (5) = 21 ms
// Full Timing (5) = 103 ms
// Trimmed Delta Timing (5) = 243 ms
// Trimmed Timing (6) = 7 ms
// Full Timing (6) = 50 ms
// Trimmed Delta Timing (6) = 86 ms
// Trimmed Timing (7) = 20 ms
// Full Timing (7) = 91 ms
// Trimmed Delta Timing (7) = 213 ms
// Trimmed Timing (8) = 239 ms
// Full Timing (8) = 1760 ms
// Trimmed Delta Timing (8) = 2837 ms

// FastIndex
// Full Timing(2) = 102 ms
// Timing(3) = 14 ms
// Trimmed Timing (4) = 21 ms
// Full Timing (4) = 98 ms
// Trimmed Delta Timing (4) = 196 ms
// Trimmed Timing (5) = 16 ms
// Full Timing (5) = 120 ms
// Trimmed Delta Timing (5) = 209 ms
// Trimmed Timing (6) = 6 ms
// Full Timing (6) = 44 ms
// Trimmed Delta Timing (6) = 57 ms
// Trimmed Timing (7) = 20 ms
// Full Timing (7) = 80 ms
// Trimmed Delta Timing (7) = 176 ms
// Trimmed Timing (8) = 222 ms
// Full Timing (8) = 1726 ms
// Trimmed Delta Timing (8) = 2407 ms
