package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import static org.junit.Assert.assertEquals;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.Collections;
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
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class XPPMatrixCorrectionTest {

	private static final double SIGMA = 3.0;
	private static final double DELTA_JAC = 0.0001;
	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);
	public static boolean DUMP = true;
	public static int MC_ITERATIONS = 20000; // 16*8000;

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
				std, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(stdMcd, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		final Report r = new Report("XPP Report - test1");
		UncertainValues resultsD = null;
		try {
			{
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
				r.addHeader("test1()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(unkMcd, Collections.singleton(stdMcd));
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
				Object tagFChiFu = XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, unkMcd, cxr);
				assertEquals(results.getEntry(tagFChiFu), 0.635, 0.001);
				Object tagFChiFs = XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, stdMcd, cxr);
				assertEquals(results.getEntry(tagFChiFs), 0.822, 0.001);
				Object tagZA = XPPMatrixCorrection.zaTag(unkMcd, stdMcd, cxr);
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

				final Object unkCompTag = new XPPMatrixCorrection.CompositionTag("Jc", unk);
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
				final UncertainValues inputs = xpp.buildInputs(unkMcd, Collections.singleton(stdMcd));
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
				assertEquals(
						results.getEntry(
								XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, unkMcd, cxr)),
						0.376, 0.001);
				assertEquals(
						results.getEntry(
								XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, stdMcd, cxr)),
						0.353, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.zaTag(unkMcd, stdMcd, cxr)), 1.078, 0.001);

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
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ZATag) {
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

	public Map<Element, Number> buildK411() {
		Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Silicon, new UncertainValue(0.25190067871134, 0.00448737523776));
		res.put(Element.Iron, new UncertainValue(0.11255374113608, 0.00209872307367));
		res.put(Element.Magnesium, new UncertainValue(0.09117902759616, 0.0012060717936));
		res.put(Element.Calcium, new UncertainValue(0.11070559976183, 0.00107203615005));
		res.put(Element.Oxygen, new UncertainValue(0.42346095279459, 0.00693579374492));
		return res;
	}

	public Map<Element, Number> buildK412() {
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
				std, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
		scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1));
		scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(stdMcd, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);
		final Report r = new Report("XPP Report - test2");
		UncertainValues resultsD = null;
		try {
			{
				r.addHeader("test2()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(unkMcd, stds.keySet());
				r.add(inputs);

				final long start = System.currentTimeMillis();
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				System.out.println("Full Timing = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof XPPMatrixCorrection.ZATag) {
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
				final UncertainValues inputs = xpp.buildInputs(unkMcd, Collections.singleton(stdMcd));
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
					if (tag instanceof XPPMatrixCorrection.ZATag) {
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
						if (tag instanceof XPPMatrixCorrection.ZATag) {
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
				std, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);
		final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		// scxr.add(CharacteristicXRay.create(Element.Silicon,
		// XRayTransition.KA1));
		// scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
		// scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
		scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2));
		// scxr.add(CharacteristicXRay.create(Element.Calcium,
		// XRayTransition.KA1));
		// scxr.add(CharacteristicXRay.create(Element.Calcium,
		// XRayTransition.L3M1));
		// scxr.add(CharacteristicXRay.create(Element.Oxygen,
		// XRayTransition.KA1));

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(stdMcd, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, stds);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<MatrixCorrectionDatum, CharacteristicXRaySet> me : stds.entrySet()) {
			final MatrixCorrectionDatum meStd = me.getKey();
			for (final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zaTag(unkMcd, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, unkMcd, cxr));
				outputs.add(XPPMatrixCorrection.tagCharacterisitic(XPPMatrixCorrection.F_CHI_F, meStd, cxr));
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
				final UncertainValues inputs = xpp.buildInputs(unkMcd, Collections.singleton(stdMcd));
				r.add(inputs);
				final long start = System.currentTimeMillis();
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				System.out.println("Timing = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(results);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<? extends Object, UncertainValue> outVals = xpp.getOutputValues(inputs);
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final Object outTag : xpp.getOutputTags()) {
					if (outTag instanceof XPPMatrixCorrection.ZATag) {
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
				final UncertainValues inputs = xpp.buildInputs(unkMcd, Collections.singleton(stdMcd));
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
					if (tag instanceof XPPMatrixCorrection.ZATag) {
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
						if (tag instanceof XPPMatrixCorrection.ZATag) {
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
				std0, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		{
			final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
			scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
			scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
			scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
			scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
			scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2));
			scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1));
			scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1));
			scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));

			stds.put(std0Mcd, scxr);
			stds.put(std1Mcd,
					CharacteristicXRaySet.build(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)));
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
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInputs(unkMcd, stds.keySet());
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds);
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs).sort();
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
					if (outTag instanceof XPPMatrixCorrection.ZATag) {
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
						"Trimmed Delta Timing = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
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
					if (tag instanceof XPPMatrixCorrection.ZATag) {
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
						if (tag instanceof XPPMatrixCorrection.ZATag) {
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
				std0, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std1Mcd = new MatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.9)) //
		);

		MatrixCorrectionDatum std2Mcd = new MatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std3Mcd = new MatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum std4Mcd = new MatrixCorrectionDatum( //
				std4, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		MatrixCorrectionDatum unkMcd = new MatrixCorrectionDatum( //
				unk, //
				new UncertainValue(15.0, 0.12), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)) //
		);

		final Map<MatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		{
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));
				stds.put(std0Mcd, scxr);
			}

			stds.put(std1Mcd,
					CharacteristicXRaySet.build(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)));
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2));
				stds.put(std2Mcd, scxr);
			}
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1));
				stds.put(std3Mcd, scxr);
			}
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
				stds.put(std4Mcd, scxr);
			}
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
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInputs(unkMcd, stds.keySet());
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unkMcd, stds);
		final UncertainValues results2 = UncertainValues.propagate(xpp2, inputs).sort();
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
					if (outTag instanceof XPPMatrixCorrection.ZATag) {
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
						"Trimmed Delta Timing = " + Long.toString(System.currentTimeMillis() - start3) + " ms");
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
					if (tag instanceof XPPMatrixCorrection.ZATag) {
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
						if (tag instanceof XPPMatrixCorrection.ZATag) {
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

}
