package gov.nist.microanalysis.roentgen.tests.quant;

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
import gov.nist.microanalysis.roentgen.microanalysis.XPPMatrixCorrection;
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
	private static final double DELTA_JAC = 0.001;
	private static final ValueToLog3 V2L3 = new ValueToLog3(1.0);
	private static final LinearToColor L2C = new LinearToColor(1.0, Color.blue, Color.red);
	public static boolean DUMP = true;
	public static int MC_ITERATIONS = 0; // 16*8000;

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
		final Object[] tags = new Object[] { XPPMatrixCorrection.beamEnergyTag(unk),
				XPPMatrixCorrection.beamEnergyTag(std), XPPMatrixCorrection.takeOffAngleTag(std),
				XPPMatrixCorrection.takeOffAngleTag(unk) };
		final RealVector vals = new ArrayRealVector(tags.length);
		vals.setEntry(0, 15.0);
		vals.setEntry(1, 15.0);
		vals.setEntry(2, Math.toRadians(40.0));
		vals.setEntry(3, Math.toRadians(40.0));
		final RealVector var = new ArrayRealVector(tags.length);
		var.setEntry(0, Math.pow(0.1, 2.0));
		var.setEntry(1, Math.pow(0.12, 2.0));
		var.setEntry(2, Math.pow(Math.toRadians(0.9), 2.0));
		var.setEntry(3, Math.pow(Math.toRadians(0.7), 2.0));

		final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));

		final UncertainValues conditions = new UncertainValues(Arrays.asList(tags), vals, var);
		final Map<Composition, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(std, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unk, stds);
		final Report r = new Report("XPP Report - test1");
		try {
			{
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
				r.addHeader("test1()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
				r.add(inputs);
				final NamedMultivariateJacobian xppI = new NamedMultivariateJacobian(xpp, inputs.getValues());
				final UncertainValues results = UncertainValues.propagate(xppI, inputs).sort();
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("A", unk, cxr.getInner())), 401.654, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("a", unk, cxr.getInner())), 11255.385, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("B", unk, cxr.getInner())), -529730.331, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("b", unk, cxr.getInner())), 12643.340, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.PHI0, unk, cxr.getInner())),
						1.252, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("A", std, cxr.getInner())), 396.744, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("a", std, cxr.getInner())), 11382.116, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("B", std, cxr.getInner())), -532506.458, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("b", std, cxr.getInner())), 12795.314, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.PHI0, std, cxr.getInner())),
						1.254, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.CHI, unk, cxr)), 2542.429,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.CHI, std, cxr)), 1038.418,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, unk, cxr)), 0.635,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, std, cxr)), 0.822,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.zaTag(unk, std, cxr)), 0.781, 0.001);

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

				final Object unkCompTag = new XPPMatrixCorrection.CompositionTag("J", unk);
				assertEquals(jac.getEntry(unkCompTag, Composition.buildMassFractionTag(unk, Element.Oxygen)), -0.027565,
						0.00001);
				assertEquals(jac.getEntry(unkCompTag, XPPMatrixCorrection.meanIonizationTag(Element.Oxygen)), 0.609601,
						0.00001);

				final UncertainValues resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1);
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();

				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("A", unk, cxr.getInner())), 2366.373, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("a", unk, cxr.getInner())), 11402.291, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("B", unk, cxr.getInner())), -1506725.664, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("b", unk, cxr.getInner())), 12050.502, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.PHI0, unk, cxr.getInner())),
						1.258, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("A", std, cxr.getInner())), 2307.215, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("a", std, cxr.getInner())), 11531.967, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("B", std, cxr.getInner())), -1505332.755, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2("b", std, cxr.getInner())), 12196.382, 0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.PHI0, std, cxr.getInner())),
						1.26, 0.001);

				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.CHI, unk, cxr)), 5836.018,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.CHI, std, cxr)), 6414.025,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, unk, cxr)), 0.376,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, std, cxr)), 0.353,
						0.001);
				assertEquals(results.getEntry(XPPMatrixCorrection.zaTag(unk, std, cxr)), 1.078, 0.001);

				r.addHeader("Analyic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances(), 1.0e-9));
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
					t.addRow(Table.th("Tag"), Table.th("V(MonteCarlo)"), Table.th("U(Monte Carlo)"),
							Table.th("V(Analytical)"), Table.th("U(Analytic)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ResultTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.NORMAL)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf));
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
	 * Computes K412 vs K411 as in SP 260-74
	 * 
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP2() throws ArgumentException, ParseException, IOException {
		// K411 and K412 as in SP 260-74
		final Composition std = Composition.combine("K411", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015)));

		final Composition unk = Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));

		final Object[] tags = new Object[] { XPPMatrixCorrection.beamEnergyTag(unk),
				XPPMatrixCorrection.beamEnergyTag(std), XPPMatrixCorrection.takeOffAngleTag(std),
				XPPMatrixCorrection.takeOffAngleTag(unk) };
		final RealVector vals = new ArrayRealVector(tags.length);
		vals.setEntry(0, 15.0);
		vals.setEntry(1, 15.0);
		vals.setEntry(2, Math.toRadians(40.0));
		vals.setEntry(3, Math.toRadians(40.0));
		final RealVector var = new ArrayRealVector(tags.length);
		var.setEntry(0, Math.pow(0.1, 2.0));
		var.setEntry(1, Math.pow(0.12, 2.0));
		var.setEntry(2, Math.pow(Math.toRadians(0.9), 2.0));
		var.setEntry(3, Math.pow(Math.toRadians(0.7), 2.0));

		final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
		scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
		scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1));
		scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1));
		scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));

		final UncertainValues conditions = new UncertainValues(Arrays.asList(tags), vals, var);
		final Map<Composition, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(std, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unk, stds);
		final Report r = new Report("XPP Report - test2");
		try {
			{
				r.addHeader("test2()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
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
					if (outTag instanceof XPPMatrixCorrection.ResultTag) {
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

				final UncertainValues resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();
				r.addHeader("Analytic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances(), 1.0e-9));
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
					if (tag instanceof XPPMatrixCorrection.ResultTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), Table.th("V(MonteCarlo)"), Table.th("U(Monte Carlo)"),
							Table.th("V(Analytical)"), Table.th("U(Analytic)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ResultTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.NORMAL)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf));
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

		final Composition unk = Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));

		final Object[] tags = new Object[] { XPPMatrixCorrection.beamEnergyTag(unk),
				XPPMatrixCorrection.beamEnergyTag(std), XPPMatrixCorrection.takeOffAngleTag(std),
				XPPMatrixCorrection.takeOffAngleTag(unk) };
		final RealVector vals = new ArrayRealVector(tags.length);
		vals.setEntry(0, 15.0);
		vals.setEntry(1, 15.0);
		vals.setEntry(2, Math.toRadians(40.0));
		vals.setEntry(3, Math.toRadians(40.0));
		final RealVector var = new ArrayRealVector(tags.length);
		var.setEntry(0, Math.pow(0.1, 2.0));
		var.setEntry(1, Math.pow(0.12, 2.0));
		var.setEntry(2, Math.pow(Math.toRadians(0.9), 2.0));
		var.setEntry(3, Math.pow(Math.toRadians(0.7), 2.0));

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

		final UncertainValues conditions = new UncertainValues(Arrays.asList(tags), vals, var);
		final Map<Composition, CharacteristicXRaySet> stds = new HashMap<>();
		stds.put(std, scxr);
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unk, stds);

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet()) {
			final Composition meStd = me.getKey();
			for (final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zaTag(unk, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, unk, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, meStd, cxr));
			}
		}

		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final Report r = new Report("XPP Report - test3");
		try {
			{
				r.addHeader("test3()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
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
					if (outTag instanceof XPPMatrixCorrection.ResultTag) {
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
				final UncertainValues resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues inputs = xpp.buildInputs(conditions);
				r.add(inputs);
				final UncertainValues results = UncertainValues.propagate(xpp, inputs).sort();
				r.addHeader("Analytic Results");
				r.add(results);

				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances(), 1.0e-8));
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
					if (tag instanceof XPPMatrixCorrection.ResultTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), Table.th("V(MonteCarlo)"), Table.th("U(Monte Carlo)"),
							Table.th("V(Analytical)"), Table.th("U(Analytic)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ResultTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.NORMAL)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf));
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
		final Composition std0 = Composition.combine("K411", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015)));

		final Composition std1 = Composition.parse("Al");

		final Composition unk = Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));

		final Object[] tags = new Object[] { XPPMatrixCorrection.beamEnergyTag(unk),
				XPPMatrixCorrection.beamEnergyTag(std0), XPPMatrixCorrection.beamEnergyTag(std1),
				XPPMatrixCorrection.takeOffAngleTag(std0), XPPMatrixCorrection.takeOffAngleTag(std1),
				XPPMatrixCorrection.takeOffAngleTag(unk) };
		final RealVector vals = new ArrayRealVector(tags.length);
		vals.setEntry(0, 15.0);
		vals.setEntry(1, 15.0);
		vals.setEntry(2, 15.0);
		vals.setEntry(3, Math.toRadians(40.0));
		vals.setEntry(4, Math.toRadians(40.0));
		vals.setEntry(5, Math.toRadians(40.0));
		final RealVector var = new ArrayRealVector(tags.length);
		var.setEntry(0, Math.pow(0.1, 2.0));
		var.setEntry(1, Math.pow(0.12, 2.0));
		var.setEntry(2, Math.pow(0.12, 2.0));
		var.setEntry(3, Math.pow(Math.toRadians(0.9), 2.0));
		var.setEntry(4, Math.pow(Math.toRadians(0.7), 2.0));
		var.setEntry(5, Math.pow(Math.toRadians(0.7), 2.0));

		final UncertainValues conditions = new UncertainValues(Arrays.asList(tags), vals, var);

		final Map<Composition, CharacteristicXRaySet> stds = new HashMap<>();
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

			stds.put(std0, scxr);
			stds.put(std1,
					CharacteristicXRaySet.build(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)));
		}

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet()) {
			final Composition meStd = me.getKey();
			for (final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zaTag(unk, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, unk, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, meStd, cxr));
			}
		}

		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unk, stds);
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInputs(conditions);
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unk, stds);
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
					if (outTag instanceof XPPMatrixCorrection.ResultTag) {
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
				final UncertainValues resultsD = UncertainValues.propagate(djac, inputs).sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValues.compareAsBitmap(results, resultsD, L2C, 8), "Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				final MCPropagator mcp = new MCPropagator(xpp, inputs, SIGMA,
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances(), 1.0e-9));
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
					if (tag instanceof XPPMatrixCorrection.ResultTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), Table.th("V(MonteCarlo)"), Table.th("U(Monte Carlo)"),
							Table.th("V(Analytical)"), Table.th("U(Analytic)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ResultTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.NORMAL)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf));
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

		final Composition unk = Composition.combine("K412", //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));

		final Object[] tags = new Object[] { XPPMatrixCorrection.beamEnergyTag(unk),
				XPPMatrixCorrection.beamEnergyTag(std0), XPPMatrixCorrection.beamEnergyTag(std1),
				XPPMatrixCorrection.beamEnergyTag(std2), XPPMatrixCorrection.beamEnergyTag(std3),
				XPPMatrixCorrection.beamEnergyTag(std4), XPPMatrixCorrection.takeOffAngleTag(unk),
				XPPMatrixCorrection.takeOffAngleTag(std0), XPPMatrixCorrection.takeOffAngleTag(std1),
				XPPMatrixCorrection.takeOffAngleTag(std2), XPPMatrixCorrection.takeOffAngleTag(std3),
				XPPMatrixCorrection.takeOffAngleTag(std4) };
		final RealVector vals = new ArrayRealVector(tags.length);
		vals.setEntry(0, 15.0);
		vals.setEntry(1, 15.0);
		vals.setEntry(2, 15.0);
		vals.setEntry(3, 15.0);
		vals.setEntry(4, 15.0);
		vals.setEntry(5, 15.0);
		vals.setEntry(6, Math.toRadians(40.0));
		vals.setEntry(7, Math.toRadians(40.0));
		vals.setEntry(8, Math.toRadians(40.0));
		vals.setEntry(9, Math.toRadians(40.0));
		vals.setEntry(10, Math.toRadians(40.0));
		vals.setEntry(11, Math.toRadians(40.0));
		final RealVector var = new ArrayRealVector(tags.length);
		var.setEntry(0, Math.pow(0.1, 2.0));
		var.setEntry(1, Math.pow(0.12, 2.0));
		var.setEntry(2, Math.pow(0.12, 2.0));
		var.setEntry(3, Math.pow(0.12, 2.0));
		var.setEntry(4, Math.pow(0.12, 2.0));
		var.setEntry(5, Math.pow(0.12, 2.0));
		var.setEntry(6, Math.pow(Math.toRadians(0.7), 2.0));
		var.setEntry(7, Math.pow(Math.toRadians(0.9), 2.0));
		var.setEntry(8, Math.pow(Math.toRadians(0.7), 2.0));
		var.setEntry(9, Math.pow(Math.toRadians(0.7), 2.0));
		var.setEntry(10, Math.pow(Math.toRadians(0.7), 2.0));
		var.setEntry(11, Math.pow(Math.toRadians(0.7), 2.0));

		final UncertainValues conditions = new UncertainValues(Arrays.asList(tags), vals, var);

		final Map<Composition, CharacteristicXRaySet> stds = new HashMap<>();
		{
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1));
				stds.put(std0, scxr);
			}

			stds.put(std1,
					CharacteristicXRaySet.build(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)));
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2));
				stds.put(std2, scxr);
			}
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1));
				stds.put(std3, scxr);
			}
			{
				final CharacteristicXRaySet scxr = new CharacteristicXRaySet();
				scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1));
				scxr.add(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1));
				stds.put(std4, scxr);
			}
		}

		final Set<Object> outputs = new HashSet<>();
		for (final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet()) {
			final Composition meStd = me.getKey();
			for (final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay()) {
				outputs.add(XPPMatrixCorrection.zaTag(unk, meStd, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, unk, cxr));
				outputs.add(XPPMatrixCorrection.tag2(XPPMatrixCorrection.F_CHI_F, meStd, cxr));
			}
		}

		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unk, stds);
		xpp.trimOutputs(outputs);
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValues inputs = xpp.buildInputs(conditions);
		final long start = System.currentTimeMillis();
		final NamedMultivariateJacobian jac = new NamedMultivariateJacobian(xpp, inputs.getValues());
		final UncertainValues results = UncertainValues.propagate(jac, inputs).sort();
		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection xpp2 = new XPPMatrixCorrection(unk, stds);
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
					if (outTag instanceof XPPMatrixCorrection.ResultTag) {
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

				final UncertainValues resultsD = UncertainValues.propagate(djac, inputs).sort();
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
						new SafeMultivariateNormalDistribution(inputs.getValues(), inputs.getCovariances(), 1.0e-9));
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
					if (tag instanceof XPPMatrixCorrection.ResultTag) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getOutputStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					Table t = new Table();
					t.addRow(Table.th("Tag"), Table.th("V(MonteCarlo)"), Table.th("U(Monte Carlo)"),
							Table.th("V(Analytical)"), Table.th("U(Analytic)"));
					for (final Object tag : xpp.getOutputTags())
						if (tag instanceof XPPMatrixCorrection.ResultTag) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.NORMAL)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(results.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(results.getUncertainty(tag), bnf));
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
