package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.awt.Color;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.MathUtilities;
import gov.nist.microanalysis.roentgen.math.uncertainty.EstimateUncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.Jacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValueEx;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.Layer;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.ZAFMultiLineLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.Shell.Principle;
import gov.nist.microanalysis.roentgen.physics.XRaySet;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat.OutputMode;
import joinery.DataFrame;

public class XPPMatrixCorrection2Test {

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
		// final Composition unk = Composition.parse("Al2SiO5");
		final List<Element> elmsU = Arrays.asList(Element.Aluminum, Element.Silicon, Element.Oxygen);
		final RealVector valsU = new ArrayRealVector(new double[] { 0.3330, 0.1733, 0.4937 });
		final RealVector varsU = new ArrayRealVector(new double[] { 1.0e-6, 0.4e-6, 4.0e-6 });
		final Composition unk = Composition.massFraction("Al<sub>2</sub>SiO<sub>5</sub>", elmsU, valsU, varsU);
		final List<Element> elmsS = Arrays.asList(Element.Silicon, Element.Oxygen);
		final RealVector valsS = new ArrayRealVector(new double[] { 0.4674, 0.5326 });
		final RealVector varsS = new ArrayRealVector(new double[] { 2.0e-6, 0.9e-6 });
		final Composition std = Composition.massFraction("SiO<sub>2</sub>", elmsS, valsS, varsS);

		final StandardMatrixCorrectionDatum stdMcd = new StandardMatrixCorrectionDatum( //
				std, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final Material unkMat = unk.getMaterial();
		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unkMat, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final KRatioLabel krlSi = new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), Method.Measured);
		final KRatioLabel krlO = new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), Method.Measured);

		final Set<KRatioLabel> skrl = new HashSet<>();
		skrl.add(krlSi);
		skrl.add(krlO);

		final UncertainValuesCalculator<EPMALabel> xppI = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);
		
		final Report r = new Report("XPP Report - test1");
		UncertainValuesBase<EPMALabel> fdResults = null;
		try {
			{
				r.addHeader("Compositions");
				r.addSubHeader("Standard");
				r.add(std, Mode.VERBOSE);
				r.addSubHeader("Unknown");
				r.add(unk, Mode.VERBOSE);

				r.addHeader("test1()");
				r.addHTML(xppI.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues<? extends EPMALabel> inputs = xppI.getInputs();
				r.add(inputs);
				final UncertainValuesBase<EPMALabel> aResults = xppI.force().sort();

				
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
				final EPMALabel tagAu = MatrixCorrectionModel2.shellLabel("A", unkMcd, cxr.getInner());
				final EPMALabel tagau = MatrixCorrectionModel2.shellLabel("a", unkMcd, cxr.getInner());
				final EPMALabel tagBu = MatrixCorrectionModel2.shellLabel("B", unkMcd, cxr.getInner());
				final EPMALabel tagbu = MatrixCorrectionModel2.shellLabel("b", unkMcd, cxr.getInner());
				final EPMALabel tagPhi0u = MatrixCorrectionModel2.phi0Label(unkMcd, cxr.getInner());

				final EPMALabel tagAs = MatrixCorrectionModel2.shellLabel("A", stdMcd, cxr.getInner());
				final EPMALabel tagas = MatrixCorrectionModel2.shellLabel("a", stdMcd, cxr.getInner());
				final EPMALabel tagBs = MatrixCorrectionModel2.shellLabel("B", stdMcd, cxr.getInner());
				final EPMALabel tagbs = MatrixCorrectionModel2.shellLabel("b", stdMcd, cxr.getInner());
				final EPMALabel tagPhi0s = MatrixCorrectionModel2.phi0Label(stdMcd, cxr.getInner());

				final EPMALabel tagChiu = MatrixCorrectionModel2.chiLabel(unkMcd, cxr);
				final EPMALabel tagChis = MatrixCorrectionModel2.chiLabel(stdMcd, cxr);
				final EPMALabel tagFChiFu = MatrixCorrectionModel2.FxFLabel(unkMcd, cxr);
				final EPMALabel tagFChiFs = MatrixCorrectionModel2.FxFLabel(stdMcd, cxr);
				final EPMALabel tagZA = MatrixCorrectionModel2.zafLabel(unkMcd, stdMcd, cxr);				
				
				assertEquals(aResults.getEntry(tagAu), 401.654, 0.001);
				assertEquals(aResults.getEntry(tagau), 11189.1680, 0.001);
				assertEquals(aResults.getEntry(tagBu), -526613.8448, 0.001);
				assertEquals(aResults.getEntry(tagbu), 12568.9576, 0.001);
				assertEquals(aResults.getEntry(tagPhi0u), 1.252, 0.001);
				assertEquals(aResults.getEntry(tagAs), 396.744, 0.001);
				assertEquals(aResults.getEntry(tagas), 11314.847, 0.001);
				assertEquals(aResults.getEntry(tagBs), -529359.3292, 0.001);
				assertEquals(aResults.getEntry(tagbs), 12719.694, 0.001);
				assertEquals(aResults.getEntry(tagPhi0s), 1.254, 0.001);
				assertEquals(aResults.getEntry(tagChiu), 2542.429, 0.001);
				assertEquals(aResults.getEntry(tagChis), 1038.418, 0.001);
				assertEquals(aResults.getEntry(tagFChiFu), 0.6331, 0.001);
				assertEquals(aResults.getEntry(tagFChiFs), 0.8206, 0.001);
				assertEquals(aResults.getEntry(tagZA), 0.7797, 0.001);

				final UncertainValuesCalculator<EPMALabel> xppFd = XPPMatrixCorrection2.buildFiniteDifference(skrl,
						unk.getValueMap(MassFraction.class), 0.001, true);

				assertEquals(aResults.getEntry(tagAu), xppFd.getEntry(tagAu), 0.001);
				assertEquals(aResults.getEntry(tagau), xppFd.getEntry(tagau), 0.001);
				assertEquals(aResults.getEntry(tagBu), xppFd.getEntry(tagBu), 0.001);
				assertEquals(aResults.getEntry(tagbu), xppFd.getEntry(tagbu), 0.001);
				assertEquals(aResults.getEntry(tagPhi0u), xppFd.getEntry(tagPhi0u), 0.001);

				assertEquals(aResults.getEntry(tagAs), xppFd.getEntry(tagAs), 0.001);
				assertEquals(aResults.getEntry(tagas), xppFd.getEntry(tagas), 0.001);
				assertEquals(aResults.getEntry(tagBs), xppFd.getEntry(tagBs), 0.001);
				assertEquals(aResults.getEntry(tagbs), xppFd.getEntry(tagbs), 0.001);
				assertEquals(aResults.getEntry(tagPhi0s), xppFd.getEntry(tagPhi0s), 0.001);

				assertEquals(aResults.getEntry(tagChiu), xppFd.getEntry(tagChiu), 0.001);
				assertEquals(aResults.getEntry(tagChis), xppFd.getEntry(tagChis), 0.001);
				assertEquals(aResults.getEntry(tagFChiFu), xppFd.getEntry(tagFChiFu), 0.001);
				assertEquals(aResults.getEntry(tagFChiFs), xppFd.getEntry(tagFChiFs), 0.001);
				assertEquals(aResults.getEntry(tagZA), xppFd.getEntry(tagZA), 0.001);

				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = xppI.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xppI.getOutputLabels()) {
					final UncertainValue uv = outVals.get(outTag);
					valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
							Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
							Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				assertTrue(xppI.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = xppI.getJacobian().get();
				assertTrue(xppFd.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> fdJac = xppFd.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aJac.toHTML(Mode.NORMAL));
					System.out.println("Jacobian(Finite Difference)");
					System.out.println(fdJac.toHTML(Mode.NORMAL));
				}

				// final Object unkCompTag = new MatrixCorrectionModel2.MaterialBasedLabel("J",
				// unkMat);
				// assertEquals(aCalc.getJacobianEntry(unkCompTag,
				// MaterialLabel.buildMassFractionTag(unkMat, Element.Oxygen)), -0.027565,
				// 0.00001);
				// assertEquals(aCalc.getJacobianEntry(unkCompTag,
				// MatrixCorrectionModel2.meanIonizationLabel(Element.Oxygen)), 0.609601,
				// 0.00001);

				fdResults = xppFd.sort();
				r.addImage(fdResults.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, fdResults, L2C, 8),
						"Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1);
				final UncertainValuesCalculator<EPMALabel> xppMc = XPPMatrixCorrection2
						.buildMonteCarlo(skrl, unk.getValueMap(MassFraction.class), MC_ITERATIONS, true).force();
				r.addHeader("Monte Carlo Results");
				r.add(xppMc);
				r.addHeader("Inputs");
				r.add(xppMc.getInputs());
				final UncertainValuesBase<EPMALabel> aResults = xppI.sort();

				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("A", unkMcd, cxr.getInner())), 2366.373,
						0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("a", unkMcd, cxr.getInner())),
						11052.873, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("B", unkMcd, cxr.getInner())),
						-1460552.6291, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("b", unkMcd, cxr.getInner())),
						11681.220, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.phi0Label(unkMcd, cxr.getInner())), 1.258, 0.001);

				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("A", stdMcd, cxr.getInner())), 2307.215,
						0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("a", stdMcd, cxr.getInner())),
						11178.191, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("B", stdMcd, cxr.getInner())),
						-1459152.312, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.shellLabel("b", stdMcd, cxr.getInner())),
						11822.223, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.phi0Label(stdMcd, cxr.getInner())), 1.26, 0.001);

				assertEquals(aResults.getEntry(MatrixCorrectionModel2.chiLabel(unkMcd, cxr)), 5836.018, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.chiLabel(stdMcd, cxr)), 6414.025, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr)), 0.367, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.FxFLabel(stdMcd, cxr)), 0.344, 0.001);
				assertEquals(aResults.getEntry(MatrixCorrectionModel2.zafLabel(unkMcd, stdMcd, cxr)), 1.0792, 0.001);

				r.addHeader("Analyic Results");
				r.add(aResults);

				final UncertainValuesBase<EPMALabel> resultsMc = xppMc.sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) xppMc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("Monte Carlo Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : resultsMc.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				r.addSubHeader("Phi0");
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel tag : xppI.getOutputLabels()) {
					r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
					r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
				}

				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					final BasicNumberFormat bnf2 = new BasicNumberFormat("0.0000");
					for (final EPMALabel tag : xppI.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(aResults.getUncertainty(tag), bnf2),
									MathUtilities.td(fdResults.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(fdResults.getUncertainty(tag), bnf2));
						}
					for (final EPMALabel tag : xppI.getOutputLabels())
						if (tag instanceof KRatioLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(aResults.getUncertainty(tag), bnf2),
									MathUtilities.td(fdResults.getValue(tag).doubleValue(), bnf2),
									MathUtilities.td(fdResults.getUncertainty(tag), bnf2));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (final Exception e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
	}

	public static Composition buildK412(
			final boolean combine
			) //
			throws ArgumentException, ParseException {
		if (combine)
			return Composition.combine("K412", false, //
					Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
					Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
					Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
					Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
					Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));
		else {
			final Map<Element, Number> res = new HashMap<Element, Number>();
			res.put(Element.Silicon, new UncertainValue(0.21226219744446, 0.00359924888862));
			res.put(Element.Iron, new UncertainValue(0.07726410130474, 0.00139914871578));
			res.put(Element.Magnesium, new UncertainValue(0.11855685731088, 0.001507589742));
			res.put(Element.Calcium, new UncertainValue(0.11034825437848, 0.00107203615005));
			res.put(Element.Aluminum, new UncertainValue(0.0494320434, 0.0015348279));
			res.put(Element.Oxygen, new UncertainValue(0.43003654616144, 0.00728714860355));
			return Composition.massFraction("K412", res);
		}
	}

	public static Composition buildK411(
			final boolean combine
			) //
			throws ArgumentException, ParseException {
		if (combine)
			return Composition.combine("K411", false, //
					Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
					Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
					Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
					Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015)));
		else {
			final Map<Element, Number> res = new HashMap<Element, Number>();
			res.put(Element.Silicon, new UncertainValue(0.25190067871134, 0.00448737523776));
			res.put(Element.Iron, new UncertainValue(0.11255374113608, 0.00209872307367));
			res.put(Element.Magnesium, new UncertainValue(0.09117902759616, 0.0012060717936));
			res.put(Element.Calcium, new UncertainValue(0.11070559976183, 0.00107203615005));
			res.put(Element.Oxygen, new UncertainValue(0.42346095279459, 0.00693579374492));
			return Composition.massFraction("K411", res);
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
		final boolean combined = false;
		final Composition std = buildK411(combined), unk = buildK412(combined);

		final StandardMatrixCorrectionDatum stdMcd = new StandardMatrixCorrectionDatum( //
				std, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), Method.Measured));

		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final Report r = new Report("XPP Report - test2");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test2()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				r.add(inputs);

				final long start = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> xppI = new UncertainValuesCalculator<>(xpp, inputs);
				final UncertainValuesBase<EPMALabel> aResults = xppI.sort();
				System.out.println("Full Timing(2) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
						unk.getValueMap(MassFraction.class), true);
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}

				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				// Can't seem to diagnose problem with input to
				// SafeMultivariateNormalDistribution
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				r.add(inputs);
				final UncertainValuesBase<EPMALabel> aResults = UncertainValuesBase.propagateAnalytical(xpp, inputs)
						.sort();
				r.addHeader("Analytic Results");
				r.add(aResults);

				final UncertainValuesBase<EPMALabel> ordered = inputs.reorder(xpp.getInputLabels());

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp, ordered);
				uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
				final UncertainValuesBase<EPMALabel> resultsMc = uvc.sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final EPMALabel tag : xpp.getOutputLabels()) {
					if (tag instanceof ZAFMultiLineLabel) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}
				}

				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
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

		final Composition std = Composition.parse("Mg");

		final boolean combined = false;
		final Composition unk = buildK412(combined);

		final StandardMatrixCorrectionDatum stdMcd = new StandardMatrixCorrectionDatum( //
				std, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, stdMcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2)), Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));

		assertEquals(xpp.getOutputDimension(), outputs.size());
		final Report r = new Report("XPP Report - test3");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test3()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				r.add(inputs);
				final long start = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> xppI = new UncertainValuesCalculator<>(xpp, inputs);
				final UncertainValuesBase<EPMALabel> aResults = xppI.sort();
				System.out.println("Timing(3) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();

				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

				final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
						unk.getValueMap(MassFraction.class), true);
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD.extract(aResults.getLabels()), L2C, 8),
						"Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				r.add(xpp);
				r.addHeader("Inputs");
				final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				r.add(inputs);
				final UncertainValuesBase<EPMALabel> aResults = UncertainValuesBase.propagateAnalytical(xpp, inputs)
						.sort();
				r.addHeader("Analytic Results");
				r.add(aResults);

				final UncertainValuesBase<EPMALabel> ordered = inputs.reorder(xpp.getInputLabels());

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp, ordered);
				uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
				final UncertainValuesBase<EPMALabel> resultsMc = UncertainValues.asUncertainValues(uvc).sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("MC Results");
				r.add(resultsMc);
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsMc, L2C, 8),
						"Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final EPMALabel tag : xpp.getOutputLabels()) {
					if (tag instanceof ZAFMultiLineLabel) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}
				}

				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (final Exception e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
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
		final boolean combined = true;
		final Composition std0 = buildK411(combined), unk = buildK412(false);
		final Composition std1 = Composition.parse("Al");
		Report.dump(std0, Mode.VERBOSE);

		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Iron, XRayTransition.LA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA2)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Calcium, XRayTransition.L3M1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1)), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std1Mcd,
				new ElementXRaySet(CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1)), Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl) {
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}

		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final long start = System.currentTimeMillis();

		final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);
		final UncertainValuesBase<EPMALabel> aResults = aCalc.sort();

		System.out.println("Trimmed Timing (4) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection2 xpp2 = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValuesBase<EPMALabel> inputs2 = xpp2.buildInput(unk.getMaterial());
		xpp2.addConstraints(xpp2.buildConstraints(inputs));
		xpp2.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesCalculator<EPMALabel> fullResults = UncertainValuesBase.propagateAnalytical(xpp2, inputs2);
		
		System.out.println("Full Timing (4) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		Set<EPMALabel> both = new HashSet<>(aResults.getLabels());
		both.retainAll(fullResults.getLabels());
		for (final EPMALabel outTag : both) {
			for (final EPMALabel outTag2 : both)
				if (outTag == outTag2) {
					final double v1 = aResults.getEntry(outTag);
					final double v2 = fullResults.getEntry(outTag);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (aResults.getVariance(outTag) > 1.0e-8)
						assertEquals(aResults.getVariance(outTag), fullResults.getVariance(outTag),
								0.01 * Math.abs(aResults.getVariance(outTag)));
				} else {
					final double cov1 = aResults.getCovariance(outTag, outTag2);
					final double cov2 = fullResults.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}
		
		final Report r = new Report("XPP Report - test4()");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test4()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs.sort());
				// r.add(inputs);
				r.addHeader("Results");
				r.add(aResults.sort());
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();

				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Uncertainty matrix");

				final long start3 = System.currentTimeMillis();

				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();
				System.out.println(
						"Trimmed Delta Timing (4) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}
				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

			}
			if (MC_ITERATIONS > 0) {
				// Can't seem to diagnose problem with input to
				// SafeMultivariateNormalDistribution
				inputs.validateCovariance();
				Report.dump(inputs.blockDiagnonalize().toSimpleHTML(new BasicNumberFormat("0.00E0")));

				final UncertainValuesBase<EPMALabel> resultsMc = UncertainValuesBase
						.propagateMonteCarlo(xpp, inputs, MC_ITERATIONS).sort();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}

				r.addHeader("Monte Carlo Results");
				r.add(resultsMc);
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsMc, L2C, 8),
						"Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				/*
				 * for (final Object tag : xpp.getOutputLabels()) { if (tag instanceof
				 * MatrixCorrectionLabel) { r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
				 * r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf)); } }
				 */
				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");
			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
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
		final Composition unk = buildK412(combined);

		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std4Mcd = new StandardMatrixCorrectionDatum( //
				std4, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std2Mcd, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std2Mcd, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA2),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Calcium, XRayTransition.L3M1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.LA1),
				Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final long start = System.currentTimeMillis();

		final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);
		final UncertainValuesBase<EPMALabel> aResults = aCalc.sort();
		System.out.println("Trimmed Timing (5) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();
		final XPPMatrixCorrection2 xpp2 = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValuesBase<EPMALabel> inputs2 = xpp2.buildInput(unk.getMaterial());
		xpp2.addConstraints(xpp2.buildConstraints(inputs));
		xpp2.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesCalculator<EPMALabel> fullResults = UncertainValuesBase.propagateAnalytical(xpp2, inputs2);

		System.out.println("Full Timing (5) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		Set<EPMALabel> both = new HashSet<>(aResults.getLabels());
		both.retainAll(fullResults.getLabels());
		for (final EPMALabel outTag : both) {
			for (final EPMALabel outTag2 : both)
				if (outTag == outTag2) {
					final double v1 = aResults.getEntry(outTag);
					final double v2 = fullResults.getEntry(outTag);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (aResults.getVariance(outTag) > 1.0e-8)
						assertEquals(aResults.getVariance(outTag), fullResults.getVariance(outTag),
								0.01 * Math.abs(aResults.getVariance(outTag)));
				} else {
					final double cov1 = aResults.getCovariance(outTag, outTag2);
					final double cov2 = fullResults.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test5()");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test5()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();
				System.out.println(
						"Trimmed Delta Timing (5) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp, inputs);
				uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
				final UncertainValuesBase<EPMALabel> resultsMc = UncertainValues.asUncertainValues(uvc).sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsMc, L2C, 8),
						"Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final EPMALabel tag : xpp.getOutputLabels()) {
					if (tag instanceof ZAFMultiLineLabel) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
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

		final Composition unk = buildK412(false);

		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std4Mcd = new StandardMatrixCorrectionDatum( //
				std4, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std2Mcd, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std2Mcd, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA2),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Calcium, XRayTransition.L3M1),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.LA1),
				Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
			}
		}

		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final long start = System.currentTimeMillis();
		final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);
		final UncertainValuesBase<EPMALabel> aResults = aCalc.sort();
		System.out.println("Trimmed Timing (6) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();

		final XPPMatrixCorrection2 xpp2 = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValuesBase<EPMALabel> inputs2 = xpp2.buildInput(unk.getMaterial());
		xpp2.addConstraints(xpp2.buildConstraints(inputs));
		xpp2.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesCalculator<EPMALabel> fullResults = UncertainValuesBase.propagateAnalytical(xpp2, inputs2);

		System.out.println("Full Timing (6) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		Set<EPMALabel> both = new HashSet<>(aResults.getLabels());
		both.retainAll(fullResults.getLabels());
		for (final EPMALabel outTag : both) {
			for (final EPMALabel outTag2 : both)
				if (outTag == outTag2) {
					final double v1 = aResults.getEntry(outTag);
					final double v2 = fullResults.getEntry(outTag);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (aResults.getVariance(outTag) > 1.0e-8)
						assertEquals(aResults.getVariance(outTag), fullResults.getVariance(outTag),
								0.01 * Math.abs(aResults.getVariance(outTag)));
				} else {
					final double cov1 = aResults.getCovariance(outTag, outTag2);
					final double cov2 = fullResults.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}
		
		final Report r = new Report("XPP Report - Test6()");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test6()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();
				System.out.println(
						"Trimmed Delta Timing (6) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp, inputs);
				uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
				final UncertainValuesBase<EPMALabel> resultsMc = UncertainValues.asUncertainValues(uvc).sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsMc, L2C, 8),
						"Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final EPMALabel tag : xpp.getOutputLabels()) {
					if (tag instanceof ZAFMultiLineLabel) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
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
		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.7), //
				MatrixCorrectionDatum.roughness(10.0, 3.2));
		// Ba, Ti, Si, O
		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9), //
				MatrixCorrectionDatum.roughness(10.0, 3.6));
		// Zn
		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7), //
				MatrixCorrectionDatum.roughness(10.0, 7.14));
		// Zr
		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7), //
				MatrixCorrectionDatum.roughness(10.0, 5.62));

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7), //
				MatrixCorrectionDatum.roughness(10.0, 3.5));

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1),
				Method.Measured));
		// Ba, Ti, Si, O
		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Barium, XRayTransition.LA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Titanium, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
				Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std1Mcd, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
				Method.Measured));
		// Zn
		skrl.add(new KRatioLabel(unkMcd, std2Mcd, ElementXRaySet.singleton(Element.Zinc, XRayTransition.KA1),
				Method.Measured));
		// Zr
		skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Zirconium, XRayTransition.LA1),
				Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl) {
			final StandardMatrixCorrectionDatum meStd = krl.getStandard();
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.zLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.aLabel(unkMcd, meStd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
				outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
				outputs.add(new KRatioLabel(unkMcd, meStd, cxr, Method.Calculated));
			}
		}

		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final long start = System.currentTimeMillis();
		final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);
		final UncertainValuesBase<EPMALabel> aResults = aCalc.sort();
		System.out.println("Trimmed Timing (7) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();

		final XPPMatrixCorrection2 xpp2 = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValuesBase<EPMALabel> inputs2 = xpp2.buildInput(unk.getMaterial());
		xpp2.addConstraints(xpp2.buildConstraints(inputs));
		xpp2.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesCalculator<EPMALabel> fullResults = UncertainValuesBase.propagateAnalytical(xpp2, inputs2);

		System.out.println("Full Timing (7) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		Set<EPMALabel> both = new HashSet<>(aResults.getLabels());
		both.retainAll(fullResults.getLabels());
		for (final EPMALabel outTag : both) {
			for (final EPMALabel outTag2 : both)
				if (outTag == outTag2) {
					final double v1 = aResults.getEntry(outTag);
					final double v2 = fullResults.getEntry(outTag);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (aResults.getVariance(outTag) > 1.0e-8)
						assertEquals(aResults.getVariance(outTag), fullResults.getVariance(outTag),
								0.01 * Math.abs(aResults.getVariance(outTag)));
				} else {
					final double cov1 = aResults.getCovariance(outTag, outTag2);
					final double cov2 = fullResults.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test7()");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test7()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();
				System.out.println(
						"Trimmed Delta Timing (7) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

				final Table res = new Table();
				res.addRow(Table.th("Line"), Table.th("ZA"), Table.th("Z"), Table.th("A"));
				final BasicNumberFormat bnf2 = new BasicNumberFormat("0.000");

				for (final KRatioLabel krl : skrl) {
					final StandardMatrixCorrectionDatum meStd = krl.getStandard();
					outputs.add(MatrixCorrectionModel2.zafLabel(krl));
					for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
						final EPMALabel zaTag = MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr);
						final EPMALabel zTag = MatrixCorrectionModel2.zLabel(unkMcd, meStd, cxr);
						final EPMALabel aTag = MatrixCorrectionModel2.aLabel(unkMcd, meStd, cxr);
						final UncertainValue za = aResults.getUncertainValue(zaTag);
						final UncertainValue a = aResults.getUncertainValue(aTag);
						final UncertainValue z = aResults.getUncertainValue(zTag);
						res.addRow(Table.td(cxr), //
								Table.td(bnf2.formatHTML(za, OutputMode.ValuePlusUncertainty)), //
								Table.td(bnf2.formatHTML(z, OutputMode.ValuePlusUncertainty)), //
								Table.td(bnf2.formatHTML(a, OutputMode.ValuePlusUncertainty)));
					}
				}
				r.addHeader("ZAF results");
				r.add(res, Mode.NORMAL);

				final Table tk = new Table();
				for (final EPMALabel tag : xpp.getOutputLabels())
					if (tag instanceof KRatioLabel) {
						tk.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
								MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf2),
								MathUtilities.td(aResults.getUncertainty(tag), bnf2),
								MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
								MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
					}
				r.addHeader("K-ratio");
				r.add(tk, Mode.NORMAL);

				for (final KRatioLabel krl : skrl) {
					final StandardMatrixCorrectionDatum meStd = krl.getStandard();
					for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
						final EPMALabel zaTag = MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr);
						final EPMALabel fUnkTag = MatrixCorrectionModel2.FxFLabel(unkMcd, cxr);
						final EPMALabel fStdTag = MatrixCorrectionModel2.FxFLabel(meStd, cxr);
						final double za = aResults.getEntry(zaTag);
						final double a = aResults.getEntry(fUnkTag) / aResults.getEntry(fStdTag);
						final double z = za / a;
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
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
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

		final Composition unk = buildK412(false);

		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std4Mcd = new StandardMatrixCorrectionDatum( //
				std4, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, std0Mcd, XRaySet.build(Element.Silicon, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd, XRaySet.build(Element.Oxygen, Principle.K, 0.001), Method.Measured));

		skrl.add(
				new KRatioLabel(unkMcd, std1Mcd, XRaySet.build(Element.Aluminum, Principle.K, 0.001), Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std2Mcd, XRaySet.build(Element.Magnesium, Principle.K, 0.001),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std3Mcd, XRaySet.build(Element.Calcium, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std3Mcd, XRaySet.build(Element.Calcium, Principle.L, 0.001), Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std4Mcd, XRaySet.build(Element.Iron, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std4Mcd, XRaySet.build(Element.Iron, Principle.L, 0.001), Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl)
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
		assertEquals(xpp.getOutputDimension(), outputs.size());
		final UncertainValuesBase<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final long start = System.currentTimeMillis();

		final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
				unk.getValueMap(MassFraction.class), true);

		final UncertainValuesBase<EPMALabel> aResults = aCalc.sort();
		System.out.println("Trimmed Timing (8) = " + Long.toString(System.currentTimeMillis() - start) + " ms");

		final long start2 = System.currentTimeMillis();

		final XPPMatrixCorrection2 xpp2 = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValuesBase<EPMALabel> inputs2 = xpp2.buildInput(unk.getMaterial());
		xpp2.addConstraints(xpp2.buildConstraints(inputs));
		xpp2.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesBase<EPMALabel> fullResults = UncertainValuesBase.propagateAnalytical(xpp2, inputs2).sort();
		System.out.println("Full Timing (8) = " + Long.toString(System.currentTimeMillis() - start2) + " ms");

		// Test untrimmed vs trimmed
		Set<EPMALabel> both = new HashSet<>(aResults.getLabels());
		both.retainAll(fullResults.getLabels());
		for (final EPMALabel outTag : both) {
			for (final EPMALabel outTag2 : both)
				if (outTag == outTag2) {
					final double v1 = aResults.getEntry(outTag);
					final double v2 = fullResults.getEntry(outTag);
					assertEquals(v1, v2, 0.01 * Math.max(Math.abs(v1), Math.abs(v2)));
					if (aResults.getVariance(outTag) > 1.0e-8)
						assertEquals(aResults.getVariance(outTag), fullResults.getVariance(outTag),
								0.01 * Math.abs(aResults.getVariance(outTag)));
				} else {
					final double cov1 = aResults.getCovariance(outTag, outTag2);
					final double cov2 = fullResults.getCovariance(outTag, outTag2);
					assertEquals(cov1, cov2, 0.01 * Math.max(Math.abs(cov1), Math.abs(cov2)));
				}
		}

		final Report r = new Report("XPP Report - Test8()");
		UncertainValuesBase<EPMALabel> resultsD = null;
		try {
			{
				r.addHeader("test8()");
				r.addHTML(xpp.toHTML(Mode.NORMAL));
				r.addHeader("Inputs");
				r.add(inputs);
				r.addHeader("Results");
				r.add(aResults);
				r.addHeader("Uncertain Values (relative to inputs, trimmed)");
				final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp, inputs)
						.getOutputValues();
				final Table valTable = new Table();
				valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
						Table.td("Value (Verbose)"));
				final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
				for (final EPMALabel outTag : xpp.getOutputLabels()) {
					if (outTag instanceof ZAFMultiLineLabel) {
						final UncertainValue uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
				}
				r.addHTML(valTable.toHTML(Mode.NORMAL));

				r.addHeader("Covariance matrix");
				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));
				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Results uncertainty matrix");

				final long start3 = System.currentTimeMillis();
				final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
						.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();
				System.out.println(
						"Trimmed Delta Timing (8) = " + Long.toString(System.currentTimeMillis() - start3) + " ms");

				assertTrue(aCalc.getJacobian().isPresent());
				final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
				final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

				for (final EPMALabel output : aJac.getOutputLabels())
					for (final EPMALabel input : aJac.getInputLabels())
						if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
							if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
							0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
									Math.abs(fdJac.getEntry(input, output)))) {
								checkEquals(input, output, aJac.getEntry(input, output),
										fdJac.getEntry(input, output),
										0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
												Math.abs(fdJac.getEntry(input, output))));
							}
						}

				resultsD = fdCalc.sort();
				r.addImage(resultsD.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsD, L2C, 8),
						"Comparing uncertainty matrix");

				if (DUMP) {
					System.out.println("Results");
					System.out.println(aResults.toCSV());

					System.out.println("Jacobian");
					System.out.println(aCalc.toCSV());
					System.out.println("Jacobian(estimated)");
					System.out.println(fdCalc.toCSV());
				}
			}
			if (MC_ITERATIONS > 0) {
				r.addHeader("Monte Carlo Results");
				Report.dump(inputs.blockDiagnonalize().toSimpleHTML(new BasicNumberFormat("0.00E0")));

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp, inputs);
				uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
				final UncertainValuesBase<EPMALabel> resultsMc = UncertainValues.asUncertainValues(uvc).sort();
				final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
						.getUncertainValues();

				if (DUMP) {
					System.out.println("MC Results");
					System.out.println(resultsMc.toCSV());
				}
				r.add(resultsMc);

				final StringBuffer sb = new StringBuffer();
				for (final EPMALabel tag : aResults.getLabels()) {
					if (sb.length() > 0)
						sb.append(",");
					sb.append(HTML.toHTML(tag, Mode.TERSE));
				}
				r.addHTML(HTML.p(sb.toString()));

				r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Analytical result matrix");
				r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "MC result matrix");
				r.addImage(UncertainValuesBase.compareAsBitmap(aResults, resultsMc, L2C, 8),
						"Comparing analytical with MC");

				final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
				for (final EPMALabel tag : xpp.getOutputLabels()) {
					if (tag instanceof ZAFMultiLineLabel) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}
				}
				{
					r.addHeader("Compare MC to Analytical");
					final Table t = new Table();
					t.addRow(Table.th("Tag"), //
							Table.th("V(MonteCarlo)"), //
							Table.th("U(Monte Carlo)"), //
							Table.th("V(Analytical)"), //
							Table.th("U(Analytic)"), //
							Table.th("V(Delta)"), //
							Table.th("U(Delta)"));
					for (final EPMALabel tag : xpp.getOutputLabels())
						if (tag instanceof ZAFMultiLineLabel) {
							t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
									MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf), //
									MathUtilities.td(resultsMc.getUncertainty(tag), bnf), //
									MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(aResults.getUncertainty(tag), bnf),
									MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf),
									MathUtilities.td(resultsD.getUncertainty(tag), bnf));
						}
					r.add(t);
				}
				r.addHeader("Done!");

			}
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);
	}

	@Test
	public void testXPP9() throws ArgumentException, ParseException, IOException {
		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al");
		final Composition std2 = Composition.parse("Mg");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final boolean combined = false;

		final Composition unk = buildK412(combined);

		final StandardMatrixCorrectionDatum std0Mcd = new StandardMatrixCorrectionDatum( //
				std0, //
				new UncertainValue(15.0, 0.1), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
				std1, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.9) //
		);

		final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
				std2, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std3Mcd = new StandardMatrixCorrectionDatum( //
				std3, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final StandardMatrixCorrectionDatum std4Mcd = new StandardMatrixCorrectionDatum( //
				std4, //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
				unk.getMaterial(), //
				new UncertainValue(15.0, 0.12), //
				UncertainValue.toRadians(40.0, 0.7) //
		);

		final Set<KRatioLabel> skrl = new HashSet<>();

		skrl.add(new KRatioLabel(unkMcd, std0Mcd, XRaySet.build(Element.Silicon, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std0Mcd, XRaySet.build(Element.Oxygen, Principle.K, 0.001), Method.Measured));

		skrl.add(
				new KRatioLabel(unkMcd, std1Mcd, XRaySet.build(Element.Aluminum, Principle.K, 0.001), Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std2Mcd, XRaySet.build(Element.Magnesium, Principle.K, 0.001),
				Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std3Mcd, XRaySet.build(Element.Calcium, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std3Mcd, XRaySet.build(Element.Calcium, Principle.L, 0.001), Method.Measured));

		skrl.add(new KRatioLabel(unkMcd, std4Mcd, XRaySet.build(Element.Iron, Principle.K, 0.001), Method.Measured));
		skrl.add(new KRatioLabel(unkMcd, std4Mcd, XRaySet.build(Element.Iron, Principle.L, 0.001), Method.Measured));

		final Set<EPMALabel> outputs = new HashSet<>();
		for (final KRatioLabel krl : skrl)
			outputs.add(MatrixCorrectionModel2.zafLabel(krl));
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, Collections.emptyList());
		final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

		final UncertainValuesBase<EPMALabel> aResults = UncertainValuesBase.propagateAnalytical(xpp,
				inputs);

		final DataFrame<Double> df = xpp.computePhiRhoZCurve(aResults.getValueMap(), 1.201e-3, 2.0e-5, 0.9);
		df.writeCsv("C:\\Users\\nicho\\OneDrive\\Desktop\\prz412.csv");
	}

	public void checkEquals(
			final EPMALabel object, final EPMALabel object2, final double v1, final double v2, final double dv
			) {
		if (Math.abs(v2 - v1) > Math.abs(dv)) {
			System.err.println(object.toString() + " and " + object2.toString() + " at (" + v1 + "," + v2 + ")");
			// assertEquals(v1, v2, dv);
		}
	}

	/**
	 * Computes Si and O in Al2SiO5 using SiO2
	 *
	 * @throws ArgumentException
	 * @throws ParseException
	 * @throws IOException
	 */
	@Test
	public void testXPP10() //
			throws ArgumentException, ParseException, IOException {
		// final Composition unk = Composition.parse("Al2SiO5");
		ExplicitMeasurementModel.sDump = null;
		try {
			final List<Element> elmsU = Arrays.asList(Element.Aluminum, Element.Silicon, Element.Oxygen);
			final RealVector valsU = new ArrayRealVector(new double[] { 0.3330, 0.1733, 0.4937 });
			final RealVector varsU = new ArrayRealVector(new double[] { 1.0e-6, 0.4e-6, 4.0e-6 });
			final Composition unk = Composition.massFraction("Al<sub>2</sub>SiO<sub>5</sub>", elmsU, valsU, varsU);
			final List<Element> elmsS = Arrays.asList(Element.Silicon, Element.Oxygen);
			final RealVector valsS = new ArrayRealVector(new double[] { 0.4674, 0.5326 });
			final RealVector varsS = new ArrayRealVector(new double[] { 2.0e-6, 0.9e-6 });
			final Composition std1 = Composition.massFraction("SiO<sub>2</sub>", elmsS, valsS, varsS);
			final Composition std2 = Composition.pureElement(Element.Aluminum);

			final StandardMatrixCorrectionDatum std1Mcd = new StandardMatrixCorrectionDatum( //
					std1, //
					new UncertainValue(15.0, 0.1), //
					UncertainValue.toRadians(40.0, 0.9), MatrixCorrectionDatum.roughness(10.0, 2.5),
					Layer.carbonCoating(new UncertainValue(10.0, 3.0))//
			);

			final StandardMatrixCorrectionDatum std2Mcd = new StandardMatrixCorrectionDatum( //
					std2, //
					new UncertainValue(15.0, 0.1), //
					UncertainValue.toRadians(40.0, 0.9), //
					MatrixCorrectionDatum.roughness(20.0, 2.5), Layer.carbonCoating(new UncertainValue(10.0, 3.0))//
			);

			final Material unkMat = unk.getMaterial();
			final UnknownMatrixCorrectionDatum unkMcd = new UnknownMatrixCorrectionDatum( //
					unkMat, //
					new UncertainValue(15.0, 0.12), //
					UncertainValue.toRadians(40.0, 0.7), //
					MatrixCorrectionDatum.roughness(20.0, 2.5), Layer.carbonCoating(new UncertainValue(12.0, 2.0)));

			final ElementXRaySet exrsSi = XRaySet.build(Element.Silicon, Principle.K, 0.01);
			final KRatioLabel krlSi = new KRatioLabel(unkMcd, std1Mcd, exrsSi, Method.Measured);
			final ElementXRaySet exrsO = XRaySet.build(Element.Oxygen, Principle.K, 0.01);
			final KRatioLabel krlO = new KRatioLabel(unkMcd, std1Mcd, exrsO, Method.Measured);
			final ElementXRaySet exrsAl = XRaySet.build(Element.Aluminum, Principle.K, 0.01);
			final KRatioLabel krlAl = new KRatioLabel(unkMcd, std2Mcd, exrsAl, Method.Measured);

			final Set<KRatioLabel> skrl = new HashSet<>();
			skrl.add(krlSi);
			skrl.add(krlO);
			skrl.add(krlAl);

			final List<EPMALabel> outputs = XPPMatrixCorrection2.buildDefaultOutputs(skrl);

			final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);

			final EPMALabel tagAu = MatrixCorrectionModel2.shellLabel("A", unkMcd, cxr.getInner());
			outputs.add(tagAu);
			final EPMALabel tagau = MatrixCorrectionModel2.shellLabel("a", unkMcd, cxr.getInner());
			outputs.add(tagau);
			final EPMALabel tagBu = MatrixCorrectionModel2.shellLabel("B", unkMcd, cxr.getInner());
			outputs.add(tagBu);
			final EPMALabel tagbu = MatrixCorrectionModel2.shellLabel("b", unkMcd, cxr.getInner());
			outputs.add(tagbu);
			final EPMALabel tagPhi0u = MatrixCorrectionModel2.phi0Label(unkMcd, cxr.getInner());
			outputs.add(tagPhi0u);

			final EPMALabel tagAs = MatrixCorrectionModel2.shellLabel("A", std1Mcd, cxr.getInner());
			outputs.add(tagAs);
			final EPMALabel tagas = MatrixCorrectionModel2.shellLabel("a", std1Mcd, cxr.getInner());
			outputs.add(tagas);
			final EPMALabel tagBs = MatrixCorrectionModel2.shellLabel("B", std1Mcd, cxr.getInner());
			outputs.add(tagBs);
			final EPMALabel tagbs = MatrixCorrectionModel2.shellLabel("b", std1Mcd, cxr.getInner());
			outputs.add(tagbs);
			final EPMALabel tagPhi0s = MatrixCorrectionModel2.phi0Label(std1Mcd, cxr.getInner());
			outputs.add(tagPhi0s);

			final EPMALabel tagChiu = MatrixCorrectionModel2.chiLabel(unkMcd, cxr);
			outputs.add(tagChiu);
			final EPMALabel tagChis = MatrixCorrectionModel2.chiLabel(std1Mcd, cxr);
			outputs.add(tagChis);
			final EPMALabel tagFChiFu = MatrixCorrectionModel2.FxFLabel(unkMcd, cxr);
			outputs.add(tagFChiFu);
			final EPMALabel tagFChiFs = MatrixCorrectionModel2.FxFLabel(std1Mcd, cxr);
			outputs.add(tagFChiFs);
			final ZAFMultiLineLabel zafLabel = MatrixCorrectionModel2.zafLabel(unkMcd, std1Mcd, exrsSi);
			// outputs.add(zafLabel); Already in...
			final EPMALabel tagZAO = MatrixCorrectionModel2.zafLabel(unkMcd, std1Mcd, exrsO);
			// outputs.add(tagZAO); Already in...

			final CharacteristicXRay cxr2 = CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1);
			final EPMALabel tagZAAl = MatrixCorrectionModel2.zafLabel(unkMcd, std2Mcd, exrsAl);
			// outputs.add(tagZAAl);
			final EPMALabel shellLabelA = MatrixCorrectionModel2.shellLabel("A", unkMcd, cxr2.getInner());
			outputs.add(shellLabelA);
			final EPMALabel shellLabela = MatrixCorrectionModel2.shellLabel("a", unkMcd, cxr2.getInner());
			outputs.add(shellLabela);
			final EPMALabel shellLabelB = MatrixCorrectionModel2.shellLabel("B", unkMcd, cxr2.getInner());
			outputs.add(shellLabelB);
			final EPMALabel shellLabelb = MatrixCorrectionModel2.shellLabel("b", unkMcd, cxr2.getInner());
			outputs.add(shellLabelb);
			final EPMALabel phi0Label = MatrixCorrectionModel2.phi0Label(unkMcd, cxr2.getInner());
			outputs.add(phi0Label);

			final EPMALabel shellLabelA2 = MatrixCorrectionModel2.shellLabel("A", std1Mcd, cxr2.getInner());
			outputs.add(shellLabelA2);
			final EPMALabel shellLabela2 = MatrixCorrectionModel2.shellLabel("a", std1Mcd, cxr2.getInner());
			outputs.add(shellLabela2);
			final EPMALabel shellLabelB2 = MatrixCorrectionModel2.shellLabel("B", std1Mcd, cxr2.getInner());
			outputs.add(shellLabelB2);
			final EPMALabel shellLabelb2 = MatrixCorrectionModel2.shellLabel("b", std1Mcd, cxr2.getInner());
			outputs.add(shellLabelb2);
			final EPMALabel phi0Label2 = MatrixCorrectionModel2.phi0Label(std1Mcd, cxr2.getInner());
			outputs.add(phi0Label2);

			final EPMALabel chiLabelu = MatrixCorrectionModel2.chiLabel(unkMcd, cxr2);
			outputs.add(chiLabelu);
			final EPMALabel chiLabels = MatrixCorrectionModel2.chiLabel(std1Mcd, cxr2);
			outputs.add(chiLabels);
			final EPMALabel fxFLabelu = MatrixCorrectionModel2.FxFLabel(unkMcd, cxr2);
			outputs.add(fxFLabelu);
			final EPMALabel fxFLabels = MatrixCorrectionModel2.FxFLabel(std1Mcd, cxr2);
			outputs.add(fxFLabels);

			final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, outputs);
			final Report r = new Report("XPP Report - test10");
			UncertainValuesBase<EPMALabel> resultsD = null;
			try {
				{
					r.addHeader("test1()");
					r.addHTML(xpp.toHTML(Mode.NORMAL));
					r.add(xpp.getOutputLabels(), Mode.VERBOSE);
					r.addHeader("Inputs");

					final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
					xpp.addConstraints(xpp.buildConstraints(inputs));
					xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

					assertTrue(UncertainValuesBase.testEquality(inputs, inputs.blockDiagnonalize()));
					r.add(inputs.blockDiagnonalize());

					final UncertainValuesCalculator<EPMALabel> xppI = new UncertainValuesCalculator<>(xpp, inputs);

					assertTrue(UncertainValuesBase.testEquality(xppI, xppI.blockDiagnonalize()));
					assertTrue(UncertainValuesBase.testEquality(xppI, xppI.sort()));
					final UncertainValuesBase<EPMALabel> aResults = xppI.sort();
					assertEquals(tagAu.toString(), aResults.getEntry(tagAu), 401.654, 0.001);
					assertEquals(aResults.getEntry(tagau), 11189.1680, 0.001);
					assertEquals(aResults.getEntry(tagBu), -526613.8449, 0.001);
					assertEquals(aResults.getEntry(tagbu), 12568.9576, 0.001);
					assertEquals(aResults.getEntry(tagPhi0u), 1.252, 0.001);

					assertEquals(aResults.getEntry(tagAs), 396.744, 0.001);
					assertEquals(aResults.getEntry(tagas), 11314.8474, 0.001);
					assertEquals(aResults.getEntry(tagBs), -529359.3293, 0.001);
					assertEquals(aResults.getEntry(tagbs), 12719.6941, 0.001);
					assertEquals(aResults.getEntry(tagPhi0s), 1.254, 0.001);

					assertEquals(aResults.getEntry(tagChiu), 2542.429, 0.001);
					assertEquals(aResults.getEntry(tagChis), 1038.418, 0.001);
					assertEquals(aResults.getEntry(tagFChiFu), 0.6324, 0.001);
					assertEquals(aResults.getEntry(tagFChiFs), 0.8198, 0.001);
					assertEquals(aResults.getEntry(zafLabel), 0.7805, 0.001);
					assertEquals(aResults.getEntry(tagZAO), 1.0731, 0.001);
					assertEquals(aResults.getEntry(tagZAAl), 0.796, 0.002);

					final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<>(xpp, inputs);

					assertEquals(aResults.getEntry(tagAu), uvc.getEntry(tagAu), 0.001);
					assertEquals(aResults.getEntry(tagau), uvc.getEntry(tagau), 0.001);
					assertEquals(aResults.getEntry(tagBu), uvc.getEntry(tagBu), 0.001);
					assertEquals(aResults.getEntry(tagbu), uvc.getEntry(tagbu), 0.001);
					assertEquals(aResults.getEntry(tagPhi0u), uvc.getEntry(tagPhi0u), 0.001);

					assertEquals(aResults.getEntry(tagAs), uvc.getEntry(tagAs), 0.001);
					assertEquals(aResults.getEntry(tagas), uvc.getEntry(tagas), 0.001);
					assertEquals(aResults.getEntry(tagBs), uvc.getEntry(tagBs), 0.001);
					assertEquals(aResults.getEntry(tagbs), uvc.getEntry(tagbs), 0.001);
					assertEquals(aResults.getEntry(tagPhi0s), uvc.getEntry(tagPhi0s), 0.001);

					assertEquals(aResults.getEntry(tagChiu), uvc.getEntry(tagChiu), 0.001);
					assertEquals(aResults.getEntry(tagChis), uvc.getEntry(tagChis), 0.001);
					assertEquals(aResults.getEntry(tagFChiFu), uvc.getEntry(tagFChiFu), 0.001);
					assertEquals(aResults.getEntry(tagFChiFs), uvc.getEntry(tagFChiFs), 0.001);
					assertEquals(aResults.getEntry(zafLabel), uvc.getEntry(zafLabel), 0.001);
					assertEquals(aResults.getEntry(tagZAO), uvc.getEntry(tagZAO), 0.001);
					assertEquals(aResults.getEntry(tagZAAl), uvc.getEntry(tagZAAl), 0.001);

					r.addHeader("Results");
					assertTrue(UncertainValuesBase.testEquality(aResults, aResults.blockDiagnonalize()));
					r.add(aResults.blockDiagnonalize());
					r.addHeader("Uncertain Values (relative to inputs)");
					final Map<EPMALabel, UncertainValueEx<EPMALabel>> outVals = new UncertainValuesCalculator<>(xpp,
							inputs).getOutputValues();
					final Table valTable = new Table();
					valTable.addRow(Table.td("Name"), Table.td("Value"), Table.td("Value (Normal)"),
							Table.td("Value (Verbose)"));
					final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
					for (final EPMALabel outTag : xpp.getOutputLabels()) {
						final UncertainValueEx<EPMALabel> uv = outVals.get(outTag);
						valTable.addRow(Table.td(HTML.toHTML(outTag, Mode.TERSE)),
								Table.td(aResults.getUncertainValue(outTag).toHTML(Mode.TERSE, bnf)),
								Table.td(uv.toHTML(Mode.TERSE, bnf)), Table.td(uv.toHTML(Mode.VERBOSE, bnf)));
					}
					r.addHTML(valTable.toHTML(Mode.NORMAL));

					r.addHeader("Covariance matrix");
					final StringBuffer sb = new StringBuffer();
					for (final EPMALabel tag : aResults.getLabels()) {
						if (sb.length() > 0)
							sb.append(",");
						sb.append(HTML.toHTML(tag, Mode.TERSE));
					}
					r.addHTML(HTML.p(sb.toString()));

					r.addImage(aResults.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

					final UncertainValuesCalculator<EPMALabel> aCalc = XPPMatrixCorrection2.buildAnalytical(skrl,
							unk.getValueMap(MassFraction.class), true);
					final UncertainValuesCalculator<EPMALabel> fdCalc = XPPMatrixCorrection2
							.buildFiniteDifference(skrl, unk.getValueMap(MassFraction.class), 0.001, true).force();

					assertTrue(aCalc.getJacobian().isPresent());
					final Jacobian<EPMALabel, EPMALabel> aJac = aCalc.getJacobian().get();
					final Jacobian<EPMALabel, EPMALabel> fdJac = fdCalc.getJacobian().get();

					for (final EPMALabel output : aJac.getOutputLabels())
						for (final EPMALabel input : aJac.getInputLabels())
							if (Math.abs(aJac.getEntry(input, output)) > 1.0e-8) {
								if (Math.abs(aJac.getEntry(input, output) - fdJac.getEntry(input, output)) > //
								0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
										Math.abs(fdJac.getEntry(input, output)))) {
									checkEquals(input, output, aJac.getEntry(input, output),
											fdJac.getEntry(input, output),
											0.01 * Math.max(Math.abs(aJac.getEntry(input, output)),
													Math.abs(fdJac.getEntry(input, output))));
								}
							}

					if (DUMP) {
						System.out.println("Results");
						System.out.println(aResults.toCSV());

						System.out.println("Jacobian");
						System.out.println(aCalc.toCSV());
						System.out.println("Jacobian(estimated)");
						System.out.println(fdCalc.toCSV());
					}

					// final Object unkCompTag = new MatrixCorrectionModel2.MaterialBasedLabel("J",
					// unkMat);
					// assertEquals(aCalc.getJacobianEntry(unkCompTag,
					// MaterialLabel.buildMassFractionTag(unkMat, Element.Oxygen)), -0.027565,
					// 0.00001);
					// assertEquals(aCalc.getJacobianEntry(unkCompTag,
					// MatrixCorrectionModel2.meanIonizationLabel(Element.Oxygen)), 0.609601,
					// 0.00001);
					resultsD = fdCalc;
					r.addImage(fdCalc.asCovarianceBitmap(8, V2L3, L2C), "Delta uncertainty matrix");
					r.addImage(UncertainValuesBase.compareAsBitmap(aCalc, fdCalc, L2C, 8), "Comparing uncertainty matrix");

				}
				if (MC_ITERATIONS > 0) {
					r.addHeader("Monte Carlo Results");
					r.add(xpp);
					r.addHeader("Inputs");
					final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
					xpp.addConstraints(xpp.buildConstraints(inputs));
					xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

					r.add(inputs);
					final UncertainValuesBase<EPMALabel> aResults = UncertainValuesBase.propagateAnalytical(xpp, inputs)
							.sort();

					assertEquals(aResults.getEntry(shellLabelA), 2366.373, 0.001);
					assertEquals(aResults.getEntry(shellLabela), 11052.8730, 0.001);
					assertEquals(aResults.getEntry(shellLabelB), -1460552.629, 0.001);
					assertEquals(aResults.getEntry(shellLabelb), 11681.219, 0.001);
					assertEquals(aResults.getEntry(phi0Label), 1.258, 0.001);

					assertEquals(aResults.getEntry(shellLabelA2), 2307.215, 0.001);
					assertEquals(aResults.getEntry(shellLabela2), 11178.191, 0.001);
					assertEquals(aResults.getEntry(shellLabelB2), -1459152.313, 0.001);
					assertEquals(aResults.getEntry(shellLabelb2), 11822.223, 0.001);
					assertEquals(aResults.getEntry(phi0Label2), 1.26, 0.001);

					assertEquals(aResults.getEntry(chiLabelu), 5836.018, 0.001);
					assertEquals(aResults.getEntry(chiLabels), 6414.025, 0.001);
					assertEquals(aResults.getEntry(fxFLabelu), 0.3554, 0.001);
					assertEquals(aResults.getEntry(fxFLabels), 0.3346, 0.001);
					assertEquals(aResults.getEntry(zafLabel), 0.7805, 0.001);

					r.addHeader("Analytic Results");
					r.add(aResults);

					final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<EPMALabel>(xpp,
							inputs);
					uvc.setCalculator(uvc.new MonteCarlo(MC_ITERATIONS));
					final UncertainValuesBase<EPMALabel> resultsMc = UncertainValues.asUncertainValues(uvc).sort();
					final EstimateUncertainValues<EPMALabel> mcp = (EstimateUncertainValues<EPMALabel>) uvc
							.getUncertainValues();
					if (DUMP) {
						System.out.println("Monte Carlo Results");
						System.out.println(resultsMc.toCSV());
					}

					r.addHeader("MC Results");
					r.add(resultsMc);

					final StringBuffer sb = new StringBuffer();
					for (final EPMALabel tag : resultsMc.getLabels()) {
						if (sb.length() > 0)
							sb.append(",");
						sb.append(HTML.toHTML(tag, Mode.TERSE));
					}
					r.addHTML(HTML.p(sb.toString()));
					// r.addImage(resultsMc.asCovarianceBitmap(8, V2L3, L2C), "Correlation matrix");

					r.addSubHeader("Phi0");
					final BasicNumberFormat bnf = new BasicNumberFormat("0.000E0");
					for (final EPMALabel tag : xpp.getOutputLabels()) {
						r.addSubHeader(HTML.toHTML(tag, Mode.NORMAL));
						r.add(MathUtilities.toHTML(mcp.getDescriptiveStatistics(tag), bnf));
					}

					{
						r.addHeader("Compare MC to Analytical");
						final Table t = new Table();
						t.addRow(Table.th("Tag"), //
								Table.th("V(MonteCarlo)"), //
								Table.th("U(Monte Carlo)"), //
								Table.th("V(Analytical)"), //
								Table.th("U(Analytic)"), //
								Table.th("V(Delta)"), //
								Table.th("U(Delta)"));
						final BasicNumberFormat bnf2 = new BasicNumberFormat("0.0000");
						for (final EPMALabel tag : xpp.getOutputLabels())
							if (tag instanceof ZAFMultiLineLabel) {
								t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
										MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
										MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
										MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf2),
										MathUtilities.td(aResults.getUncertainty(tag), bnf2),
										MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
										MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
							}
						for (final EPMALabel tag : xpp.getOutputLabels())
							if (tag instanceof KRatioLabel) {
								t.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), //
										MathUtilities.td(resultsMc.getValue(tag).doubleValue(), bnf2), //
										MathUtilities.td(resultsMc.getUncertainty(tag), bnf2), //
										MathUtilities.td(aResults.getValue(tag).doubleValue(), bnf2),
										MathUtilities.td(aResults.getUncertainty(tag), bnf2),
										MathUtilities.td(resultsD.getValue(tag).doubleValue(), bnf2),
										MathUtilities.td(resultsD.getUncertainty(tag), bnf2));
							}
						r.add(t);
					}
					r.addHeader("Done!");
				}
			} catch (final Throwable e) {
				r.addThrowable(e);
				r.inBrowser(Mode.VERBOSE);
				throw e;
			}
			r.inBrowser(Mode.VERBOSE);
		} finally {
			ExplicitMeasurementModel.sDump = null;
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
	public void testXPP11() throws ArgumentException, ParseException, IOException {
		// K412 as in SP 260-74 using elements and simple compounds

		final Composition std0 = Composition.parse("SiO2");
		final Composition std1 = Composition.parse("Al2O3");
		final Composition std2 = Composition.parse("MgO");
		final Composition std3 = Composition.parse("CaF2");
		final Composition std4 = Composition.parse("Fe");

		final Composition unk = buildK412(false);
		final UncertainValue toa = UncertainValue.toRadians(40.0, 0.5);
		final Layer coating = Layer.carbonCoating(new UncertainValue(10.0, 2.0));
		final double roughness = MatrixCorrectionDatum.roughness(1.0e-8, 3.0);

		final Report r = new Report("XPP Report - Test11()");
		final int MIN_E = 9;
		final int MAX_E = 31;
		final Map<Integer, Map<EPMALabel, UncertainValueEx<EPMALabel>>> outVals = new TreeMap<>();
		List<? extends EPMALabel> outLabels = null, inLabels = null;

		try {
			final Table valTable = new Table();
			for (int i = MIN_E; i < MAX_E; ++i) {
				final double e0 = i;
				final UncertainValue e0u = new UncertainValue(e0, 0.1);
				final StandardMatrixCorrectionDatum std0Mcd = //
						new StandardMatrixCorrectionDatum(std0, e0u, toa, roughness, coating);
				final StandardMatrixCorrectionDatum std1Mcd = //
						new StandardMatrixCorrectionDatum(std1, e0u, toa, roughness, coating);
				final StandardMatrixCorrectionDatum std2Mcd = //
						new StandardMatrixCorrectionDatum(std2, e0u, toa, roughness, coating);
				final StandardMatrixCorrectionDatum std3Mcd = //
						new StandardMatrixCorrectionDatum(std3, e0u, toa, roughness, coating);
				final StandardMatrixCorrectionDatum std4Mcd = //
						new StandardMatrixCorrectionDatum(std4, e0u, toa, roughness, coating);
				final UnknownMatrixCorrectionDatum unkMcd = //
						new UnknownMatrixCorrectionDatum(unk.getMaterial(), e0u, toa, roughness, coating);

				final Set<KRatioLabel> skrl = new HashSet<>();
				skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Silicon, XRayTransition.KA1),
						Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std0Mcd, ElementXRaySet.singleton(Element.Oxygen, XRayTransition.KA1),
						Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std1Mcd,
						ElementXRaySet.singleton(Element.Aluminum, XRayTransition.KA1), Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std2Mcd,
						ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA1), Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std2Mcd,
						ElementXRaySet.singleton(Element.Magnesium, XRayTransition.KA2), Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std3Mcd, ElementXRaySet.singleton(Element.Calcium, XRayTransition.KA1),
						Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std3Mcd,
						ElementXRaySet.singleton(Element.Calcium, XRayTransition.L3M1), Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.KA1),
						Method.Measured));
				skrl.add(new KRatioLabel(unkMcd, std4Mcd, ElementXRaySet.singleton(Element.Iron, XRayTransition.LA1),
						Method.Measured));

				final Set<EPMALabel> outputs = new HashSet<>();
				for (final KRatioLabel krl : skrl) {
					outputs.add(MatrixCorrectionModel2.zafLabel(krl));
					final StandardMatrixCorrectionDatum meStd = krl.getStandard();
					for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
						outputs.add(MatrixCorrectionModel2.zafLabel(unkMcd, meStd, cxr));
						outputs.add(MatrixCorrectionModel2.FxFLabel(unkMcd, cxr));
						outputs.add(MatrixCorrectionModel2.FxFLabel(meStd, cxr));
					}
				}
				final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(skrl, new ArrayList<>(outputs));
				final UncertainValues<EPMALabel> inputs = xpp.buildInput(unk.getMaterial());
				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(unk.getValueMap(MassFraction.class));

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<>(xpp, inputs);

				outVals.put(i, uvc.getOutputValues(0.0));
				if (i == MIN_E) {
					outLabels = uvc.getOutputLabels();
					inLabels = uvc.getInputLabels();
				}
			}
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
				if (outTag instanceof ZAFMultiLineLabel) {
					final ZAFMultiLineLabel mcl = (ZAFMultiLineLabel) outTag;
					for (int ii = 0; ii < inLabels.size(); ++ii) {
						final EPMALabel inTag = inLabels.get(ii);
						final List<Item> row = new ArrayList<>();
						row.add(Table.td(mcl));
						row.add(Table.td(mcl.getElementXRaySet()));
						row.add(Table.td(inTag));
						boolean addRow = false;
						for (int i = MIN_E; i < MAX_E; ++i) {
							final Map<EPMALabel, UncertainValueEx<EPMALabel>> oVals = outVals.get(i);
							final UncertainValueEx<EPMALabel> tmp = getByXRT(oVals, mcl);
							if (tmp != null)
								row.add(Table.td(bnf.format(tmp.doubleValue())));
							else
								row.add(Table.td("---"));
						}
						for (int i = MIN_E; i < MAX_E; ++i) {
							final Map<EPMALabel, UncertainValueEx<EPMALabel>> oVals = outVals.get(i);
							final UncertainValueEx<EPMALabel> tmp = getByXRT(oVals, mcl);
							if (tmp != null) {
								final double cbn = getComponentByName(tmp, inTag);
								if (cbn != 0.0)
									addRow = true;
								row.add(Table.td(bnf.format(cbn)));
							} else
								row.add(Table.td("---"));
						}
						for (int i = MIN_E; i < MAX_E; ++i) {
							final Map<EPMALabel, UncertainValueEx<EPMALabel>> oVals = outVals.get(i);
							final UncertainValueEx<EPMALabel> tmp = getByXRT(oVals, mcl);
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
			r.addHTML(valTable.toHTML(Mode.NORMAL));
		} catch (final Throwable e) {
			r.addThrowable(e);
			r.inBrowser(Mode.VERBOSE);
			throw e;
		}
		r.inBrowser(Mode.VERBOSE);

	}

	public double getComponentByName(
			final UncertainValueEx<EPMALabel> tmp, final EPMALabel tag
		) {
		final String name = tag.toString();
		for (final Map.Entry<EPMALabel, Double> me : tmp.getComponents().entrySet())
			if (me.getKey().toString().equals(name))
				return me.getValue().doubleValue();
		return 0.0;
	}

	public UncertainValueEx<EPMALabel> getByXRT(
			final Map<EPMALabel, UncertainValueEx<EPMALabel>> oVals, final ZAFMultiLineLabel mcl
		) {
		for (final Map.Entry<EPMALabel, UncertainValueEx<EPMALabel>> me : oVals.entrySet()) {
			if (me.getKey() instanceof ZAFMultiLineLabel) {
				final ZAFMultiLineLabel mcl2 = (ZAFMultiLineLabel) me.getKey();
				if (mcl2.toString().equals(mcl.toString()))
					return me.getValue();
			}
		}
		return null;
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
