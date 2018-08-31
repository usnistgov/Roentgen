package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.io.IOException;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import com.duckandcover.html.Report;
import com.duckandcover.html.IToHTML.Mode;

import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.math.MathUtilities;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.EDSMatrixCorrection;
import gov.nist.microanalysis.roentgen.matrixcorrection.EDSMatrixCorrection.XRayWeightTag;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;

public class EDSMatrixCorrectionTest {

	private class Datum {

		private final CharacteristicXRay mXRay;
		private final Composition mStd;
		private final double mZaf;
		private final double mWeight;

		Datum(String line, String std, double zaf, double weight) {
			mXRay = CharacteristicXRay.parse(line);
			Composition tmp = null;
			try {
				tmp = Composition.parse(std);
			} catch (ParseException e) {
				e.printStackTrace();
			}
			mStd = tmp;
			mZaf = zaf;
			mWeight = weight;
		}
	}

	private final Datum[] test1Data = { //
			new Datum("O K-L3", "SiO2", 1.0628, 0.666666666666667),
			new Datum("O K-L2", "SiO2", 1.0628, 0.333333333333333),
			new Datum("Si K-L3", "SiO2", 0.8817, 0.644163875289874),
			new Datum("Si K-L2", "SiO2", 0.8817, 0.325302757021386),
			new Datum("Si K-M3", "SiO2", 0.8936, 0.02035557845916),
			new Datum("Si K-M2", "SiO2", 0.8942, 0.01017778922958) };

	private final Datum[] test2Data = { //
			new Datum("Na K-L3", "Na", 0.5494, 0.666666666666667),
			new Datum("Na K-L2", "Na", 0.5494, 0.333333333333333) };
	
	private final Datum[] test3Data = { //
			new Datum("Al K-L3", "Al", 0.7336, 0.652230628750326),
			new Datum("Al K-L2", "Al", 0.7336, 0.330680928776415),
			new Datum("Al K-M3", "Al", 0.7532, 0.017088442473259), };

	@SuppressWarnings("unchecked")
	@Test
	public void test1() throws ParseException, IOException {
		final Datum[] data = test1Data;
		MatrixCorrectionDatum unk = new MatrixCorrectionDatum( //
				Composition.parse("NaAlSi3O8"), false, new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));
		MatrixCorrectionDatum std = new MatrixCorrectionDatum( //
				data[0].mStd, //
				true, // 
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));

		List<Object> inputs = new ArrayList<Object>();
		RealVector point = new ArrayRealVector(2 * data.length);
		RealVector vars = new ArrayRealVector(2 * data.length);
		Map<Element, ElementXRaySet> mxrs = new HashMap<>();
		for (int i = 0; i < data.length; ++i) {
			Datum d = data[i];
			inputs.add(new XRayWeightTag(d.mXRay));
			point.setEntry(2 * i, d.mWeight);
			vars.setEntry(2 * i, Math.pow(d.mWeight * 0.05, 2.0));
			inputs.add(XPPMatrixCorrection.zaTag(unk, std, d.mXRay));
			point.setEntry(2 * i + 1, d.mZaf);
			vars.setEntry(2 * i + 1, Math.pow(d.mZaf * 0.01, 2.0));
			ElementXRaySet t = mxrs.get(d.mXRay.getElement());
			if (t != null)
				t.add(d.mXRay);
			else
				mxrs.put(d.mXRay.getElement(), new ElementXRaySet(d.mXRay));
		}
		List<ElementXRaySet> lexrs = new ArrayList<>(mxrs.values());
		EDSMatrixCorrection emc = new EDSMatrixCorrection(unk, std, lexrs);
		Pair<RealVector, RealMatrix> res = emc.evaluate(point);

		Report r = new Report("EDSMatrixCorrection - 1");
		r.addHeader("Values");
		r.addHTML(MathUtilities.toHTML_Vertical(res.getFirst(), (List<Object>) emc.getOutputTags(),
				NumberFormat.getInstance()));
		r.addHeader("Jacobian");
		r.addHTML(MathUtilities.toHTML(res.getSecond(), (List<Object>) emc.getOutputTags(), inputs,
				NumberFormat.getInstance()));
		UncertainValues inp = new UncertainValues(emc.getInputTags(), point, vars);
		r.addHeader("Inputs");
		r.add(inp);
		UncertainValues output = UncertainValues.propagate(emc, inp);
		r.addHeader("Results");
		r.add(output);
		r.inBrowser(Mode.NORMAL);
	}

	@SuppressWarnings("unchecked")
	@Test
	public void test2() throws ParseException, IOException {
		final Datum[] data = test2Data;
		MatrixCorrectionDatum unk = new MatrixCorrectionDatum( //
				Composition.parse("NaAlSi3O8"), false, new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));
		MatrixCorrectionDatum std = new MatrixCorrectionDatum( //
				data[0].mStd, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));

		List<Object> inputs = new ArrayList<Object>();
		RealVector point = new ArrayRealVector(2 * data.length);
		RealVector vars = new ArrayRealVector(2 * data.length);
		Map<Element, ElementXRaySet> mxrs = new HashMap<>();
		for (int i = 0; i < data.length; ++i) {
			Datum d = data[i];
			inputs.add(new XRayWeightTag(d.mXRay));
			point.setEntry(2 * i, d.mWeight);
			vars.setEntry(2 * i, Math.pow(d.mWeight * 0.05, 2.0));
			inputs.add(XPPMatrixCorrection.zaTag(unk, std, d.mXRay));
			point.setEntry(2 * i + 1, d.mZaf);
			vars.setEntry(2 * i + 1, Math.pow(d.mZaf * 0.01, 2.0));
			ElementXRaySet t = mxrs.get(d.mXRay.getElement());
			if (t != null)
				t.add(d.mXRay);
			else
				mxrs.put(d.mXRay.getElement(), new ElementXRaySet(d.mXRay));
		}
		List<ElementXRaySet> lexrs = new ArrayList<>(mxrs.values());
		EDSMatrixCorrection emc = new EDSMatrixCorrection(unk, std, lexrs);
		Pair<RealVector, RealMatrix> res = emc.evaluate(point);

		Report r = new Report("EDSMatrixCorrection - 2");
		r.addHeader("Values");
		r.addHTML(MathUtilities.toHTML_Vertical(res.getFirst(), (List<Object>) emc.getOutputTags(),
				NumberFormat.getInstance()));
		r.addHeader("Jacobian");
		r.addHTML(MathUtilities.toHTML(res.getSecond(), (List<Object>) emc.getOutputTags(), inputs,
				NumberFormat.getInstance()));
		UncertainValues inp = new UncertainValues(emc.getInputTags(), point, vars);
		r.addHeader("Inputs");
		r.add(inp);
		UncertainValues output = UncertainValues.propagate(emc, inp);
		r.addHeader("Results");
		r.add(output);
		r.inBrowser(Mode.NORMAL);
	}

	@SuppressWarnings("unchecked")
	@Test
	public void test3() throws ParseException, IOException {
		final Datum[] data = test3Data;
		MatrixCorrectionDatum unk = new MatrixCorrectionDatum( //
				Composition.parse("NaAlSi3O8"), false, new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));
		MatrixCorrectionDatum std = new MatrixCorrectionDatum( //
				data[0].mStd, true, //
				new UncertainValue(15.0, 0.1), //
				new UncertainValue(Math.toRadians(40.0), Math.toRadians(0.7)));

		List<Object> inputs = new ArrayList<Object>();
		RealVector point = new ArrayRealVector(2 * data.length);
		RealVector vars = new ArrayRealVector(2 * data.length);
		Map<Element, ElementXRaySet> mxrs = new HashMap<>();
		for (int i = 0; i < data.length; ++i) {
			Datum d = data[i];
			inputs.add(new XRayWeightTag(d.mXRay));
			point.setEntry(2 * i, d.mWeight);
			vars.setEntry(2 * i, Math.pow(d.mWeight * 0.05, 2.0));
			inputs.add(XPPMatrixCorrection.zaTag(unk, std, d.mXRay));
			point.setEntry(2 * i + 1, d.mZaf);
			vars.setEntry(2 * i + 1, Math.pow(d.mZaf * 0.01, 2.0));
			ElementXRaySet t = mxrs.get(d.mXRay.getElement());
			if (t != null)
				t.add(d.mXRay);
			else
				mxrs.put(d.mXRay.getElement(), new ElementXRaySet(d.mXRay));
		}
		List<ElementXRaySet> lexrs = new ArrayList<>(mxrs.values());
		EDSMatrixCorrection emc = new EDSMatrixCorrection(unk, std, lexrs);
		Pair<RealVector, RealMatrix> res = emc.evaluate(point);

		Report r = new Report("EDSMatrixCorrection - 3");
		r.addHeader("Values");
		r.addHTML(MathUtilities.toHTML_Vertical(res.getFirst(), (List<Object>) emc.getOutputTags(),
				NumberFormat.getInstance()));
		r.addHeader("Jacobian");
		r.addHTML(MathUtilities.toHTML(res.getSecond(), (List<Object>) emc.getOutputTags(), inputs,
				NumberFormat.getInstance()));
		UncertainValues inp = new UncertainValues(emc.getInputTags(), point, vars);
		r.addHeader("Inputs");
		r.add(inp);
		UncertainValues output = UncertainValues.propagate(emc, inp);
		r.addHeader("Results");
		r.add(output);
		r.inBrowser(Mode.NORMAL);
	}

	
}
