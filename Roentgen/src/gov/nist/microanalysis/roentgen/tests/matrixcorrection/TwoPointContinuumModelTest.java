package gov.nist.microanalysis.roentgen.tests.matrixcorrection;

import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.HTMLList;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTML.Mode;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.ModelLabels;
import gov.nist.microanalysis.roentgen.matrixcorrection.TwoPointContinuumModel;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

public class TwoPointContinuumModelTest {
	
	private static double sqr(double x) {
		return x*x;
	}

	@Test
	public void test1() throws ArgumentException, ParseException, IOException {
		Composition comp = Composition.parse("SiO2");

		MatrixCorrectionDatum mcd = new MatrixCorrectionDatum(comp, true, new UncertainValue(20.0, 0.1),
				new UncertainValue(40.0, 0.2));
		CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
		RealVector vals = new ArrayRealVector(12);
		RealVector vars = new ArrayRealVector(12);
		List<Object> labels = Arrays.asList(new Object[12]);
		// LOW
		labels.set(0, ModelLabels.buildSpectrometerPosition(mcd, cxr, TwoPointContinuumModel.LOW_BACK));
		vals.setEntry(0, 69.96629);
		vars.setEntry(0, sqr(0.01));
		labels.set(1, ModelLabels.buildRawIntensity(mcd, cxr, TwoPointContinuumModel.LOW_BACK));
		vals.setEntry(1, 242);
		vars.setEntry(1, 242);
		labels.set(2, ModelLabels.buildProbeCurrent(mcd, cxr, TwoPointContinuumModel.LOW_BACK));
		vals.setEntry(2, 10.0);
		vars.setEntry(2, sqr(0.01));
		labels.set(3, ModelLabels.buildLiveTime(mcd, cxr, TwoPointContinuumModel.LOW_BACK));
		vals.setEntry(3, 5.0);
		vars.setEntry(3, sqr(0.01));
		// PEAK
		labels.set(4, ModelLabels.buildSpectrometerPosition(mcd, cxr, TwoPointContinuumModel.ON_PEAK));
		vals.setEntry(4, 77.40314);
		vars.setEntry(4, sqr(0.01));
		labels.set(5, ModelLabels.buildRawIntensity(mcd, cxr, TwoPointContinuumModel.ON_PEAK));
		vals.setEntry(5, 57013);
		vars.setEntry(5, 57013);
		labels.set(6, ModelLabels.buildProbeCurrent(mcd, cxr, TwoPointContinuumModel.ON_PEAK));
		vals.setEntry(6, 10.0);
		vars.setEntry(6, sqr(0.01));
		labels.set(7, ModelLabels.buildLiveTime(mcd, cxr, TwoPointContinuumModel.ON_PEAK));
		vals.setEntry(7, 20.0);
		vars.setEntry(7, sqr(0.01));
		// HIGH
		labels.set(8, ModelLabels.buildSpectrometerPosition(mcd, cxr, TwoPointContinuumModel.HIGH_BACK));
		vals.setEntry(8, 84.64756);
		vars.setEntry(8, sqr(0.01));
		labels.set(9, ModelLabels.buildRawIntensity(mcd, cxr, TwoPointContinuumModel.HIGH_BACK));
		vals.setEntry(9, 148);
		vars.setEntry(9, 148);
		labels.set(10, ModelLabels.buildProbeCurrent(mcd, cxr, TwoPointContinuumModel.HIGH_BACK));
		vals.setEntry(10, 10.0);
		vars.setEntry(10, sqr(0.01));
		labels.set(11, ModelLabels.buildLiveTime(mcd, cxr, TwoPointContinuumModel.HIGH_BACK));
		vals.setEntry(11, 5.0);
		vars.setEntry(11, sqr(0.01));
		UncertainValues uvs = new UncertainValues(labels, vals, vars);

		Report r = new Report("2Point");
		try {
			r.addHeader("Inputs");
			TwoPointContinuumModel tpcm = new TwoPointContinuumModel(mcd, cxr);
			{
				HTMLList l = new HTMLList();
				l.addAll(tpcm.getInputLabels());
				r.addHTML(l.toHTML(Mode.NORMAL));
			}
			r.add(uvs);
			r.addHeader("Outputs");
			{
				HTMLList l = new HTMLList();
				l.addAll(tpcm.getOutputLabels());
				r.addHTML(l.toHTML(Mode.NORMAL));
			}
			UncertainValues res = UncertainValues.propagate(tpcm, uvs);
			r.add(res);
			
			HashMap<? extends Object, UncertainValue> ovals=tpcm.getOutputValues(uvs, 0.0);
			Table tt = new Table();
			for(Map.Entry<? extends Object, UncertainValue> me : ovals.entrySet()) {
				tt.addRow(Table.td(me.getKey()), Table.td((IToHTML)me.getValue()));
			}
			r.addHTML(tt.toHTML(Mode.VERBOSE));
		} finally {
			r.inBrowser(Mode.NORMAL);
		}

	}

}
