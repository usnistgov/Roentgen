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

import com.duckandcover.html.HTMLList;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.wds.ModelLabels;
import gov.nist.microanalysis.roentgen.wds.TwoPointContinuumModel;

public class TwoPointContinuumModelTest {

	private static double sqr(final double x) {
		return x * x;
	}

	@Test
	public void test1() throws ArgumentException, ParseException, IOException {
		final Composition comp = Composition.parse("SiO2");

		final MatrixCorrectionDatum mcd = new StandardMatrixCorrectionDatum(comp, new UncertainValue(20.0, 0.1),
				new UncertainValue(40.0, 0.2));
		final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1);
		final RealVector vals = new ArrayRealVector(12);
		final RealVector vars = new ArrayRealVector(12);
		final List<ModelLabels<?,?>> labels = Arrays.asList(new ModelLabels[12]);
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
		final UncertainValues<ModelLabels<?,?>> uvs = new UncertainValues<>(labels, vals, vars);

		final Report r = new Report("2Point");
		try {
			r.addHeader("Inputs");
			final TwoPointContinuumModel tpcm = new TwoPointContinuumModel(mcd, cxr);
			{
				final HTMLList l = new HTMLList();
				l.addAll(tpcm.getInputLabels());
				r.addHTML(l.toHTML(Mode.NORMAL));
			}
			r.add(uvs);
			r.addHeader("Outputs");
			{
				final HTMLList l = new HTMLList();
				l.addAll(tpcm.getOutputLabels());
				r.addHTML(l.toHTML(Mode.NORMAL));
			}
			final UncertainValuesBase<ModelLabels<?,?>> res = UncertainValues.propagate(tpcm, uvs);
			r.add(res);

			final HashMap<? extends Object, UncertainValue> ovals = tpcm.getOutputValues(uvs, 0.0);
			final Table tt = new Table();
			for (final Map.Entry<? extends Object, UncertainValue> me : ovals.entrySet()) {
				tt.addRow(Table.td(me.getKey()), Table.td((IToHTML) me.getValue()));
			}
			r.addHTML(tt.toHTML(Mode.VERBOSE));
		} finally {
			r.inBrowser(Mode.NORMAL);
		}

	}

}
