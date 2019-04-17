package gov.nist.microanalysis.roentgen.tests.wds;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import com.duckandcover.html.Report;
import com.duckandcover.html.IToHTML.Mode;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.wds.ModelLabels;
import gov.nist.microanalysis.roentgen.wds.TwoPointContinuumModel;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class TwoPointContinuumModelTest {

	/**
	 * @throws Exception 
	 * 
	 */
	@Test
	public void test1() throws Exception {
		
		Material unknown = new Material("Unknown", Arrays.asList( Element.Sodium, Element.Aluminum, Element.Oxygen, Element.Tellurium ));
		
		Composition albite = Composition.pureElement(Element.Tellurium);
		
		UnknownMatrixCorrectionDatum unk = new UnknownMatrixCorrectionDatum(unknown, UncertainValue.valueOf(15.0,0.1), UncertainValue.toRadians(40.0, 1.0));
		StandardMatrixCorrectionDatum std = new StandardMatrixCorrectionDatum(albite, UncertainValue.valueOf(15.0,0.1), UncertainValue.toRadians(40.0, 1.0));		
		
		final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Tellurium, XRayTransition.L3M5);
		TwoPointContinuumModel tpcm = new TwoPointContinuumModel(unk, std, cxr);
		
		Map<EPMALabel, UncertainValue> map = new HashMap<>();
		
		map.put(ModelLabels.buildLiveTime(std, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.valueOf(2.0, 0.01));
		map.put(ModelLabels.buildLiveTime(std, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.valueOf(10.0, 0.01));
		map.put(ModelLabels.buildLiveTime(std, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.valueOf(2.0, 0.01));
		map.put(ModelLabels.buildLiveTime(unk, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.valueOf(4.0, 0.01));
		map.put(ModelLabels.buildLiveTime(unk, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.valueOf(60.0, 0.01));
		map.put(ModelLabels.buildLiveTime(unk, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.valueOf(4.0, 0.01));
		
		map.put(ModelLabels.buildSpectrometerPosition(std, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.valueOf(102.0, 0.01));
		map.put(ModelLabels.buildSpectrometerPosition(std, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.valueOf(105.35, 0.01));
		map.put(ModelLabels.buildSpectrometerPosition(std, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.valueOf(108.0, 0.01));
		map.put(ModelLabels.buildSpectrometerPosition(unk, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.valueOf(102.0, 0.01));
		map.put(ModelLabels.buildSpectrometerPosition(unk, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.valueOf(105.35, 0.01));
		map.put(ModelLabels.buildSpectrometerPosition(unk, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.valueOf(108.0, 0.01));

		map.put(ModelLabels.buildProbeCurrent(std,cxr, TwoPointContinuumModel.LOW_BACK),UncertainValue.valueOf(150.0, 1.0));
		map.put(ModelLabels.buildProbeCurrent(std,cxr, TwoPointContinuumModel.ON_PEAK),UncertainValue.valueOf(151.0, 1.0));
		map.put(ModelLabels.buildProbeCurrent(std,cxr, TwoPointContinuumModel.HIGH_BACK),UncertainValue.valueOf(149.0, 1.0));
		map.put(ModelLabels.buildProbeCurrent(unk,cxr, TwoPointContinuumModel.LOW_BACK),UncertainValue.valueOf(350.0, 1.0));
		map.put(ModelLabels.buildProbeCurrent(unk,cxr, TwoPointContinuumModel.ON_PEAK),UncertainValue.valueOf(351.0, 1.0));
		map.put(ModelLabels.buildProbeCurrent(unk,cxr, TwoPointContinuumModel.HIGH_BACK),UncertainValue.valueOf(349.0, 1.0));
		
		
		map.put(ModelLabels.buildRawIntensity(std, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.normal(2.0*1.2e3));
		map.put(ModelLabels.buildRawIntensity(std, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.normal(10.0*17.5e3));
		map.put(ModelLabels.buildRawIntensity(std, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.normal(2.0*0.82e3));
		map.put(ModelLabels.buildRawIntensity(unk, cxr, TwoPointContinuumModel.LOW_BACK), UncertainValue.normal(4.0*0.8e3));
		map.put(ModelLabels.buildRawIntensity(unk, cxr, TwoPointContinuumModel.ON_PEAK), UncertainValue.normal(60.0*1.2e3));
		map.put(ModelLabels.buildRawIntensity(unk, cxr, TwoPointContinuumModel.HIGH_BACK), UncertainValue.normal(4.0*0.65e3));

		UncertainValues<EPMALabel> inputs = new UncertainValues<>(map);
		
		UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<>(tpcm, inputs);
		
		UncertainValuesCalculator<EPMALabel> uvcMc = new UncertainValuesCalculator<>(tpcm, inputs);
		uvcMc.setCalculator(uvcMc.new MonteCarlo(100000, true));
		
		
		Report r = new Report("WDS Na K&alpha;1");
		BasicNumberFormat bnf = new BasicNumberFormat();
		r.addHeader("WDS Na K&alpha;1");
		r.addSubHeader("Inputs");
		r.add(inputs, Mode.NORMAL);
		r.addHeader("Outputs");
		r.addSubHeader("Analytical");
		r.addHTML(uvc.toSimpleHTML(bnf));
		r.addSubHeader("Monte Carlo");
		r.addHTML(uvcMc.toSimpleHTML(bnf));
		r.inBrowser(Mode.NORMAL);
	
		assertTrue(UncertainValuesBase.testSimiliarity(uvc, uvcMc, 1.0e-2));
	}

}
