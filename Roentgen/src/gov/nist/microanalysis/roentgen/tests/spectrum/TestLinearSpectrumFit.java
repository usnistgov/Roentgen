package gov.nist.microanalysis.roentgen.tests.spectrum;

import java.util.HashSet;
import java.util.Set;

import org.junit.Before;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.spectrum.EDSSpectrum;
import gov.nist.microanalysis.roentgen.spectrum.FilterFit;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestLinearSpectrumFit {

	private EDSSpectrum mUnknown;
	private EDSSpectrum mSiO2Std;
	private EDSSpectrum mMgStd;
	private EDSSpectrum mAlStd;
	private EDSSpectrum mCaF2Std;
	private EDSSpectrum mFeStd;

	@Before
	public void setup() throws Exception {
		mUnknown = TestReadSpectrum.fromResource("K412 15 keV.msa");
		mSiO2Std = TestReadSpectrum.fromResource("SiO2 15 keV.msa");
		mMgStd = TestReadSpectrum.fromResource("Mg 15 keV.msa");
		mAlStd = TestReadSpectrum.fromResource("Al 15 keV.msa");
		mCaF2Std = TestReadSpectrum.fromResource("CaF2 15 keV.msa");
		mFeStd = TestReadSpectrum.fromResource("Fe 15 keV.msa");
	}

	/**
	 * Constructs a TestLinearSpectrumFit
	 */
	// @Test
	public TestLinearSpectrumFit() {
		final double e0 = 15.0e3;
		final FilterFit lsf = new FilterFit(mUnknown.size(), mUnknown.getEnergyCalibration(),
				mUnknown.getLineshapeCalibration());
		{
			final Set<Element> elms = new HashSet<>();
			elms.add(Element.Silicon);
			elms.add(Element.Oxygen);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Silicon, e0), elms, mSiO2Std);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Oxygen, e0), elms, mSiO2Std);
		}
		{
			final Set<Element> elms = new HashSet<>();
			elms.add(Element.Magnesium);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Magnesium, e0), elms, mMgStd);
		}
		{
			final Set<Element> elms = new HashSet<>();
			elms.add(Element.Aluminum);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Aluminum, e0), elms, mAlStd);
		}
		{
			final Set<Element> elms = new HashSet<>();
			elms.add(Element.Calcium);
			elms.add(Element.Fluorine);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Calcium, e0), elms, mCaF2Std);
		}
		{
			final Set<Element> elms = new HashSet<>();
			elms.add(Element.Iron);
			lsf.addReference(ElementXRaySet.buildElementXRaySet(Element.Iron, e0), elms, mFeStd);
		}
		// lsf.compute(mUnknown);
	}

}
