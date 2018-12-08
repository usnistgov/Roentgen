package gov.nist.microanalysis.roentgen.spectrum;

import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealVector;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class FilterFit {

	public enum Mode {
		AdaptiveGaussian, AdaptiveTophat
	};

	private final Mode mMode = Mode.AdaptiveTophat;
	private final int mChannelCount;
	private final EnergyCalibration mEnergy;
	private final LineshapeCalibration mLineshape;

	private class ElementData {
		private Element mElement;
		private Map<Element, ElementXRaySet> mXRays;
		private EDSSpectrum mSpectrum;
		private RealVector mFiltered;

	};

	private final SimplyLazy<EDSFittingFilter> mFilter = new SimplyLazy<EDSFittingFilter>() {

		@Override
		protected EDSFittingFilter initialize() {
			switch (mMode) {
			case AdaptiveGaussian:
				return new AdaptiveGaussianFilter(mChannelCount, mEnergy, mLineshape);
			case AdaptiveTophat:
			default:
				return new AdaptiveTophatFilter(mChannelCount, mEnergy, mLineshape);
			}
		}

	};

	/**
	 * Constructs a FilterFit
	 */
	public FilterFit(final int nChannels, final EnergyCalibration ec, final LineshapeCalibration ls) {
		mChannelCount = nChannels;
		mEnergy = ec;
		mLineshape = ls;
	}

	public FilterFit addReference(final ElementXRaySet xrays, final Set<Element> elms, final EDSSpectrum spec) {

		return this;
	}

}
