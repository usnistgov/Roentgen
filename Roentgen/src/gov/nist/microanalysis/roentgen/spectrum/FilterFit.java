package gov.nist.microanalysis.roentgen.spectrum;

import java.util.Set;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
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

	private final SimplyLazy<EDSFittingFilter> mFilter = new SimplyLazy<EDSFittingFilter>() {

		@Override
		protected EDSFittingFilter initialize() {
			try {
				switch (mMode) {
				case AdaptiveGaussian:
					return new AdaptiveGaussianFilter(mChannelCount, mEnergy, mLineshape);
				case AdaptiveTophat:
				default:
					return new AdaptiveTophatFilter(mChannelCount, mEnergy, mLineshape);
				}
			} catch (final ArgumentException e) {
				e.printStackTrace();
			}
			return null;
		}

	};

	/**
	 * Constructs a FilterFit
	 */
	public FilterFit(
			final int nChannels, final EnergyCalibration ec, final LineshapeCalibration ls
	) {
		mChannelCount = nChannels;
		mEnergy = ec;
		mLineshape = ls;
	}

	public FilterFit addReference(
			final ElementXRaySet xrays, final Set<Element> elms, final EDSSpectrum spec
	) {
		mFilter.get();
		return this;
	}

}
