package gov.nist.microanalysis.roentgen.spectrum;

import org.apache.commons.math3.linear.RealVector;

import gov.nist.microanalysis.roentgen.math.IntInterval;

/**
 * <p>
 * A calibrated reference is a filtered set of peak shape references spectrum
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class CalibratedReference {

	private final IntInterval mChannels;

	private final RealVector mFiltered;

	/**
	 * Constructs a CalibratedReference
	 */
	public CalibratedReference(final RealVector filtered, final IntInterval channels) {
		mChannels = channels;
		mFiltered = filtered;
	}

}
