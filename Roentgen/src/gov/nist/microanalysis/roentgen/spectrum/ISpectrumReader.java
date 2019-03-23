package gov.nist.microanalysis.roentgen.spectrum;

import java.io.IOException;
import java.io.InputStream;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>
 * This interface defines a mechanism for reading spectrum files from an
 * InputStream. Implement this interface to read spectrum files of various
 * formats.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 312 $
 */
public interface ISpectrumReader {

	/**
	 * Reads the i-th spectrum in the File f
	 *
	 * @param det
	 * @param f
	 * @param i
	 * @return Spectrum
	 * @throws IOException
	 * @throws ArgumentException
	 */
	public EDSSpectrum read(
			final InputStream f, final int i
	) throws IOException, ArgumentException;

	/**
	 * Returns the number of spectra contained within f.
	 *
	 * @param f
	 * @return 0 for none
	 */
	public int getSpectrumCount(
			InputStream is
	);

}