package gov.nist.microanalysis.roentgen.spectrum;

import java.io.IOException;
import java.io.InputStream;

/**
 * <p>
 * This interface defines a mechanism for reading spectrum files from an
 * InputStream. Implement this interface to read spectrum files of various
 * formats.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
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
    */
   public EDSSpectrum read(final InputStream f, final int i)
         throws IOException;

   /**
    * Returns the number of spectra contained within f.
    *
    * @param f
    * @return 0 for none
    */
   public int getSpectrumCount(InputStream is);

}