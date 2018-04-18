package gov.nist.microanalysis.roentgen;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.spectrum.EMSAReader;
import gov.nist.microanalysis.roentgen.spectrum.SpectrumFileReader;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class Globals {

   public static final String LIB_NAME = "Roentgen";

   private static final SimplyLazy<Logger> LOGGER = new SimplyLazy<Logger>() {

      @Override
      protected Logger initialize() {
         final Logger res = LogManager.getLogger("Roentgen library");
         return res;
      }
   };

   private final static SimplyLazy<SpectrumFileReader> SPECTRUM_FILE_READER = new SimplyLazy<SpectrumFileReader>() {
      @Override
      protected SpectrumFileReader initialize() {
         final SpectrumFileReader res = new SpectrumFileReader();
         res.register(new EMSAReader.EMSAFileReaderFactory());
         return res;
      }
   };

   static public Logger getLogger() {
      return LOGGER.get();
   }

   static public SpectrumFileReader getSpectrumFileReader() {
      return SPECTRUM_FILE_READER.get();
   }
}
