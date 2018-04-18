package gov.nist.microanalysis.roentgen.spectrum;

import javax.swing.filechooser.FileFilter;

/**
 * <p>
 * An interface implemented by classed that provide a mechanism to read
 * {@link EDSSpectrum} objects using the {@link ISpectrumReader} interface.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 199 $
 */
public interface ISpectrumReaderFactory {

   public ISpectrumReader get();

   /**
    * Gets a FileFilter to associate with this spectrum file type.
    *
    * @return
    */
   public FileFilter getFileFilter();

}
