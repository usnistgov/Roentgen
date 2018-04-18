package com.duckandcover.html;

import java.io.File;
import java.io.IOException;

/**
 * <p>
 * An extension to IToHTML that adds the ability to include links to images or
 * other files within the HTML content.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public interface IToHTMLExt
   extends
   IToHTML {

   /**
    * Specifies the location to place the output files. The HTML file is assumed
    * to be stored in 'base' and so the output files are referenced relative to
    * 'base' in 'dir'.
    *
    * @param mode VERBOSE, NORMAL or TERSE
    * @param base The directory into which to write the images / other data
    *           items.
    */
   public String toHTML(Mode mode, File base, String dir)
         throws IOException;
}
