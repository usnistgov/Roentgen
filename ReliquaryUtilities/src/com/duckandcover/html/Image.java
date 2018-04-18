package com.duckandcover.html;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * <p>
 * Description
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class Image
   implements
   IToHTMLExt {

   private final BufferedImage mImage;
   private final String mCaption;

   /**
    * Constructs a Image
    */
   public Image(final BufferedImage bi, final String caption) {
      mImage = bi;
      mCaption = caption;
   }

   /**
    * @param mode
    * @return
    * @see com.duckandcover.roentgen.html.IToHTML#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      if(mode == Mode.VERBOSE)
         return HTML.b("Image[" + mImage.getWidth() + " pixels &times; " + mImage.getHeight() + " pixels]: ") + mCaption;
      else
         return "&nbsp;";
   }

   /**
    * @param mode
    * @param base
    * @param dir
    * @return
    * @see com.duckandcover.roentgen.html.IToHTMLExt#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode,
    *      java.io.File, java.lang.String)
    */
   @Override
   public String toHTML(final Mode mode, final File base, final String dir)
         throws IOException {
      return HTML.image(mImage, mode != Mode.TERSE ? mCaption : null, base, dir);
   }

}
