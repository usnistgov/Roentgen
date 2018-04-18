package com.duckandcover.html;

import java.awt.image.BufferedImage;
import java.awt.image.BufferedImageOp;
import java.io.File;
import java.io.IOException;

import org.imgscalr.Scalr;
import org.imgscalr.Scalr.Method;

/**
 * <p>
 * LinkedReport is a mechanism to embed a link to a new HTML page into the
 * current report.
 * </p>
 * <ul>
 * <li>link - An IToHTML object or HTML string with the HTML to be placed within
 * the link</li>
 * <li>linkedItem - An IToHTML or IToHTMLExt object with HTML which is placed in
 * a new file and is the target of the link.</li>
 * <li>dir - The directory into which to write the linked item.</li>
 * </ul>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class LinkedReport
   implements
   IToHTMLExt {

   private static Method METHOD = Method.BALANCED;

   private final IToHTML mLink;
   private final Report mLinkedItem;

   public static void setDefaultMethod(final Method meth) {
      METHOD = meth;
   }

   /**
    * Constructs a LinkedItem
    *
    * @param link HTML to display in the link in the main report
    * @param linkedItem HTML to display in a separate HTML page
    * @param dir The directory into which to write the separate HTML page
    */
   public LinkedReport(final IToHTML link, final IToHTML linkedItem) {
      mLink = link;
      mLinkedItem = new Report("Linked");
      mLinkedItem.add(linkedItem);
   }

   /**
    * Constructs a LinkedItem
    *
    * @param link HTML to display in the link in the main report
    * @param linkedItem HTML to display in a separate HTML page
    * @param dir The directory into which to write the separate HTML page
    */
   public LinkedReport(final IToHTML link, final Report linkedItem) {
      mLink = link;
      mLinkedItem = linkedItem;
   }

   /**
    * Constructs a LinkedReport
    *
    * @param link HTML to display in the link in the main report
    * @param linkedItem HTML to display in a separate HTML page
    * @param dir The directory into which to write the separate HTML page
    */
   public LinkedReport(final String link, final IToHTML linkedItem) {
      mLink = Transforms.createHTML(link);
      mLinkedItem = new Report("Linked");
      mLinkedItem.add(linkedItem);
   }

   /**
    * Creates an embedded reduced size image and a linked file with the full
    * size image.
    *
    * @param bi BufferedImage
    * @param caption Image caption
    * @param width Size of image in main file
    * @param dir
    * @return LinkedReport
    */
   public static LinkedReport image(final BufferedImage bi, final String caption, final int width) {
      final BufferedImage sbi = Scalr.resize(bi, METHOD, width, new BufferedImageOp[0]);
      final Report rep = new Report("Image");
      rep.addImage(bi, caption);
      return new LinkedReport(new Image(sbi, caption), rep);
   }

   public Report getReport() {
      return mLinkedItem;
   }

   /**
    * @see com.duckandcover.roentgen.html.IToHTML#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      final Report report = new Report("Linked");
      report.add(Transforms.demote(mLink));
      report.add(Transforms.demote(mLinkedItem));
      return report.toHTML(mode);
   }

   /**
    * Creates a link in this HTML document to a second distinct HTML document
    * which is also generated.
    *
    * @param mode
    * @param base
    * @param dir
    * @return String
    * @see com.duckandcover.roentgen.html.IToHTMLExt#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode,
    *      java.io.File, java.lang.String)
    */
   @Override
   public String toHTML(final Mode mode, final File base, final String dir)
         throws IOException {
      final File outdir = new File(base, dir);
      final File outfile = File.createTempFile("index", ".html", outdir);
      mLinkedItem.toFile(outfile, mode);
      return HTML.p(HTML.link(dir, outfile.getName(), Transforms.promote(mLink).toHTML(mode, base, dir)));
   }

}
