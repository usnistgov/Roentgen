package com.duckandcover.html;

import java.awt.Desktop;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * <p>
 * Provides a mechanism to combine objects implementing IToHTML and IToHTMLExt
 * into a single enlongated report.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class Report
   implements
   IToHTMLExt {

   private final String mName;
   private final List<IToHTML> mItems;
   private final Map<IToHTML, IToHTML.Mode> mCustomMode;

   /**
    * Constructs a HTMLReport
    */
   public Report(final String name) {
      mName = name;
      mItems = new ArrayList<>();
      mCustomMode = new HashMap<>();
   }

   public void addHeader(final IToHTML item) {
      mItems.add(Transforms.h3(item));
   }

   public void addHeader(final String html) {
      addHeader(Transforms.createHTML(html));
   }

   public void addSubHeader(final IToHTML item) {
      mItems.add(Transforms.h4(item));
   }

   public void addSubHeader(final String html) {
      addSubHeader(Transforms.createHTML(html));
   }

   public void add(final List<Object> items, final IToHTML.Mode mode) {
      add(Transforms.createList(items), mode);
   }

   public void add(final List<Object> items) {
      mItems.add(Transforms.createList(items));
   }

   public void addHTML(final String html) {
      mItems.add(Transforms.createHTML(html));
   }

   public void addTable(final List<String> colHeaders, final List<String> rowHeaders, final List<List<String>> items) {
      addHTML(HTML.asTable(colHeaders, rowHeaders, items));
   }

   public void addImage(final BufferedImage img, final String caption) {
      add(Transforms.scrollPane(new Image(img, caption)));
   }

   /**
    * Adds a block of text. Special characters are escaped to make the text HTML
    * friendly so no HTML code should be present.
    *
    * @param note
    */
   public void addNote(final String note) {
      addHTML(HTML.p(HTML.escape(note)));
   }

   public void addParagraph(final IToHTML item) {
      add(Transforms.createParagraph(item));
   }

   public void add(final IToHTML item, final IToHTML.Mode mode) {
      mItems.add(item);
      mCustomMode.put(item, mode);
   }

   public void add(final IToHTML item) {
      mItems.add(item);
   }

   public void add(final Table table) {
      mItems.add(Transforms.scrollPane(table));
   }

   @Override
   public String toHTML(final Mode mode) {
      final StringBuffer sb = new StringBuffer();
      for(final IToHTML item : mItems)
         sb.append(item.toHTML(mCustomMode.getOrDefault(item, mode)));
      return sb.toString();
   }

   public String highest(final IToHTML item, final Mode mode, final File base, final String dir)
         throws IOException {
      if(item instanceof IToHTMLExt)
         return ((IToHTMLExt) item).toHTML(mCustomMode.getOrDefault(item, mode), base, dir);
      else
         return item.toHTML(mCustomMode.getOrDefault(item, mode));
   }

   @Override
   public String toHTML(final Mode mode, final File base, final String dir)
         throws IOException {
      final File fd = new File(base, dir);
      final StringBuffer sb = new StringBuffer();
      if(fd.isDirectory() || fd.mkdirs()) {
         for(final IToHTML item : mItems)
            sb.append(highest(item, mode, base, dir));
      } else {
         sb.append(HTML.error("Unable to create the report directory. Defaulting to the base report."));
         for(final IToHTML item : mItems)
            sb.append(item.toHTML(mode));

      }
      return sb.toString();
   }

   /**
    * Writes the report to the specified file.
    *
    * @param f {@link File}
    * @param mode {@link Mode}
    * @throws FileNotFoundException
    */
   public void toFile(final File f, final Mode mode)
         throws IOException {
      HTML.toFile(this, f, mode, mName);
   }

   /**
    * Writes the report to a temporary file and then opens the file in the
    * browser.
    *
    * @param mode
    * @throws IOException
    */
   public void inBrowser(final Mode mode)
         throws IOException {
      final File f = File.createTempFile("temp", ".html");
      HTML.toFile(this, f, mode, mName);
      Desktop.getDesktop().browse(f.toURI());
   }

   @Override
   public String toString() {
      return mName;
   }
}
