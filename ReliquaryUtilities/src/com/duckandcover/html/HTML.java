package com.duckandcover.html;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URI;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.commons.text.StringEscapeUtils;


import com.duckandcover.html.IToHTML.Mode;

/**
 * <p>
 * These static functions make it easier to write consistent accurate HTML.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2018
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class HTML {

   public static String escape(final String html) {
      if(html.startsWith("<html>") || html.startsWith("<HTML>"))
         return html.substring(6);
      else
         return StringEscapeUtils.escapeHtml4(html);
   }

   public static String error(final String html) {
      return "<p class=\"error\"><b>ERROR:</b> " + html + "</p>";
   }

   public static String warning(final String html) {
      return "<p class=\"warn\"><b>WARNING:</b> " + html + "</p>";
   }

   public static String result(final String string) {
      return "<p class=\"result\">" + escape(string) + "</p>";
   }

   public static String link(final String dir, final String filename, final String html) {
      return "<a href=\"" + dir + "/" + filename + "\">" + html + "</a>";
   }

   public static String link(final URI uri, final String html)
         throws MalformedURLException {
      return "<a href=\"" + uri.toURL().toString() + "\">" + html + "</a>";
   }

   public static String header(final String html) {
      return "<h3>" + html + "</h3>";
   }

   public static String subHeader(final String html) {
      return "<h4>" + html + "</h4>";
   }

   public static String tiny(final String html) {
      return "<font size=\"1\">" + html + "</font>";
   }

   public static String small(final String html) {
      return "<font size=\"2\">" + html + "</font>";
   }

   public static String normal(final String html) {
      return "<font size=\"3\">" + html + "</font>";
   }

   public static String large(final String html) {
      return "<font size=\"4\">" + html + "</font>";
   }

   public static String Large(final String html) {
      return "<font size=\"5\">" + html + "</font>";
   }

   public static String huge(final String html) {
      return "<font size=\"6\">" + html + "</font>";
   }

   public static String Huge(final String html) {
      return "<font size=\"7\">" + html + "</font>";
   }

   public static String scrollPane(final String html) {
      return "<div style=\\\"overflow-x:auto;\">" + html + "</div>";
   }

   public static String image(final BufferedImage img, final String caption, final File base, final String dir)
         throws IOException {
      final File outpath = new File(base, dir);
      if(!outpath.exists())
         if(!outpath.mkdirs())
            return HTML.error("Unable to create output path for " + caption);
      final StringBuffer sb = new StringBuffer();
      File pngFile;
      pngFile = File.createTempFile("img", ".png", outpath);
      if(writePNG(img, pngFile)) {
         sb.append("<img ");
         sb.append("width=\"" + Integer.toString(img.getWidth()) + "\"");
         sb.append(" height=\"" + Integer.toString(img.getHeight()) + "\" ");
         sb.append("src=\"" + dir + "/" + pngFile.getName() + "\"");
         sb.append(" alt=\"" + caption + "\" /></a>");
         if((caption != null) && (caption.length() > 0))
            sb.append(HTML.p(HTML.b("Figure: ") + caption));
      } else
         sb.append(HTML.error("Error writing image " + caption));
      return sb.toString();
   }

   /**
    * Converts any object to HTML using the IToHTML interface when available or
    * escaped toString() otherwise.
    *
    * @param obj
    * @param mode
    * @return String in HTML
    */
   public static String toHTML(final Object obj, final Mode mode) {
      if(obj == null)
         return HTML.error("NULL");
      if(obj instanceof IToHTML)
         return ((IToHTML) obj).toHTML(mode);
      return HTML.escape(obj.toString());
   }
   
   
   /**
    * Converts any object to HTML using the IToHTML interface when available or
    * escaped toString() otherwise.
    *
    * @param obj
    * @param mode
    * @return String in HTML
 * @throws IOException 
    */
   public static String toHTMLExt(final Object obj, final Mode mode, File base, String dir) throws IOException {
      if(obj == null)
         return HTML.error("NULL");
      if(obj instanceof IToHTMLExt)
    	  return ((IToHTMLExt)obj).toHTML(mode, base, dir);
      if(obj instanceof IToHTML)
         return ((IToHTML) obj).toHTML(mode);
      return HTML.escape(obj.toString());
   }

   /**
    * Occasionally writing to pngFile fails. Delete it and try a second time.
    *
    * @param img
    * @param pngFile
    * @return
    * @throws IOException
    */
   private static boolean writePNG(final BufferedImage img, final File pngFile)
         throws IOException {
      pngFile.delete();
      return ImageIO.write(img, "png", pngFile);
   }

   /**
    * Returns the specified integer as an ordinal string of the form 1st, 2nd,
    * 3rd etc.
    *
    * @param n
    * @return String in HTML
    */
   public static String asOrdinalHTML(final int n) {
      if((n >= 11) && (n <= 13))
         return Integer.toString(n) + HTML.sup("th");
      else
         switch(n % 10) {
            case 1:
               return Integer.toString(n) + HTML.sup("st");
            case 2:
               return Integer.toString(n) + HTML.sup("nd");
            case 3:
               return Integer.toString(n) + HTML.sup("rd");
            default:
               return Integer.toString(n) + HTML.sup("th");
         }
   }

   private static String td(final String html) {
      return "<td>" + html + "</td>";
   }

   private static String th(final String html) {
      return "<th>" + html + "</th>";
   }

   /**
    * Builds a simple table with optional column and row headers.
    *
    * @param colHeaders A list of column headers (in HTML) (maybe null)
    * @param rowHeaders A list of column headers (in HTML) (maybe null)
    * @param items A list of lists of table items (in HTML)
    * @return String as HTML
    */
   public static String asTable(final List<String> colHeaders, final List<String> rowHeaders, final List<List<String>> items) {
      final StringBuffer sb = new StringBuffer();
      sb.append("<table>");
      if(colHeaders != null) {
         sb.append("<tr>");
         if(rowHeaders != null)
            sb.append(td("&nbsp;"));
         for(final String item : colHeaders)
            sb.append(th(item));
         sb.append("</tr>");
      }
      for(int r = 0; r < items.size(); ++r) {
         sb.append("<tr>");
         if(rowHeaders != null)
            sb.append(th(r < rowHeaders.size() ? rowHeaders.get(r) : "&nbsp;"));
         for(final String item : items.get(r))
            sb.append(td(item));
         sb.append("</tr>");
      }
      sb.append("</table>");
      return sb.toString();
   }

   public static String alignLeft(final String html) {
      return "<div align=\"left\">" + html + "</div>";
   }

   public static String alignCenter(final String html) {
      return "<div align=\"center\">" + html + "</div>";
   }

   public static String alignRight(final String html) {
      return "<div align=\"right\">" + html + "</div>";

   }

   public static String p(final String html) {
      return "<p>" + html + "</p>";
   }

   public static String i(final String html) {
      return "<i>" + html + "</i>";
   }

   public static String b(final String html) {
      return "<b>" + html + "</b>";
   }

   public static String br() {
      return "<br/>";
   }

   public static String sup(final String html) {
      return "<sup>" + html + "</sup>";
   }

   public static String sub(final String html) {
      return "<sub>" + html + "</sub>";
   }

   public static String createPageEnd() {
      final StringBuffer sb = new StringBuffer();
      sb.append("</body>\n");
      sb.append("</html>");
      return sb.toString();
   }

   public static String createPageHeader(final String reportName, final String cssName) {
      final StringBuffer sb = new StringBuffer();
      sb.append("<html>\n");
      sb.append("<head>\n");
      if(cssName != null)
         sb.append(" <link rel=stylesheet type=\"text/css\" href=\"" + cssName + "\" />\n");
      sb.append(" <title>" + HTML.escape(reportName) + "</title>\n");
      sb.append("</head>\n");
      sb.append("<body>");
      return sb.toString();
   }

   public static void generateDefaultStyleSheet(final File reportDir, final String cssName)
         throws IOException {
      if(!(new File(reportDir, cssName)).exists())
         try (final InputStream cssStream = HTML.class.getResourceAsStream("style.css")) {
            reportDir.mkdirs();
            final Path path = Paths.get(reportDir.getPath(), cssName);
            Files.copy(cssStream, path, new CopyOption[] {});
         }
   }

   /**
    * Wraps the item in page header / tailer and writes the contents to the
    * specified file.
    *
    * @param item
    * @param f
    * @param mode
    * @param dir
    * @throws FileNotFoundException
    */
   public static void toFile(final IToHTML item, final File f, final Mode mode, final String dir)
         throws IOException {
      try (final PrintWriter pw = new PrintWriter(f)) {
         String cssName = "report.css";
         try {
            HTML.generateDefaultStyleSheet(f.getParentFile(), cssName);
         }
         catch(final Exception e) {
            cssName = null;
         }
         pw.println(HTML.createPageHeader("Report", cssName));
         pw.println(Transforms.promote(item).toHTML(mode, f.getParentFile(), dir));
         pw.println(HTML.createPageEnd());
      }
   }

   /**
    * Removes all HTML tags and replaces all HTML letter encodings with the
    * unicode equivalents.
    *
    * @param html
    * @return String
    */
   public static String stripTags(final String html) {
      final StringBuffer sb = new StringBuffer();
      int tagCx = 0;
      for(int i = 0; i < html.length(); ++i) {
         final char c = html.charAt(i);
         if(c == '<')
            ++tagCx;
         else if(c == '>')
            --tagCx;
         else if(tagCx == 0)
            sb.append(c);
      }
      return StringEscapeUtils.unescapeHtml4(sb.toString());
   }
}
