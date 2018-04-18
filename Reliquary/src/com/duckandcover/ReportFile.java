package com.duckandcover;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

/**
 * <p>
 * A helper class to build HTML files
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class ReportFile {

   public static File create(final String base)
         throws IOException {
      final Calendar c = Calendar.getInstance();
      c.setTime(new Date(System.currentTimeMillis()));
      final File baseDir = new File(base);
      final File year = new File(baseDir, Integer.toString(c.get(Calendar.YEAR)));
      final File month = new File(year, Integer.toString(c.get(Calendar.MONTH) + 1));
      final File day = new File(month, Integer.toString(c.get(Calendar.DAY_OF_MONTH)));
      final File report = new File(day, "index" + Integer.toString(c.get(Calendar.HOUR_OF_DAY)) + "_"
            + Integer.toString(c.get(Calendar.MINUTE)) + "_" + Integer.toString(c.get(Calendar.SECOND)) + "."
            + Long.toString(System.currentTimeMillis() % 1000) + ".html");
      buildBaseReport(c, report);
      return report;
   }

   /**
    * @param c
    * @param report
    * @throws IOException
    */
   private static void buildBaseReport(final Calendar c, final File report)
         throws IOException {
      report.mkdirs();
      final Report r = new Report(Reliquary.APP_NAME + " Report");
      final SimpleDateFormat sdf = new SimpleDateFormat();
      sdf.setCalendar(c);
      final StringBuffer html = new StringBuffer();
      html.append(sdf.format(c.getTime()));
      html.append(HTML.br());
      html.append(report.getPath());
      html.append(HTML.br());
      html.append(System.getProperty("java.vendor") + " Java " + System.getProperty("java.version"));
      html.append(HTML.br());
      html.append(System.getProperty("java.vm.vendor") + " " + System.getProperty("java.vm.name") + " "
            + System.getProperty("java.vm.version"));
      r.addHTML(HTML.alignRight(html.toString()));
      // r.addHeader(Reliquary.APP_NAME + HTML.normal(" - " + Reliquary.SLOGAN));
      r.toFile(report, Mode.NORMAL);
   }

   public static File backup()
         throws IOException {
      final File f = File.createTempFile("index", ".html");
      final Calendar c = Calendar.getInstance();
      c.setTime(new Date(System.currentTimeMillis()));
      buildBaseReport(c, f);
      return f;
   }

}
