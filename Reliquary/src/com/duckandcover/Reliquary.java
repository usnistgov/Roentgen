package com.duckandcover;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.logging.StreamHandler;
import java.util.prefs.Preferences;

import javax.swing.JFrame;
import javax.swing.UIManager;

import com.duckandcover.lazy.SimplyLazy;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class Reliquary {

   public final static String APP_PATH = "com.duckandcover.Relinquary";
   public final static String APP_NAME = "Relinquary";
   public final static String SLOGAN = "Scriptable Application Container";

   private final static SimplyLazy<Preferences> PREFERENCES = new SimplyLazy<Preferences>() {
      @Override
      protected Preferences initialize() {
         return Preferences.userNodeForPackage(Reliquary.class);
      }
   };

   final static SimplyLazy<File> REPORT = new SimplyLazy<File>() {

      @Override
      protected File initialize() {
         try {
            return ReportFile.create(System.getProperty("user.home"));
         }
         catch(final IOException e) {
            try {
               return ReportFile.backup();
            }
            catch(final IOException e1) {
               e1.printStackTrace();
               Reliquary.getLogger().log(Level.SEVERE, "Can't create HTML report file", e1);
               System.exit(1);
            }
         }
         return null;
      }
   };

   private final static SimplyLazy<Logger> LOGGER = new SimplyLazy<Logger>() {

      @Override
      protected Logger initialize() {
         final Logger res = Logger.getLogger(APP_PATH);
         try {
            final File f = new File(System.getProperty("user.dir"), APP_NAME + ".log");
            final FileOutputStream fos = new FileOutputStream(f);
            final StreamHandler streamHandler = new StreamHandler(fos, new SimpleFormatter());
            streamHandler.setLevel(Level.WARNING);
            res.addHandler(streamHandler);
         }
         catch(final Exception fe) {
            final StreamHandler streamHandler = new StreamHandler(System.err, new SimpleFormatter());
            streamHandler.setLevel(Level.WARNING);
            res.addHandler(streamHandler);
            res.log(Level.WARNING, "Unable to create file logging stream", fe);
         }
         return res;
      }
   };

   private static final SimplyLazy<Reliquary> mInstance = new SimplyLazy<Reliquary>() {

      @Override
      protected Reliquary initialize() {
         return new Reliquary();
      }

   };

   private final JReliquaryFrame mMainFrame;

   public static int[] getJavaVersion() {
      final String tmp = System.getProperty("java.version");
      final String[] javaVer = tmp.split("\\.");
      final String[] last = javaVer[javaVer.length - 1].split("_");
      final int len = javaVer.length - 1 + last.length;
      final int[] res = new int[len];
      for(int i = 0; i < javaVer.length - 1; ++i)
         res[i] = Integer.parseInt(javaVer[i]);
      for(int i = 0; i < last.length; ++i)
         res[i + javaVer.length - 1] = Integer.parseInt(last[i]);
      return res;
   }

   public static boolean versionHigherThanOrEqualTo(final int[] current, final int[] test) {
      for(int i = 0; i < current.length; ++i)
         if(current[i] < test[i])
            return false;
      return true;
   }

   public static boolean argExists(final String[] args, final String argItem) {
      for(final String arg : args)
         if(arg.startsWith(argItem))
            return true;
      return false;
   }

   private Reliquary() {
      mMainFrame = new JReliquaryFrame();
      mMainFrame.setVisible(true);
   }

   public JReliquaryFrame getMainFrame() {
      return mMainFrame;
   }

   public static Reliquary instance() {
      return mInstance.get();
   }
   
   public static boolean inEclipse() {
       String inEclipseStr = System.getProperty("runInEclipse");
       return "true".equalsIgnoreCase(inEclipseStr);
   }

   /**
    * @param args
    */
   public static void main(final String[] args) {
      final String os = System.getProperty("os.name").toLowerCase();
      if((os.indexOf("mac") >= 0) || (os.indexOf("os x") >= 0)) {
         System.setProperty("apple.laf.useScreenMenuBar", "true");
         System.setProperty("apple.laf.smallTabs", "true");
      }
      // This eliminates an exception in Jython 2.7.0
      // "console: Failed to install '':
      // java.nio.charset.UnsupportedCharsetException: cp0"
      System.setProperty("python.console.encoding", "UTF-8");
      try {
         final String osName = System.getProperty("os.name");
         // Windows 10 Creator's Update introduced a bug that causes the system
         // look-and-feel to crash for JRE versions less than 1.8.0_144. Use
         // Nimbus look-and-feel instead.
         final int[] testVer = new int[] {
            1,
            8,
            0,
            144
         };
         if((osName.endsWith("10") && (!versionHigherThanOrEqualTo(getJavaVersion(), testVer))) || argExists(args, "-nimbus"))
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
         else
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      }
      catch(final Exception e) {
         try {
            e.printStackTrace();
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
         }
         catch(final Exception e1) {
            e1.printStackTrace();
         }
      }
      JFrame.setDefaultLookAndFeelDecorated(false);
      Reliquary.instance();
   }

   public static Logger getLogger() {
      return LOGGER.get();
   }

   public static Preferences getPreferences() {
      return PREFERENCES.get();
   }

   public static File getReport() {
      return REPORT.get();
   }
}
