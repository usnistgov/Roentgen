package com.duckandcover.scripting;

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.prefs.Preferences;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextPane;
import javax.swing.KeyStroke;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkEvent.EventType;
import javax.swing.event.HyperlinkListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.BadLocationException;
import javax.swing.text.EditorKit;
import javax.swing.text.StyleConstants;
import javax.swing.text.html.HTML;
import javax.swing.text.html.HTMLDocument;

import com.duckandcover.Reliquary;
import com.jgoodies.forms.builder.ButtonBarBuilder;

import jsyntaxpane.DefaultSyntaxKit;

/**
 * <p>
 * This panel represents a self-contained mechanism for editing, running and
 * viewing the output of Jythons commands and scripts.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 309 $
 */
public class JPythonPanel
   extends
   JPanel {

   /**
    *
    */
   private static final int SAVED_HIST_SIZE = 100;

   private static final long serialVersionUID = -65363407319317390L;

   private final JEditorPane jEditorPane_Python = new JEditorPane();
   private final JTextPane jTextPane_Report = new JTextPane();
   private final JSplitPane jSplitPane = new JSplitPane();

   private File mReport;
   private final ScriptingWorker mScripter;
   private boolean mModified = true;
   private final ArrayList<String> mHistory = new ArrayList<>();
   private int mHistoryPosition = 0;
   private int mLastMark = 0;

   private static final class OpenPyScriptListener
      implements
      HyperlinkListener {
      @Override
      public void hyperlinkUpdate(final HyperlinkEvent e) {
         if(e.getEventType() == EventType.ACTIVATED) {
            final String fn = e.getURL().getFile().replace("%20", " ");
            if(fn.endsWith(".py"))
               try {
                  Desktop.getDesktop().open(new File(fn));
               }
               catch(final Exception e1) {
                  Reliquary.getLogger().log(Level.SEVERE, "Unable to open script: " + fn, e1.getMessage());
               }
         }
      }
   }

   private static final class OpenHTMLListener
      implements
      HyperlinkListener {
      @Override
      public void hyperlinkUpdate(final HyperlinkEvent e) {
         if(e.getEventType() == EventType.ACTIVATED) {
            final String fn = e.getURL().getFile();
            if(fn.endsWith(".html") || fn.endsWith(".htm"))
               try {
                  Desktop.getDesktop().browse(e.getURL().toURI());
               }
               catch(final Exception e1) {
                  Reliquary.getLogger().log(Level.SEVERE, "Unable to open HTML file: " + fn, e1);
               }
         }
      }
   }

   public class OpenAction
      extends
      AbstractAction {

      private static final long serialVersionUID = 102210209668238962L;

      public OpenAction(final String title) {
         super(title);
      }

      @Override
      public void actionPerformed(final ActionEvent e) {
         final JFileChooser jfc = new JFileChooser();
         jfc.addChoosableFileFilter(new FileNameExtensionFilter("Python script", "py"));
         jfc.setAcceptAllFileFilterUsed(true);
         jfc.setMultiSelectionEnabled(false);
         final String path = Reliquary.getPreferences().get("Script path", System.getProperty("user.home"));
         jfc.setCurrentDirectory(new File(path));
         final int res = jfc.showOpenDialog(JPythonPanel.this);
         if(res == JFileChooser.APPROVE_OPTION) {
            JPythonPanel.this.execute(jfc.getSelectedFile());
            Reliquary.getPreferences().put("Script path", jfc.getSelectedFile().getPath());
         }
      }
   }

   private class TerminateAction
      extends
      AbstractAction {

      private static final long serialVersionUID = 102210209668238962L;

      private TerminateAction() {
         super("Terminate");
      }

      @Override
      public void actionPerformed(final ActionEvent e) {
         mScripter.terminate();

      }
   }

   private class ExecuteAction
      extends
      AbstractAction {

      private static final long serialVersionUID = 102210209668238962L;

      private ExecuteAction() {
         super("Execute");
      }

      @Override
      public void actionPerformed(final ActionEvent e) {
         final String str = jEditorPane_Python.getText();
         jEditorPane_Python.setText("");
         mHistory.add(str);
         mHistoryPosition = 0;
         savePrefs();
         mScripter.execute(str);
      }
   }

   private class HistoryAction
      extends
      AbstractAction {

      private static final long serialVersionUID = 4245294612740716783L;
      private final int mStep;

      private HistoryAction(final int step) {
         super(step > 0 ? "Previous command" : "Next command");
         mStep = step;
      }

      @Override
      public void actionPerformed(final ActionEvent e) {
         final int newPos = bound(mHistoryPosition - mStep, 1, mHistory.size() + 1);
         if((newPos != mHistoryPosition) && (mHistory.size() > 0)) {
            mHistoryPosition = newPos;
            jEditorPane_Python.setText(mHistory.get(mHistory.size() - mHistoryPosition));
         }
      }

   }

   private int bound(final int val, final int min, final int max) {
      return val < min ? min : (val >= max ? max - 1 : val);
   }

   public void loadHTML(final File html) {
      try {
         jTextPane_Report.setPage(html.toURI().toURL());
         mReport = html;
      }
      catch(final Exception e2) {
         Reliquary.getLogger().log(Level.SEVERE, "Error in loadHTML", e2);
      }
   }

   public JPythonPanel() {
      super(new BorderLayout());
      init();
      restorePrefs();
      mScripter = new ScriptingWorker(this);
      jTextPane_Report.setEditable(false);
      jTextPane_Report.addHyperlinkListener(new OpenPyScriptListener());
      jTextPane_Report.addHyperlinkListener(new OpenHTMLListener());
   }

   public void start() {
      mScripter.execute();
   }

   public void addHyperlinkListener(final HyperlinkListener hl) {
      jTextPane_Report.addHyperlinkListener(hl);
   }

   public void removeHyperlinkListener(final HyperlinkListener hl) {
      jTextPane_Report.removeHyperlinkListener(hl);
   }

   private void savePrefs() {
      final Preferences pref = Reliquary.getPreferences();
      pref.putInt("jpythonpanel.splitter.location", jSplitPane.getDividerLocation());
      if(mHistory.size() > 0) {
         final int i = 0;
         for(final String hist : mHistory.subList(Math.max(0, mHistory.size() - SAVED_HIST_SIZE), mHistory.size()))
            pref.put("jpythonpanel.history" + Integer.toString(i), hist);
      }
   }

   private void restorePrefs() {
      final Preferences up = Preferences.userNodeForPackage(JPythonPanel.class);
      final int pos = up.getInt("jpythonpanel.splitter.location", -1);
      if(pos > 10)
         jSplitPane.setDividerLocation(pos);
      for(int i = 0; i < SAVED_HIST_SIZE; ++i) {
         final String hist = up.get("jpythonpanel.history" + Integer.toString(i), null);
         if(hist != null)
            mHistory.add(hist);
      }
   }

   private void init() {
      final ButtonBarBuilder bb = new ButtonBarBuilder();
      bb.addButton(new ExecuteAction());
      bb.addRelatedGap();
      bb.addButton(new TerminateAction());
      bb.addUnrelatedGap();
      bb.addButton(new OpenAction("Open"));
      final JPanel buttons = bb.getPanel();
      buttons.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
      add(buttons, BorderLayout.NORTH);
      jSplitPane.setOrientation(JSplitPane.VERTICAL_SPLIT);
      DefaultSyntaxKit.initKit();
      /*
       * final Configuration c =
       * DefaultSyntaxKit.getConfig(DefaultSyntaxKit.class); for (final
       * Map.Entry<String, String> me : c.entrySet())
       * Roentgen.getLogger().info(me.getKey() + " = " + me.getValue());
       */
      final JScrollPane jScrollPane = new JScrollPane(jEditorPane_Python);
      Dimension minSize=jScrollPane.getMinimumSize();
      minSize.height*=4;
      jScrollPane.setMinimumSize(minSize);
      jSplitPane.setTopComponent(jScrollPane);
      jEditorPane_Python.setContentType("text/python");
      jEditorPane_Python.addComponentListener(new ComponentAdapter() {
         @Override
         public void componentResized(final ComponentEvent arg0) {
            savePrefs();
         }

         @Override
         public void componentHidden(final ComponentEvent e) {
            savePrefs();
         }
      });
      addEditorKeyEvents();
      this.jSplitPane.setBottomComponent(new JScrollPane(jTextPane_Report));
      add(jSplitPane, BorderLayout.CENTER);
   }

   private void addEditorKeyEvents() {
      final ActionMap map = jEditorPane_Python.getActionMap();
      map.put("prev-history", new HistoryAction(-1));
      map.put("next-history", new HistoryAction(+1));
      map.put("execute-script", new ExecuteAction());

      final InputMap inputMap = jEditorPane_Python.getInputMap(JComponent.WHEN_FOCUSED);
      inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.CTRL_MASK), "prev-history");
      inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.CTRL_MASK), "next-history");
      inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, InputEvent.CTRL_MASK), "execute-script");
   }

   public void execute(final File f) {
      mScripter.execute(f);
   }

   public void execute(final String code) {
      mScripter.execute(code);
   }

   /**
    * Appends a block of HTML code to the end of this document.
    *
    * @param html
    * @throws BadLocationException
    * @throws IOException
    */
   public void append(final String html) {
      append(html, false);
   }

   /**
    * Appends a block of HTML code to the end of this document.
    *
    * @param html
    * @param scrollTo - Scroll to the beginning?
    * @throws BadLocationException
    * @throws IOException
    */
   public void append(final String html, final boolean scrollTo) {
      try {
         // Find the document body and insert just before this...
         final HTMLDocument document = (HTMLDocument) jTextPane_Report.getDocument();
         javax.swing.text.Element body = null;
         {

            final javax.swing.text.Element root = document.getDefaultRootElement();
            for(int i = 0; i < root.getElementCount(); i++) {
               final javax.swing.text.Element element = root.getElement(i);
               if(element.getAttributes().getAttribute(StyleConstants.NameAttribute) == HTML.Tag.BODY) {
                  body = element;
                  break;
               }
            }
            assert body != null;
         }
         if(body != null) {
            if(scrollTo) {
               jTextPane_Report.scrollToReference("ITEM_" + Long.toHexString(mLastMark));
               final String mark = "<A NAME=\"ITEM_" + Integer.toHexString(++mLastMark) + "\" />\n";
               document.insertBeforeEnd(body, mark + html);
               jTextPane_Report.scrollToReference("ITEM_" + Long.toHexString(mLastMark));
            } else
               document.insertBeforeEnd(body, html);
            mModified = true;
         } else
            Reliquary.getLogger().log(Level.WARNING, "Error finding end: " + html);
      }
      catch(final Exception e) {
         Reliquary.getLogger().log(Level.WARNING, "Error appending: " + html, e);
      }
   }

   public void flush() {
      if(mModified && (mReport != null)) {
         if(mLastMark > 0)
            jTextPane_Report.scrollToReference("ITEM_" + Long.toHexString(mLastMark));
         try (final Writer out = new FileWriter(mReport)) {
            final EditorKit editor = jTextPane_Report.getEditorKit();
            final HTMLDocument document = (HTMLDocument) jTextPane_Report.getDocument();
            editor.write(out, document, 0, document.getLength());
            append("<A NAME=\"ITEM_" + Integer.toHexString(++mLastMark) + "\" />\n");
            mModified = false;
         }
         catch(final Exception e) {
            Reliquary.getLogger().log(Level.WARNING, "Error flushing HTML", e);
         }
      }
   }

   /**
    * Writes an image to the report directory and then returns the name of the
    * file in the report path. Add this name to the user generated HTML code.
    *
    * @param img
    * @return String
    * @throws IOException
    */
   public String writeImage(final RenderedImage img)
         throws IOException {
      final File f = File.createTempFile("img", ".png", mReport.getParentFile());
      ImageIO.write(img, "PNG", f);
      return f.getName();
   }

   public void setStartupScript(final File startup) {
      mScripter.setStartupScript(startup);
   }

}
