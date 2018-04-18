package com.duckandcover.swing;

import java.awt.Color;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Image;
import java.awt.LayoutManager;
import java.awt.SystemColor;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.WindowConstants;
import javax.swing.border.Border;
import javax.swing.border.CompoundBorder;

import com.jgoodies.forms.builder.ButtonBarBuilder;
import com.jgoodies.forms.factories.CC;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

/**
 * <p>
 * JWizardDialog<W> is a step-by-step GUI mechanism to edit the contents of an
 * object of type W. Each step is implemented by an instance of
 * {@link JWizardPanel}. {@link JWizardPanel} initializes itself based on the
 * output of the previous step and when the panel is hidden it updates its
 * internal copy of the W object which then serves to initialize the next panel
 * (or act as the final result.)
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class JWizardDialog<W>
   extends
   JDialog {
   /**
    *
    */
   private static final long serialVersionUID = 5516307745683470014L;
   // Controls
   private final JLabel jLabel_Banner = new JLabel();
   private final JLabel jLabel_PrevLabel = new JLabel();
   private final JLabel jLabel_NextLabel = new JLabel();
   private final JPanel jPanel_Main = new JPanel();
   private final JPanel jPanel_Contents;
   private final JButton jButton_Back = new JButton("Back");
   private final JButton jButton_Next = new JButton("Next");
   private final JButton jButton_Cancel = new JButton("Cancel");
   private final JButton jButton_Finish = new JButton("Finish");
   private final JLabel jLabel_ErrorLabel = new JLabel();
   private final JButton jButton_ErrorBtn = new JButton("More...");
   private JLabel jLabel_Icon = new JLabel();
   // Status information
   private final ArrayList<JWizardPanel<W>> mPreviousPanels = new ArrayList<>();
   private final ArrayList<String> mPreviousBanners = new ArrayList<>();
   private JWizardPanel<W> jWizardPanel_Active = null;
   private JWizardPanel<W> jWizardPanel_Next = null;
   private String mNextBanner = null;
   private boolean mCancelled = false;
   private String mDialogMsg = "";
   private final ArrayList<String> mLongErrorMsg = new ArrayList<>();

   private final W mInitial;
   private W mResult;

   public static final int FINISHED = 1;
   public static final int CANCELLED = 2;

   final AbstractAction mCancelAction = new AbstractAction() {
      private static final long serialVersionUID = 585668420885106696L;

      @Override
      public void actionPerformed(final ActionEvent e) {
         mCancelled = true;
         mResult = null;
         dispose();
      }
   };

   public JWizardDialog(final W initial, final Frame owner, final String title) {
      this(initial, owner, title, new Dimension(400, 170));
   }

   public JWizardDialog(final W initial, final Dialog owner, final String title) {
      this(initial, owner, title, new Dimension(400, 170));
   }

   public JWizardDialog(final W initial, final Frame owner, final String title, final Dimension dim) {
      super(owner, title, true);
      mInitial = initial;
      jPanel_Contents = new JPanel(new FormLayout("center:" + dim.width + "dlu", "center:" + dim.width + "dlu"));

      try {
         wizInitialize(dim);
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   public JWizardDialog(final W initial, final Dialog owner, final String title, final Dimension dim) {
      super(owner, title, true);
      mInitial = initial;
      jPanel_Contents = new JPanel(new FormLayout("center:" + dim.width + "dlu", "center:" + dim.width + "dlu"));
      try {
         wizInitialize(dim);
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   private void wizInitialize(final Dimension dim)
         throws Exception {
      setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
      final LayoutManager layout = new FormLayout( //
            "8dlu, " + Integer.toString(dim.width + 20) + "dlu, 8dlu", //
            "8dlu, pref, 8dlu, " + Integer.toString(dim.height + 20) + "dlu, 8dlu, pref, 8dlu, pref, 8dlu" //
      );
      jPanel_Main.setLayout(layout);
      {
         createTop();
      }

      jButton_Back.addActionListener(new ActionListener() {
         @Override
         public void actionPerformed(final ActionEvent e) {
            assert (jWizardPanel_Active != null);
            assert (mPreviousPanels.size() > 0);
            assert (mPreviousPanels.size() == mPreviousBanners.size());
            if(jWizardPanel_Active != null) {
               jPanel_Contents.remove(jWizardPanel_Active);
               setNextPanel(jWizardPanel_Active, jLabel_Banner.getText());
            }
            final JWizardPanel<W> panel = mPreviousPanels.remove(mPreviousPanels.size() - 1);
            final String banner = mPreviousBanners.remove(mPreviousBanners.size() - 1);
            setActivePanel(panel, banner);
            jButton_Back.setEnabled(mPreviousPanels.size() > 0);
         }
      });

      jButton_Next.addActionListener(new ActionListener() {
         @Override
         public void actionPerformed(final ActionEvent e) {
            if(jWizardPanel_Active.permitNext()) {
               assert (jWizardPanel_Active != null);
               assert (jWizardPanel_Next != null);
               mPreviousPanels.add(jWizardPanel_Active);
               mPreviousBanners.add(jLabel_Banner.getText());
               jPanel_Contents.remove(jWizardPanel_Active);
               final JWizardPanel<W> next = jWizardPanel_Next;
               final String banner = mNextBanner;
               setNextPanel(null, null);
               next.initialize(jWizardPanel_Active.getData());
               setActivePanel(next, banner);
               jButton_Back.setEnabled(mPreviousPanels.size() > 0);
            } else
               Toolkit.getDefaultToolkit().beep();
         }
      });

      jButton_Finish.addActionListener(new ActionListener() {
         @Override
         public void actionPerformed(final ActionEvent e) {
            if(jWizardPanel_Active.permitNext()) {
               jWizardPanel_Active.onHide();
               mResult = jWizardPanel_Active.getData();
               setVisible(false);
            }
         }
      });
      jButton_Cancel.addActionListener(mCancelAction);

      final ButtonBarBuilder bbb = new ButtonBarBuilder();
      bbb.addGlue();
      bbb.addButton(jButton_Back, jButton_Next);
      bbb.addUnrelatedGap();
      bbb.addButton(jButton_Finish);
      bbb.addUnrelatedGap();
      bbb.addButton(jButton_Cancel);
      final JPanel btns = bbb.build();
      jButton_Next.setMnemonic(KeyEvent.VK_N);
      jButton_Back.setMnemonic(KeyEvent.VK_B);
      jButton_Cancel.setMnemonic(KeyEvent.VK_ESCAPE);

      jPanel_Main.add(btns, CC.xy(2, 8));
      setContentPane(jPanel_Main);
      setResizable(false);
      addCancelByEscapeKey();
      setForeground(SystemColor.controlText);
      setBackground(SystemColor.control);
   }

   private void addCancelByEscapeKey() {
      final String CANCEL_ACTION_KEY = "CANCEL_ACTION_KEY";
      final int noModifiers = 0;
      final KeyStroke escapeKey = KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, noModifiers, false);
      final InputMap inputMap = getRootPane().getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
      inputMap.put(escapeKey, CANCEL_ACTION_KEY);
      getRootPane().getActionMap().put(CANCEL_ACTION_KEY, mCancelAction);
   }

   public boolean isFinished() {
      return (mResult != null) && (!mCancelled);
   }

   public boolean isCancelled() {
      return mCancelled;
   }

   /**
    * Show the wizard dialog and wait for the user to complete it.
    *
    * @return FINISHED or CANCELLED depending upon whether the Wizard was exited
    *         using the Finish button
    */
   public int showWizard() {
      setVisible(true);
      return (mResult != null) ? FINISHED : CANCELLED;
   }

   public void setMessageText(final String msg) {
      if(mLongErrorMsg.size() == 0) {
         jLabel_ErrorLabel.setForeground(SystemColor.textText);
         jLabel_ErrorLabel.setText(msg);
         jButton_ErrorBtn.setEnabled(false);
      }
   }

   public void setErrorText(final String error) {
      jLabel_ErrorLabel.setForeground(Color.red);
      jLabel_ErrorLabel.setText(error);
      jButton_ErrorBtn.setEnabled(false);
   }

   public void setExceptionText(final Throwable th) {
      setExceptionText(th.toString(), th);
   }

   public void setExceptionText(final String shortMsg, final Throwable th) {
      setExceptionText(shortMsg, shortMsg, th);
   }

   public void setExceptionText(final String shortMsg, final String dialogMsg, final Throwable th) {
      th.printStackTrace();
      String longMsg = th.getMessage();
      if(longMsg == null) {
         final StringWriter sw = new StringWriter();
         final PrintWriter pw = new PrintWriter(sw);
         th.printStackTrace(pw);
         longMsg = sw.toString();
      }
      setExtendedError(shortMsg, dialogMsg, longMsg != null ? longMsg : "No additional information available.");
   }

   public void setExtendedError(final String shortMsg, final String longMsg) {
      setExtendedError(shortMsg, shortMsg, longMsg);
   }

   public void setExtendedError(final String errorMsg, final String dialogMsg, final String longMsg) {
      mDialogMsg = dialogMsg != null ? dialogMsg : "Click details for more information";
      mLongErrorMsg.add(mDialogMsg + ": " + (longMsg != null ? longMsg : "No addional information available."));
      if(mLongErrorMsg.size() == 1)
         setErrorText(errorMsg);
      else
         setErrorText(mLongErrorMsg.size() + " errors...");
      jButton_ErrorBtn.setEnabled(true);
      // If this method is called more than once, we must ensure that the Error
      // Button will not be called more than once with older error messages
      final ActionListener[] al = jButton_ErrorBtn.getActionListeners();
      for(int index = 0; jButton_ErrorBtn.getActionListeners().length > 0; index++)
         jButton_ErrorBtn.removeActionListener(al[index]);
      jButton_ErrorBtn.addActionListener(new ActionListener() {
         @Override
         public void actionPerformed(final ActionEvent e) {
            final StringBuffer sb = new StringBuffer();
            sb.append(mLongErrorMsg.size() + " errors:");
            for(int i = 0; i < mLongErrorMsg.size(); ++i) {
               sb.append("\nError " + (i + 1) + "\n");
               sb.append(mLongErrorMsg.get(i));
            }
            JErrorDialog.createErrorMessage(JWizardDialog.this, "Extended error message", mLongErrorMsg.size() == 1 ? mDialogMsg
                  : mLongErrorMsg.size() + " errors...", sb.toString());
         }
      });
   }

   public void clearMessageText() {
      jLabel_ErrorLabel.setText("");
      jButton_ErrorBtn.setEnabled(false);
      mLongErrorMsg.clear();
   }

   public void centerDialog(final Dialog d) {
      d.setLocationRelativeTo(this);
   }

   /**
    * Give ourselves a fresh start.
    */
   public void clearWizard() {
      jWizardPanel_Active = null;
      jLabel_Banner.setText("Nothing");
      jLabel_PrevLabel.setText("Start");
      jLabel_NextLabel.setText("Finish");
      mPreviousPanels.clear();
      mPreviousBanners.clear();
      jWizardPanel_Next = null;
      mNextBanner = "";
      jButton_Back.setEnabled(false);
      jButton_Next.setEnabled(false);
      jButton_Finish.setEnabled(true);
   }

   public void setActivePanel(final JWizardPanel<W> panel, final String banner) {
      jLabel_ErrorLabel.setText("");
      if(jWizardPanel_Active != null) {
         jWizardPanel_Active.onHide();
         jPanel_Contents.remove(jWizardPanel_Active);
         jWizardPanel_Active = null;
      }
      if(panel != null) {
         jPanel_Contents.add(panel, CC.xy(1, 1));
         clearMessageText();
         panel.onShow();
         jPanel_Contents.repaint();
      }
      jWizardPanel_Active = panel;
      jLabel_Banner.setText(banner != null ? banner : "None");
      if(jWizardPanel_Next != null)
         jLabel_NextLabel.setText("<html>Next: <em>" + mNextBanner + "</em>");
      if(mPreviousBanners.size() > 0) {
         final String str = mPreviousBanners.get(mPreviousBanners.size() - 1);
         jLabel_PrevLabel.setText("<html>Previous: <em>" + str + "</em>");
      } else
         jLabel_PrevLabel.setText("First page");
      jButton_Back.setEnabled(mPreviousPanels.size() > 0);
      jButton_Next.setEnabled(jWizardPanel_Next != null);
   }

   public void setNextPanel(final JWizardPanel<W> panel, final String banner) {
      jWizardPanel_Next = panel;
      mNextBanner = banner;
      if(jWizardPanel_Next != null)
         jLabel_NextLabel.setText("<html>Next: <em>" + mNextBanner + "</em>");
      else
         jLabel_NextLabel.setText("Finish");
      jButton_Next.setEnabled(jWizardPanel_Next != null);
   }

   public void enableFinish(final boolean b) {
      if(jButton_Finish.isEnabled() != b) {
         jButton_Finish.setEnabled(b);
         getRootPane().setDefaultButton(jButton_Finish.isEnabled() ? jButton_Finish : null);
      }
   }

   public void enableNext(final boolean b) {
      if(jButton_Next.isEnabled() != b) {
         jButton_Next.setEnabled(b);
         getRootPane().setDefaultButton(jButton_Next.isEnabled() ? jButton_Next : null);
      }
   }

   public void setBackEnabled(final boolean b) {
      if(b != jButton_Back.isEnabled())
         jButton_Back.setEnabled(b);
   }

   public boolean isBackEnabled() {
      return jButton_Back.isEnabled();
   }

   private JPanel createIconPanel() {
      final JPanel header = new JPanel(new FormLayout("pref, 3dlu, default", "pref"));

      jLabel_Icon = new JLabel();
      jLabel_Icon.setIcon(new ImageIcon(JWizardDialog.class.getResource("ClipArt/alien_sm.png")));
      header.add(jLabel_Icon, CC.xy(1, 1));
      return header;
   }

   public Icon getIcon() {
      return jLabel_Icon.getIcon();
   }

   public void setIcon(final Icon icon) {
      if(icon instanceof ImageIcon) {
         final ImageIcon imgIcon = (ImageIcon) icon;
         jLabel_Icon.setIcon(new ImageIcon(imgIcon.getImage().getScaledInstance(64, 59, Image.SCALE_AREA_AVERAGING)));
      } else
         jLabel_Icon.setIcon(icon);
   }

   private void createTop() {
      final JPanel header = createIconPanel();
      {
         final JPanel banner = new JPanel(new FormLayout("320dlu", "12dlu, 1dlu, 20dlu, 1dlu, 12dlu"));
         jLabel_Banner.setFont(new Font("Dialog", Font.BOLD, 24));
         jLabel_Banner.setText("Test banner");
         banner.add(jLabel_Banner, CC.xy(1, 3));
         jLabel_PrevLabel.setFont(new Font("Dialog", Font.PLAIN, 12));
         jLabel_PrevLabel.setText("Start ");
         banner.add(jLabel_PrevLabel, CC.xy(1, 1, CellConstraints.RIGHT, CellConstraints.CENTER));
         jLabel_NextLabel.setFont(new Font("Dialog", Font.PLAIN, 12));
         jLabel_NextLabel.setText("Finish ");
         banner.add(jLabel_NextLabel, CC.xy(1, 5, CellConstraints.RIGHT, CellConstraints.CENTER));
         header.add(banner, CC.xy(3, 1));
      }
      jPanel_Contents.setBorder(createPanelBorder());
      jPanel_Main.add(jPanel_Contents, CC.xy(2, 4));
      header.setBorder(createPanelBorder());
      jPanel_Main.add(header, CC.xy(2, 2));
      {
         final JPanel msgPanel = new JPanel(new FormLayout("5dlu, right:40dlu, 5dlu, 260dlu, 5dlu, pref, 5dlu", "pref"));
         msgPanel.setBorder(createPanelBorder());
         msgPanel.add(new JLabel("Message:"), CC.xy(2, 1));
         msgPanel.add(jLabel_ErrorLabel, CC.xy(4, 1));
         msgPanel.add(jButton_ErrorBtn, CC.xy(6, 1));
         jPanel_Main.add(msgPanel, CC.xy(2, 6));
      }
   }

   private CompoundBorder createPanelBorder() {
      final Border line = SwingUtils.createDefaultBorder();
      final Border empty = BorderFactory.createEmptyBorder(8, 8, 8, 8);
      return BorderFactory.createCompoundBorder(line, empty);
   }

   /**
    * Gets an instance of the data items used to initialize the configuration of
    * this wizard dialog.
    *
    * @return An instance of the templated class
    */
   public W getInitial() {
      return mInitial;
   }

   /**
    * Gets an instance of the data items used to initialize the configuration of
    * this wizard dialog.
    *
    * @return An instance of the templated class W or null
    */
   public W getResult() {
      return mResult;
   }
}
