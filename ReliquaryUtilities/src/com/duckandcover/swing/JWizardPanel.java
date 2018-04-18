package com.duckandcover.swing;

import java.awt.LayoutManager;

import javax.swing.JPanel;

import com.thoughtworks.xstream.XStream;

public class JWizardPanel<W>
   extends
   JPanel {

   /**
    *
    */
   private static final long serialVersionUID = 6629485617850376788L;
   private final JWizardDialog<W> mWizard;
   /**
    * The data to initialize the panel in onShow().
    */
   protected W mCurrent;

   /**
    * Constructs a JWizardPanel associated with the specified WizardDialog
    *
    * @param wiz
    */
   public JWizardPanel(final JWizardDialog<W> wiz) {
      super();
      assert (wiz != null);
      mWizard = wiz;
   }

   /**
    * Constructs a JWizardPanel associated with the specified WizardDialog and
    * using the specified LayoutManager
    *
    * @param lo
    * @param wiz
    */
   public JWizardPanel(final JWizardDialog<W> wiz, final LayoutManager lo) {
      super(lo);
      assert (wiz != null);
      mWizard = wiz;
   }

   /**
    * Returns the instance of WizardDialog with which this JWizardPanel is
    * associated.
    *
    * @return JWizardDialog
    */
   public JWizardDialog<W> getWizardDialog() {
      assert (mWizard != null);
      return mWizard;
   }

   /**
    * Creates an internal clone of the data which is used to initialize the
    * panel in onShow().
    *
    * @param data
    */
   @SuppressWarnings("unchecked")
   public void initialize(final W data) {
      // Creates a clone...
      final XStream xs = new XStream();
      final String xml = xs.toXML(data);
      mCurrent = (W) xs.fromXML(xml);
   }

   /**
    * Implement the onShow to perform some action when this panel is displayed
    * in the WizardDialog. Should take the data in mCurrent and initialize the
    * panel contents.
    */
   public void onShow() {
      // Don't do anything by default
   }

   /**
    * Implement the onHide to perform some action when this panel is hidden by
    * the WizardDialog. Should take the data from the panel and update mCurrent
    * to reflect the updated data
    */
   public void onHide() {
      // Don't do anything by default
   }

   /**
    * Called when the user presses the next or finish button. If the panel does
    * not want to loose focus then this function should be overriden to return
    * false. This mechanism allows the panel to veto the next/finished button to
    * enforce data verification.
    *
    * @return true to permit moving to the next panel or false otherwise.
    */
   public boolean permitNext() {
      return true;
   }

   /**
    * Returns an updated copy of the
    *
    * @return
    */
   public W getData() {
      return mCurrent;
   }
}
