package com.duckandcover.swing;

import java.util.ArrayList;

import javax.swing.JPanel;

public abstract class JPreferencePanel
   extends
   JPanel {

   private static final long serialVersionUID = -2201836961819522100L;
   private final String mDescription;
   private final String mLabel;
   private final ArrayList<JPreferencePanel> mChildren = new ArrayList<>();

   public JPreferencePanel(final String label, final String desc) {
      mLabel = label;
      mDescription = desc;
   }

   public void commitAll() {
      commit();
      for(final JPreferencePanel pp : mChildren)
         pp.commitAll();
   }

   /**
    * Override to apply and commit all modifications to the preferences on this
    * page.
    */
   public abstract void commit();

   /**
    * Initialize the GUI fields with preference data
    */
   public abstract void initialize();

   /**
    * Is the data on this page valid? Should the PreferenceDialog allow the user
    * to switch to another page?
    */
   public boolean permitExit() {
      return true;
   }

   /**
    * Returns any child preference panels associated with this preference panel.
    * The child panels may too have child panels.
    *
    * @return PreferencePanel[]
    */
   public JPreferencePanel[] getChildren() {
      final JPreferencePanel[] res = new JPreferencePanel[mChildren.size()];
      return mChildren.toArray(res);
   }

   /**
    * Add a child PreferencePanel to this panel. Eventually this mechanism will
    * be made dynamic to allow children to be added while the user modifies the
    * PreferenceDialog.
    *
    * @param child
    */
   public void addChildPanel(final JPreferencePanel child) {
      mChildren.add(child);
   }

   /**
    * Gets the current value assigned to description
    *
    * @return Returns the description.
    */
   public String getDescription() {
      return mDescription;
   }

   /**
    * Gets the current value assigned to name
    *
    * @return Returns the name.
    */
   @Override
   public String toString() {
      return mLabel;
   }
}