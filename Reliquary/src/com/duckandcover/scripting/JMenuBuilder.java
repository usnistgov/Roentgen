package com.duckandcover.scripting;

import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.JMenuItem;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

public class JMenuBuilder
   implements
   IToHTML {

   static private class Indexer
      implements
      Comparable<Indexer> {
      private final int mPriority;
      final String mName;

      private Indexer(final int p, final String name) {
         mPriority = p;
         mName = name;
      }

      /**
       * @see java.lang.Comparable#compareTo(java.lang.Object)
       */
      @Override
      public int compareTo(final Indexer o) {
         final int res = Integer.compare(mPriority, o.mPriority);
         return res == 0 ? mName.compareTo(o.mName) : res;
      }

   }

   private final String mName;
   private final SortedMap<Indexer, JMenuBuilder> jMenuBuilder_SubMenus;
   private final SortedMap<Indexer, JMenuItem> jMenuItems;

   JMenuBuilder(final String name) {
      mName = name;
      jMenuBuilder_SubMenus = new TreeMap<>();
      jMenuItems = new TreeMap<>();
   }

   public void add(final int priority, final Action action) {
      jMenuItems.put(new Indexer(priority, action.toString()), new JMenuItem(action));
   }

   public void reprioritize(final int newPriority, final String name) {
      for(final Map.Entry<Indexer, JMenuBuilder> me : jMenuBuilder_SubMenus.entrySet())
         if(me.getKey().mName.equals(name)) {
            final JMenuBuilder mb = me.getValue();
            jMenuBuilder_SubMenus.remove(me.getKey());
            jMenuBuilder_SubMenus.put(new Indexer(newPriority, name), mb);
            return;
         }
      jMenuBuilder_SubMenus.put(new Indexer(newPriority, name), new JMenuBuilder(name));
   }

   private JMenuBuilder find(final String name) {
      for(final Map.Entry<Indexer, JMenuBuilder> me : jMenuBuilder_SubMenus.entrySet())
         if(me.getKey().mName.equals(name))
            return me.getValue();
      return null;
   }

   /**
    * Returns the JMenuBuilder associated with the sub-menu with the specified
    * name or creates one with the specified name and priority if one does not
    * exist. The priority of an existing menu will not be changed.
    *
    * @param priority
    * @param name
    * @return JMenuBuilder
    */
   public JMenuBuilder findOrAddSubMenu(final int priority, final String name) {
      JMenuBuilder res = find(name);
      if(name == null) {
         res = new JMenuBuilder(name);
         jMenuBuilder_SubMenus.put(new Indexer(priority, name), res);
      }
      return res;
   }

   public JMenu build() {
      final JMenu res = new JMenu(mName);
      for(final JMenuBuilder mb : jMenuBuilder_SubMenus.values())
         res.add(mb.build());
      for(final JMenuItem mi : jMenuItems.values())
         res.add(mi);
      return res;
   }

   public String getName() {
      return mName;
   }

   public void setEnabled(final boolean value) {
      for(final JMenuBuilder jmb : jMenuBuilder_SubMenus.values())
         jmb.setEnabled(value);
      for(final JMenuItem mi : jMenuItems.values())
         mi.setEnabled(value);
   }

   public void setEnabled(final String name, final boolean value) {
      for(final JMenuBuilder jmb : jMenuBuilder_SubMenus.values())
         if(jmb.getName().equals(name))
            jmb.setEnabled(value);
         else
            jmb.setEnabled(name, value);
      for(final JMenuItem mi : jMenuItems.values())
         mi.setEnabled(value);
   }

   /**
    * @param mode
    * @return String as HTML
    * @see com.duckandcover.html.IToHTML#toHTML(com.duckandcover.html.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      final Table table = new Table();
      for(final JMenuBuilder jmb : jMenuBuilder_SubMenus.values())
         table.addRow(Table.td(jmb.mName), Table.td(jmb.toHTML(mode)));
      for(final JMenuItem mi : jMenuItems.values())
         table.addRow(Table.td("Menu item"), Table.td(mi.getText()));
      return table.toHTML(mode);
   }

}