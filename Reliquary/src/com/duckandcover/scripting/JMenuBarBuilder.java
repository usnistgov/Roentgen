package com.duckandcover.scripting;

import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.swing.JMenuBar;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

/**
 * <p>
 * A simple mechanism to build a modifiable JMenuBar.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class JMenuBarBuilder
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

   private final SortedMap<Indexer, JMenuBuilder> jMenu_Top;

   /**
    * Constructs a JMenuBarBuilder
    */
   public JMenuBarBuilder() {
      jMenu_Top = new TreeMap<>();
   }

   public void setEnabled(final String name, final boolean value) {
      for(final JMenuBuilder jmb : jMenu_Top.values())
         if(jmb.getName().equals(name))
            jmb.setEnabled(value);
         else
            jmb.setEnabled(name, value);
   }

   private JMenuBuilder find(final String name) {
      for(final Map.Entry<Indexer, JMenuBuilder> me : jMenu_Top.entrySet())
         if(me.getKey().mName.equals(name))
            return me.getValue();
      return null;
   }

   public void reprioritize(final int newPriority, final String name) {
      for(final Map.Entry<Indexer, JMenuBuilder> me : jMenu_Top.entrySet())
         if(me.getKey().mName.equals(name)) {
            final JMenuBuilder mb = me.getValue();
            jMenu_Top.remove(me.getKey());
            jMenu_Top.put(new Indexer(newPriority, name), mb);
            return;
         }
      jMenu_Top.put(new Indexer(newPriority, name), new JMenuBuilder(name));
   }

   public JMenuBuilder findOrAdd(final String name, final int priority) {
      JMenuBuilder res = find(name);
      if(res == null) {
         res = new JMenuBuilder(name);
         jMenu_Top.put(new Indexer(priority, name), res);
      }
      return res;
   }

   public void setEnabled(final boolean value) {
      for(final JMenuBuilder jmb : jMenu_Top.values())
         jmb.setEnabled(value);
   }

   public JMenuBar build() {
      final JMenuBar res = new JMenuBar();
      for(final JMenuBuilder mb : jMenu_Top.values())
         res.add(mb.build());
      return res;
   }

   /**
    * @param mode
    * @return
    * @see com.duckandcover.html.IToHTML#toHTML(com.duckandcover.html.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      final Table table = new Table();
      final JMenuBar res = new JMenuBar();
      for(final JMenuBuilder mb : jMenu_Top.values())
         table.addRow(new Table.TD(mb.getName(), 1, 1, false), new Table.TD(mb.toHTML(mode), 1, 1, false));
      return table.toHTML(mode);
   }

}
