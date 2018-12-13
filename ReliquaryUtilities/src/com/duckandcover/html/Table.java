package com.duckandcover.html;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;

/**
 * <p>
 * A more sophisticated method for creating tables that may include header and
 * table items anywhere and the items may also have "rowspan" and "colspan" that
 * are non-unity. The contents of the table can vary depending upon the mode
 * argument to toHTML(...)
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2016
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
public class Table
   implements
   IToHTML,
   IToHTMLExt {

   private class Index
      implements
      Comparable<Index> {
      private final int mRow;
      private final int mCol;

      private Index(final int row, final int col) {
         assert row >= 0;
         assert col >= 0;
         mRow = row;
         mCol = col;
      }

      @Override
      public int compareTo(final Index o) {
         final int res = Integer.compare(mRow, o.mRow);
         return res == 0 ? Integer.compare(mCol, o.mCol) : res;
      }

   }

   private final TreeMap<Index, Item> mItems;
   private final TreeMap<Index, Mode> mMode;
   private boolean mMatrix;

   public static abstract class Item
      implements
      IToHTMLExt {
      protected final IToHTML mHTML;
      protected final int mRowSpan;
      protected final int mColSpan;
      protected final boolean mCenter;

      protected Item(final IToHTML html, final int rowSpan, final int colSpan, final boolean center) {
         mHTML = html;
         mRowSpan = rowSpan;
         mColSpan = colSpan;
         mCenter = center;
      }

      protected String span() {
         final StringBuffer sb = new StringBuffer();
         if(mColSpan > 1)
            sb.append(" colspan=\"" + mColSpan + "\"");
         if(mRowSpan > 1)
            sb.append(" rowspan=\"" + mRowSpan + "\"");
         return sb.toString();
      }

      protected String center() {
         return mCenter ? " align = \"center\"" : "";
      }
   }

   public static class TD
      extends
      Item {
      public TD(final IToHTML html, final int rowSpan, final int colSpan, final boolean center) {
         super(html, rowSpan, colSpan, center);
      }

      public TD(final String html, final int rowSpan, final int colSpan, final boolean center) {
         super(Transforms.createHTML(html), rowSpan, colSpan, center);
      }

      @Override
      public String toHTML(final Mode mode, final File base, final String dir)
            throws IOException {
    	  if(mHTML instanceof IToHTMLExt)
    		  return "<td" + span() + center() + ">" + ((IToHTMLExt)mHTML).toHTML(mode, base, dir) + "</td>";
    	  else
    		  return "<td" + span() + center() + ">" + mHTML.toHTML(mode) + "</td>";
      }

      @Override
      public String toHTML(final Mode mode) {
         return "<td" + span() + center() + ">" + mHTML.toHTML(mode) + "</td>";
      }
   }

   public static class TH
      extends
      Item {
      public TH(final IToHTML html, final int rowSpan, final int colSpan, final boolean center) {
         super(html, rowSpan, colSpan, center);
      }

      public TH(final String html, final int rowSpan, final int colSpan, final boolean center) {
         super(Transforms.createHTML(html), rowSpan, colSpan, center);
      }

      @Override
      public String toHTML(final Mode mode, final File base, final String dir)
            throws IOException {
    	  if(mHTML instanceof IToHTMLExt)
    		  return "<th" + span() + center() + ">" + ((IToHTMLExt)mHTML).toHTML(mode, base, dir) + "</th>";
    	  else
    		  return "<th" + span() + center() + ">" + mHTML.toHTML(mode) + "</th>";
      }

      @Override
      public String toHTML(final Mode mode) {
         return "<th" + span() + center() + ">" + mHTML.toHTML(mode) + "</th>";
      }
   }

   public static class NullItem
      extends
      Item {
      NullItem() {
         super(null, 1, 1, false);
      }

      @Override
      public String toHTML(final Mode mode, final File base, final String dir) {
         return "";
      }

      @Override
      public String toHTML(final Mode mode) {
         return "";
      }
   }

   /**
    * Constructs a HTMLTable
    */
   public Table() {
      mItems = new TreeMap<>();
      mMode = new TreeMap<>();
      mMatrix = false;
   }

   public Table asMatrix() {
      mMatrix = true;
      return this;
   }

   public void add(final int row, final int col, final Item item) {
      final Index index = new Index(row, col);
      mItems.put(index, item);
      mMode.remove(index);
   }

   public void add(final int row, final int col, final Item item, final Mode mode) {
      final Index index = new Index(row, col);
      mItems.put(index, item);
      mMode.put(index, mode);
   }

   public int rowCount() {
      int maxRow = -1;
      for(final Index idx : mItems.keySet())
         maxRow = Math.max(maxRow, idx.mRow);
      return maxRow;
   }

   public int colCount() {
      int maxCol = -1;
      for(final Index idx : mItems.keySet())
         maxCol = Math.max(maxCol, idx.mCol);
      return maxCol;
   }

   public void addRow(final List<Item> row) {
      final int lastRow = rowCount();
      int col = 0;
      for(final Item item : row) {
         add(lastRow + 1, col, item);
         col += item.mColSpan;
      }
   }

   public void addRow(final Item... items) {
      addRow(Arrays.asList(items));
   }

   public void addCol(final List<Item> col) {
      final int lastCol = rowCount();
      int row = 0;
      for(final Item item : col) {
         add(row, lastCol + 1, item);
         row += item.mRowSpan;
      }
   }

   public void addCol(final Item... items) {
      addCol(Arrays.asList(items));
   }

   public static Item th(final String html, final int rowSpan, final int colSpan) {
      return new TH(html, rowSpan, colSpan, false);
   }

   public static Item thc(final String html, final int rowSpan, final int colSpan) {
      return new TH(html, rowSpan, colSpan, true);
   }

   public static Item th(final String html, final int colSpan) {
      return new TH(html, 1, colSpan, false);
   }

   public static Item thc(final String html, final int colSpan) {
      return new TH(html, 1, colSpan, true);
   }

   public static Item th(final String html) {
      return new TH(html, 1, 1, false);
   }

   public static Item thc(final String html) {
      return new TH(html, 1, 1, true);
   }

   public static Item th(final IToHTML html, final int rowSpan, final int colSpan) {
      return new TH(html, rowSpan, colSpan, false);
   }

   public static Item thc(final IToHTML html, final int rowSpan, final int colSpan) {
      return new TH(html, rowSpan, colSpan, true);
   }

   public static Item th(final IToHTML html, final int colSpan) {
      return new TH(html, 1, colSpan, false);
   }

   public static Item thc(final IToHTML html, final int colSpan) {
      return new TH(html, 1, colSpan, true);
   }

   public static Item th(final IToHTML html) {
      return new TH(html, 1, 1, false);
   }

   public static Item thc(final IToHTML html) {
      return new TH(html, 1, 1, true);
   }

   public static Item td(final String html, final int rowSpan, final int colSpan) {
      return new TD(html, rowSpan, colSpan, false);
   }

   public static Item tdc(final String html, final int rowSpan, final int colSpan) {
      return new TD(html, rowSpan, colSpan, true);
   }

   public static Item td(final String html, final int colSpan) {
      return new TD(html, 1, colSpan, false);
   }
   
   public static Item tdc(final String html, final int colSpan) {
      return new TD(html, 1, colSpan, true);
   }

   public static Item td(final String html) {
      return new TD(html, 1, 1, false);
   }

   public static Item tdc(final String html) {
      return new TD(html, 1, 1, true);
   }

   public static Item td(final Number number) {
      return new TD(number.toString(), 1, 1, false);
   }

   public static Item tdc(final Number number) {
      return new TD(number.toString(), 1, 1, true);
   }

   public static Item td(final double number) {
      return new TD(Double.toString(number), 1, 1, false);
   }

   public static Item tdc(final double number) {
      return new TD(Double.toString(number), 1, 1, true);
   }

   public static Item td() {
      return new TD("&nbsp;", 1, 1, false);
   }

   public static Item td(final IToHTML html, final int rowSpan, final int colSpan) {
      return new TD(html, rowSpan, colSpan, false);
   }

   public static Item tdc(final IToHTML html, final int rowSpan, final int colSpan) {
      return new TD(html, rowSpan, colSpan, true);
   }

   public static Item td(final IToHTML html, final int colSpan) {
      return new TD(html, 1, colSpan, false);
   }

   public static Item tdc(final IToHTML html, final int colSpan) {
      return new TD(html, 1, colSpan, true);
   }

   public static Item td(final IToHTML html) {
      return new TD(html, 1, 1, false);
   }

   public static Item tdc(final IToHTML html) {
      return new TD(html, 1, 1, true);
   }
   
   public static Item td(Object key) {
		return td(HTML.toHTML(key, Mode.NORMAL));
	}


   /**
    * This function implements rowspan and colspan.
    *
    * @return Item[][]
    */
   private Item[][] buildItems() {
      final int maxRow = rowCount() + 1, maxCol = colCount() + 1;
      final Item[][] items = new Item[maxRow][maxCol];
      final NullItem ni = new NullItem();
      final Item empty = new TD("&nbsp;", 1, 1, false);
      for(int ri = 0; ri < maxRow; ri++)
         for(int ci = 0; ci < maxCol; ci++)
            if(items[ri][ci] == null) {
               final Index idx = new Index(ri, ci);
               final Item item = mItems.getOrDefault(idx, empty);
               items[ri][ci] = item;
               for(int rr = ri; rr < (ri + item.mRowSpan); rr++)
                  if(rr != ri)
                     items[rr][ci] = ni;
               for(int cc = ci; cc < (ci + item.mColSpan); cc++)
                  if(cc != ci)
                     items[ri][cc] = ni;
            }
      return items;
   }

   /**
    * @see com.duckandcover.roentgen.html.IToHTMLExt#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode,
    *      java.io.File, java.lang.String)
    */
   @Override
   public String toHTML(final Mode mode, final File base, final String dir)
         throws IOException {
      final Item[][] items = buildItems();
      final StringBuffer sb = new StringBuffer();
      if(mMatrix)
         sb.append("<table class = \"matrix\">");
      else
         sb.append("<table>");
      int r = 0;
      for(final Item[] row : items) {
         sb.append("<tr>");
         int c = 0;
         for(final Item item : row) {
            final Index idx = new Index(r, c);
            if(item != null)
               sb.append(item.toHTML(mMode.getOrDefault(idx, mode), base, dir));
            ++c;
         }
         sb.append("</tr>");
         ++r;
      }
      sb.append("</table>");
      sb.append("");
      return sb.toString();
   }

   /**
    * @see com.duckandcover.roentgen.html.IToHTML#toHTML(com.duckandcover.roentgen.html.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      final Item[][] items = buildItems();
      final StringBuffer sb = new StringBuffer();
      if(mMatrix)
         sb.append("<table class = \"matrix\">");
      else
         sb.append("<table>");
      int r = 0;
      for(final Item[] row : items) {
         sb.append("<tr>");
         int c = 0;
         for(final Item item : row) {
            final Index idx = new Index(r, c);
            sb.append(item.toHTML(mMode.getOrDefault(idx, mode)));
            ++c;
         }
         sb.append("</tr>");
         ++r;
      }
      sb.append("</table>");
      return sb.toString();
   }
}
