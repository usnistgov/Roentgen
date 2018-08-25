package gov.nist.microanalysis.roentgen.math;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * Static functions for performing math on either Interval or UncertainValue
 * numbers as appropriate for each. Mixing the types converts Interval objects
 * to UncertainValue objects for the calculation.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 199 $
 */
public class MathUtilities {

   /**
    * Hidden
    */
   private MathUtilities() {
   }

   /**
    * Parses "number" to a Double, "(number&plusmn;number)" / "number+-number"
    * to an UncertainValue and "[number,number]" / "[number-number]" to a
    * Interval.
    *
    * @param number
    * @return Number
    */
   public static Number parseNumber(final String number, final String dname)
         throws NumberFormatException {
      String tmp = number.trim();
      if(tmp.contains("+-") || tmp.contains("\u00B1")) {
         if(tmp.startsWith("(") && tmp.endsWith(")"))
            tmp = tmp.substring(1, tmp.length() - 1);
         // The +- symbol or "+-"
         final String[] terms = tmp.split("\\xB1|\\+\\-");
         assert terms.length > 1;
         return new UncertainValue(Double.parseDouble(terms[0].trim()), "dname", Double.parseDouble(terms[1].trim()));
      } else
         return new Double(Double.parseDouble(tmp));
   }

   public static String toHTML(final RealMatrix rm, final NumberFormat nf) {
      return toHTML(rm, null, null, nf);
   }

   public static String toHTML(final RealMatrix rm, final List<Object> rowLabel, final List<Object> colLabel, final NumberFormat nf) {
      final Table t = new Table();
      if(colLabel != null) {
         final List<Item> header = new ArrayList<>();
         header.add(Table.td());
         for(int c = 0; c < rm.getColumnDimension(); ++c)
            header.add(Table.th(HTML.toHTML(colLabel.get(c),Mode.TERSE)));
         t.addRow(header);
      }
      for(int r = 0; r < rm.getRowDimension(); ++r) {
         final List<Item> items = new ArrayList<>();
         if(rowLabel != null)
            items.add(Table.th(HTML.toHTML(rowLabel.get(r),Mode.TERSE)));
         for(int c = 0; c < rm.getColumnDimension(); ++c)
            items.add(Table.td(nf.format(rm.getEntry(r, c))));
         t.addRow(items);
      }
      return t.toHTML(Mode.NORMAL);
   }

   public static String toHTML_Vertical(final RealVector rv, final List<Object> label, final NumberFormat nf) {
      final Table t = new Table();
      for(int i = 0; i < rv.getDimension(); ++i) {
         final List<Item> items = new ArrayList<>();
         if(label != null)
            items.add(Table.td(HTML.toHTML(label.get(i),Mode.TERSE)));
         items.add(Table.td(nf.format(rv.getEntry(i))));
         t.addRow(items);
      }
      return t.toHTML(Mode.NORMAL);
   }

   public static String toHTML_Horizontal(final RealVector rv, final List<Object> label, final NumberFormat nf) {
      final Table t = new Table();
      final List<Item> hdr = new ArrayList<>();
      final List<Item> data = new ArrayList<>();
      for(int i = 0; i < rv.getDimension(); ++i) {
         hdr.add(Table.td(HTML.toHTML(label.get(i),Mode.TERSE)));
         data.add(Table.td(nf.format(rv.getEntry(i))));
      }
      t.addRow(hdr);
      t.addRow(data);
      return t.toHTML(Mode.NORMAL);
   }

   public static String toHTML_Terse(final UncertainValue uv, final NumberFormat nf) {
      return nf.format(uv.doubleValue());
   }

   public static String toHTML_Normal(final UncertainValue uv, final NumberFormat nf) {
      return nf.format(uv.doubleValue()) + "&pm;" + nf.format(uv.uncertainty());
   }

   public static String toHTML_Verbose(final UncertainValue uv, final NumberFormat nf) {
      final StringBuffer sb = new StringBuffer();
      sb.append(nf.format(uv.doubleValue()));
      sb.append("&pm;");
      boolean first = true;
      for(final Object obj : uv.getComponentNames()) {
         if(!first) {
            sb.append(",");
            first = false;
         }
         sb.append("(");
         sb.append(obj.toString());
         sb.append(":");
         sb.append(nf.format(uv.getComponent(obj)));
         sb.append(")");
      }
      return sb.toString();
   }

   public static String toHTML_Image(final RealMatrix rm, final String caption, final File base, final String dir)
         throws IOException {
      double min = 1.0e300, max = -1.0e-300;
      for(int r = 0; r < rm.getRowDimension(); ++r)
         for(int c = 0; c < rm.getColumnDimension(); ++c) {
            final double val = rm.getEntry(r, c);
            min = Math.min(min, val);
            max = Math.max(max, val);
         }
      final IDoubleAsColor colorize = new PositiveNegativeAsColor(min, Color.red, max, Color.green, Color.yellow, 0.5);
      final BufferedImage img = RealMatrixAsBitmap(rm, colorize);
      return HTML.image(img, caption, base, dir);
   }

   public static double min(final RealMatrix rm) {
      double min = Double.MAX_VALUE;
      for(int r = 0; r < rm.getRowDimension(); ++r)
         for(int c = 0; c < rm.getColumnDimension(); ++c)
            if(rm.getEntry(r, c) < min)
               min = rm.getEntry(r, c);
      return min;
   }

   public static double max(final RealMatrix rm) {
      double max = -Double.MAX_VALUE;
      for(int r = 0; r < rm.getRowDimension(); ++r)
         for(int c = 0; c < rm.getColumnDimension(); ++c)
            if(rm.getEntry(r, c) > max)
               max = rm.getEntry(r, c);
      return max;
   }

   public interface IDoubleAsColor {

      /**
       * Converts a double value to a Color
       *
       * @param v
       * @return Color
       */
      public Color compute(double v);
   }

   public static class PositiveNegativeAsColor
      implements
      IDoubleAsColor {
      final double mRange;
      final int[] mMinRgb;
      final int[] mMaxRgb;
      final double mGamma;
      final Color mOutOfRange;
      final Color mTransparent = new Color(0,0,0,0xFF);

      public PositiveNegativeAsColor(final double min, final Color minColor, final double max, final Color maxColor, final Color outOfRange, final double gamma) {
         mRange = Math.max(Math.abs(min), Math.abs(max));
         mMinRgb = new int[] {
            minColor.getRed(),
            minColor.getGreen(),
            minColor.getBlue()
         };
         mMaxRgb = new int[] {
            maxColor.getRed(),
            maxColor.getGreen(),
            maxColor.getBlue()
         };
         mGamma = gamma;
         mOutOfRange = outOfRange;
      }

      /**
       * A linear mapping from one color to the other.
       *
       * @see gov.nist.microanalysis.roentgen.math.MathUtilities.IDoubleAsColor#compute(double)
       */
      @Override
      public Color compute(final double val) {
         if((val >= -mRange) && (val <= mRange) && (!Double.isNaN(val))) {
            final double sc = Math.pow(Utility.bound(Math.abs(val) / mRange, 0.0, 1.0), mGamma);
            if(val==0.0)
            	return Color.white;
            else if(val < 0.0) {
               return new Color((int) (mMinRgb[0] * sc), (int) (mMinRgb[1] * sc), (int) (mMinRgb[2] * sc));
            } else {
               return new Color((int) (mMaxRgb[0] * sc), (int) (mMaxRgb[1] * sc), (int) (mMaxRgb[2] * sc));
            }
         } else
            return mOutOfRange;
      }
   }

   public static class LinearDoubleAsColor
      implements
      IDoubleAsColor {
      final int[] mMinRgb;
      final int[] mMaxRgb;
      final double mMin;
      final double mMax;
      final double mGamma;
      final Color mOutOfRange;

      public LinearDoubleAsColor(final double min, final Color minColor, final double max, final Color maxColor, final Color outOfRange, final double gamma) {
         mMin = min;
         mMinRgb = new int[] {
            minColor.getRed(),
            minColor.getGreen(),
            minColor.getBlue()
         };
         mMax = max;
         mMaxRgb = new int[] {
            maxColor.getRed(),
            maxColor.getGreen(),
            maxColor.getBlue()
         };
         mGamma = gamma;
         mOutOfRange = outOfRange;
      }

      public static LinearDoubleAsColor blackAndWhite(final double min, final double max, final double gamma) {
         return new LinearDoubleAsColor(min, Color.BLACK, max, Color.WHITE, Color.YELLOW, gamma);
      }

      /**
       * A linear mapping from one color to the other.
       *
       * @see gov.nist.microanalysis.roentgen.math.MathUtilities.IDoubleAsColor#compute(double)
       */
      @Override
      public Color compute(final double val) {
         if((val >= mMin) && (val <= mMax) && (!Double.isNaN(val))) {
            final double sc = 255.0 * Math.pow(Utility.bound((val - mMin) / (mMax - mMin), 0.0, 1.0), mGamma);
            return new Color( //
                  mMinRgb[0] + (int) ((mMaxRgb[0] - mMinRgb[0]) * sc), //
                  mMinRgb[1] + (int) ((mMaxRgb[1] - mMinRgb[1]) * sc), //
                  mMinRgb[2] + (int) ((mMaxRgb[2] - mMinRgb[2]) * sc) //
            );
         } else
            return mOutOfRange;
      }
   }

   public static BufferedImage RealMatrixAsBitmap(final RealMatrix rm, final IDoubleAsColor colorize) {
      final BufferedImage res = new BufferedImage(rm.getColumnDimension(), rm.getRowDimension(), BufferedImage.TYPE_3BYTE_BGR);
      for(int r = 0; r < rm.getRowDimension(); ++r)
         for(int c = 0; c < rm.getColumnDimension(); ++c)
            res.setRGB(c, r, colorize.compute(rm.getEntry(r, c)).getRGB());
      return res;
   }

   public static String getKurtosisDescription(final double kurtosis) {
      if(Math.abs(kurtosis) < 0.3)
         return "Similar to normal";
      else if(Math.abs(kurtosis) < 1000.0)
         return Math.abs(kurtosis) > 0.0 ? "Tail heavy" : "Peaky";
      else
         return Math.abs(kurtosis) > 0.0 ? "Very tail heavy" : "Very peaky";

   }

   public static String getSkewnessDescription(final double skewness) {
      if(Math.abs(skewness) < 0.1)
         return "Highly symmetric";
      else if(Math.abs(skewness) < 1.0)
         return skewness > 0.0 ? "Moderately skew positive" : "Moderately skew negative";
      else if(Math.abs(skewness) < 2.0)
         return skewness > 0.0 ? "Skew positive" : "Skew negative";
      else
         return skewness > 0.0 ? "Highly skew positive" : "Highly skew negative";
   }

   public static IToHTML toHTML(final DescriptiveStatistics desc, final NumberFormat bnf) {

      class DescStatTable
         implements
         IToHTML {

         final DescriptiveStatistics mDescStat;
         final NumberFormat mFormat;

         DescStatTable(final DescriptiveStatistics d, final NumberFormat b) {
            mDescStat = d;
            mFormat = b;
         }

         @Override
         public String toHTML(final Mode mode) {
            final List<List<String>> items = new ArrayList<>();
            final List<String> colHeaders = Arrays.asList("Statistic", "Value", "Note");
            items.add(Arrays.asList("Mean (&mu;)", mFormat.format(mDescStat.getMean()), ""));
            items.add(Arrays.asList("Standard Deviation (&sigma;)", mFormat.format(mDescStat.getStandardDeviation()), mFormat.format((100.0
                  * mDescStat.getStandardDeviation()) / mDescStat.getMean()) + " %"));
            if((mode == Mode.NORMAL) || (mode == Mode.VERBOSE)) {
               final double skewness = mDescStat.getSkewness();
               items.add(Arrays.asList("Skewness", mFormat.format(skewness), getSkewnessDescription(skewness)));
               final double kurtosis = mDescStat.getKurtosis();
               items.add(Arrays.asList("Kurtosis", mFormat.format(kurtosis), getKurtosisDescription(kurtosis)));
               if(mode == Mode.VERBOSE) {
                  final int[] pcts = new int[] {
                     1,
                     5,
                     10,
                     25,
                     50,
                     75,
                     90,
                     95,
                     99
                  };
                  for(final int pct : pcts)
                     items.add(Arrays.asList(HTML.asOrdinalHTML(pct)
                           + " percentile", mFormat.format(mDescStat.getPercentile(pct)), ""));
               }
            }
            return HTML.asTable(colHeaders, null, items);
         }

      }
      return new DescStatTable(desc, bnf);
   }
   
   public static final Table.Item td(Number n, BasicNumberFormat bnf){
	   return Table.td(bnf.formatHTML(n));
   }
   
   public static final Table.Item td(double n, BasicNumberFormat bnf){
	   return Table.td(bnf.formatHTML(n));
   }


}
