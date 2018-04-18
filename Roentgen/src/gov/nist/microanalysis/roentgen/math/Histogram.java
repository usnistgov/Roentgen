package gov.nist.microanalysis.roentgen.math;

import java.text.MessageFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * A class for creating histograms of series of data points. It is possible to
 * create histograms of equally spaced bins, ratio spaced bins or arbitrarily
 * spaced bins.
 * </p>
 * <p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class Histogram
   implements
   IToHTML {
   private double[] mBinMin;
   private int[] mCounts;

   public class BinName
      implements
      Comparable<BinName> {
      private final int mBin;
      private final String mString;

      private BinName(final int bin, final String format) {
         assert (bin >= -1);
         assert (bin < (binCount() + 1));
         mBin = bin;
         mString = "[" + (mBin < 0 ? "-inf"
               : MessageFormat.format(format, new Object[] {
                  new Double(minValue(mBin))
               })) + "-" + (mBin >= binCount() ? "inf"
                     : MessageFormat.format(format, new Object[] {
                        new Double(maxValue(mBin))
                     })) + ")";
      }

      private BinName(final int bin, final NumberFormat nf) {
         assert (bin >= -1);
         assert (bin < (binCount() + 1));
         mBin = bin;
         mString = "[" + (mBin < 0 ? "-inf" : nf.format(minValue(mBin))) + "-"
               + (mBin >= binCount() ? "inf" : nf.format(maxValue(mBin))) + ")";

      }

      @Override
      public String toString() {
         return mString;
      }

      @Override
      public int compareTo(final BinName o) {
         return Integer.compare(mBin, o.mBin);
      }
   }

   /**
    * <p>
    * Histogram - Creates a histogram representing the specified range with bins
    * for over-range and under-range.
    * </p>
    * <p>
    * Bins are numbered such that 0..(nBins-1) correspond to the full range.
    * Bin[0] -&gt; [min,min+delta) where delta = (max-min)/nBins, Bin[nBins-1]
    * -&gt; [max-delta,max). There are two bins for under and overrange values.
    * Bin[-1]-&gt;[-inf,min), Bin[nBins] -&gt; [max,inf) where inf is infinity.
    * </p>
    *
    * @param min double
    * @param max double
    * @param nBins int - The number of bins (excluding over- and under-range)
    */
   public Histogram(final double min, final double max, final int nBins) {
      Preconditions.checkArgument(max > min);
      Preconditions.checkArgument(nBins > 0);
      mCounts = new int[nBins + 2]; // nBins + under + over
      mBinMin = new double[nBins + 1];
      final double delta = (max - min) / nBins;
      mBinMin[0] = min;
      for(int i = 1; i < mBinMin.length; ++i)
         mBinMin[i] = min + (i * delta);
      Arrays.sort(mBinMin);
   }

   /**
    * Constructs a Histogram object to represent bins between min (inclusive)
    * and max (excluded)
    *
    * @param min
    * @param max
    * @param ratio
    * @throws EPQException
    */
   public Histogram(double min, final double max, final double ratio) {
      Preconditions.checkArgument(ratio > 1.0);
      Preconditions.checkArgument(max > min);
      final int nBins = (int) (Math.log(max / min) / Math.log(ratio));
      mCounts = new int[nBins + 2];
      mBinMin = new double[nBins + 1];
      for(int i = 0; i < mBinMin.length; ++i, min *= ratio)
         mBinMin[i] = min;
   }

   /**
    * Constructs a Histogram object to bin into the arbitrary bins defined by
    * lower ends define by binMins. The max end of the top bin is defined by
    * max.
    *
    * @param binMins
    * @param max
    */
   public Histogram(final double[] binMins, final double max) {
      mBinMin = new double[binMins.length + 1];
      for(int i = 0; i < binMins.length; ++i)
         mBinMin[i] = binMins[i];
      mBinMin[binMins.length] = max;
      Arrays.sort(mBinMin);
      Preconditions.checkArgument(max > mBinMin[mBinMin.length - 1]);
      mCounts = new int[binMins.length + 2];
   }

   public Histogram(final Histogram hist) {
      mBinMin = new double[hist.mBinMin.length];
      mCounts = new int[hist.mCounts.length];
      System.arraycopy(hist.mBinMin, 0, mBinMin, 0, mBinMin.length);
      System.arraycopy(hist.mCounts, 0, mCounts, 0, mCounts.length);
   }

   /**
    * Build a Histogram suitable for displaying the data in the
    * DescriptiveStatistics object with approximately the specified number of
    * bins. The algorithm attempts to chose reasonable min and max limits based
    * on the data.
    *
    * @param ds A {@link DescriptiveStatistics} object
    * @param approxNBins The approximate number of bins
    * @return Histogram
    */
   public static Histogram build(final DescriptiveStatistics ds, final int approxNBins) {
      final double min = ds.getMin(), max = ds.getMax();
      final double sc = Math.pow(10.0, Math.round(Math.log10((max - min) / approxNBins)));
      final int xtra = Math.max(1, approxNBins / 40);
      final double minR = Math.floor(min / sc) * sc, maxR = Math.ceil(max / sc) * sc;
      final int nBins = (int) Math.round((maxR - minR) / sc) + (2 * xtra);
      return new Histogram(minR - (xtra * sc), maxR + (xtra * sc), nBins);
   }

   public int maxBin() {
      int maxBin = binCount() - 1;
      int max = counts(maxBin);
      for(int i = binCount() - 2; i >= 0; --i)
         if(counts(i) > max) {
            max = counts(i);
            maxBin = i;
         }
      return maxBin;
   }

   public void addBin(final double binMin) {
      final double[] newBinMin = new double[mBinMin.length + 1];
      System.arraycopy(mBinMin, 0, newBinMin, 0, mBinMin.length);
      newBinMin[mBinMin.length] = binMin;
      Arrays.sort(newBinMin);
      mCounts = new int[newBinMin.length + 2];
      mBinMin = newBinMin;
   }

   /**
    * bin - Returns the bin into which the val fits.
    *
    * @param val double
    * @return int
    */
   public int bin(final double val) {
      int i = Arrays.binarySearch(mBinMin, val);
      i = (i >= 0 ? i : -i - 2);
      assert i >= -1 : "index is " + Integer.toString(i) + " for " + Double.toString(val);
      assert i < mCounts.length;
      return i;
   }

   /**
    * minValue - Returns the minimum value stored in the specified bin
    *
    * @param bin int
    * @return double
    */
   public double minValue(final int bin) {
      return bin > -1 ? mBinMin[bin] : Double.NEGATIVE_INFINITY;
   }

   /**
    * maxValue - Returns the upper limit for values stored in this bin. Actually
    * this value is excluded from the bin and included in the next larger bin.
    *
    * @param bin int
    * @return double
    */
   public double maxValue(final int bin) {
      return (bin + 1) < mBinMin.length ? mBinMin[bin + 1] : Double.POSITIVE_INFINITY;
   }

   /**
    * add - Add the specified value to the histogram.
    *
    * @param val double
    */
   public void add(final double val) {
      ++mCounts[bin(val) + 1];
   }

   /**
    * add - Add the specified array of values to the histogram.
    *
    * @param vals double[]
    */
   public void add(final double[] vals) {
      for(final double v : vals)
         add(v);
   }

   /**
    * binCount - Returns the number of bins (not counting over-range and
    * under-range bins)
    *
    * @return int
    */
   public int binCount() {
      return mBinMin.length - 1;
   }

   /**
    * binName - Returns a string that describes the specified bin. The String
    * format is used to format numbers using MessageFormat.format(...). An
    * common usage might be <code>binName(i,"{0,number,#.##}")</code>.
    *
    * @param bin int
    * @param format String
    * @return String
    */
   public String binName(final int bin, final String format) {
      return (new BinName(bin, format)).toString();
   }

   /**
    * binName - Returns a string that describes the specified bin. The String
    * format is used to format numbers using MessageFormat.format(...). An
    * common usage might be <code>binName(i,"{0,number,#.##}")</code>.
    *
    * @param bin int
    * @param nf NumberFormat
    * @return String
    */
   public String binName(final int bin, final NumberFormat nf) {
      return (new BinName(bin, nf)).toString();
   }

   /**
    * counts - Returns the number of counts in the specified bin.
    *
    * @param bin int - -1 to binCount()
    * @return int
    */
   public int counts(final int bin) {
      return mCounts[bin + 1];
   }

   /**
    * Reset the counts to 0.
    */
   public void clear() {
      Arrays.fill(mCounts, 0);
   }

   /**
    * overrange - Returns the number of events that fell into the overrange bin
    * (larger than the max argument of the constructor.)
    *
    * @return int
    */
   public int overrange() {
      return mCounts[mCounts.length - 1];
   }

   /**
    * underrange -Returns the number of events that fell into the underrange bin
    * (less than the min argument of the constructor.)
    *
    * @return int
    */
   public int underrange() {
      return mCounts[0];
   }

   /**
    * Returns the total number of events recorded in the histogram.
    *
    * @return int
    */
   public int totalCounts() {
      int res = 0;
      for(final int c : mCounts)
         res += c;
      return res;
   }

   public boolean isBinMin(final double binMin) {
      final int i = Arrays.binarySearch(mBinMin, binMin);
      return i >= 0;
   }

   public void removeBin(final int binNum) {
      Preconditions.checkArgument((binNum < 0) || (binNum > (mCounts.length - 2)));
      final double[] newBinMin = new double[mBinMin.length - 1];
      for(int index = 0; index < mBinMin.length; index++)
         if(index < binNum)
            newBinMin[index] = mBinMin[index];
         else if(index > binNum)
            newBinMin[index - 1] = mBinMin[index];
      mBinMin = newBinMin;
      Arrays.sort(mBinMin);

      final int[] newCounts = new int[mCounts.length - 1];
      // Since we're deleting the bin, move the counts into
      // the bin before the one being removed
      mCounts[binNum] += mCounts[binNum + 1];
      for(int index = 0; index < mCounts.length; index++)
         if(index < (binNum + 1))
            newCounts[index] = mCounts[index];
         else if(index > (binNum + 1))
            newCounts[index - 1] = mCounts[index];
      mCounts = newCounts;
   }

   public TreeMap<BinName, Integer> getResultMap(final String format) {
      final TreeMap<BinName, Integer> res = new TreeMap<>();
      for(int i = -1; i < (binCount() + 1); ++i)
         res.put(new BinName(i, format), counts(i));
      return res;
   }

   @Override
   public String toHTML(final Mode mode) {
      final Report report = new Report("Histogram");
      if((mode == Mode.NORMAL) || (mode == Mode.VERBOSE)) {
         if(mode == Mode.VERBOSE)
            report.addHeader("Histogram table");
         final BasicNumberFormat bnf = new BasicNumberFormat();
         final List<List<String>> items = new ArrayList<>();
         final List<String> colHeaders = Arrays.asList("Index", "Bin", "Min", "Max", "Count");

         for(int i = -1; i < (binCount() + 1); ++i) {
            final List<String> row = Arrays.asList(Integer.toString(i), HTML.escape(binName(i, bnf)), i >= 0
                  ? bnf.format(minValue(i))
                  : "-&infin;", i < binCount() ? bnf.format(maxValue(i)) : "+&infin;", Integer.toString(counts(i)));
            items.add(row);
         }
         report.addTable(colHeaders, null, items);
      } else
         report.addNote("Histogram with " + binCount() + " bins containing " + totalCounts() + " events.");
      return report.toHTML(mode);
   }
}
