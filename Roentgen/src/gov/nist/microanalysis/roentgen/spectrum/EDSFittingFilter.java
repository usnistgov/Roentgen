package gov.nist.microanalysis.roentgen.spectrum;

import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.IntInterval;
import gov.nist.microanalysis.roentgen.math.uncertainty.MultiLinearJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.spectrum.LineshapeCalibration.InvertMode;

/**
 * <p>
 * A base class for EDS Schamber filter fitting style filters.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
abstract public class EDSFittingFilter
   extends
   MultiLinearJacobianFunction {

   final static protected double GCONST = 1.0 / Math.sqrt(2.0 * Math.PI);

   protected final EnergyCalibration mEnergy;
   protected final LineshapeCalibration mLineshape;
   protected double mMinE = 100.0; // eV

   /**
    * Constructs a EDSFittingFilter
    *
    * @param nInput
    * @param inputPrefix
    * @param outputPrefix
    */
   public EDSFittingFilter(final int nChannels, final EnergyCalibration ec, final LineshapeCalibration ls) {
      super(nChannels, "C", "S");
      mEnergy = ec;
      mLineshape = ls;
   }

   /**
    * Computes the extent of the convolved filter and spectrum data.
    *
    * @param e The energy of the Gaussian shaped peak
    * @param frac Fraction of the total filtered intensity to use to define the
    *           extent
    * @return int[] with [lowCh, highCh)
    */
   public IntInterval extent(final double e, final double frac) {
      final RealVector tmp = new ArrayRealVector(getInputDimension());
      final double low = e - mLineshape.invert(e, InvertMode.LOWER_HALF, 1.0e-2 * frac);
      final double high = e + mLineshape.invert(e, InvertMode.UPPER_HALF, 1.0e-2 * frac);
      final int lowCh = Math.max(0, mEnergy.channelIndex(low));
      final int highCh = Math.min(mEnergy.channelIndex(high) + 1, getInputDimension());
      for(int ch = lowCh; ch < highCh; ++ch) {
         final double val = mLineshape.compute(e, mEnergy.minEnergyForChannel(ch), mEnergy.maxEnergyForChannel(ch));
         tmp.setEntry(ch, val);
      }
      final Pair<RealVector, RealMatrix> res2 = evaluate(tmp);
      final RealVector filtered = res2.getFirst();
      double sum = 0.0;
      for(int ch = 0; ch < filtered.getDimension(); ++ch)
         sum += Math.abs(filtered.getEntry(ch));
      final int ch = mEnergy.channelIndex(e);
      double remaining = sum - Math.abs(filtered.getEntry(ch));
      int lowRes = ch, highRes = ch + 1;
      while(remaining > frac * sum) {
         final double lowVal = lowRes > 0 ? Math.abs(filtered.getEntry(lowRes - 1)) : 0.0;
         final double highVal = highRes < filtered.getDimension() ? Math.abs(filtered.getEntry(highRes)) : 0.0;
         if(lowVal > highVal) {
            remaining -= lowVal;
            lowRes--;
         } else {
            remaining -= highVal;
            highRes++;
         }
      }
      return new IntInterval(lowRes, highRes - 1);
   }

   public Set<IntInterval> extents(final Element elm, final double eMax, final double frac) {
      final TreeSet<IntInterval> res = new TreeSet<>();
      final CharacteristicXRay[] xrays = CharacteristicXRay.forElement(elm, mMinE, eMax);
      for(final CharacteristicXRay cxr : xrays) {
         final double wgt = cxr.getWeight();
         if(wgt > frac) {
            final IntInterval extent = extent(cxr.getEnergy(), Math.max(0.01, frac / wgt));
            IntInterval.add(res, extent);
         }
      }
      return res;
   }

   /**
    * <p>
    * Sets all channels <b>not</b> associated with specified intervals to zero.
    * Apply this to reference spectra to eliminate the contribution of those
    * regions that are not associate with characteristic peaks.
    * </p>
    * <p>
    * Typically used with EDSFittingFilter.extents(...).
    * </p>
    *
    * @param Set&lt;IntInterval&gt;
    * @param eMax double All channels above this energy are also zeroed.
    * @param filtered The filtered spectrum data
    * @return RealVector The filtered data with intervening regions zeroed.
    */
   public RealVector zero(final Set<IntInterval> exts, final double eMax, final RealVector filtered) {
      final RealVector res = new ArrayRealVector(filtered);
      final SortedSet<IntInterval> ext = new TreeSet<>(exts);
      final int maxCh = Math.min(filtered.getDimension(), mEnergy.channelIndex(eMax));
      int min = 0;
      for(final IntInterval ii : ext) {
         final int max = Math.min(ii.lowerLimit(), maxCh);
         for(int ch = min; ch < max; ++ch)
            res.setEntry(ch, 0.0);
         min = ii.upperLimit() + 1;
         if(max == maxCh) {
            for(int ch = max; ch < filtered.getDimension(); ++ch)
               res.setEntry(ch, 0.0);
            break;
         }
      }
      return res;
   }
}
