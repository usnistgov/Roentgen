package gov.nist.microanalysis.roentgen.math;

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * <p>
 * Replaces {@link MultivariateNormalDistribution} to address the shortcoming
 * that {@link MultivariateNormalDistribution} doesn't handle zero variances.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class FixedMultivariateNormalDistribution
   extends
   AbstractMultivariateRealDistribution {

   final double[] mValues;
   final int[] mIndex;
   final MultivariateNormalDistribution mDistribution;

   public FixedMultivariateNormalDistribution(final RandomGenerator rg, final double[] values, final double[][] covariances, final double tol) {
      super(rg, values.length);
      mIndex = new int[values.length];
      mValues = values;
      // Build the smaller MVND containing only elements with
      // non-zero variances;
      int p = 0;
      for(int i = 0; i < values.length; ++i)
         if(covariances[i][i] > 0.0) {
            mIndex[i] = p;
            ++p;
         } else
            mIndex[i] = -1;
      final double[] newVals = new double[p];
      final double[][] newCov = new double[p][p];
      for(int i = 0; i < p; ++i) {
         newVals[i] = values[mIndex[i]];
         for(int j = 0; j < p; ++j)
            newCov[i][j] = covariances[mIndex[i]][mIndex[j]];
      }
      mDistribution = new MultivariateNormalDistribution(newVals, newCov);
   }

   @Override
   public double[] sample() {
      final double[] res = mValues.clone();
      final double[] samp = mDistribution.sample();
      for(int i = 0; i < mIndex.length; ++i)
         if(mIndex[i] != -1)
            res[i] = samp[mIndex[i]];
      return res;
   }

   @Override
   public double density(final double[] x) {
      final double[] x2 = new double[mDistribution.getDimension()];
      int j = 0;
      for(int i = 0; i < mIndex.length; ++i) {
         if(mIndex[i] == -1) {
            if(mValues[i] != x[i])
               return 0.0;
         } else {
            x2[j] = x[i];
            ++j;
         }
      }
      return mDistribution.density(x2);
   }
}