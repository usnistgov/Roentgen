package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * @author Nicholas
 */
public class MCPropagator {

   public interface Constraint {

      public double limit(double val);

   }

   private final NamedMultivariateJacobianFunction mFunction;
   private final Constraint[] mConstraints;
   private final UncertainValues mValues;
   private final List<RealVector> mEvaluations = new ArrayList<>();
   private final List<RealVector> mPoints = new ArrayList<>();

   public class None
      implements
      Constraint {

      @Override
      public double limit(final double val) {
         return val;
      }

   }

   public class Range
      implements
      Constraint {

      private final double mMin;
      private final double mMax;

      public Range(final double min, final double max) {
         mMin = Math.min(min, max);
         mMax = Math.max(min, max);
      }

      @Override
      public double limit(final double val) {
         if(val < mMin)
            return mMin;
         else if(val > mMax)
            return mMax;
         else
            return val;
      }
   }

   private final class SafeMultivariateNormalDistribution
      extends
      AbstractMultivariateRealDistribution {

      final MultivariateNormalDistribution mDistribution;
      final int[] mMap;
      final RealVector mMeans;

      SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov) {
         super(RandomGeneratorFactory.createRandomGenerator(new Random()), vals.getDimension());
         int next = 0;
         mMeans = vals;
         mMap = new int[vals.getDimension()];
         for(int i = 0; i < mMap.length; ++i) {
            if(cov.getEntry(i, i) > 1.0e-10) {
               mMap[i] = next;
               ++next;
            } else
               mMap[i] = -1;
         }
         final double[] mval = new double[next];
         final double[][] mvar = new double[next][next];
         int nextR = 0;
         for(int r = 0; r < mMap.length; ++r) {
            if(mMap[r] != -1) {
               mval[nextR] = vals.getEntry(mMap[r]);
               int nextC = 0;
               for(int c = 0; c < mMap.length; ++c) {
                  if(mMap[c] != -1) {
                     mvar[nextR][nextC] = cov.getEntry(mMap[r], mMap[c]);
                     ++nextC;
                  }
               }
               ++nextR;
            }
         }
         mDistribution = new MultivariateNormalDistribution(mval, mvar);
      }

      @Override
      public double density(final double[] x) {
         final double[] mns = new double[mDistribution.getDimension()];
         for(int i = 0; i < mMap.length; ++i) {
            if(mMap[i] == -1) {
               if(mMeans.getEntry(i) != x[i])
                  return 0.0;
            } else {
               mns[mMap[i]] = x[i];
            }
         }
         return mDistribution.density(mns);
      }

      @Override
      public double[] sample() {
         final double[] res = mMeans.toArray();
         final double[] dist = mDistribution.sample();
         for(int r = 0; r < mMap.length; ++r)
            if(mMap[r] != -1)
               res[r] = dist[mMap[r]];
         return res;
      }
   }

   public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv) {
      mFunction = nmvjf;
      mConstraints = new Constraint[mFunction.getInputDimension()];
      for(int i = 0; i < mConstraints.length; ++i)
         mConstraints[i] = new None();
      mValues = uv;
   }

   public MCPropagator(final NamedMultivariateJacobianFunction nmvjf, final UncertainValues uv, final double nSigmaConstraint) {
      mFunction = nmvjf;
      mConstraints = new Constraint[mFunction.getInputDimension()];
      mValues = uv;
      for(int i = 0; i < mConstraints.length; ++i) {
         final double sigma = Math.max(0.001 * mValues.getEntry(i), Math.sqrt(mValues.getVariance(i)));
         mConstraints[i] = new Range(mValues.getEntry(i) - nSigmaConstraint * sigma, mValues.getEntry(i)
               + nSigmaConstraint * sigma);
      }
   }

   public void setConstraint(final Object tag, final Constraint con) {
      final int idx = mFunction.inputIndex(tag);
      mConstraints[idx] = con;
   }

   public UncertainValues compute(final int nEvals) {
      final SafeMultivariateNormalDistribution smnd = new SafeMultivariateNormalDistribution(mValues.getValues(), mValues.getCovariances());
      for(int eval = 0; eval < nEvals; ++eval) {
         final RealVector pt = MatrixUtils.createRealVector(smnd.sample());
         for(int i = 0; i < mConstraints.length; ++i)
            pt.setEntry(i, mConstraints[i].limit(pt.getEntry(i)));
         mPoints.add(pt);
         mEvaluations.add(mFunction.evaluate(pt).getFirst());
      }
      final int len = mFunction.getOutputDimension();
      final RealVector sum = new ArrayRealVector(len);
      for(final RealVector eval : mEvaluations)
         sum.combineToSelf(1.0, 1.0, eval);
      final RealVector avg = sum.mapDivide(mEvaluations.size());
      final RealMatrix cov = MatrixUtils.createRealMatrix(len, len);
      for(int r = 0; r < len; r++)
         for(int c = r; c < len; c++) {
            double tmp = 0.0;
            for(final RealVector eval : mEvaluations)
               tmp += (eval.getEntry(r) - avg.getEntry(r)) * (eval.getEntry(c) - avg.getEntry(c));
            tmp /= mEvaluations.size();
            cov.setEntry(r, c, tmp);
            cov.setEntry(c, r, tmp);
         }
      return new UncertainValues(mFunction.getOutputTags(), avg, cov);
   }

   public DescriptiveStatistics getOutputStatistics(final Object tag) {
      final DescriptiveStatistics ds = new DescriptiveStatistics(mEvaluations.size());
      final int idx = mFunction.outputIndex(tag);
      for(final RealVector eval : mEvaluations)
         ds.addValue(eval.getEntry(idx));
      return ds;
   }

   public DescriptiveStatistics getInputStatistics(final Object tag) {
      final DescriptiveStatistics ds = new DescriptiveStatistics(mPoints.size());
      final int idx = mFunction.inputIndex(tag);
      for(final RealVector pt : mPoints)
         ds.addValue(pt.getEntry(idx));
      return ds;
   }
}
