package gov.nist.microanalysis.roentgen.math;

import java.util.Random;

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGeneratorFactory;

public class SafeMultivariateNormalDistribution extends AbstractMultivariateRealDistribution {

	final MultivariateNormalDistribution mDistribution;
	final int[] mMap;
	final RealVector mMeans;

	public SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov) {
		this(vals, cov, 1.0e-6);
	}

	public SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov, final double tol) {
		super(RandomGeneratorFactory.createRandomGenerator(new Random()), vals.getDimension());
		mMeans = vals;
		mMap = new int[vals.getDimension()];
		assert cov.getRowDimension() == vals.getDimension();
		assert cov.getColumnDimension() == vals.getDimension();
		double maxCov = 0.0;
		for (int i = 0; i < mMap.length; ++i)
			maxCov = Math.max(maxCov, cov.getEntry(i, i));
		int count = 0;
		for (int i = 0; i < mMap.length; ++i) {
			if (cov.getEntry(i, i) > tol * maxCov) {
				mMap[i] = count;
				++count;
			} else
				mMap[i] = -1;
		}

		final double[] mval = new double[count];
		final double[][] mvar = new double[count][count];
		for (int r = 0; r < mMap.length; ++r) {
			final int mapR = mMap[r];
			if (mapR >= 0) {
				assert mapR < count;
				mval[mapR] = vals.getEntry(r);
				for (int c = 0; c < mMap.length; ++c) {
					final int mapC = mMap[c];
					if (mapC >= 0) {
						assert mapC < count;
						assert (mapR != mapC) || (cov.getEntry(r, c) > tol * maxCov) //
						: mapR + ", " + mapC + ", " + r + ", " + c + ", " + cov.getEntry(r, c) + ", " + maxCov;
						mvar[mapR][mapC] = cov.getEntry(r, c);
					}
				}
			}
		}
		mDistribution = new MultivariateNormalDistribution(mval, mvar);
	}

	@Override
	public double density(final double[] x) {
		final double[] mns = new double[mDistribution.getDimension()];
		for (int i = 0; i < mMap.length; ++i) {
			if (mMap[i] == -1) {
				if (mMeans.getEntry(i) != x[i])
					return 0.0;
			} else
				mns[mMap[i]] = x[i];
		}
		return mDistribution.density(mns);
	}

	@Override
	public double[] sample() {
		final double[] res = mMeans.toArray();
		final double[] dist = mDistribution.sample();
		for (int r = 0; r < mMap.length; ++r)
			if (mMap[r] != -1)
				res[r] = dist[mMap[r]];
		return res;
	}
}

