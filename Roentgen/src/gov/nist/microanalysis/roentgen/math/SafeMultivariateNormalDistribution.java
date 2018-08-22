package gov.nist.microanalysis.roentgen.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.AbstractMultivariateRealDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGeneratorFactory;

/**
 * This class is designed to avoid some of the problems associated with
 * generating random variants from a correlated uncertainties.
 * <ol>
 * <li>The {@link MultivariateNormalDistribution} class doesn't handle variables
 * with zero variance.
 * <li>The {@link MultivariateNormalDistribution} class doesn't handle variables
 * which range over a large range of mean values.
 * <li>Address the occasional NaN returned by
 * {@link MultivariateNormalDistribution}
 * </ol>
 * 
 * The SafeMultivariateNormalDistribution class addresses these points by
 * explicitly returning the mean value when the variance is zero and by breaking
 * the problem up into a number of smaller pieces based on mutual correlation.
 * 
 * @author nicholas
 *
 */
public class SafeMultivariateNormalDistribution extends AbstractMultivariateRealDistribution {

	final MultivariateNormalDistribution[] mDistribution;
	final NormalDistribution[] mNormal;
	final int[] mDistIndex;
	final int[] mInnerIndex;
	final double[] mMean;

	public SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov) {
		this(vals, cov, 1.0e-9);
	}

	public SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov, double tol) {
		super(RandomGeneratorFactory.createRandomGenerator(new Random()), vals.getDimension());
		mMean = vals.toArray();
		int dim = vals.getDimension();
		assert cov.getRowDimension() == dim;
		assert cov.getColumnDimension() == dim;
		List<TreeSet<Integer>> correl = new ArrayList<>();
		for (int r = 0; r < dim; ++r) {
			for (int c = r + 1; c < dim; ++c) {
				if (Math.abs(cov.getEntry(r, c)) > tol * Math.sqrt(cov.getEntry(r, r) * cov.getEntry(c, c))) {
					Integer rr = Integer.valueOf(r);
					boolean done = false;
					for (TreeSet<Integer> rows : correl)
						if (rows.contains(rr)) {
							rows.add(Integer.valueOf(c));
							done = true;
							break;
						}
					if (!done) {
						TreeSet<Integer> row = new TreeSet<>();
						row.add(rr);
						row.add(Integer.valueOf(c));
						correl.add(row);
					}
				}
			}
		}
		mDistIndex = new int[dim];
		mInnerIndex = new int[dim];
		Arrays.fill(mDistIndex, -1);
		Arrays.fill(mInnerIndex, -1);
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < correl.size(); ++j)
				if (correl.get(j).contains(Integer.valueOf(i)))
					mDistIndex[i] = j;
		}
		mDistribution = new MultivariateNormalDistribution[correl.size()];
		for (int j = 0; j < correl.size(); ++j) {
			TreeSet<Integer> itemJ = correl.get(j);
			// Build one MvND per itemJ
			final int dimJ = itemJ.size();
			assert dimJ > 1;
			final double[] means = new double[dimJ];
			final double[][] covariances = new double[dimJ][dimJ];
			int i = 0;
			for (int row : itemJ) {
				mDistIndex[row] = j;
				mInnerIndex[row] = i;
				means[i] = vals.getEntry(row);
				int k = 0;
				for (int col : itemJ) {
					covariances[i][k] = cov.getEntry(row, col);
					++k;
				}
				++i;
			}
			assert MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(covariances), 1.0e-9);
			mDistribution[j] = new MultivariateNormalDistribution(random, means, covariances);
		}
		mNormal = new NormalDistribution[mDistIndex.length];
		for (int i = 0; i < mDistIndex.length; ++i)
			if ((mDistIndex[i] == -1) && (cov.getEntry(i, i) > tol * vals.getEntry(i)))
				mNormal[i] = new NormalDistribution(vals.getEntry(i), Math.sqrt(cov.getEntry(i, i)));
	}

	/*
	 * Not yet tested... (non-Javadoc)
	 * 
	 * @see
	 * org.apache.commons.math3.distribution.MultivariateRealDistribution#density(
	 * double[])
	 */
	@Override
	public double density(final double[] x) {
		double res = 1.0;
		for (int i = 0; i < mNormal.length; ++i)
			if (mNormal[i] != null)
				res = mNormal[i].density(x[i]);
			else if ((mDistIndex[i] == -1) && (x[i] != mMean[i]))
				return 0.0;
		// For each MND compute the density
		for (int i = 0; i < mDistribution.length; ++i) {
			final MultivariateNormalDistribution mnd = mDistribution[i];
			final double[] xi = new double[mnd.getDimension()];
			for (int r = 0; r < mDistIndex.length; ++r)
				if (mDistIndex[r] == i)
					xi[mInnerIndex[i]] = x[r];
			res *= mnd.density(xi);
		}
		return res;
	}

	@Override
	public double[] sample() {
		final double[] res = mMean.clone();
		for (int i = 0; i < mNormal.length; ++i)
			if (mNormal[i] != null) {
				res[i] = Double.NaN;
				while (Double.isNaN(res[i]))
					res[i] = mNormal[i].sample();
			}
		for (int i = 0; i < mDistribution.length; ++i) {
			// Handle the occasional NaN returned by MultivariateNormalDistribution
			boolean done = false;
			while (!done) {
				final double[] samp = mDistribution[i].sample();
				for (int j = 0; j < mDistIndex.length; ++j)
					if (mDistIndex[j] == i) {
						res[j] = samp[mInnerIndex[j]];
						if (Double.isNaN(res[j]))
							continue;
					}
				done = true;
			}
		}
		return res;
	}
}
