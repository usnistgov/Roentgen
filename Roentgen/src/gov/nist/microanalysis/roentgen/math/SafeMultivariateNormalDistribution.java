package gov.nist.microanalysis.roentgen.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;
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
 * @author Nicholas W. M. Ritchie
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

	public SafeMultivariateNormalDistribution(final RealVector vals, final RealMatrix cov, final double tol) {
		super(RandomGeneratorFactory.createRandomGenerator(new Random()), vals.getDimension());
		mMean = vals.toArray();
		final int dim = vals.getDimension();
		assert cov.getRowDimension() == dim;
		assert cov.getColumnDimension() == dim;
		final List<Set<Integer>> groups = new ArrayList<>();
		for (int r = 0; r < dim; ++r)
			for (int c = r + 1; c < dim; ++c)
				if (Math.abs(cov.getEntry(r, c)) > tol * Math.sqrt(cov.getEntry(r, r) * cov.getEntry(c, c))) {
					assert Math.abs(cov.getEntry(c, r)) > tol * Math.sqrt(cov.getEntry(r, r) * cov.getEntry(c, c));
					final Integer v1 = Integer.valueOf(r);
					Set<Integer> g1 = null;
					for (final Set<Integer> group : groups)
						if (group.contains(v1)) {
							g1 = group;
							break;
						}
					final Integer v2 = Integer.valueOf(c);
					Set<Integer> g2 = null;
					for (final Set<Integer> group : groups)
						if (group.contains(v2)) {
							g2 = group;
							break;
						}
					if ((g1 != null) && (g2 != null)) {
						// Combine g1 and g2
						if (g1 != g2) {
							g1.addAll(g2);
							assert g1.contains(v1) && g1.contains(v2);
							groups.remove(g2);
						}
					} else if (g1 != null) {
						assert g2 == null;
						g1.add(v2);
					} else if (g2 != null) {
						assert g1 == null;
						g2.add(v1);
					} else {
						assert (g1 == null) && (g2 == null);
						final Set<Integer> newGroup = new TreeSet<>();
						newGroup.add(v1);
						newGroup.add(v2);
						groups.add(newGroup);
					}
				}
		// Which distribution is this row/column located?
		mDistIndex = new int[dim];
		Arrays.fill(mDistIndex, -1);
		for (int groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
			final Set<Integer> group = groups.get(groupIdx);
			for (int row = 0; row < dim; ++row)
				if (group.contains(Integer.valueOf(row))) {
					assert mDistIndex[row] == -1;
					mDistIndex[row] = groupIdx;
				}
		}
		mDistribution = new MultivariateNormalDistribution[groups.size()];
		mInnerIndex = new int[dim];
		Arrays.fill(mInnerIndex, -1);
		for (int groupIdx = 0; groupIdx < groups.size(); ++groupIdx) {
			final Set<Integer> rows = groups.get(groupIdx);
			// Build one MvND per itemJ
			final int dimRows = rows.size();
			assert dimRows > 1;
			final double[] means = new double[dimRows];
			final double[][] covariances = new double[dimRows][dimRows];
			int i = 0;
			for (final int row : rows) {
				assert mDistIndex[row] == groupIdx;
				mInnerIndex[row] = i;
				means[i] = vals.getEntry(row);
				int k = 0;
				for (final int col : rows) {
					covariances[i][k] = cov.getEntry(row, col);
					++k;
				}
				++i;
			}
			assert MatrixUtils.isSymmetric(MatrixUtils.createRealMatrix(covariances), 1.0e-9);
			assert validateCovariances(covariances);
			// System.err.println("Means[" + groupIdx + "] = " + Arrays.toString(means));
			// System.err.println("Covariances[" + groupIdx + "] = " +
			// MatrixUtils.createRealMatrix(covariances).toString());
			mDistribution[groupIdx] = new MultivariateNormalDistribution(random, means, covariances);
		}
		mNormal = new NormalDistribution[mDistIndex.length];
		for (int i = 0; i < mDistIndex.length; ++i)
			if ((mDistIndex[i] == -1) && (cov.getEntry(i, i) > tol * vals.getEntry(i)))
				mNormal[i] = new NormalDistribution(vals.getEntry(i), Math.sqrt(cov.getEntry(i, i)));
	}

	private boolean validateCovariances(final double[][] covs) {
		boolean result = true;
		for (int r = 0; r < covs.length; ++r) {
			final double rv = covs[r][r];
			if (rv < 0.0) {
				System.err.println(rv + " is less than zero at row " + r);
				result = false;
			}
			for (int c = r + 1; c < covs.length; ++c) {
				final double cv = covs[c][c];
				if (cv < 0.0) {
					System.err.println(cv + " is less than zero at column " + c);
					result = false;
				}
				if (Math.abs(covs[r][c] - covs[c][r])>1.0e-6*Math.max(Math.abs(covs[r][c]), Math.abs(covs[c][r]))) {
					System.err.println("Not symmetric at " + r + ", " + c + " -> " + covs[r][c] + " != " + covs[c][r]);
					result = false;
				}
				if (covs[r][c] * covs[r][c] > rv * cv) {
					System.err.println("Covariance[" + r + ", " + c + "] -> " + covs[r][c] + ", " + rv + ", " + cv);
					result = false;
				}
			}
		}
		return result;
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
