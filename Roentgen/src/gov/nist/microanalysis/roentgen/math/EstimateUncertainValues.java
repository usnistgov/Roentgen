/**
 * 
 */
package gov.nist.microanalysis.roentgen.math;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;

/**
 * Computes an estimated UncertainValues object give a set of samples
 * representing measurements of the values specified in 'tags'
 * 
 * 
 * @author Nicholas
 *
 */
public class EstimateUncertainValues {

	private final int mDimension;
	private final ArrayList<? extends Object> mTags;

	private final ArrayList<RealVector> mSamples = new ArrayList<RealVector>();

	/**
	 * Create an object to accumulate samples to estimate an {@link UncertainValues}
	 * object.
	 * 
	 * @param tags The names of the parameters
	 */
	public EstimateUncertainValues(List<? extends Object> tags) {
		mDimension = tags.size();
		mTags = new ArrayList<>(tags);
	}

	public void add(MultivariateRealDistribution mvd, int nSamples) {
		for (int i = 0; i < nSamples; ++i)
			add(mvd.sample());
	}

	public double[] add(double[] sample) {
		assert sample.length == mDimension;
		synchronized (mSamples) {
			mSamples.add(new ArrayRealVector(sample));
		}
		return sample;
	}

	public RealVector add(RealVector rv) {
		assert rv.getDimension() == mDimension;
		synchronized (mSamples) {
			mSamples.add(rv);
		}
		return rv;
	}

	public UncertainValues estimateDistribution() {
		final int len = mDimension;
		synchronized (mSamples) {
			final RealVector sum = new ArrayRealVector(len);
			for (final RealVector samp : mSamples)
				sum.combineToSelf(1.0, 1.0, samp);
			final RealVector avg = sum.mapDivide(mSamples.size());
			final RealMatrix cov = MatrixUtils.createRealMatrix(len, len);
			for (int r = 0; r < len; r++)
				for (int c = r; c < len; c++) {
					double tmp = 0.0;
					for (final RealVector samp : mSamples)
						tmp += (samp.getEntry(r) - avg.getEntry(r)) * (samp.getEntry(c) - avg.getEntry(c));
					tmp /= mSamples.size();
					cov.setEntry(r, c, tmp);
					cov.setEntry(c, r, tmp);
				}
			return new UncertainValues(mTags, avg, cov);
		}
	}

	public List<? extends Object> getTags() {
		return mTags;
	}

	public int getSampleCount() {
		return mSamples.size();
	}

	public List<RealVector> getSamples() {
		return new ArrayList<>(mSamples);
	}

	public DescriptiveStatistics getDescriptiveStatistics(final Object tag) {
		final DescriptiveStatistics ds = new DescriptiveStatistics(getSampleCount());
		final int idx = mTags.indexOf(tag);
		for (final RealVector eval : mSamples)
			ds.addValue(eval.getEntry(idx));
		return ds;
	}

}
