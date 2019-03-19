package gov.nist.microanalysis.roentgen.math;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;

/**
 * Computes an estimated UncertainValues object give a set of samples
 * representing measurements of the values specified in 'tags'
 *
 *
 * @author Nicholas
 *
 */
public class EstimateUncertainValues<H> extends UncertainValuesBase<H> {

	private final int mDimension;

	private final ArrayList<RealVector> mSamples = new ArrayList<RealVector>();

	private final SimplyLazy<UncertainValues<H>> mResults = new SimplyLazy<UncertainValues<H>>() {

		@Override
		protected UncertainValues<H> initialize() {
			synchronized (mSamples) {
				final int len = mDimension;
				final RealVector sum = new ArrayRealVector(len);
				for (final RealVector samp : mSamples)
					sum.combineToSelf(1.0, 1.0, samp);
				final RealVector avg = sum.mapDivide(mSamples.size());
				final RealMatrix cov = MatrixUtils.createRealMatrix(len, len);

				for (int r = 0; r < len; r++) {
					for (int c = r; c < len; c++) {
						double tmp = 0.0;
						for (final RealVector samp : mSamples)
							tmp += samp.getEntry(r) * samp.getEntry(c);
						final double coval = tmp / mSamples.size() - avg.getEntry(r) * avg.getEntry(c);
						cov.setEntry(r, c, coval);
						cov.setEntry(c, r, coval);
					}
				}
				return new UncertainValues<H>(getLabels(), avg, cov);
			}
		}
	};

	/**
	 * Create an object to accumulate samples to estimate an {@link UncertainValues}
	 * object.
	 *
	 * @param tags The names of the parameters
	 */
	public EstimateUncertainValues(final List<? extends H> tags) {
		super(new ArrayList<H>(tags));
		mDimension = tags.size();
	}

	/**
	 * Compute <code>nSamples</code> draws from the specified
	 * {@link MultivariateRealDistribution} and add them to the list of samples.
	 *
	 * @param mvd      {@link MultivariateRealDistribution
	 * @param nSamples int
	 */
	public void add(final MultivariateRealDistribution mvd, final int nSamples) {
		for (int i = 0; i < nSamples; ++i)
			add(mvd.sample());
	}

	/**
	 * Add a measurement / sample to the data set used to compute the estimated
	 * uncertainties.
	 *
	 *
	 * @param sample double[]
	 * @return {@link RealVector} The input vector
	 */
	public double[] add(final double[] sample) {
		assert sample.length == mDimension;
		synchronized (mSamples) {
			mSamples.add(new ArrayRealVector(sample));
		}
		return sample;
	}

	/**
	 * Add a measurement / sample to the data set used to compute the estimated
	 * uncertainties.
	 *
	 *
	 * @param rv {@link RealVector}
	 * @return {@link RealVector} The input vector
	 */
	public RealVector add(final RealVector rv) {
		assert rv.getDimension() == mDimension;
		synchronized (mSamples) {
			mSamples.add(rv);
		}
		return rv;
	}

	/**
	 * Number of samples / measurements which will be used to compute the estimated
	 * uncertain values object.
	 *
	 * @return int
	 */
	public int getSampleCount() {
		return mSamples.size();
	}

	/**
	 * Returns an unmodifiable list view of the set of sample points from which the
	 * estimated uncertainty will be computed.
	 *
	 * @return List&lt;RealVector&gt;
	 */
	public List<RealVector> getSamples() {
		return Collections.unmodifiableList(mSamples);
	}

	/**
	 * Returns a {@link DescriptiveStatistics} object that provides insight into the
	 * variable associated with the specified tag.
	 *
	 * @param label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getDescriptiveStatistics(final H label) {
		final DescriptiveStatistics ds = new DescriptiveStatistics(getSampleCount());
		final int idx = indexOf(label);
		for (final RealVector eval : mSamples)
			ds.addValue(eval.getEntry(idx));
		return ds;
	}

	@Override
	public RealVector getValues() {
		return mResults.get().getValues();
	}

	@Override
	public RealMatrix getCovariances() {
		return mResults.get().getCovariances();
	}

}
