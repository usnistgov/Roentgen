package gov.nist.juncertainty;

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

/**
 * <p>
 * An implementation of the {@link UncertainValuesBase} class that computes an
 * estimate of the values and associated uncertainty matrix for a set of samples
 * representing 'measurements.'
 * </p>
 * <p>
 * Computes the expectation values E[y<sub>i</sub>] and
 * E[(y<sub>i</sub>-E[y<sub>i</sub>])(y<sub>j</sub>-E[y<sub>j</sub>])]. The
 * value E[y<sub>i</sub>] represents the best estimate of the value of
 * Y<sub>i</sub> and the value
 * E[(y<sub>i</sub>-E[y<sub>i</sub>])(y<sub>j</sub>-E[y<sub>j</sub>])]
 * represents the best estimate of the covariance matrix element
 * C<sub>i,j</sub>.
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class EstimateUncertainValues<H> //
		extends UncertainValuesBase<H> {

	/**
	 * A list of the sampled values. The dimension and ordering of the values in
	 * each RealVector is defined by the labels in getLabels().
	 */
	private final ArrayList<RealVector> mSamples = new ArrayList<RealVector>();

	private final SimplyLazy<UncertainValues<H>> mResults = new SimplyLazy<UncertainValues<H>>() {

		@Override
		protected UncertainValues<H> initialize() {
			synchronized (mSamples) {
				final int len = getDimension();
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
	 * Create an object to accumulate samples to estimate the values and covariance
	 * matrix associated with the specified set of values identified by the labels.
	 * object.
	 *
	 * @param labels A list of parameter labels
	 */
	public EstimateUncertainValues(
			final List<? extends H> labels
	) {
		super(new ArrayList<H>(labels));
	}

	/**
	 * Compute <code>nSamples</code> draws from the specified
	 * {@link MultivariateRealDistribution} and add them to the list of samples.
	 *
	 * @param mvd      {@link MultivariateRealDistribution}
	 * @param nSamples The number of samples to compute
	 */
	public void add(
			final MultivariateRealDistribution mvd, //
			final int nSamples
	) {
		for (int i = 0; i < nSamples; ++i)
			add(mvd.sample());
	}

	/**
	 * Add a measurement / sample to the data set used to compute the estimated
	 * uncertainties. The length of sample must match the size and order of the
	 * labels List passed to the constructor.
	 *
	 *
	 * @param sample double[]
	 * @return {@link RealVector} The input vector
	 */
	public double[] add(
			final double[] sample
	) {
		add(new ArrayRealVector(sample));
		return sample;
	}

	/**
	 * Add a measurement / sample to the data set used to compute the estimated
	 * uncertainties. The length of rv must match the size and order of the labels
	 * List passed to the constructor.
	 *
	 *
	 * @param rv {@link RealVector}
	 * @return {@link RealVector} The input vector
	 */
	public RealVector add(
			final RealVector rv
	) {
		assert rv.getDimension() == getDimension();
		if (rv.getDimension() == getDimension())
			synchronized (mSamples) {
				mSamples.add(rv);
				mResults.reset();
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
	 * @param label A label
	 * @return {@link DescriptiveStatistics}
	 */
	public DescriptiveStatistics getDescriptiveStatistics(
			final H label
	) {
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
