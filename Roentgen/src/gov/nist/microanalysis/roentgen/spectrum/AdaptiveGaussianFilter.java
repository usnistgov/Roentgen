package gov.nist.microanalysis.roentgen.spectrum;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.Utility;

/**
 * <p>
 * A filter constructed from the second derivative of the PDF for the normal
 * distibution. The width adapts to match the nominal resolution of the detector
 * at the specified channel.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class AdaptiveGaussianFilter extends EDSFittingFilter {

	final private static double FILTER_WIDTH = 1.0 / 2.35;
	final private static double EXTRA = 3.0;

	/**
	 * Constructs a AdaptiveTophatFilter
	 *
	 * @param nChannels Number of channels in the spectrum data.
	 * @param ec        EnergyCalibration
	 * @param ls        LineshapeCalibration
	 * @throws ArgumentException
	 */
	public AdaptiveGaussianFilter(
			final int nChannels, //
			final EnergyCalibration ec, //
			final LineshapeCalibration ls
	) throws ArgumentException {
		super(nChannels, ec, ls);
	}

	/**
	 * Construct the Jacobian matrix for the transform from spectrum channels into
	 * filtered spectrum channels.
	 *
	 * @see gov.nist.microanalysis.roentgen.math.uncertainty.MultiLinearMeasurementModel#buildLinearTransform(int)
	 */
	@Override
	public RealMatrix buildLinearTransform(
			final int nCh, final int ignored
	) {
		final RealMatrix res = MatrixUtils.createRealMatrix(nCh, nCh);
		for (int fCh = 0; fCh < nCh; ++fCh) {
			final double eCh = mEnergy.averageEnergyForChannel(fCh);
			final double sigma = FILTER_WIDTH * mLineshape.getFWHM(eCh);
			final int minO = Math.max(0, mEnergy.channelIndex(eCh - EXTRA * sigma));
			final int maxO = Math.min(nCh - 1, mEnergy.channelIndex(eCh + EXTRA * sigma));
			double sum = 0.0;
			for (int ch = minO; ch < maxO; ++ch) {
				final double de = eCh - mEnergy.averageEnergyForChannel(ch);
				final double g = Math.exp(-0.5 * Math.pow(de / sigma, 2.0));
				final double entry = g * (1.0 - Math.pow(de / sigma, 2.0));
				res.setEntry(fCh, ch, entry);
				sum += entry;
			}
			sum /= (maxO - minO);
			for (int i = minO; i < maxO; ++i)
				res.setEntry(fCh, i, res.getEntry(fCh, i) - sum);
			assert Math.abs(Utility.sum(res.getRow(fCh))) < 1.0e-6;
		}
		return res;
	}
}
