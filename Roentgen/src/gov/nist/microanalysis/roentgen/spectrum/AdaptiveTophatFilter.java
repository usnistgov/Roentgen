package gov.nist.microanalysis.roentgen.spectrum;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import gov.nist.microanalysis.roentgen.math.IntInterval;
import gov.nist.microanalysis.roentgen.math.Utility;
import gov.nist.microanalysis.roentgen.spectrum.LineshapeCalibration.InvertMode;

/**
 * <p>
 * A filter constructed from the class Schamber-style top-hat filter. The width
 * adapts to match the nominal FWHM at the specified channel.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class AdaptiveTophatFilter extends EDSFittingFilter {

	final static double FILTER_WIDTH = 1.0;
	final static double EXTRA = 1.0;
	private static boolean FAST_EXTENT = false;

	/**
	 * Constructs a AdaptiveTophatFilter
	 *
	 * @param nChannels Number of channels in the spectrum data.
	 * @param ec        EnergyCalibration
	 * @param ls        LineshapeCalibration
	 */
	public AdaptiveTophatFilter(final int nChannels, final EnergyCalibration ec, final LineshapeCalibration ls) {
		super(nChannels, ec, ls);
	}

	public static void setUseFastExtent(final boolean fast) {
		FAST_EXTENT = fast;
	}

	/**
	 * Construct the Jacobian matrix for the transform from spectrum channels into
	 * filtered spectrum channels.
	 *
	 * @see gov.nist.microanalysis.roentgen.math.uncertainty.MultiLinearJacobianFunction#buildLinearTransform(int)
	 */
	@Override
	public RealMatrix buildLinearTransform(final int nCh, final int ignored) {
		final RealMatrix res = MatrixUtils.createRealMatrix(nCh, nCh);
		for (int fCh = 0; fCh < nCh; ++fCh) {
			final double eCh = mEnergy.averageEnergyForChannel(fCh);
			final double width = FILTER_WIDTH * mLineshape.getFWHM(eCh);
			final int minO = Math.max(0, mEnergy.channelIndex(eCh - EXTRA * width));
			final int minI = Math.max(0, mEnergy.channelIndex(eCh - 0.5 * EXTRA * width));
			assert minO <= minI;
			final int maxI = Math.min(nCh - 1, mEnergy.channelIndex(eCh + 0.5 * EXTRA * width));
			final int maxO = Math.min(nCh - 1, mEnergy.channelIndex(eCh + EXTRA * width));
			assert maxO >= maxI : maxO + ">" + maxI;
			final double nNorm = -1.0 / ((minI - minO) + (maxO - maxI));
			final double pNorm = 1.0 / (maxI - minI);
			for (int i = minO; i < minI; ++i)
				res.setEntry(fCh, i, nNorm);
			for (int i = minI; i < maxI; ++i)
				res.setEntry(fCh, i, pNorm);
			for (int i = maxI; i < maxO; ++i)
				res.setEntry(fCh, i, nNorm);
			assert Math.abs(Utility.sum(res.getRow(fCh))) < pNorm * 1.0e-6;
		}
		return res;
	}

	@Override
	public IntInterval extent(final double e, final double frac) {
		if (FAST_EXTENT) {
			final double filtW = EXTRA * FILTER_WIDTH * mLineshape.getFWHM(e);
			final double f = 3.0e1 * frac;
			if (f < 1.0) {
				final double low = mLineshape.invert(e, InvertMode.LOWER_HALF, 3.0e1 * frac);
				final double high = mLineshape.invert(e, InvertMode.UPPER_HALF, 3.0e1 * frac);
				final int minCh = mEnergy.channelIndex(e - (filtW + low));
				final int maxCh = mEnergy.channelIndex(e + (filtW + high));
				// return new IntInterval(minCh, maxCh);
				return new IntInterval(Math.max(0, minCh), Math.min(maxCh, getInputDimension()));
			} else
				return IntInterval.NULL_INTERVAL;
		} else
			return super.extent(e, frac);
	}

}
