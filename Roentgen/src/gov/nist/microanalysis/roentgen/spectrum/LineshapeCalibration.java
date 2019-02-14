package gov.nist.microanalysis.roentgen.spectrum;

import java.util.Objects;

import org.apache.commons.math3.special.Erf;

/**
 * <p>
 * A lineshape calibration abstracts the notion of resolution. Given an incident
 * X-ray at energy E, the LineshapeCalibration gives the distribution of
 * measured energies.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 312 $
 */
abstract public class LineshapeCalibration {

	static final public double MN_KA = 5898.7; // eV

	protected LineshapeCalibration() {
		super();
	}

	/**
	 * Returns the x-ray intensity that would be measured on average between the
	 * energies minE and maxE for an x-ray with energy xRayEnergy. Note that
	 * compute(xRayEnergy, -infinity, infinity) == 1.0 for normalization.
	 *
	 * @param xRayEnergy
	 * @param minE
	 * @param maxE
	 * @return double on range (0.0, 1.0)
	 */
	abstract public double compute(double xRayEnergy, double minE, double maxE);

	/**
	 * Returns the x-ray intensity that would be measured at energy E for an x-ray
	 * with energy xRayEnergy.
	 *
	 * @param xRayEnergy
	 * @param e
	 * @return double on range (0.0, 1.0)
	 */
	abstract public double compute(double xRayEnergy, double e);

	public static enum InvertMode {
		LOWER_HALF, FULL, UPPER_HALF
	};

	/**
	 * For a line centered at xRayEnergy, this function returns the distance until
	 * the initial intensity has dropped to 'fraction' of its intial value.
	 *
	 * @param xRayEnergy
	 * @param mode
	 * @param fraction   fraction on (0,1.0]
	 * @return double
	 */
	abstract public double invert(double xRayEnergy, InvertMode mode, double fraction);

	/**
	 * Returns the FWHM in eV at Mn Ka.
	 *
	 * @return double
	 */
	abstract public double getFWHM(double e);

	/**
	 * Takes the raw counts in <code>data</code> and convolves them with the
	 * LineshapeCalibration to produce a spectrum with the equivalent integral.
	 *
	 * @param ec
	 * @param data
	 * @return double[] the data convolved with the lineshape function.
	 */
	public double[] convolve(final EnergyCalibration ec, final double[] data) {
		final double[] res = new double[data.length];
		for (int ch = 0; ch < res.length; ch++) {
			final double e = ec.compute(ch + 0.5);
			final double minE = e - invert(e, InvertMode.LOWER_HALF, 0.001);
			final double maxE = e + invert(e, InvertMode.UPPER_HALF, 0.001);
			final int minCh = Math.max(0, ec.channelIndex(minE));
			final int maxCh = Math.min(data.length, ec.channelIndex(maxE));
			for (int ch2 = minCh; ch2 < maxCh; ++ch2)
				res[ch] += data[ch2] * compute(e, 0.5 * (ec.minEnergyForChannel(ch2) + ec.maxEnergyForChannel(ch2)));
		}
		return res;
	}

	static public class Gaussian extends LineshapeCalibration {
		private final double SQRT_2PI_INV = 1.0 / Math.sqrt(2.0 * Math.PI);

		final double mFano;
		final double mNoise;
		final double mEnergyPerEHPair;

		/**
		 * The average energy required to create one electron-hole pair for Si at
		 * typical SDD temperatures.
		 */
		final static public double SDD_EV_PER_EH = 3.64;
		/**
		 * The average energy required to create one electron-hole pair for Si at
		 * typical Si(Li) temperatures.
		 */
		final static public double SILI_EV_PER_EH = 3.71;

		/**
		 * The energy of the Mn Ka in eV.
		 */
		public final static double MN_KA = 5898.7;

		/**
		 * Compute the FWHM from the Gaussian width
		 *
		 * @param gaussianWidth
		 * @return double
		 */
		static public double computeFWHM(final double gaussianWidth) {
			return 2.354820045030949382023138652918 * gaussianWidth;
		}

		/**
		 * Compute the Gaussian width from the FWHM
		 *
		 * @param fwhm
		 * @return double
		 */
		static public double computeGaussianWidth(final double fwhm) {
			return fwhm / 2.354820045030949382023138652918;
		}

		public Gaussian(final double fano, final double noise, final double eVpEH) {
			super();
			mFano = fano;
			mNoise = noise;
			mEnergyPerEHPair = eVpEH;
		}

		@Override
		public double invert(final double xRayEnergy, final InvertMode mode, final double fraction) {
			if (fraction >= 1.0)
				return 0.0;
			if (fraction <= 0.0)
				return Double.POSITIVE_INFINITY;
			final double g = gaussianWidth(xRayEnergy);
			final double hw = Math.sqrt(-2.0 * Math.log(fraction) * g * g);
			return mode == InvertMode.FULL ? 2.0 * hw : hw;
		}

		/**
		 * Constructs a Gaussian
		 *
		 * @param fwhm The full-width half-max at Mn Ka (5898.7 eV)
		 * @param sdd  true for an SDD, false otherwise
		 */
		public Gaussian(final double fwhm, final double eVpEH) {
			super();
			mFano = 0.12;
			mEnergyPerEHPair = eVpEH;
			final double gRes = computeGaussianWidth(fwhm);
			mNoise = Math.sqrt(Math.pow(gRes / mEnergyPerEHPair, 2.0) - ((5898.7 * mFano) / mEnergyPerEHPair));
		}

		public double gaussianWidth(final double energy) {
			return mEnergyPerEHPair * Math.sqrt((mNoise * mNoise) + ((energy * mFano) / mEnergyPerEHPair));
		}

		@Override
		public double compute(final double xRayEnergy, final double minE, final double maxE) {
			final double w = gaussianWidth(xRayEnergy);
			return 0.5 * Math.abs(Erf.erf((minE - xRayEnergy) / w, (maxE - xRayEnergy) / w));
		}

		@Override
		public double compute(final double xRayEnergy, final double e) {
			final double w = gaussianWidth(xRayEnergy);
			final double xx = (xRayEnergy - e) / w;
			return (SQRT_2PI_INV / w) * Math.exp(-0.5 * xx * xx);
		}

		@Override
		public double getFWHM(final double e) {
			return computeFWHM(gaussianWidth(e));
		}

		@Override
		public int hashCode() {
			return super.hashCode() ^ Objects.hash(mEnergyPerEHPair, mFano, mNoise);
		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (!super.equals(obj))
				return false;
			if (getClass() != obj.getClass())
				return false;
			final Gaussian other = (Gaussian) obj;
			return super.equals(obj) && Objects.equals(mEnergyPerEHPair, other.mEnergyPerEHPair)
					&& Objects.equals(mFano, other.mFano) && Objects.equals(mNoise, other.mNoise);
		}

	}

}
