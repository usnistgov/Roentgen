package gov.nist.microanalysis.roentgen.math;

/**
 * An interface for limiting the range of values accessible to a random
 * variable.
 *
 * @author Nicholas
 *
 */
public interface Constraint {

	/**
	 * Unconstrained.
	 *
	 * @author Nicholas
	 *
	 */
	public static class None //
			implements Constraint {
		@Override
		public double limit(final double val) {
			return val;
		}
	}

	/**
	 * Constrains a value to within a [min,max] range. Values less than min are set
	 * to min and values larger than max are set to max.
	 *
	 * @author Nicholas
	 *
	 */
	public static class Range //
			implements Constraint {

		private final double mMin;
		private final double mMax;

		public Range(final double min, final double max) {
			mMin = Math.min(min, max);
			mMax = Math.max(min, max);
		}

		@Override
		public double limit(final double val) {
			if (val < mMin)
				return mMin;
			else if (val > mMax)
				return mMax;
			else
				return val;
		}
	}
	
	
	/**
	 * Takes the value <code>val</code> and ensures that it is within the bounds of
	 * the constraint.
	 *
	 * @param val
	 * @return double val constrained
	 */
	public double limit(double val);
}