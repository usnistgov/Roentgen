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
		public double limit(
				final double val
		) {
			return val;
		}

		@Override
		public String toString() {
			return "None";
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

		public Range(
				final double min, final double max
		) {
			mMin = Math.min(min, max);
			mMax = Math.max(min, max);
		}

		public Range(
				final double nominal, final double lowPct, final double highPct
		) {
			this(nominal * lowPct, nominal * highPct);
		}

		@Override
		public double limit(
				final double val
		) {
			if (val < mMin)
				return mMin;
			else if (val > mMax)
				return mMax;
			else
				return val;
		}

		@Override
		public String toString() {
			return "Range[" + mMin + "," + mMax + "]";
		}
	}

	/**
	 * Takes the value <code>val</code> and ensures that it is within the bounds of
	 * the constraint.
	 *
	 * @param val
	 * @return double val constrained
	 */
	public double limit(
			double val
	);
}