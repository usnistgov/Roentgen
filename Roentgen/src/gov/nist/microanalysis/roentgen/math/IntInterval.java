package gov.nist.microanalysis.roentgen.math;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import com.google.common.base.Objects;

/**
 * <p>
 * The Interval object encapsulates the idea of an interval on the real number
 * line. It provides methods for determining the intersection
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 199 $
 */
public class IntInterval implements Comparable<IntInterval> {

	static public enum Location {
		/**
		 * Strictly below the Interval
		 */
		BELOW,
		/**
		 * On the lowerLimit() boundary
		 */
		LOWER_BOUNDARY,
		/**
		 * Inside the interval
		 */
		INSIDE,
		/**
		 * On the upperLimit() boundary
		 */
		UPPER_BOUNDARY,
		/**
		 * Strictly above the Interval
		 */
		ABOVE,
		/**
		 * These relationships aren't defined relative to a null interval
		 * (NULL_INTERVAL)
		 */
		UNDEFINED
	};

	private final int mLower;
	private final int mUpper;

	public static final IntInterval NULL_INTERVAL = new IntInterval(Integer.MIN_VALUE, Integer.MIN_VALUE);

	/**
	 * Constructs an Interval encompassing both l1 and l2.
	 */
	public IntInterval(final int l1, final int l2) {
		mLower = Math.min(l1, l2);
		mUpper = Math.max(l1, l2);
	}

	public boolean isNull() {
		return (mLower == Integer.MIN_VALUE) && (mUpper == Integer.MIN_VALUE);
	}

	/**
	 * Get the Location of the value n with respect to this Interval.
	 *
	 * @param val
	 * @return {@link Location}
	 */
	public Location getLocation(final int val) {
		if (isNull())
			return Location.UNDEFINED;
		else if (val < mLower)
			return Location.BELOW;
		else if (val == mLower)
			return Location.LOWER_BOUNDARY;
		else if (val < mUpper)
			return Location.INSIDE;
		else if (val == mUpper)
			return Location.UPPER_BOUNDARY;
		else
			return Location.ABOVE;
	}

	/**
	 * The distance from val to the closest point on the Interval.
	 *
	 * @param val
	 * @return int 0 if inside, otherwise distance to closer of lower or upper
	 *         bound, MAX_VALUE if interval is null.
	 */
	public int distance(final int val) {
		if (isNull())
			return Integer.MAX_VALUE;
		if (val < mLower)
			return mLower - val;
		else if (val > mUpper)
			return val - mUpper;
		else
			return 0;
	}

	/**
	 * Is val inside (or on a boundary) of this Interval?
	 *
	 * @param val The value to test
	 * @return boolean true if inside, false otherwise
	 */
	public boolean contains(final int val) {
		return (!isNull()) && (val >= mLower) && (val <= mUpper);
	}

	/**
	 * Does this Interval intersect the Interval r2? An interval intersects if even
	 * one point touches.
	 *
	 * @param r2
	 * @return boolean
	 */
	public boolean intersects(final IntInterval r2) {
		return (!isNull()) && (mLower <= r2.mUpper) && (mUpper >= r2.mLower);
	}

	/**
	 * Returns the Interval representing the overlap of r1 and r2 if
	 * r1.intersects(r2) or NULL_INTERVAL otherwise.
	 *
	 * @param r1
	 * @param r2
	 * @return Interval or null
	 */
	public static IntInterval intersection(final IntInterval r1, final IntInterval r2) {
		if (r1.isNull() || r2.isNull())
			return IntInterval.NULL_INTERVAL;
		final int min = Math.max(r1.mLower, r2.mLower);
		final int max = Math.min(r1.mUpper, r2.mUpper);
		return min <= max ? new IntInterval(min, max) : NULL_INTERVAL;
	}

	/**
	 * Returns the Interval which fully contains r1 and r2. r1 and r2 don't need to
	 * intersect.
	 *
	 * @param r1 IntInterval (may be NULL_INTERVAL to facilitate building new
	 *           extents)
	 * @param r2 A IntInterval by which to extend r1
	 * @return Interval
	 */
	public static IntInterval extent(final IntInterval r1, final IntInterval r2) {
		assert r1 != null;
		if (r1.isNull())
			return r2;
		else {
			final int min = Math.min(r1.mLower, r2.mLower);
			final int max = Math.max(r1.mUpper, r2.mUpper);
			return new IntInterval(min, max);
		}
	}

	/**
	 * Creates an Interval offset by the specified amount from this one.
	 *
	 * @param offset
	 * @return Interval
	 */
	public IntInterval offset(final int offset) {
		return !isNull() ? new IntInterval(mLower + offset, mUpper + offset) : NULL_INTERVAL;
	}

	/**
	 * Returns the Interval which fully contains r1 and val. r1 and val don't need
	 * to intersect.
	 *
	 * @param r1  Interval (may be NULL_INTERVAL)
	 * @param val A point to include
	 * @return Interval
	 */
	public static IntInterval extent(final IntInterval r1, final int val) {
		if (r1.isNull()) {
			return new IntInterval(val, val);
		} else {
			final int min = Math.min(r1.mLower, val);
			final int max = Math.max(r1.mUpper, val);
			return new IntInterval(min, max);
		}
	}

	/**
	 * Adds r to the set of Range objects. If r intersects a Range in rs then the r
	 * and the intersecting Range object are combined into a single Range.
	 *
	 * @param rs
	 * @param r  An Interval to add(). May be NULL_INTERVAL.
	 * @return Set&lt;Range&gt;
	 */
	public static Set<IntInterval> add(final Set<IntInterval> rs, final IntInterval r) {
		assert r != null;
		if (!r.isNull()) {
			IntInterval tmp = r;
			for (final Iterator<IntInterval> i = rs.iterator(); i.hasNext();) {
				final IntInterval rr = i.next();
				if (rr.intersects(r)) {
					i.remove();
					tmp = extent(rr, tmp);
				}
			}
			rs.add(tmp);
		}
		return rs;
	}

	/**
	 * Returns all intervals in the set of IntIntervals rs that are covered by r.
	 *
	 * @param rs Set&lt;IntInterval&gt;
	 * @param r  IntInterval
	 * @return Set&lt;IntInterval&gt;
	 */
	public static Set<IntInterval> intersection(final Set<IntInterval> rs, final IntInterval r) {
		final Set<IntInterval> res = new HashSet<>();
		for (final IntInterval ii : rs)
			IntInterval.add(res, IntInterval.intersection(ii, r));
		return res;

	}

	/**
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hashCode(mLower, mUpper);
	}

	/**
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof IntInterval))
			return false;
		final IntInterval other = (IntInterval) obj;
		return (mUpper == other.mUpper) && (mLower == other.mLower);
	}

	@Override
	public int compareTo(final IntInterval r2) {
		int res = Integer.compare(mLower, r2.mLower);
		if (res == 0)
			res = Integer.compare(mUpper, r2.mUpper);
		return res;
	}

	@Override
	public String toString() {
		return "[" + Integer.toString(mLower) + "," + Integer.toString(mUpper) + "]";
	}

	/**
	 * Lowest value in contiguous interval (inclusive)
	 *
	 * @return int
	 */
	public int lowerLimit() {
		return mLower;
	}

	/**
	 * Highest value in contiguous interval (inclusive)
	 *
	 * @return int
	 */
	public int upperLimit() {
		return mUpper;
	}

	/**
	 * upper - lower. The width of NULL_INTERVAL = 0
	 *
	 * @return int upperLimit() - lowerLimit() + 1 unless NULL_INTERVAL
	 */
	public int width() {
		return isNull() ? 0 : mUpper - mLower + 1;
	}

	/**
	 * Merges the sets of IntInterval objects into the minimal set of IntInterval
	 * objects that cover the same set of values.
	 *
	 * @param s1
	 * @param s2
	 * @return
	 */
	public Set<IntInterval> merge(final Set<IntInterval> s1, final Set<IntInterval> s2) {
		final Set<IntInterval> res = new HashSet<>(s1);
		for (final IntInterval ii2 : s2)
			add(res, ii2);
		return res;
	}

}
