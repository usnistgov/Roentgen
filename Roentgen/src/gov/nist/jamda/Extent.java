package gov.nist.jamda;

public class Extent {

	private final Index mMin;
	private final Index mMax;

	public Extent(final Index minI, final Index maxI) {
		assert minI.size() == maxI.size();
		final int[] min = minI.indices();
		final int[] max = maxI.indices();
		for (int i = 0; i < min.length; ++i) {
			if (min[i] > max[i]) {
				final int tmp = min[i];
				min[i] = max[i];
				max[i] = tmp;
			}
		}
		mMin = new Index(min);
		mMax = new Index(max);
	}

	public Index getLower() {
		return mMin;
	}

	public Index getUpper() {
		return mMax;
	}

	/**
	 * Returns the number of dimensions enclosed by this Extent. This number will be
	 * less than or equal to the length of the size() of the Index objects.
	 *
	 * <ul>
	 * <li>A 1D extent is a portion of a column</li>
	 * <li>A 2D extent is a portion of a table</li>
	 * <li>A 3D extent is a portion of a 3D block</li>
	 * </ul>
	 *
	 * @return int
	 */
	public int getDimensionality() {
		int dim = 0;
		for (int i = 0; i < mMin.size(); ++i)
			if (mMin.index(i) != mMax.index(i))
				++dim;
		return dim;
	}

}
