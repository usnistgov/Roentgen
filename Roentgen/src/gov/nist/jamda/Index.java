package gov.nist.jamda;

import java.util.Arrays;

public class Index implements Comparable<Index> {

	private final int[] mIndices;

	Index(final int[] indices) {
		mIndices = indices.clone();
	}

	public int[] indices() {
		return mIndices.clone();
	}

	public int index(final int dim) {
		return mIndices[dim];
	}

	@Override
	public String toString() {
		return "Index" + Arrays.toString(mIndices);
	}

	public int size() {
		return mIndices.length;
	}

	@Override
	public int compareTo(final Index o) {
		assert mIndices.length == o.mIndices.length;
		int res = 0;
		for (int i = 0; (res == 0) && (i < mIndices.length); ++i) {
			if (o.mIndices.length > i)
				res = Integer.compare(mIndices[i], o.mIndices[i]);
			else
				return 1;
		}
		return res;
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(mIndices);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Index other = (Index) obj;
		return Arrays.equals(mIndices, other.mIndices);
	}
}