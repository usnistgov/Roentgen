package gov.nist.microanalysis.roentgen.DataStore;

import java.util.Objects;
import java.util.UUID;

/**
 * @author nicho
 *
 */
public class UniqueString implements Comparable<UniqueString>, CharSequence {

	private final String mString;
	private final UUID mUUID;

	/**
	 *
	 */
	public UniqueString(final String str) {
		mString = str;
		mUUID = UUID.randomUUID();
	}

	@Override
	public int hashCode() {
		return Objects.hash(mString, mUUID);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UniqueString other = (UniqueString) obj;
		return Objects.equals(mString, other.mString) && Objects.equals(mUUID, other.mUUID);
	}

	@Override
	public int compareTo(final UniqueString o) {
		int res = mString.compareTo(o.mString);
		if (res == 0)
			res = mUUID.compareTo(o.mUUID);
		return res;
	}

	@Override
	public String toString() {
		return mString;
	}

	@Override
	public char charAt(final int index) {
		return mString.charAt(index);
	}

	@Override
	public int length() {
		return mString.length();
	}

	@Override
	public CharSequence subSequence(final int start, final int end) {
		return mString.subSequence(start, end);
	}
}
