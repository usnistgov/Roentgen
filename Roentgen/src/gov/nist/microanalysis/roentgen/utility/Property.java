package gov.nist.microanalysis.roentgen.utility;

import java.util.Objects;
import java.util.UUID;

import com.thoughtworks.xstream.annotations.XStreamAlias;

/**
 * <p>
 * This is the base class for objects implementing properties or objects that
 * describe aspects of an object.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
public class Property {

	@XStreamAlias("propertyID")
	private final UUID mUUID;
	@XStreamAlias("propertyName")
	private final String mFriendlyName;

	protected Property(final String name) {
		mUUID = UUID.randomUUID();
		mFriendlyName = name;
	}

	public UUID getUUID() {
		assert mUUID != null;
		return mUUID;
	}

	@Override
	public String toString() {
		assert mFriendlyName != null;
		return mFriendlyName;
	}

	/**
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hash(mUUID, mFriendlyName);
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
		if (!(obj instanceof Property))
			return false;
		return mUUID.equals(((Property) obj).mUUID);
	}

}
