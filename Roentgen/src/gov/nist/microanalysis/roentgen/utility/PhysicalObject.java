package gov.nist.microanalysis.roentgen.utility;

import java.util.UUID;

import com.google.common.base.Objects;
import com.thoughtworks.xstream.annotations.XStreamAlias;

/**
 * <p>
 * This class represents the base class for classes representing physical
 * objects. Every "thing" should be represented by one-and-only-one
 * PhysicalObject derived instance. PhysicalObjects may have properties and the
 * properties may change over time, but the actual object itself either exists
 * or has been destroyed. For example, a sample is a PhysicalObject. It's owner
 * may change, its description may change but the sample remains fundamentally
 * the same object. If the sample is sub-divided, the original sample is
 * destroyed but now there are a series of new PhysicalObject representing the
 * parts of the original object.
 * </p>
 * <p>
 * Physical objects are labeled with a UUID. There is no way to enforce it but a
 * PhysicalObject should be referenced by the same UUID by everyone
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class PhysicalObject {

	@XStreamAlias("physicalObjectID")
	private final UUID mUUID;
	@XStreamAlias("physicalObjectName")
	private final String mFriendlyName;

	protected PhysicalObject(final String name) {
		mUUID = UUID.randomUUID();
		mFriendlyName = name;
	}

	public UUID getUUID() {
		return mUUID;
	}

	@Override
	public String toString() {
		return mFriendlyName;
	}

	/**
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hashCode(mUUID, mFriendlyName);
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
		if (!(obj instanceof PhysicalObject))
			return false;
		return mUUID.equals(((PhysicalObject) obj).mUUID);
	}

}
