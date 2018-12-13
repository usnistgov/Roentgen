package gov.nist.microanalysis.roentgen.DataStore;

import java.net.URI;

import com.thoughtworks.xstream.annotations.XStreamAlias;
import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * <p>
 * A wrapper around objects that have been curated to allow them to:
 * </p>
 * <ol>
 * <li>Referenced by their location.</li>
 * <li>Kept in a cache by location.</li>
 * <li>Marked when modified and requiring versioning in the archive.</li>
 * </ol>
 * <p>
 * Objects once curated must not be modified.
 * </p>
 * <p>
 * If modification is necessary, a copy should be constructed and the copy
 * modified and re-curated.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
@XStreamAlias("curated")
public class Curated<T> {

	@XStreamAlias("address")
	@XStreamAsAttribute
	private String mAddress;

	@XStreamOmitField
	private final T mObject;

	public Curated(final URI uri, final T obj) {
		this(uri.toString(), obj);
	}

	public Curated(final String address, final T obj) {
		mAddress = address;
		mObject = obj;
	}

	public T getObject() {
		return mObject;
	}

	public String getLocation() {
		return mAddress;
	}

	@Override
	public String toString() {
		return "Curated[" + getObject().toString() + "," + getLocation().toString() + "]";
	}
}
