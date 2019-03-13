package gov.nist.microanalysis.roentgen.utility;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Seems to produce a 10 % - 20 % improvement in overall evaluation speed.
 *
 * @author Nicholas W. M. Ritchie
 *
 * @param <H> H must implement hashCode() and equals()
 */
public class FastIndex<H> extends ArrayList<H> {

	private static final long serialVersionUID = 2429500537433349495L;

	private final Map<H, Integer> mIndex;

	public FastIndex() {
		this(Collections.emptyList());
	}

	public FastIndex(final Collection<H> list) {
		super();
		mIndex = new HashMap<>();
		for (final H h : list)
			add(h);
	}

	@Override
	public int indexOf(final Object h) {
		final Integer res = mIndex.get(h);
		assert (res == null) || (get(res.intValue()).equals(h));
		return res == null ? -1 : res.intValue();
	}

	@Override
	public boolean add(final H obj) {
		if (indexOf(obj) == -1) {
			final int len = super.size();
			if (super.add(obj)) {
				assert super.indexOf(obj) == len : obj + " " + super.indexOf(obj) + "!=" + len;
				mIndex.put(obj, len);
				return true;
			} else
				return false;
		} else
			return false;
	}

	@Override
	public boolean addAll(final Collection<? extends H> hs) {
		boolean added = false;
		for (final H h : hs)
			added |= add(h);
		return added;
	}

	@Override
	public H set(final int index, final H element) {
		// It is not possible to support this in a self-consistent manner
		throw new UnsupportedOperationException("This method is not supported");
	}

	@Override
	public void add(final int index, final H element) {
		// It is not possible to support this in a self-consistent manner
		throw new UnsupportedOperationException("This method is not supported");
	}
}