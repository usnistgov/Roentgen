package gov.nist.microanalysis.roentgen.utility;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
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

	public FastIndex(final List<H> list) {
		super(list);
		mIndex = new HashMap<>();
		for (int i = 0; i < list.size(); ++i)
			mIndex.put(list.get(i), i);
	}

	@Override
	public int indexOf(final Object h) {
		final Integer res = mIndex.get(h);
		assert (res == null) || (get(res.intValue()).equals(h));
		return res == null ? -1 : res.intValue();
	}

	public boolean add(H obj) {
		throw new UnsupportedOperationException("add is not supported by FastIndex.");
	}

	public boolean addAll(Collection<? extends H> hs) {
		throw new UnsupportedOperationException("addAll is not supported by FastIndex.");
	}
	
	public H set(int index, H element) {
		throw new UnsupportedOperationException("set is not supported by FastIndex.");
	}
	
	public void add(int index, H element) {
		throw new UnsupportedOperationException("add is not supported by FastIndex.");
	}
}