package gov.nist.juncertainty.utility;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A {@link List} with an associated index to speed the indexOf(...) function.
 * Does not allow adding an item twice
 * 
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

	public FastIndex(final List<? extends H> outputs) {
		super(outputs.size() == 0 ? 256 : outputs.size());
		mIndex = new HashMap<>(outputs.size() == 0 ? 256 : outputs.size());
		for (H h : outputs)
			add(h);
	}

	public FastIndex<H> copy() {
		return new FastIndex<H>(this);
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
			throw new UnsupportedOperationException("The item " + obj + " already exists in the FastIndex.");
	}

	public boolean addIfMissing(final H obj) {
		if (indexOf(obj) == -1) {
			final int len = super.size();
			if (super.add(obj)) {
				assert super.indexOf(obj) == len : obj + " " + super.indexOf(obj) + "!=" + len;
				mIndex.put(obj, len);
				return true;
			} else {
				assert false : "Failed to add " + obj;
				return false;
			}
		} else
			return false;
	}

	public boolean addMissing(final Collection<? extends H> hs) {
		boolean added = false;
		for (final H h : hs)
			added |= addIfMissing(h);
		return added;
	}

	@Override
	public boolean addAll(final Collection<? extends H> hs) {
		boolean added = false;
		for (final H h : hs)
			added |= add(h);
		return added;
	}

	@Override
	public H remove(int i) {
		throw new UnsupportedOperationException("This method is not supported");
	}

	@Override
	public boolean remove(Object o) {
		throw new UnsupportedOperationException("This method is not supported");
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