package gov.nist.jamda;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.util.Pair;

public class Series<H, J> {

	private final List<H> mData;
	private final Map<J, Integer> mIndex;

	protected Series(final H[] h, final J[] index) {
		mData = new ArrayList<>(Arrays.asList(h));
		mIndex = new HashMap<J, Integer>();
		if (index != null)
			for (int i = 0; i < Math.min(h.length, index.length); ++i)
				if (index[i] != null)
					mIndex.put(index[i], Integer.valueOf(i));
	}

	protected Series(final List<H> h, final List<J> index) {
		mData = new ArrayList<>(h);
		mIndex = new HashMap<J, Integer>();
		if (index != null)
			for (int i = 0; i < Math.min(h.size(), index.size()); ++i)
				if (index.get(i) != null)
					mIndex.put(index.get(i), Integer.valueOf(i));
	}

	protected Series(final Map<J, H> vals) {
		mData = new ArrayList<>();
		mIndex = new HashMap<J, Integer>();
		for (final Map.Entry<J, H> val : vals.entrySet()) {
			mData.add(val.getValue());
			mIndex.put(val.getKey(), mData.size() - 1);
		}
	}

	protected Series(final H[] h) {
		this(h, null);
	}

	public int indexOf(final J obj) {
		final Integer i = mIndex.get(obj);
		return i != null ? i.intValue() : -1;
	}

	public H get(final int i) {
		return mData.get(i);
	}

	public H get(final J obj) {
		return mData.get(indexOf(obj));
	}

	public int size() {
		return mData.size();
	}

	public void put(final J key, final H value) {
		final Integer idx = mIndex.get(key);
		if (idx == null) {
			mData.add(value);
			mIndex.put(key, mData.size() - 1);
		} else
			mData.set(idx.intValue(), value);
	}

	public void add(final H value) {
		mData.add(value);
	}

	public J key(final int idx) {
		for (final Map.Entry<J, Integer> me : mIndex.entrySet())
			if (me.getValue().intValue() == idx)
				return me.getKey();
		return null;
	}

	public Series<H, J> slice(final int minI, final int maxE) {
		final List<H> vals = new ArrayList<>();
		final List<J> keys = new ArrayList<>();
		for (int i = minI; i < maxE; ++i) {
			vals.add(get(i));
			keys.add(key(i));
		}
		return new Series<H, J>(vals, keys);
	}

	public Series<H, J> slice(final J minI, final J maxI) {
		return slice(indexOf(minI), indexOf(maxI) + 1);
	}

	public Series<H, J> sortKeys(final Comparator<J> compareKeys) {
		final SortedSet<J> sortedKeys = new TreeSet<>(compareKeys);
		sortedKeys.addAll(mIndex.keySet());
		final List<H> vals = new ArrayList<>();
		final List<J> keys = new ArrayList<>();
		for (final J key : sortedKeys) {
			vals.add(mData.get(mIndex.get(key).intValue()));
			keys.add(key);
		}
		return new Series<H, J>(vals, keys);
	}

	public Series<H, J> sortValues(final Comparator<H> compareValues) {
		final ArrayList<Pair<H, J>> items = new ArrayList<>();
		for (int i = 0; i < mData.size(); ++i) {
			final H value = mData.get(i);
			final J key = key(i);
			items.add(Pair.create(value, key));
		}
		items.sort((p, q) -> compareValues.compare(p.getFirst(), q.getFirst()));
		final List<H> vals = new ArrayList<>();
		final List<J> keys = new ArrayList<>();
		for (final Pair<H, J> pr : items) {
			vals.add(pr.getFirst());
			keys.add(pr.getSecond());
		}
		return new Series<H, J>(vals, keys);
	}

	public int[] indexes(final Set<J> keys) {
		final int[] res = new int[keys.size()];
		int i = 0;
		for (final J key : keys) {
			res[i] = indexOf(key);
			++i;
		}
		return res;
	}

	public int[] matchesKey(final IPredicate<J> pred) {
		final ArrayList<Integer> res = new ArrayList<>();
		for (final Map.Entry<J, Integer> me : mIndex.entrySet())
			if (pred.evaluate(me.getKey()))
				res.add(me.getValue());
		final int[] resa = new int[res.size()];
		for (int i = 0; i < resa.length; ++i)
			resa[i] = res.get(i);
		return resa;
	}

	public int[] matchesValue(final IPredicate<H> pred) {
		final ArrayList<Integer> res = new ArrayList<>();
		for (int i = 0; i < mData.size(); ++i)
			if (pred.evaluate(mData.get(i)))
				res.add(i);
		final int[] resa = new int[res.size()];
		for (int i = 0; i < resa.length; ++i)
			resa[i] = res.get(i);
		return resa;
	}

	public Series<H, J> subSeries(final Set<J> keys) {
		return subSeries(indexes(keys));
	}

	public Series<H, J> subSeriesByKeys(final IPredicate<J> pred) {
		return subSeries(matchesKey(pred));
	}

	public Series<H, J> subSeriesByValues(final IPredicate<H> pred) {
		return subSeries(matchesValue(pred));
	}

	public Series<H, J> subSeries(final int[] indexes) {
		final List<H> vals = new ArrayList<>();
		final List<J> keys = new ArrayList<>();
		for (final int indexe : indexes) {
			vals.add(mData.get(indexe));
			keys.add(key(indexe));
		}
		return new Series<H, J>(vals, keys);
	}

	public Series<H, J> subSeries(final Collection<J> keys) {
		final List<H> vals = new ArrayList<>();
		final List<J> keyl = new ArrayList<>(keys);
		for (final J key : keyl)
			vals.add(get(key));
		return new Series<H, J>(vals, keyl);
	}
}
