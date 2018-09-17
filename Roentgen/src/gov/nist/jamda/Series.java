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

	protected Series(H[] h, J[] index) {
		mData = new ArrayList<>(Arrays.asList(h));
		mIndex = new HashMap<J, Integer>();
		if (index != null)
			for (int i = 0; i < Math.min(h.length, index.length); ++i)
				if (index[i] != null)
					mIndex.put(index[i], Integer.valueOf(i));
	}

	protected Series(List<H> h, List<J> index) {
		mData = new ArrayList<>(h);
		mIndex = new HashMap<J, Integer>();
		if (index != null)
			for (int i = 0; i < Math.min(h.size(), index.size()); ++i)
				if (index.get(i) != null)
					mIndex.put(index.get(i), Integer.valueOf(i));
	}

	protected Series(Map<J, H> vals) {
		mData = new ArrayList<>();
		mIndex = new HashMap<J, Integer>();
		for (Map.Entry<J, H> val : vals.entrySet()) {
			mData.add(val.getValue());
			mIndex.put(val.getKey(), mData.size() - 1);
		}
	}

	protected Series(H[] h) {
		this(h, null);
	}

	public int indexOf(J obj) {
		Integer i = mIndex.get(obj);
		return i != null ? i.intValue() : -1;
	}

	public H get(int i) {
		return mData.get(i);
	}

	public H get(J obj) {
		return mData.get(indexOf(obj));
	}

	public int size() {
		return mData.size();
	}

	public void put(J key, H value) {
		Integer idx = mIndex.get(key);
		if (idx == null) {
			mData.add(value);
			mIndex.put(key, mData.size() - 1);
		} else
			mData.set(idx.intValue(), value);
	}

	public void add(H value) {
		mData.add(value);
	}
	
	public J key(int idx) {
		for(Map.Entry<J, Integer> me: mIndex.entrySet())
			if(me.getValue().intValue()==idx)
				return me.getKey();
		return null;
	}
	
	public Series<H,J> slice(int minI, int maxE){
		List<H> vals = new ArrayList<>();
		List<J> keys = new ArrayList<>();
		for(int i=minI;i<maxE;++i) {
			vals.add(get(i));
			keys.add(key(i));
		}
		return new  Series<H,J>(vals, keys);
	}
	
	public Series<H,J> slice(J minI, J maxI){
		return slice(indexOf(minI),indexOf(maxI)+1);
	}
	
	public Series<H,J> sortKeys(Comparator<J> compareKeys){
		SortedSet<J> sortedKeys = new TreeSet<>(compareKeys);
		sortedKeys.addAll(mIndex.keySet());
		List<H> vals = new ArrayList<>();
		List<J> keys = new ArrayList<>();
		for(J key : sortedKeys) {
			vals.add(mData.get(mIndex.get(key).intValue()));
			keys.add(key);
		}
		return new Series<H,J>(vals,keys);
	}
	
	public Series<H,J> sortValues(Comparator<H> compareValues){
		ArrayList<Pair<H,J>> items = new ArrayList<>();
		for(int i=0;i<mData.size();++i) {
			H value = mData.get(i);
			J key = key(i);
			items.add(Pair.create(value, key));
		}
		items.sort(( p, q ) -> compareValues.compare(p.getFirst(), q.getFirst()) );
		List<H> vals = new ArrayList<>();
		List<J> keys = new ArrayList<>();
		for(Pair<H,J> pr : items) {
			vals.add(pr.getFirst());
			keys.add(pr.getSecond());
		}
		return new Series<H,J>(vals,keys);
	}
	
	public int[] indexes(Set<J> keys) {
		int[] res = new int[keys.size()];
		int i=0;
		for(J key : keys) {
			res[i]=indexOf(key);
			++i;
		}
		return res;
	}
	
	public int[] matchesKey(IPredicate<J> pred) {
		ArrayList<Integer> res = new ArrayList<>();
		for(Map.Entry<J, Integer>  me : mIndex.entrySet())
			if(pred.evaluate(me.getKey()))
				res.add(me.getValue());
		int[] resa = new int[res.size()];
		for(int i=0;i<resa.length;++i)
			resa[i]=res.get(i);
		return resa;
	}
	
	public int[] matchesValue(IPredicate<H> pred) {
		ArrayList<Integer> res = new ArrayList<>();
		for(int i=0;i<mData.size();++i)
			if(pred.evaluate(mData.get(i)))
				res.add(i);
		int[] resa = new int[res.size()];
		for(int i=0;i<resa.length;++i)
			resa[i]=res.get(i);
		return resa;
	}

	
	
	public Series<H,J> subSeries(Set<J> keys) {
		return subSeries(indexes(keys));
	}
	
	
	public Series<H,J> subSeriesByKeys(IPredicate<J> pred) {
		return subSeries(matchesKey(pred));
	}
	
	public Series<H,J> subSeriesByValues(IPredicate<H> pred) {
		return subSeries(matchesValue(pred));
	}
	
	
	public Series<H,J> subSeries(int[] indexes) {
		List<H> vals = new ArrayList<>();
		List<J> keys = new ArrayList<>();
		for(int i=0;i<indexes.length;++i) {
			vals.add(mData.get(indexes[i]));
			keys.add(key(indexes[i]));
		}
		return new Series<H,J>(vals,keys);
	}
	
	
	public Series<H,J> subSeries(Collection<J> keys) {
		List<H> vals = new ArrayList<>();
		List<J> keyl = new ArrayList<>(keys);
		for(J key : keyl)
			vals.add(get(key));
		return new Series<H,J>(vals, keyl);
	}
}
