package gov.nist.jamda;

import java.util.Map;
import java.util.TreeMap;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>
 * DataFrame is a mechanism to store multi-dimensional data in a
 * sparse-compatible format. Unassigned values will by default have the value
 * NULL_DATA with representation "N/A".
 * </p>
 *
 * <p>
 * Data items are index
 *
 *
 * @author Nicholas
 *
 */
public class DataFrame {

	public static final class NullDatum {

		@Override
		public String toString() {
			return "---";
		}
	};

	public final static Object NULL_DATUM = new NullDatum();

	private final Map<Index, Object> mData;
	@SuppressWarnings("rawtypes")
	private final Dimension[] mDimensions;
	private Object mDefaultDatum = NULL_DATUM;

	public DataFrame(final Dimension<?, ?>... dims) {
		mData = new TreeMap<>();
		mDimensions = new Dimension[dims.length];
		for (int i = 0; i < mDimensions.length; ++i)
			mDimensions[i] = dims[i];
	}

	public Dimension<?, ?> getDimension(final int i) {
		return mDimensions[i];
	}

	public int getDimension() {
		return mDimensions.length;
	}

	public int find(final int dim, final Object label) {
		return mDimensions[dim].find(label);
	}

	@SuppressWarnings("unchecked")
	public int findOrAdd(final int dim, final Object label) //
			throws ArgumentException {
		return mDimensions[dim].findOrAdd(label);
	}

	public Object[] getLabels(final Index idx) {
		return getLabels(idx.indices());
	}

	private Object[] getLabels(final int[] indices) {
		final Object[] res = new Object[indices.length];
		for (int i = 0; i < indices.length; ++i)
			res[i] = mDimensions[i].getLabel(indices[i]);
		return res;
	}

	public void put(final Object[] labels, final IGenerate gen) {
		put(find(labels), gen);
	}

	public void put(final Index index, final IGenerate gen) {
		put(index, gen.value(getLabels(index)));
	}

	private void generate(final int dim, final int[] index, final int[] lower, final int[] upper, final IGenerate gen) {
		assert dim < lower.length;
		final int[] tmp = lower.clone();
		if (dim == lower.length - 1) {
			for (int i = lower[dim]; i < upper[dim]; ++i) {
				tmp[i] = i;
				final Index ti = new Index(tmp);
				put(ti, gen.value(getLabels(ti)));
			}
		} else {
			for (int i = lower[dim]; i < upper[dim]; ++i) {
				tmp[dim] = i;
				generate(dim + 1, tmp, lower, upper, gen);
			}
		}
	}

	public void put(final Extent extent, final IGenerate gen) {
		final int[] lower = extent.getLower().indices();
		final int[] upper = extent.getUpper().indices();
		assert lower.length == upper.length;
		final int[] tmp = lower.clone();
		if (lower.length == 1) {
			for (int i = lower[0]; i < upper[0]; ++i) {
				tmp[i] = i;
				final Index ti = new Index(tmp);
				put(ti, gen.value(getLabels(ti)));
			}
		} else {
			for (int i = lower[0]; i < upper[0]; ++i) {
				tmp[0] = i;
				generate(1, tmp, lower, upper, gen);
			}
		}
	}

	public void put(final Extent extent, final Object datum) {
		put(extent, (labels) -> datum);
	}

	public void clear(final Extent extent) {
		put(extent, NULL_DATUM);
	}

	public Index find(final Object[] tags) {
		assert tags.length == mDimensions.length;
		final int[] idx = new int[mDimensions.length];
		for (int i = 0; i < mDimensions.length; ++i)
			idx[i] = mDimensions[i].find(tags[i]);
		return new Index(idx);
	}

	@SuppressWarnings("unchecked")
	public void put(final Object[] tags, final Object datum) //
			throws ArgumentException {
		assert tags.length == mDimensions.length;
		final int[] idx = new int[mDimensions.length];
		for (int i = 0; i < mDimensions.length; ++i)
			idx[i] = mDimensions[i].findOrAdd(tags[i]);
		mData.put(new Index(idx), datum);
	}

	public void setDefaultDatum(final Object defaultDatum) {
		mDefaultDatum = defaultDatum;
	}

	public Object get(final Index idx) {
		final Object res = mData.get(idx);
		return res != null ? res : mDefaultDatum;
	}

	public Object put(final Index idx, final Object datum) {
		return mData.put(idx, datum);
	}
}
