package gov.nist.jamda;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import gov.nist.microanalysis.roentgen.ArgumentException;

public class Dimension<I, J> {

	final private I mName;
	final private Class<?> mJClass;
	final private List<J> mLabels = new ArrayList<>();

	public Dimension(final I name, final J[] labels) throws ArgumentException {
		mName = name;
		mJClass = labels[0].getClass();
		for(int i=0;i<labels.length;++i)
			if (!mJClass.isInstance(labels[i]))
				throw new ArgumentException(labels[i] + " is not derived from " + mJClass);
		mLabels.addAll(Arrays.asList(labels));
	}
	
	public Dimension(final I name, Class<?> clss, final J[] labels) throws ArgumentException {
		mName = name;
		mJClass = clss;
		for(int i=0;i<labels.length;++i)
			if (!mJClass.isInstance(labels[i]))
				throw new ArgumentException(labels[i] + " is not derived from " + mJClass);
		mLabels.addAll(Arrays.asList(labels));
	}

	public Dimension(final I name, Class<?> clss) {
		mName = name;
		mJClass = clss;
	}

	@Override
	public String toString() {
		return mName.toString();
	}

	public void add(final J label) {
		if (!mLabels.contains(label))
			mLabels.add(label);
	}

	public int find(final Object label) {
		return mLabels.indexOf(label);
	}

	public J getLabel(int i) {
		return mLabels.get(i);
	}

	public int findOrAdd(final J label) throws ArgumentException {
		int res = mLabels.indexOf(label);
		if (res == -1) {
			assert mJClass.getClass().isInstance(label);
			if (!mJClass.isInstance(label))
				throw new ArgumentException(label + " is not derived from " + mJClass);
			mLabels.add(label);
			res = mLabels.size() - 1;
		}
		return res;
	}
}