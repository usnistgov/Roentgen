package gov.nist.jamda;

public class Extent {
	
	private Index mMin;
	private Index mMax;
	
	public Extent(Index minI, Index maxI) {
		assert minI.size()==maxI.size();
		int[] min = minI.indices();
		int[] max = maxI.indices();
		for(int i=0;i<min.length;++i) {
			if(min[i]>max[i]) {
				int tmp=min[i];
				min[i]=max[i];
				max[i]=tmp;
			}
		}
		mMin = new Index(min);
		mMax = new Index(max);
	}
	
	public Index getLower() {
		return mMin;
	}
	
	public Index getUpper() {
		return mMax;
	}
}
