package gov.nist.jamda;

public class Datum<H> {
	private final H mDatum;

	protected Datum(H h) {
		mDatum = h;
	}

	public H get() {
		return mDatum;
	}

	public String toString() {
		return mDatum != null ? mDatum.toString() : "N/A";
	}
	
	public static class NumericDatum extends Datum<Double> {

		protected NumericDatum(double h) {
			super(h);
		}

		public String toString() {
			return get().toString();
		}
	}

	public static class TextDatum extends Datum<String> {

		protected TextDatum(String h) {
			super(h);
			assert h != null;
		}

		public String toString() {
			return get();
		}
	}

}