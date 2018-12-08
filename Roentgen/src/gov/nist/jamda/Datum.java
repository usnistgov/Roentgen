package gov.nist.jamda;

public class Datum<H> {
	private final H mDatum;

	protected Datum(final H h) {
		mDatum = h;
	}

	public H get() {
		return mDatum;
	}

	@Override
	public String toString() {
		return mDatum != null ? mDatum.toString() : "N/A";
	}

	public static class NumericDatum extends Datum<Double> {

		protected NumericDatum(final double h) {
			super(h);
		}

		@Override
		public String toString() {
			return get().toString();
		}
	}

	public static class TextDatum extends Datum<String> {

		protected TextDatum(final String h) {
			super(h);
			assert h != null;
		}

		@Override
		public String toString() {
			return get();
		}
	}

}