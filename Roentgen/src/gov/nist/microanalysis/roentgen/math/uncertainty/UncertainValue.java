package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collection;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class UncertainValue
		extends Number //
		implements Comparable<UncertainValue>, IToHTML {

	public static final UncertainValue ONE = new UncertainValue(1.0);
	public static final UncertainValue ZERO = new UncertainValue(0.0);
	public static final UncertainValue NaN = new UncertainValue(Double.NaN);
	public static final UncertainValue POSITIVE_INFINITY = new UncertainValue(Double.POSITIVE_INFINITY);
	public static final UncertainValue NEGATIVE_INFINITY = new UncertainValue(Double.NEGATIVE_INFINITY);

	private static final long serialVersionUID = -7284125207920225793L;
	final Double mValue;
	final double mSigma;

	public UncertainValue(final double value) {
		this(value, 0.0);
	}

	public UncertainValue(final Number n) {
		this(n.doubleValue(), n instanceof UncertainValue ? ((UncertainValue) n).uncertainty() : 0.0);
	}

	/**
	 *
	 */
	public UncertainValue(final double value, final double sigma) {
		mValue = Double.valueOf(value);
		mSigma = sigma;
	}

	@Override
	public double doubleValue() {
		return mValue.doubleValue();
	}

	@Override
	public float floatValue() {
		return mValue.floatValue();
	}

	@Override
	public int intValue() {
		return mValue.intValue();
	}

	@Override
	public long longValue() {
		return mValue.longValue();
	}

	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mode, new BasicNumberFormat());
	}

	public String toHTML(final Mode mode, final BasicNumberFormat bnf) {
		switch (mode) {
		case TERSE:
			return bnf.formatHTML(mValue);
		case NORMAL:
		case VERBOSE:
		default:
			return bnf.formatHTML(mValue) + "&nbsp;&#177;&nbsp;" + bnf.formatHTML(uncertainty());
		}
	}

	/**
	 * Returns the one-&sigma; uncertainty.
	 *
	 * @return double
	 */
	public double uncertainty() {
		return mSigma;
	}

	/**
	 * Returns the variance = sigma<sup>2</sup>
	 *
	 * @return variance
	 */
	public double variance() {
		return mSigma * mSigma;
	}

	public static double mean(final Number n) {
		return n.doubleValue();
	}

	public static UncertainValue normal(final double v) {
		return new UncertainValue(v, Math.sqrt(v));
	}

	public static UncertainValue toRadians(final double degrees, final double ddegrees) {
		return new UncertainValue(Math.toRadians(degrees), Math.toRadians(ddegrees));
	}

	public static double uncertainty(final Number n) {
		return n instanceof UncertainValue ? ((UncertainValue) n).uncertainty() : 0.0;
	}

	public static double fractionalUncertainty(final Number n) {
		return uncertainty(n) / n.doubleValue();
	}

	public UncertainValue multiply(double k) {
		return new UncertainValue(k * mValue, k * mSigma);
	}

	/**
	 * <p>
	 * Computes the variance weighted mean - the maximum likelyhood estimator of the
	 * mean under the assumption that the samples are independent and normally
	 * distributed.
	 * </p>
	 *
	 * @param cuv
	 * @return UncertainValue
	 */
	static public Number weightedMean(final Collection<? extends UncertainValue> cuv) throws Exception {
		double varSum = 0.0, sum = 0.0;
		for (final UncertainValue uv : cuv) {
			final double ivar = (isSpecialNumber(uv) ? 1.0 / uv.variance() : Double.NaN);
			if (Double.isNaN(ivar))
				throw new Exception(
						"Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
			varSum += ivar;
			sum += ivar * uv.doubleValue();
		}
		final double iVarSum = 1.0 / varSum;
		return Double.isNaN(iVarSum) ? UncertainValue.NaN : new UncertainValue(sum / varSum, Math.sqrt(1.0 / varSum));
	}

	static final private boolean isSpecialNumber(final Number n) {
		return n instanceof UncertainValue;
	}

	static final public Number unwrap(final Number n) {
		if (n instanceof UncertainValue) {
			final UncertainValue uv = (UncertainValue) n;
			if (!uv.isUncertain())
				return Double.valueOf(n.doubleValue());
		}
		return n;
	}

	@Override
	public int compareTo(final UncertainValue arg0) {
		int res = Double.compare(mValue, arg0.mValue);
		if (res == 0)
			res = Double.compare(mSigma, arg0.mSigma);
		return res;
	}

	public String formatLong(BasicNumberFormat bnf) {
		return bnf.format(mValue) + "\u00B1" + bnf.format(mSigma);
	}

	public String format(BasicNumberFormat bnf) {
		return bnf.format(mValue) + "\u00B1" + bnf.format(mSigma);
	}

	public double fractionalUncertainty() {
		return mValue / mSigma;
	}

	public static boolean isUncertain(Number n) {
		return (n instanceof UncertainValue) && ((UncertainValue)n).isUncertain();
	}

	public boolean isUncertain() {
		return mSigma > 0.0;
	}

}
