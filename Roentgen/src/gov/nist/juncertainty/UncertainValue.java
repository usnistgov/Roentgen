package gov.nist.juncertainty;

import java.util.Collection;
import java.util.Objects;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * UncertainValue represents a number and an associated uncertainty (which may
 * be zero.) UncertainValue is derived from {@link Number}.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class UncertainValue extends Number //
		implements Comparable<UncertainValue>, IToHTML {

	public static final UncertainValue ONE = new UncertainValue(1.0);
	public static final UncertainValue ZERO = new UncertainValue(0.0);
	public static final UncertainValue NaN = new UncertainValue(Double.NaN);
	public static final UncertainValue POSITIVE_INFINITY = new UncertainValue(Double.POSITIVE_INFINITY);
	public static final UncertainValue NEGATIVE_INFINITY = new UncertainValue(Double.NEGATIVE_INFINITY);

	private static final long serialVersionUID = -7284125207920225793L;

	public static double fractionalUncertainty(
			final Number n
	) {
		return uncertainty(n) / n.doubleValue();
	}

	static final public boolean isSpecialNumber(
			final Number n
	) {
		return n instanceof UncertainValue;
	}

	public static boolean isUncertain(
			final Number n
	) {
		return (n instanceof UncertainValue) && ((UncertainValue) n).isUncertain();
	}

	public static double mean(
			final Number n
	) {
		return n.doubleValue();
	}

	public static UncertainValue normal(
			final double v
	) {
		return new UncertainValue(v, Math.sqrt(v));
	}

	/**
	 * Parses a string on the form "Double ± Double" or "Double" returning the
	 * result as an UncertainValue. The "±" character can be replaced with "+-" or
	 * "-+".
	 *
	 * @param str A string containing a text representation of an {@link UncertainValue}
	 * @return UncertainValue
	 */
	public static UncertainValue parse(
			final String str
	) {
		final String[] pms = { "\u00B1", "+-", "-+" };
		for (final String pm : pms) {
			final int idx = str.indexOf(pm);
			if (idx != -1) {
				final double value = Double.parseDouble(str.substring(0, idx).trim());
				final double sigma = Double.parseDouble(str.substring(idx + pm.length()).trim());
				return new UncertainValue(value, sigma);
			}
		}
		final double value = Double.parseDouble(str.trim());
		return new UncertainValue(value);
	}

	public static UncertainValue toRadians(
			final double degrees, final double ddegrees
	) {
		return new UncertainValue(Math.toRadians(degrees), Math.toRadians(ddegrees));
	}

	public static double uncertainty(
			final Number n
	) {
		return n instanceof UncertainValue ? ((UncertainValue) n).uncertainty() : 0.0;
	}
	
	public static double variance(
			final Number n
	) {
		return n instanceof UncertainValue ? ((UncertainValue) n).variance() : 0.0;
	}


	static final public Number unwrap(
			final Number n
	) {
		if (n instanceof UncertainValue) {
			final UncertainValue uv = (UncertainValue) n;
			if (!uv.isUncertain())
				return Double.valueOf(n.doubleValue());
		}
		return n;
	}

	public static UncertainValue valueOf(
			final double val, final double unc
	) {
		return new UncertainValue(val, unc);
	}

	public static UncertainValue valueOf(
			final double val
	) {
		return new UncertainValue(val);
	}

	/**
	 * <p>
	 * Computes the variance weighted mean - the maximum likelyhood estimator of the
	 * mean under the assumption that the samples are independent and normally
	 * distributed.
	 * </p>
	 *
	 * @param cuv A {@link Collection} of UncertainValue objects
	 * @return UncertainValue
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	static public Number weightedMean(
			final Collection<? extends UncertainValue> cuv
	) throws ArgumentException {
		double varSum = 0.0, sum = 0.0;
		for (final UncertainValue uv : cuv) {
			final double ivar = (isSpecialNumber(uv) ? 1.0 / uv.variance() : Double.NaN);
			if (Double.isNaN(ivar))
				throw new ArgumentException(
						"Unable to compute the weighted mean when one or more datapoints have zero uncertainty.");
			varSum += ivar;
			sum += ivar * uv.doubleValue();
		}
		final double iVarSum = 1.0 / varSum;
		return Double.isNaN(iVarSum) ? UncertainValue.NaN : new UncertainValue(sum / varSum, Math.sqrt(1.0 / varSum));
	}

	final Double mValue;

	final double mSigma;

	/**
	 * Constructs an UncertainValue with uncertainty equal to 0.0.
	 * 
	 * @param value The value
	 */
	public UncertainValue(
			final double value
	) {
		this(value, 0.0);
	}

	/**
	 * <p>
	 * Constructs an UncertainValue equal to "value ±  sigma".
	 * </p>
	 * 
	 * @param value The value
	 * @param sigma The uncertainty (should be &gt;= 0.0)
	 * 
	 */
	public UncertainValue(
			final double value, //
			final double sigma
	) {
		mValue = Double.valueOf(value);
		mSigma = Math.abs(sigma);
	}

	/**
	 * Constructs an UncertainValue from an instance of any Number derived class.
	 * If the Number is an {@link UncertainValue} then the result is a copy.
	 * 
	 * @param n Number
	 */
	public UncertainValue(
			final Number n
	) {
		this(n.doubleValue(), uncertainty(n));
	}

	/**
	 * First compares the values and then compares the uncertainties using
	 * Double.compare(...).
	 *
	 * @param uv1 The {@link UncertainValue} against which to compare
	 * @return int
	 */
	@Override
	public int compareTo(
			final UncertainValue uv1
	) {
		int res = mValue.compareTo(uv1.mValue);
		if (res == 0)
			res = Double.compare(mSigma, uv1.mSigma);
		return res;
	}

	@Override
	public double doubleValue() {
		return mValue.doubleValue();
	}

	@Override
	public boolean equals(
			final Object obj
	) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValue other = (UncertainValue) obj;
		return Double.doubleToLongBits(mSigma) == Double.doubleToLongBits(other.mSigma)
				&& Objects.equals(mValue, other.mValue);
	}

	@Override
	public float floatValue() {
		return mValue.floatValue();
	}

	public String format(
			final BasicNumberFormat bnf
	) {
		if (mSigma == 0.0)
			return bnf.format(mValue);
		else
			return bnf.format(mValue) + "\u00B1" + bnf.format(mSigma);
	}

	public String formatLong(
			final BasicNumberFormat bnf
	) {
		return bnf.format(mValue) + "\u00B1" + bnf.format(mSigma);
	}

	/**
	 * Returns the fractional uncertainty.
	 * 
	 * @return sigma/value
	 */
	public double fractionalUncertainty() {
		return mSigma / mValue;
	}

	@Override
	public int hashCode() {
		return Objects.hash(mSigma, mValue);
	}

	@Override
	public int intValue() {
		return mValue.intValue();
	}

	/**
	 * @return true if the uncertainty component is non-zero
	 */
	public boolean isUncertain() {
		return mSigma > 0.0;
	}

	@Override
	public long longValue() {
		return mValue.longValue();
	}

	public UncertainValue multiply(
			final double k
	) {
		return new UncertainValue(k * mValue, k * mSigma);
	}

	@Override
	public String toHTML(
			final Mode mode
	) {
		return toHTML(mode, new BasicNumberFormat());
	}

	public String toHTML(
			final Mode mode, final BasicNumberFormat bnf
	) {
		switch (mode) {
		case TERSE:
			if (uncertainty() != 0)
				return bnf.formatHTML(mValue) + "&#177;" + bnf.formatHTML(uncertainty());
			else
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

	public String toString() {
		return mValue + "\u00B1" + mSigma;
	}
}
