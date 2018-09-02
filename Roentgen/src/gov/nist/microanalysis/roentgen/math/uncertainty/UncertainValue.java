package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.thoughtworks.xstream.annotations.XStreamAlias;
import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamImplicit;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * The UncertainValue class implements a class for handling values with zero or
 * more component normally distributed uncertainties. The class implements a
 * number of static methods for performing basic mathematical operations on
 * numbers while propagating the component uncertainties.
 * </p>
 * <p>
 * Each uncertain value is represented by a value and a series of named
 * component uncertainties. The component names identify the source of the
 * uncertainty. Uncertainty components associated with the same name are assumed
 * to be 100% correlated (r=1.0) and are accumulated each time an operation is
 * performed. Uncertainties associated with different names are accumulated
 * separately. The named uncertainties can be reduced to a single uncertainty as
 * a final step. This step may either assume the components are independent or
 * correlated using the Correlation imbedded class.
 * </p>
 * <p>
 * Each component of uncertainty propagates through the mathematical operations
 * as though it was the only source of uncertainty.
 * </p>
 * <p>
 * XML looks like:
 * </p>
 * <p>
 * &lt;UncertainValue&gt;</br>
 * &nbsp;&nbsp;&nbsp;&lt;value&gt;1.0&lt;/value&gt;</br>
 * &nbsp;&nbsp;&nbsp;&lt;sigma name="dx1"&gt;0.1&lt;/sigma&gt;</br>
 * &nbsp;&nbsp;&nbsp;&lt;sigma name="dx2"&gt;0.3&lt;/sigma&gt;</br>
 * &lt;/UncertainValue&gt;</br>
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author nritchie
 * @version $Rev: 199 $
 */
@XStreamAlias("UncertainValue")
final public class UncertainValue extends Number implements Comparable<UncertainValue>, IToHTML {

	private static final long serialVersionUID = 119495064970078787L;

	/**
	 * The value
	 */
	@XStreamAlias("sigma")
	// XML looks like <sigma name="dx">0.1</sigma>
	private class Sigma {
		@XStreamAsAttribute
		@XStreamAlias("name")
		private final Object mName;
		@XStreamAlias("value")
		private final double mValue;

		private Sigma(final Object name, final double val) {
			mName = name;
			mValue = Math.abs(val);
		}
	}

	@XStreamAlias("value")
	private final Double mValue;

	/**
	 * A map of one-sigma width uncertainty components.
	 */
	@XStreamImplicit
	@XStreamAlias("sigmas")
	private final List<Sigma> mSigmas = new ArrayList<>();;

	public static final UncertainValue ONE = new UncertainValue(1.0);
	public static final UncertainValue ZERO = new UncertainValue(0.0);
	public static final UncertainValue NaN = new UncertainValue(Double.NaN);
	public static final UncertainValue POSITIVE_INFINITY = new UncertainValue(Double.POSITIVE_INFINITY);
	public static final UncertainValue NEGATIVE_INFINITY = new UncertainValue(Double.NEGATIVE_INFINITY);

	private static final Random sRandom = new Random(System.currentTimeMillis());

	static final private boolean isSpecialNumber(final Number n) {
		return n instanceof UncertainValue;
	}

	static final private Number unwrap(final Number n) {
		if (n instanceof UncertainValue) {
			final UncertainValue uv = (UncertainValue) n;
			if (!uv.isUncertain())
				return Double.valueOf(n.doubleValue());
		}
		return n;
	}

	/**
	 * Constructs a UncertainValue with value <code>v</code> and uncertainty
	 * <code>dv</code>. The uncertainty is given a default name "DEFAULTXXX"= where
	 * XXX is an auto-incrementing index to ensure every unnamed uncertainty is
	 * considered independent.
	 *
	 * @param v  The value
	 * @param dv The associated uncertainty
	 */
	public UncertainValue(final double v, final double dv) {
		this(v, "dx", dv);
	}

	private UncertainValue(final double v, final List<Sigma> sigmas) {
		mValue = v;
		if (sigmas != null)
			mSigmas.addAll(sigmas);
	}

	public UncertainValue(final Number n) {
		this(n.doubleValue(), n instanceof UncertainValue ? ((UncertainValue) n).mSigmas : null);
	}

	/**
	 * Constructs a UncertainValue with value <code>v</code> and no uncertainty.
	 *
	 * @param v The value
	 */
	public UncertainValue(final double v) {
		this(v, "dv", 0.0);
	}

	/**
	 * Constructs a UncertainValue
	 *
	 * @param v      double
	 * @param source String The name of the source of the uncertainty
	 * @param dv     The associate uncertainty
	 */
	public UncertainValue(final double v, final Object source, final double dv) {
		mValue = v;
		assignComponent(source, dv);
	}

	public static UncertainValue normal(final double v, final Object src) {
		return new UncertainValue(v, src, Math.sqrt(v));
	}

	/**
	 * Constructs a UncertainValue with value <code>v</code> and uncertainties
	 * <code>sigmas</code>.
	 *
	 * @param v      The value
	 * @param sigmas The associated uncertainties
	 */
	public UncertainValue(final double v, final Map<Object, Double> sigmas) {
		mValue = v;
		if (sigmas != null)
			for (final Map.Entry<Object, Double> me : sigmas.entrySet())
				assignComponent(me.getKey(), me.getValue());
	}

	/**
	 * @see java.lang.Number#doubleValue()
	 */
	@Override
	public double doubleValue() {
		return mValue.doubleValue();
	}

	/**
	 * @see java.lang.Number#floatValue()
	 */
	@Override
	public float floatValue() {
		return mValue.floatValue();
	}

	/**
	 * @see java.lang.Number#intValue()
	 */
	@Override
	public int intValue() {
		return mValue.intValue();
	}

	/**
	 * @see java.lang.Number#longValue()
	 */
	@Override
	public long longValue() {
		return mValue.longValue();
	}

	/**
	 * Assigns the magnitude (abs(sigma)) of the specified source of uncertainty.
	 * Any previous value assigned to this source is replaced. If sigma is zero then
	 * any earlier value is erased.
	 *
	 * @param name  The name of the source of the uncertainty
	 * @param sigma The magnitude of the uncertainty
	 */
	public void assignComponent(final Object name, final double sigma) {
		if (sigma != 0.0) {
			final int idx = indexOf(name);
			if (idx >= 0)
				mSigmas.set(idx, new Sigma(name, sigma));
			else
				mSigmas.add(new Sigma(name, Math.abs(sigma)));
		} else
			mSigmas.remove(name);
	}

	private int indexOf(final Object name) {
		for (int i = 0; i < mSigmas.size(); ++i)
			if (mSigmas.get(i).mName.equals(name))
				return i;
		return -1;
	}

	/**
	 * Flattens the various different uncertainties associated with this
	 * UncertainValue down to a single value associated with a single source.
	 *
	 * @param name
	 * @returns UncertainValue
	 */
	public Number flatten(final Object name) {
		return unwrap(new UncertainValue(mValue, name, uncertainty()));
	}

	/**
	 * Flattens the various different uncertainties associated with this
	 * UncertainValue down to a single value associated with a single source.
	 *
	 * @param name
	 * @returns UncertainValue
	 */
	static public Number flatten(final Number n, final Object name) {
		if (n instanceof UncertainValue)
			return ((UncertainValue) n).flatten(name);
		else
			return n;
	}

	@Override
	public String toString() {
		if (mSigmas.size() > 0)
			return Double.toString(mValue) + " \u00B1 " + Double.toString(uncertainty());
		else
			return Double.toString(mValue);
	}

	/**
	 * Formats the {@link UncertainValue} as a val +- uncertainty using the
	 * specified NumberFormat for both.
	 *
	 * @param nf
	 * @return String
	 */
	public String format(final NumberFormat nf) {
		if (mSigmas.size() > 0)
			return nf.format(mValue) + "\u00B1" + nf.format(uncertainty());
		else
			return nf.format(mValue);
	}

	public static String format(final NumberFormat nf, final Number n) {
		if (n instanceof UncertainValue)
			return ((UncertainValue) n).format(nf);
		else
			return nf.format(n);
	}

	/**
	 * Formats the {@link UncertainValue} as a val +- dv1(src1) +- dv2(src2)...
	 * using the specified NumberFormat.
	 *
	 * @param nf
	 * @return String
	 */
	public String formatLong(final NumberFormat nf) {
		final StringBuffer sb = new StringBuffer();
		sb.append(nf.format(mValue));
		for (final Sigma sigma : mSigmas) {
			sb.append("\u00B1");
			sb.append(nf.format(sigma.mValue));
			sb.append("(");
			sb.append(sigma.mName);
			sb.append(")");
		}
		return sb.toString();
	}

	/**
	 * Format the specified uncertainty source in the format "U(src)=df".
	 *
	 * @param src
	 * @param nf
	 * @return String
	 */
	public String format(final Object src, final NumberFormat nf) {
		return "U(" + src + ")=" + nf.format(getComponent(src));
	}

	/**
	 * Returns the specified source of uncertainty.
	 *
	 * @param src
	 * @return The specified uncertainty or 0.0 if src not defined.
	 */
	public double getComponent(final Object src) {
		final int idx = indexOf(src);
		final Double dv = (idx >= 0 ? mSigmas.get(idx).mValue : Double.valueOf(0.0));
		return dv.doubleValue();
	}

	public String formatComponent(final String comp, final NumberFormat nf) {
		return nf.format(getComponent(comp)) + "(" + comp + ")";
	}

	/**
	 * Is this component defined?
	 *
	 * @param src
	 * @return true or false
	 */
	public boolean hasComponent(final Object src) {
		return indexOf(src) != -1;
	}

	public Map<Object, Double> getComponents() {
		final Map<Object, Double> res = new HashMap<>();
		for (final Sigma s : mSigmas)
			res.put(s.mName, s.mValue);
		return Collections.unmodifiableMap(res);
	}

	public Set<Object> getComponentNames() {
		final Set<Object> res = new HashSet<>();
		for (final Sigma s : mSigmas)
			res.add(s.mName);
		return Collections.unmodifiableSet(res);
	}

	/**
	 * Rename an uncertainty component. Fails with an EPQException if a component
	 * with the new name already exists.
	 *
	 * @param oldName
	 * @param newName
	 */
	public void renameComponent(final String oldName, final String newName) {
		final int idx = indexOf(oldName);
		if (idx != -1)
			mSigmas.set(idx, new Sigma(newName, mSigmas.get(idx).mValue));
	}

	static private UncertainValue asUncertainValue(final Number n) {
		if (n instanceof UncertainValue)
			return (UncertainValue) n;
		else
			return new UncertainValue(n.doubleValue());
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
		return Double.isNaN(iVarSum) ? UncertainValue.NaN
				: new UncertainValue(sum / varSum, "WM", Math.sqrt(1.0 / varSum));
	}

	public static double uncertainty(final Number n) {
		return isSpecialNumber(n) ? asUncertainValue(n).uncertainty() : 0.0;
	}

	public static double mean(final Number n) {
		return n.doubleValue();
	}

	/**
	 * True if the uncertainty associated with this item is non-zero.
	 *
	 * @return boolean
	 */
	public boolean isUncertain() {
		return mSigmas.size() > 0;
	}

	/**
	 * The uncertainty in the UncertainValue.
	 *
	 * @return A double
	 */
	public double uncertainty() {
		return Math.sqrt(variance());
	}

	/**
	 * The standard deviation in the UncertainValue (equal to uncertainty()).
	 *
	 * @return A double
	 */
	public double getStandardDeviation() {
		return Math.sqrt(variance());
	}

	/**
	 * The variance associated with this UncertainValue
	 *
	 * @return double
	 */
	public double variance() {
		double sigma2 = 0.0;
		for (final Sigma s : mSigmas)
			sigma2 += s.mValue * s.mValue;
		return sigma2;
	}

	/**
	 * The fractional uncertainty = Math.abs(uncertainty()/doubleValue())
	 *
	 * @return A double
	 */
	public double fractionalUncertainty() {
		return Double.isNaN(1.0 / mValue) ? Double.NaN : Math.abs(uncertainty() / mValue);
	}

	public boolean isNaN() {
		return Double.isNaN(mValue);
	}

	public static boolean isNaN(final Number n) {
		return Double.isNaN(n.doubleValue());
	}

	/**
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hash(mValue, mSigmas);
	}

	/**
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValue other = (UncertainValue) obj;
		return mSigmas.equals(other.mSigmas) && (mValue == other.mValue);
	}

	/**
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(final UncertainValue o) {
		final int res = Double.compare(doubleValue(), o.doubleValue());
		return res != 0 ? res : Double.compare(uncertainty(this), uncertainty(o));
	}

	public UncertainValue positiveDefinite() {
		return doubleValue() >= 0.0 ? this : new UncertainValue(0.0, this.mSigmas);
	}

	/**
	 * A random variate distributed such that it has a mean equal to the nominal
	 * value of this UncertainValue and a variance equal to the variance of this
	 * UncertainValue.
	 *
	 * @return double
	 */
	public double randomVariate() {
		return mValue + (sRandom.nextGaussian() * uncertainty());
	}

	public String toHTML(final Mode mode, final BasicNumberFormat bnf) {
		switch (mode) {
		case TERSE: {
			return bnf.formatHTML(mValue) + "&nbsp;&pm;&nbsp;" + bnf.formatHTML(uncertainty());
		}
		case NORMAL: {
			final StringBuffer sb = new StringBuffer();
			sb.append(bnf.formatHTML(mValue));
			if (!mSigmas.isEmpty()) {
				boolean first = true;
				sb.append("&pm;(");
				List<Sigma> sigmas = new ArrayList<>(mSigmas);
				sigmas.sort(new Comparator<Sigma>() {
					@Override
					public int compare(Sigma o1, Sigma o2) {
						return -Double.compare(o1.mValue, o2.mValue);
					}
				});
				for (final Sigma sigma : sigmas) {
					if (!first)
						sb.append(",");
					sb.append(HTML.toHTML(sigma.mName, Mode.TERSE));
					sb.append(":");
					sb.append(bnf.formatHTML(sigma.mValue));
					first = false;
				}
				sb.append(")");
			}
			return sb.toString();
		}
		case VERBOSE:
		default: {
			final StringBuffer sb = new StringBuffer();
			BasicNumberFormat pct = new BasicNumberFormat("0.0");
			sb.append("<table><tr><td>");
			sb.append(bnf.formatHTML(mValue));
			sb.append("</td>");
			if (!mSigmas.isEmpty()) {
				sb.append("<td>&pm;</td>");
				sb.append("<td><table class=\"matrix\">");
				List<Sigma> sigmas = new ArrayList<>(mSigmas);
				sigmas.sort(new Comparator<Sigma>() {
					@Override
					public int compare(Sigma o1, Sigma o2) {
						return -Double.compare(o1.mValue, o2.mValue);
					}
				});
				for (final Sigma sigma : sigmas) {
					sb.append("<tr>");
					sb.append("<td>U(" + HTML.toHTML(sigma.mName, Mode.TERSE) + ")</td>");
					sb.append("<td>" + bnf.formatHTML(sigma.mValue) + "</td>");
					sb.append("<td>&nbsp;&nbsp;" + pct.formatHTML(100.0 * sigma.mValue / mValue) + " %</td>");
					sb.append("</tr>");
				}
				sb.append("</table></td>");
			}
			sb.append("</tr></table>");
			return sb.toString();
		}
		}

	}

	/**
	 * @param mode
	 * @return
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mode, new BasicNumberFormat());
	}
}
