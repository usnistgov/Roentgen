package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class UncertainValue2<H> //
		extends Number //
		implements Comparable<UncertainValue2<H>>, IToHTML {

	private class Sigma {
		private final H mLabel;
		private final double mOneSigma;

		private Sigma(final H name, final double val) {
			mLabel = name;
			mOneSigma = Math.abs(val);
		}

		@Override
		public String toString() {
			return mLabel + " = " + mOneSigma;
		}
	}
	public static final UncertainValue2<?> ONE = new UncertainValue2<Object>(1.0);
	public static final UncertainValue2<?> ZERO = new UncertainValue2<Object>(0.0);
	public static final UncertainValue2<?> NaN = new UncertainValue2<Object>(Double.NaN);
	public static final UncertainValue2<?> POSITIVE_INFINITY = new UncertainValue2<Object>(Double.POSITIVE_INFINITY);
	public static final UncertainValue2<?> NEGATIVE_INFINITY = new UncertainValue2<Object>(Double.NEGATIVE_INFINITY);

	private static final long serialVersionUID = -7284125207920225793L;
	final Double mValue;
	final double mSigma;

	final List<Sigma> mSigmas;

	public UncertainValue2(final double value) {
		this(value, 0.0, Collections.emptyMap());
	}

	/**
	 *
	 */
	public UncertainValue2(final double value, final double sigma) {
		this(value, sigma, Collections.emptyMap());
	}

	/**
	 *
	 */
	public UncertainValue2(final double value, final double sigma, final Map<? extends H, Double> variances) {
		mValue = Double.valueOf(value);
		mSigma = sigma;
		mSigmas = new ArrayList<>();
		for (final Map.Entry<? extends H, Double> me : variances.entrySet())
			mSigmas.add(new Sigma(me.getKey(), me.getValue()));
	}

	@Override
	public int compareTo(final UncertainValue2<H> arg0) {
		int res = Double.compare(mValue, arg0.mValue);
		if (res == 0)
			res = Double.compare(mSigma, arg0.mSigma);
		if (res == 0) {
			final List<Sigma> ss = sortedSigmas();
			final List<Sigma> ss0 = arg0.sortedSigmas();
			for (int i = 0; i < Math.min(ss.size(), ss0.size()); ++i) {
				final Sigma s = ss.get(i), s0 = ss0.get(i);
				res = Double.compare(s.mOneSigma, s0.mOneSigma);
				if (res != 0)
					return res;
			}
			if (ss.size() > ss0.size())
				return 1;
			else if (ss0.size() > ss.size())
				return -1;
		}
		return res;
	}

	@Override
	public double doubleValue() {
		return mValue.doubleValue();
	}

	@Override
	public float floatValue() {
		return mValue.floatValue();
	}

	public Map<H, Double> getComponents() {
		final Map<H, Double> res = new HashMap<>();
		for (final Sigma s : mSigmas)
			res.put(s.mLabel, s.mOneSigma);
		return res;
	}

	public double getCovariance(final H label) {
		for (final Sigma s : mSigmas)
			if (s.mLabel.equals(label))
				return s.mOneSigma;
		return 0.0;
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
		case TERSE: {
			return bnf.formatHTML(mValue) + "&nbsp;&#177;&nbsp;" + bnf.formatHTML(uncertainty());
		}
		case NORMAL: {
			final StringBuffer sb = new StringBuffer();
			sb.append(bnf.formatHTML(mValue));
			if (!mSigmas.isEmpty()) {
				boolean first = true;
				sb.append("&#177;(");
				for (final Sigma sigma : sortedSigmas()) {
					if (!first)
						sb.append(",");
					sb.append(HTML.toHTML(sigma.mLabel, Mode.TERSE));
					sb.append(":");
					sb.append(bnf.formatHTML(sigma.mOneSigma));
					first = false;
				}
				sb.append(")");
			}
			return sb.toString();
		}
		case VERBOSE:
		default: {
			final StringBuffer sb = new StringBuffer();
			final BasicNumberFormat pct = new BasicNumberFormat("0.0");
			sb.append("<table><tr><td>");
			sb.append(bnf.formatHTML(mValue));
			sb.append("</td>");
			if (!mSigmas.isEmpty()) {
				sb.append("<td>&#177;</td>");
				sb.append("<td><table class=\"matrix\">");
				for (final Sigma sigma : sortedSigmas()) {
					sb.append("<tr>");
					sb.append("<td>U(" + HTML.toHTML(sigma.mLabel, Mode.TERSE) + ")</td>");
					sb.append("<td>" + bnf.formatHTML(sigma.mOneSigma) + "</td>");
					sb.append("<td>&nbsp;&nbsp;" + pct.formatHTML(100.0 * sigma.mOneSigma / mValue) + " %</td>");
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

	private List<Sigma> sortedSigmas() {
		final List<Sigma> sigmas = new ArrayList<>(mSigmas);
		sigmas.sort(new Comparator<Sigma>() {
			@Override
			public int compare(final Sigma o1, final Sigma o2) {
				return -Double.compare(o1.mOneSigma, o2.mOneSigma);
			}
		});
		return sigmas;
	}

}
