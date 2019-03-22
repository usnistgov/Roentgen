/**
 *
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.duckandcover.html.HTML;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * @author nicho
 *
 */
public class UncertainValueEx<H> extends UncertainValue {

	private static final long serialVersionUID = -5637027270542233886L;

	final List<Sigma> mSigmas;

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

	/**
	 * @param value
	 * @param sigma
	 * @param variances
	 */
	public UncertainValueEx(final double value, final double sigma, final Map<? extends H, Double> variances) {
		super(value, sigma);
		mSigmas = new ArrayList<>();
		for (final Map.Entry<? extends H, Double> me : variances.entrySet())
			mSigmas.add(new Sigma(me.getKey(), me.getValue()));
	}

	/**
	 * @param value
	 * @param sigma
	 */
	public UncertainValueEx(final double value, final double sigma) {
		this(value, sigma, Collections.emptyMap());
	}

	/**
	 * @param value
	 * @param sigma
	 * @param variances
	 */
	public UncertainValueEx(final double value) {
		this(value, 0.0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(final UncertainValue arg0) {
		int res = Double.compare(mValue, arg0.mValue);
		if (res == 0)
			res = Double.compare(mSigma, arg0.mSigma);
		if (arg0 instanceof UncertainValueEx<?>)
			if (res == 0) {
				final List<Sigma> ss = sortedSigmas();
				@SuppressWarnings({ "unchecked", "rawtypes" })
				final List<Sigma> ss0 = ((UncertainValueEx) arg0).sortedSigmas();
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
	
	public List<H> getComponentNames(){
		final List<H> res = new ArrayList<>();
		for (final Sigma s : mSigmas)
			res.add(s.mLabel);
		return res;
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

	public void assignComponent(final H label, final double val) {
		for (final Sigma s : mSigmas)
			if (s.mLabel.equals(label)) {
				mSigmas.remove(s);
				break;
			}
		mSigmas.add(new Sigma(label, val));
	}

	public double getComponent(final Object label) {
		for (final Sigma s : mSigmas)
			if (s.mLabel.equals(label))
				return s.mOneSigma;
		return 0.0;
	}

}
