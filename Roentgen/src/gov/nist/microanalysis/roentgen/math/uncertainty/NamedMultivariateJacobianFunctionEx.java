package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.duckandcover.html.HTML;
import com.duckandcover.html.Table;

/**
 * Adds support for passing constant values to a
 * {@link NamedMultivariateJacobianFunction} instance. This can be useful for
 * situations in which some input values are not being treated as sources of
 * uncertainty. The functions <code>getValue(...)</code> and
 * <code>writeCovariance(...)</code> function can be used as a mechanism to
 * switch between using a variable as a source of uncertainty or not.
 *
 *
 * @author nicholas
 *
 */
abstract public class NamedMultivariateJacobianFunctionEx extends NamedMultivariateJacobianFunction {

	private final Map<Object, Double> mConstants = new HashMap<>();

	public NamedMultivariateJacobianFunctionEx(final List<? extends Object> inputTags,
			final List<? extends Object> outputTags) {
		super(inputTags, outputTags);
	}

	public NamedMultivariateJacobianFunctionEx( //
			final List<? extends Object> inputTags, //
			final List<? extends Object> outputTags, //
			Map<Object, Double> constants//
	) {
		super(inputTags, outputTags);
		for (Map.Entry<Object, Double> me : constants.entrySet()) {
			assert inputIndex(me.getKey()) == -1;
			mConstants.put(me.getKey(), me.getValue());
		}
	}

	/**
	 * Initializes the constant tags with the specified values. The tag for value
	 * <code>vals.getEntry(i)</code> is <code>list.get(i)</code>. Checks first to
	 * see is the tag is being used as an input variable and won't define it as a
	 * constant if it is an input item.
	 * 
	 * @param list
	 * @param vals
	 */
	public void initializeConstants(final List<? extends Object> list, final RealVector vals) {
		for (int i = 0; i < list.size(); ++i)
			if (inputIndex(list.get(i)) == -1)
				mConstants.put(list.get(i), vals.getEntry(i));
	}

	/**
	 * Initializes the constant tags with the specified values. The tag for value
	 * <code>vals.getEntry(i)</code> is <code>list.get(i)</code>. 
	 * 
	 * @param list
	 * @param vals
	 */
	public void initializeConstants(final Map<Object, Double> mod) {
		for (Map.Entry<Object, Double> me : mod.entrySet())
			mConstants.put(me.getKey(), me.getValue());
	}

	/**
	 * Returns the constant value associated with <code>tag</code>.
	 * 
	 * @param tag
	 * @return double
	 */
	public double getConstant(final Object tag) {
		assert inputIndex(tag) == -1 : "Tag " + tag + " is not a constant.";
		return mConstants.get(tag).doubleValue();
	}

	/**
	 * Check whether <code>tag</code> is defined as a constant value.
	 * 
	 * @param tag A tag
	 * @return true if <code>tag</code> is initialized as a constant, false
	 *         otherwise.
	 */
	public boolean isConstant(final Object tag) {
		return mConstants.containsKey(tag);
	}

	/**
	 * First checks to see if <code>tag</code> is an input variable in which case it
	 * gets the associated value from <code>point</code>. Otherwise, it returns the
	 * constant value associate with <code>tag</code>.
	 *
	 * @param tag
	 * @param point
	 * @return double
	 */
	public double getValue(final Object tag, final RealVector point) {
		final int p = inputIndex(tag);
		if (p != -1) {
			assert !isConstant(tag);
			return point.getEntry(p);
		} else {
			assert isConstant(tag);
			return getConstant(tag);
		}
	}

	/**
	 * Only writes <code>value</code> to the covariance matrix <code>cov</code> if
	 * <code>tag</code> is an input variable (not a constant.)
	 *
	 *
	 * @param row
	 * @param tag
	 * @param value
	 * @param cov
	 */
	public void writeCovariance(final int row, final Object tag, final double value, final RealMatrix cov) {
		final int p = inputIndex(tag);
		if (p == -1) {
			assert !isConstant(tag);
			cov.setEntry(row, p, value);
		}
	}

	/**
	 * Returns an unmodifiable view of the map of tags and the associated constant
	 * values.
	 * 
	 * @return Map<Object, Double>
	 */
	public Map<Object, Double> getConstants() {
		return Collections.unmodifiableMap(mConstants);
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE: {
			return HTML.escape("V[" + getOutputTags().size() + " values]=F(" + getInputTags().size() + " arguments, "
					+ getConstants().size() + " constants)");
		}
		case NORMAL: {
			final StringBuffer sb = new StringBuffer();
			for (final Object tag : getOutputTags())
				if (sb.length() == 0)
					sb.append("<" + tag.toString());
				else
					sb.append("," + tag.toString());
			sb.append(">=F(");
			boolean first = true;
			for (final Object tag : getOutputTags()) {
				if (!first)
					sb.append(",");
				sb.append(tag.toString());
				first = false;
			}
			first = true;
			for (final Object tag : getConstants().keySet()) {
				sb.append(first ? ";" : ",");
				sb.append(tag.toString());
				first = false;
			}
			sb.append(")");
			return HTML.escape(sb.toString());
		}
		default:
		case VERBOSE: {
			final Table res = new Table();
			final Table outs = new Table();
			for (final Object tag : getOutputTags())
				outs.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			final Table args = new Table();
			for (final Object tag : getInputTags())
				args.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			final Table consts = new Table();
			for (final Object tag : getConstants().keySet())
				consts.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)));
			res.addRow(//
					Table.td(outs.toHTML(Mode.NORMAL)), //
					Table.td(" = F("), //
					Table.td(args.toHTML(Mode.NORMAL)), //
					Table.td(";"), //
					Table.td(consts.toHTML(Mode.NORMAL)), //
					Table.td(")")); //
			return res.toHTML(Mode.NORMAL);
		}
		}
	}

}