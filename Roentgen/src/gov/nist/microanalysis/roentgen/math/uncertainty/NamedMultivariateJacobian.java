package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A Jacobian function evaluated at a specific point in input tag space. This is
 * particularly useful for complex multistep calculations in which the early
 * part of the calculation is reused many times in the later part of the
 * calculation.
 *
 * @author Nicholas
 */
public class NamedMultivariateJacobian extends NamedMultivariateJacobianFunctionEx implements IToHTML {

	private final RealVector mPoint;
	private final Pair<RealVector, RealMatrix> mResult;

	public NamedMultivariateJacobian(final List<? extends Object> inp, final List<? extends Object> out,
			final RealVector pt, final Pair<RealVector, RealMatrix> result) {
		super(inp, out);
		mPoint = pt;
		mResult = result;
	}

	public static NamedMultivariateJacobian compute(final NamedMultivariateJacobianFunction nmjf, final RealVector pt) {
		return new NamedMultivariateJacobian(nmjf, pt);
	}

	public static NamedMultivariateJacobian computeDelta(final NamedMultivariateJacobianFunction nmjf,
			final UncertainValues uv, final double sc) {
		assert uv.getDimension() == nmjf.getInputDimension();
		final RealMatrix rm = new NullableRealMatrix(nmjf.getOutputDimension(), nmjf.getInputDimension());
		final RealVector inp = new ArrayRealVector(uv.getValues());
		final RealVector vals = nmjf.compute(inp);
		for (int c = 0; c < nmjf.getInputDimension(); ++c) {
			final RealVector pt0 = new ArrayRealVector(inp), pt1 = new ArrayRealVector(inp);
			final double deltaX = sc * Math.max(uv.getUncertainty(c), Math.abs(0.01 * pt0.getEntry(c)));
			assert (deltaX > 0.0);
			pt0.setEntry(c, pt0.getEntry(c) + 0.5 * deltaX);
			pt1.setEntry(c, pt1.getEntry(c) - 0.5 * deltaX);
			final RealVector output0 = nmjf.compute(pt0), output1 = nmjf.compute(pt1);
			for (int r = 0; r < nmjf.getOutputDimension(); ++r) {
				rm.setEntry(r, c, (output0.getEntry(r) - output1.getEntry(r)) / deltaX);
			}
		}
		return new NamedMultivariateJacobian(nmjf.getInputTags(), nmjf.getOutputTags(), inp, Pair.create(vals, rm));
	}

	public NamedMultivariateJacobian(final NamedMultivariateJacobianFunction nmjf, final RealVector pt) {
		this(nmjf.getInputTags(), nmjf.getOutputTags(), pt, nmjf.evaluate(pt));
	}

	public RealVector getInputValues() {
		return mPoint;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		if (point.equals(mPoint))
			return mResult;
		else
			return null;
	}

	/**
	 * Extract a subset of the output values and the associated Jacobian elements.
	 *
	 * @param output
	 * @return {@link NamedMultivariateJacobian}
	 */
	public NamedMultivariateJacobian extract(final List<? extends Object> output) {
		final RealVector rv = new ArrayRealVector(output.size());
		final RealMatrix rm = MatrixUtils.createRealMatrix(output.size(), getInputDimension());
		final RealVector vals = mResult.getFirst();
		final RealMatrix jac = mResult.getSecond();
		for (int r = 0; r < output.size(); ++r) {
			final int idx = inputIndex(output.get(r));
			rv.setEntry(r, vals.getEntry(idx));
			for (int c = 0; c < getInputDimension(); ++c)
				rm.setEntry(r, c, jac.getEntry(idx, c));
		}
		return new NamedMultivariateJacobian(getInputTags(), output, mPoint, Pair.create(rv, rm));
	}

	public String toCSV() {
		final RealVector vals = mResult.getFirst();
		final RealMatrix jac = mResult.getSecond();
		final StringBuffer sb = new StringBuffer(4096);
		final List<? extends Object> inputTags = getInputTags();
		final List<? extends Object> outputTags = getOutputTags();
		sb.append("\"Name\",\"Value\"");
		for (int c = 0; c < jac.getColumnDimension(); ++c) {
			sb.append(",\"");
			sb.append(inputTags.get(c));
			sb.append("\"");
		}
		sb.append("\n");
		for (int r = 0; r < vals.getDimension(); ++r) {
			sb.append("\"");
			sb.append(outputTags.get(r));
			sb.append("\",");
			sb.append(vals.getEntry(r));
			for (int c = 0; c < jac.getColumnDimension(); ++c) {
				sb.append(",");
				sb.append(jac.getEntry(r, c));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	public double getEntry(final int outIdx, final int inIdx) {
		return mResult.getSecond().getEntry(outIdx, inIdx);
	}

	public double getEntry(final Object tag1, final Object tag2) {
		final int outputIndex = outputIndex(tag1);
		final int inputIndex = inputIndex(tag2);
		assert outputIndex != -1 : tag1;
		assert inputIndex != -1 : tag2;
		return mResult.getSecond().getEntry(outputIndex, inputIndex);
	}

	/**
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mode, new BasicNumberFormat());
	}

	public String toHTML(final Mode mode, final BasicNumberFormat nf) {
		final RealVector val = mResult.getFirst();
		final RealMatrix jac = mResult.getSecond();
		switch (mode) {
		case TERSE: {
			final Table table = new Table();
			final List<Item> header = new ArrayList<>();
			final List<Item> vals = new ArrayList<>();
			for (int i = 0; i < getOutputDimension(); ++i) {
				header.add(Table.th(HTML.toHTML(getOutputTags().get(i), Mode.TERSE)));
				vals.add(Table.tdc(nf.formatHTML(val.getEntry(i))));
			}
			table.addRow(header);
			table.addRow(vals);
			return table.toHTML(Mode.NORMAL);
		}
		case NORMAL:
		case VERBOSE:
		default: {
			final Table valt = new Table();
			final Table jact = new Table();
			final List<Item> items = new ArrayList<>();
			for (int r = 0; r < getInputDimension(); ++r)
				items.add(Table.td(HTML.toHTML(getInputTags().get(r), Mode.TERSE)));
			jact.addRow(items);
			valt.addRow(Table.td("Output"), Table.td("Values"));
			for (int r = 0; r < getOutputDimension(); ++r) {
				final Object tag = getOutputTags().get(r);
				valt.addRow(Table.td(HTML.toHTML(tag, Mode.TERSE)), Table.td(nf.formatHTML(val.getEntry(r))));
				items.clear();
				// items.add(Table.td(HTML.toHTML(getInputTags().get(r),
				// Mode.TERSE)));
				for (int c = 0; c < getInputDimension(); ++c)
					items.add(Table.td(nf.formatHTML(jac.getEntry(r, c))));
				jact.addRow(items);
			}

			final Table res = new Table();
			// res.addRow(Table.td("Values"), Table.td(), Table.td("Jacobian"));
			res.addRow(Table.td(valt.toHTML(Mode.NORMAL)), Table.tdc("<bold>J</bold> ="),
					Table.td(jact.toHTML(Mode.NORMAL)));
			return res.toHTML(Mode.NORMAL);
		}
		}
	}
}
