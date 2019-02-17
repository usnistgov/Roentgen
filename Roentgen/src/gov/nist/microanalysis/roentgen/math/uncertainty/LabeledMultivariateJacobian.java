package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
 * A Jacobian function evaluated at a specific point in input label space. This
 * is particularly useful to optimize a complex multistep calculations in which
 * result from an early part of the calculation is reused many times in the
 * later part of the calculation.
 *
 * @author Nicholas
 */
public class LabeledMultivariateJacobian //
		extends LabeledMultivariateJacobianFunction //
		implements IToHTML {

	/**
	 * The evaluation point
	 */
	private final RealVector mPoint;
	/**
	 * The values at the evaluation point
	 */
	private final RealVector mValues;
	/**
	 * The Jacobian at the evaluation point
	 */
	private final RealMatrix mJacobian;

	private final Map<Object, Double> mConstants;

	/**
	 * @param inpLabels List of input labels
	 * @param outLabels List of output labels
	 * @param pt        The point at which the Jacobian is evaluated
	 * @param result    The quantities and associated covariances at the evaluation
	 *                  point
	 */
	private LabeledMultivariateJacobian( //
			final List<? extends Object> inpLabels, //
			final List<? extends Object> outLabels, //
			final RealVector pt, //
			final Pair<RealVector, //
					RealMatrix> result //
	) {
		super(inpLabels, outLabels);
		mPoint = pt;
		mValues = result.getFirst();
		mJacobian = result.getSecond();
		mConstants = new HashMap<>();
	}

	/**
	 * Construct a {@link LabeledMultivariateJacobian} from a
	 * {@link LabeledMultivariateJacobianFunction} evaluated at a specific
	 * {@link RealVector} point.
	 *
	 * @param nmjf A {@link LabeledMultivariateJacobianFunction}
	 * @param pt   The evaluation point
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public static LabeledMultivariateJacobian compute( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final RealVector pt //
	) {
		return new LabeledMultivariateJacobian(nmjf, pt);
	}

	/**
	 * Construct a {@link LabeledMultivariateJacobian} from a
	 * {@link LabeledMultivariateJacobianFunction} evaluated at a specific
	 * {@link UncertainValues} object
	 *
	 * @param nmjf A {@link LabeledMultivariateJacobianFunction}
	 * @param uvs  The evaluation point represented as an {@link UncertainValues}
	 *             object
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public static LabeledMultivariateJacobian compute( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues uvs //
	) {
		return new LabeledMultivariateJacobian(nmjf, uvs);
	}

	/**
	 * Computes the values and an estimate of the Jacobian using a finite difference
	 * algorithm.
	 *
	 * @param nmjf {@link LabeledMultivariateJacobianFunction}
	 * @param uv   {@link UncertainValues} representing the point at which to
	 *             evaluate nmjf
	 * @param sc   A scale factor which together with the uncertainties in uv
	 *             determines the size of the finite difference in each dimension.
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public static LabeledMultivariateJacobian computeDelta( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues uv, //
			final double sc //
	) {
		assert uv.getDimension() == nmjf.getInputDimension();
		final RealMatrix rm = NullableRealMatrix.build(nmjf.getOutputDimension(), nmjf.getInputDimension());
		final RealVector inp = new ArrayRealVector(uv.getValues());
		final RealVector vals = nmjf.compute(inp);
		for (int c = 0; c < nmjf.getInputDimension(); ++c) {
			final RealVector pt0 = new ArrayRealVector(inp), pt1 = new ArrayRealVector(inp);
			final double deltaX = sc * Math.max(uv.getUncertainty(c), Math.abs(0.01 * pt0.getEntry(c)));
			assert deltaX > 0.0 : //
			nmjf.getInputLabels().get(c) + ", sc=" + sc + ", uv.unc=" + uv.getUncertainty(c) + ", 0.01*pt0="
					+ Math.abs(0.01 * pt0.getEntry(c));
			pt0.setEntry(c, pt0.getEntry(c) + 0.5 * deltaX);
			pt1.setEntry(c, pt1.getEntry(c) - 0.5 * deltaX);
			final RealVector output0 = nmjf.compute(pt0), output1 = nmjf.compute(pt1);
			for (int r = 0; r < nmjf.getOutputDimension(); ++r)
				rm.setEntry(r, c, (output0.getEntry(r) - output1.getEntry(r)) / deltaX);
		}
		return new LabeledMultivariateJacobian(nmjf.getInputLabels(), nmjf.getOutputLabels(), inp,
				Pair.create(vals, rm));
	}

	/**
	 * Construct a {@link LabeledMultivariateJacobian} from a
	 * {@link LabeledMultivariateJacobianFunction} evaluated at a specific
	 * {@link RealVector} point.
	 *
	 * @param nmjf A {@link LabeledMultivariateJacobianFunction}
	 * @param pt   The evaluation point
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public LabeledMultivariateJacobian( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final RealVector pt //
	) {
		this(nmjf.getInputLabels(), nmjf.getOutputLabels(), pt, nmjf.evaluate(pt));
	}

	private static final RealVector extractPoint(final LabeledMultivariateJacobianFunction func,
			final UncertainValues args) {
		return args.extractValues(func.getInputLabels());
	}

	/**
	 * Construct a {@link LabeledMultivariateJacobian} from a
	 * {@link LabeledMultivariateJacobianFunction} evaluated at a specific
	 * {@link RealVector} point.
	 *
	 * @param nmjf A {@link LabeledMultivariateJacobianFunction}
	 * @param pt   The evaluation point
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public LabeledMultivariateJacobian( //
			final LabeledMultivariateJacobianFunction nmjf, //
			final UncertainValues args //
	) {
		this(nmjf.getInputLabels(), //
				nmjf.getOutputLabels(), //
				extractPoint(nmjf, args), //
				nmjf.evaluate(extractPoint(nmjf, args))//
		);
		mConstants.putAll(nmjf.getConstants());
	}

	@Override
	public Map<Object, Double> getConstants() {
		return Collections.unmodifiableMap(mConstants);
	}

	/**
	 * Returns a {@link RealVector} containing the point at which the
	 * {@link LabeledMultivariateJacobianFunction} was evaluated.
	 *
	 * @return
	 */
	public RealVector getInputValues() {
		return mPoint;
	}

	/**
	 * Implements the value(...) method in the MultivariateJacobianFunction
	 * interface.
	 *
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 *      value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		assert point.equals(
				mPoint) : "The argument to NamedMultivariateJacobian.value(...) does not match the one used to create it.";
		return point.equals(mPoint) ? Pair.create(mValues, mJacobian) : null;
	}

	/**
	 * Extract a subset of the output values and the associated Jacobian elements.
	 *
	 * @param outputLabels
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public LabeledMultivariateJacobian extract(final List<? extends Object> outputLabels) {
		final RealVector rv = new ArrayRealVector(outputLabels.size());
		final RealMatrix rm = MatrixUtils.createRealMatrix(outputLabels.size(), getInputDimension());
		for (int r = 0; r < outputLabels.size(); ++r) {
			final int row = outputIndex(outputLabels.get(r));
			assert row != -1 : "Input index associated with " + outputLabels.get(r) + " is not available.";
			rv.setEntry(r, mValues.getEntry(row));
			for (int c = 0; c < getInputDimension(); ++c)
				rm.setEntry(r, c, mJacobian.getEntry(row, c));
		}
		return new LabeledMultivariateJacobian(getInputLabels(), outputLabels, mPoint, Pair.create(rv, rm));
	}

	/**
	 * Extract a subset of the input and output values and the associated Jacobian
	 * elements.
	 *
	 * @param input
	 * @param output
	 * @return {@link LabeledMultivariateJacobian}
	 */
	public LabeledMultivariateJacobian extract(final List<? extends Object> input,
			final List<? extends Object> output) {
		final RealVector rv = new ArrayRealVector(output.size());
		final RealMatrix rm = MatrixUtils.createRealMatrix(output.size(), input.size());
		final int[] colIdx = new int[input.size()];
		for (int c = 0; c < input.size(); ++c) {
			colIdx[c] = inputIndex(input.get(c));
			assert colIdx[c] != -1 : "Input index associated with " + input.get(c) + " is not available.";
		}

		for (int r = 0; r < output.size(); ++r) {
			final int row = outputIndex(output.get(r));
			assert row != -1 : "Output index associated with " + output.get(r) + " is not available.";
			rv.setEntry(r, mValues.getEntry(row));
			for (int c = 0; c < colIdx.length; ++c)
				rm.setEntry(r, c, mJacobian.getEntry(row, colIdx[c]));
		}
		return new LabeledMultivariateJacobian(input, output, mPoint, Pair.create(rv, rm));
	}

	public String toCSV() {
		final StringBuffer sb = new StringBuffer(4096);
		final List<? extends Object> inputLabels = getInputLabels();
		final List<? extends Object> outputLabels = getOutputLabels();
		sb.append("\"Name\",\"Value\"");
		for (int c = 0; c < mJacobian.getColumnDimension(); ++c) {
			sb.append(",\"");
			sb.append(inputLabels.get(c));
			sb.append("\"");
		}
		sb.append("\n");
		for (int r = 0; r < mValues.getDimension(); ++r) {
			sb.append("\"");
			sb.append(outputLabels.get(r));
			sb.append("\",");
			sb.append(mValues.getEntry(r));
			for (int c = 0; c < mJacobian.getColumnDimension(); ++c) {
				sb.append(",");
				sb.append(mJacobian.getEntry(r, c));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	public double getEntry(final int outIdx, final int inIdx) {
		return mJacobian.getEntry(outIdx, inIdx);
	}

	public double getEntry(final Object label1, final Object label2) {
		final int outputIndex = outputIndex(label1);
		final int inputIndex = inputIndex(label2);
		assert outputIndex != -1 : label1;
		assert inputIndex != -1 : label2;
		return getEntry(outputIndex, inputIndex);
	}

	/**
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mode, new BasicNumberFormat());
	}

	public String toHTML(final Mode mode, final BasicNumberFormat nf) {
		switch (mode) {
		case TERSE: {
			final Table table = new Table();
			final List<Item> header = new ArrayList<>();
			final List<Item> vals = new ArrayList<>();
			for (int i = 0; i < getOutputDimension(); ++i) {
				header.add(Table.th(HTML.toHTML(getOutputLabels().get(i), Mode.TERSE)));
				vals.add(Table.tdc(nf.formatHTML(mValues.getEntry(i))));
			}
			table.addRow(header);
			table.addRow(vals);
			return table.toHTML(Mode.NORMAL);
		}
		case NORMAL:
		default: {
			final Table jact = new Table();
			final List<Item> items = new ArrayList<>();
			items.add(Table.td("Output"));
			items.add(Table.td("Value"));
			for (int r = 0; r < getInputDimension(); ++r)
				items.add(Table.td(HTML.toHTML(getInputLabels().get(r), Mode.TERSE)));
			jact.addRow(items);
			for (int r = 0; r < getOutputDimension(); ++r) {
				final Object label = getOutputLabels().get(r);
				items.clear();
				items.add(Table.td(HTML.toHTML(label, Mode.TERSE)));
				items.add(Table.td(nf.format(mValues.getEntry(r))));
				for (int c = 0; c < getInputDimension(); ++c)
					items.add(Table.td(nf.format(mJacobian.getEntry(r, c))));
				jact.addRow(items);
			}
			return jact.toHTML(Mode.NORMAL);
		}
		case VERBOSE: {
			final Table jact = new Table();
			final List<Item> items = new ArrayList<>();
			items.add(Table.td("Output"));
			items.add(Table.td("Value"));
			for (int r = 0; r < getInputDimension(); ++r)
				items.add(Table.td(HTML.toHTML(getInputLabels().get(r), Mode.TERSE)));
			jact.addRow(items);
			for (int r = 0; r < getOutputDimension(); ++r) {
				final Object label = getOutputLabels().get(r);
				items.clear();
				items.add(Table.td(HTML.toHTML(label, Mode.TERSE)));
				items.add(Table.td(nf.formatHTML(mValues.getEntry(r))));
				for (int c = 0; c < getInputDimension(); ++c)
					items.add(Table.td(nf.formatHTML(mJacobian.getEntry(r, c))));
				jact.addRow(items);
			}
			return jact.toHTML(Mode.NORMAL);
		}
		}
	}
}
