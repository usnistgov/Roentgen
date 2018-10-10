package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

/**
 * Solve C<sub>y</sub> U<sub>y</sub> C<sub>y</sub><sup>T</sup> = C<sub>x,z</sub>
 * U<sub>x,z</sub> C<sub>x,z</sub><sup>T</sup> for
 * C<sub>y</sub><sup>-1</sup>C<sub>x,z</sub> and <b>Y</b>.
 * <ul>
 * <li><b>h</b>(<b>X</b>,<b>Y</b>,<b>Z</b>) = <b>0</b>
 * <li>C<sub>y</sub> = &delta;h/&delta;C<sub>y</sub>
 * <li>C<sub>x,z</sub> = &delta;h/&delta;C<sub>x</sub> and
 * &delta;h/&delta;C<sub>z</sub>,
 * </ul>
 *
 * <p>
 * The values associated with the output labels are passed in as constant values
 * to conform to ensure that the input and output labels are disjoint.
 * </p>
 *
 *
 * @author Nicholas
 *
 */
public class ImplicitMeasurementModel //
		extends LabeledMultivariateJacobianFunction {

	public static class HLabel extends BaseLabel<Object, Object, Object> {

		public HLabel(final Object label) {
			super("h", label);
		}
	}

	/**
	 * Builds a special set of labels to represent the function
	 * <b>h</b>(<b>X</b>,<b>Y</b>,<b>Z</b>) = <b>0</b>.
	 * 
	 * @param outlabels
	 * @return List&lt;? extends Object&gt;
	 */
	public static List<? extends Object> buildHLabels(//
			final List<? extends Object> outlabels //
	) {
		final List<Object> res = new ArrayList<>();
		for (final Object label : outlabels)
			res.add(new HLabel(label));
		return res;
	}

	private final LabeledMultivariateJacobianFunction mHFunction;

	private static List<? extends Object> buildInputs(//
			final LabeledMultivariateJacobianFunction h, //
			final List<? extends Object> outputLabels //
	) {
		final List<? extends Object> res = new ArrayList<>(h.getInputLabels());
		res.removeAll(outputLabels);
		return res;
	}

	public ImplicitMeasurementModel(//
			final LabeledMultivariateJacobianFunction h, //
			final List<? extends Object> outputLabels //
	) {
		super(buildInputs(h, outputLabels), outputLabels);
		mHFunction = h;
	}

	/**
	 * Implements value(...) for an implicit measurement model by evaluating the
	 * h-function splitting it apart into <b>J<sub>y</sub></b> and
	 * <b>J<sub>x</sub></b> parts from which is computed <b>J</b> =
	 * <b>J<sub>y</sub></b><sup>-1</sup><b>J<sub>x</sub></b>
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 *      value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector hPoint = new ArrayRealVector(mHFunction.getInputDimension());
		for (int r = 0; r < hPoint.getDimension(); ++r)
			hPoint.setEntry(r, getValue(mHFunction.getInputLabels().get(r), point));
		final Pair<RealVector, RealMatrix> hpr = mHFunction.evaluate(hPoint);
		// hpr.getFirst() values should all be close to zero...
		final RealMatrix hrm = hpr.getSecond();
		final List<? extends Object> hIn = mHFunction.getInputLabels();
		final int outDim = getOutputDimension();
		final int inDim = getInputDimension();
		final RealMatrix mcy = MatrixUtils.createRealMatrix(outDim, outDim);
		for (int c = 0; c < outDim; ++c) {
			final int cIdx = hIn.indexOf(getOutputLabels().get(c));
			for (int r = 0; r < outDim; ++r)
				mcy.setEntry(r, c, hrm.getEntry(r, cIdx));
		}
		final RealMatrix mcyInv = MatrixUtils.inverse(mcy);
		final RealMatrix mcx = MatrixUtils.createRealMatrix(outDim, inDim);
		for (int c = 0; c < inDim; ++c) {
			final int cIdx = hIn.indexOf(getInputLabels().get(c));
			for (int r = 0; r < outDim; ++r)
				mcx.setEntry(r, c, hrm.getEntry(r, cIdx));
		}
		final RealVector rv = new ArrayRealVector(outDim);
		for (int r = 0; r < rv.getDimension(); ++r)
			rv.setEntry(r, getConstant(getOutputLabels().get(r)));
		return Pair.create(rv, mcyInv.multiply(mcx));
	}
}
