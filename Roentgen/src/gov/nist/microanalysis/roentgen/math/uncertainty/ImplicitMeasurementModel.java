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
 * The values associated with the output tags are passed in as constant values
 * to conform to ensure that the input and output tags are dijoint.
 * </p>
 *
 *
 * @author Nicholas
 *
 */
public class ImplicitMeasurementModel
		extends NamedMultivariateJacobianFunction {

	public static class hTag extends BaseTag<Object, Object, Object> {

		public hTag(final Object tag) {
			super("h", tag);
		}
	}

	public static List<? extends Object> buildHTags(//
			final List<? extends Object> outtags //
	) {
		final List<Object> res = new ArrayList<>();
		for (final Object tag : outtags)
			res.add(new hTag(tag));
		return res;
	}

	private final NamedMultivariateJacobianFunction mHFunction;

	private static List<? extends Object> buildInputs(//
			final NamedMultivariateJacobianFunction h, //
			final List<? extends Object> outputTags //
	) {
		final List<? extends Object> res = new ArrayList<>(h.getInputTags());
		res.removeAll(outputTags);
		return res;
	}

	public ImplicitMeasurementModel(final NamedMultivariateJacobianFunction h,
			final List<? extends Object> outputTags) {
		super(buildInputs(h, outputTags), outputTags);
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
			hPoint.setEntry(r, getValue(mHFunction.getInputTags().get(r), point));
		final Pair<RealVector, RealMatrix> hpr = mHFunction.evaluate(hPoint);
		// hpr.getFirst() values should all be close to zero...
		final RealMatrix hrm = hpr.getSecond();
		final List<? extends Object> hIn = mHFunction.getInputTags();
		final int outDim = getOutputDimension();
		final int inDim = getInputDimension();
		final RealMatrix mcy = MatrixUtils.createRealMatrix(outDim, outDim);
		for (int c = 0; c < outDim; ++c) {
			final int cIdx = hIn.indexOf(getOutputTags().get(c));
			for (int r = 0; r < outDim; ++r)
				mcy.setEntry(r, c, hrm.getEntry(r, cIdx));
		}
		final RealMatrix mcyInv = MatrixUtils.inverse(mcy);
		final RealMatrix mcx = MatrixUtils.createRealMatrix(outDim, inDim);
		for (int c = 0; c < inDim; ++c) {
			final int cIdx = hIn.indexOf(getInputTags().get(c));
			for (int r = 0; r < outDim; ++r)
				mcx.setEntry(r, c, hrm.getEntry(r, cIdx));
		}
		final RealVector rv = new ArrayRealVector(outDim);
		for (int r = 0; r < rv.getDimension(); ++r)
			rv.setEntry(r, getConstant(getOutputTags().get(r)));
		return Pair.create(rv, mcyInv.multiply(mcx));
	}
}
