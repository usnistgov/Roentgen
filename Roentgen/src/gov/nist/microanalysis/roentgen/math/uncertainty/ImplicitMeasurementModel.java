package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
 * @author Nicholas
 *
 */
abstract public class ImplicitMeasurementModel extends NamedMultivariateJacobianFunction {

	private final NamedMultivariateJacobianFunction mCy;
	private final NamedMultivariateJacobianFunction mCx;

	private static final boolean NAIVE = true;

	private static List<?> combine(final List<?>... los) {
		final Set<Object> res = new HashSet<>();
		for (final List<? extends Object> lo : los)
			res.addAll(lo);
		return new ArrayList<>(res);
	}
	
	/**
	 * <li>C<sub>y</sub> = &delta;h/&delta;C<sub>y</sub>
	 * <li>C<sub>x</sub> = &delta;h/&delta;C<sub>x</sub>
	 * @param cy
	 * @param cx
	 * @param outputTags
	 */

	public ImplicitMeasurementModel( //
			final NamedMultivariateJacobianFunction cy, //
			final NamedMultivariateJacobianFunction cx,
			final List<? extends Object> outputTags) {
		super(combine(cx.getInputTags(), cy.getInputTags()), outputTags);
		mCy = cy;
		mCx = cx;
	}

	/**
	 * Extracts from <code>point</code> the input values required to evaluate
	 * <code>nmvj</code>.
	 *
	 * @param point A RealVector of length getInputDimension()
	 * @param nmvj
	 * @return RealVector of length nmvj.getInputDimension()
	 */
	private RealVector extract(final RealVector point, final NamedMultivariateJacobianFunction nmvj) {
		final List<? extends Object> inpCx = nmvj.getInputTags();
		final RealVector res = new ArrayRealVector(inpCx.size());
		for (int i = 0; i < inpCx.size(); ++i) {
			final Object obj = inpCx.get(i);
			final int idx = inputIndex(obj);
			assert idx >= 0;
			res.setEntry(i, point.getEntry(idx));
		}
		return res;
	}

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
	 * @param point 
	 * @return Pair&lt;RealVector, RealMatrix&gt;
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector y = extract(point, mCy);
		final Pair<RealVector, RealMatrix> pcy = mCy.evaluate(extract(point, mCy));
		final Pair<RealVector, RealMatrix> pcx = mCx.evaluate(extract(point, mCx));
		final RealVector vx = pcx.getFirst(), vy = pcy.getFirst();
		assert vx.getDistance(vy) / (vx.getNorm() + vy.getNorm()) < 1.0e-8;
		if (NAIVE) {
			final RealMatrix cyi = MatrixUtils.inverse(pcy.getSecond());
			final RealMatrix cx = pcx.getSecond();
			final RealMatrix rm = cyi.multiply(cx);
			return Pair.create(y, rm);
		} else {
			assert false;
		}
		return null;
	}
}
