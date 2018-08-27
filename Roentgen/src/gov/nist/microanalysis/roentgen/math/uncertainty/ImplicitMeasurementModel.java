package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.Report;

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
abstract public class ImplicitMeasurementModel extends NamedMultivariateJacobianFunctionEx {

	private final NamedMultivariateJacobianFunctionEx mCy;
	private final NamedMultivariateJacobianFunctionEx mCx;

	private static final boolean NAIVE = true;

	private static Map<Object, Double> buildTempConstants(List<? extends Object> inps) {
		Map<Object, Double> res = new HashMap<>();
		for (Object tag : inps)
			res.put(tag, Double.valueOf(Double.NaN));
		return res;
	}

	/**
	 * <li>C<sub>y</sub> = &delta;h/&delta;C<sub>y</sub>
	 * <li>C<sub>x</sub> = &delta;h/&delta;C<sub>x</sub>
	 * 
	 * @param cy
	 * @param cx
	 * @param outputTags
	 */

	public ImplicitMeasurementModel( //
			final NamedMultivariateJacobianFunctionEx cy, //
			final NamedMultivariateJacobianFunctionEx cx, //
			final List<? extends Object> outputTags) {
		super(cx.getInputTags(), cy.getInputTags());
		initializeConstants(buildTempConstants(cy.getInputTags()));
		mCy = cy;
		mCy.initializeConstants(buildTempConstants(cx.getInputTags()));
		mCx = cx;
		mCx.initializeConstants(buildTempConstants(cy.getInputTags()));
	}

	private RealVector constantsToPoint(NamedMultivariateJacobianFunctionEx nmjfe) {
		final List<? extends Object> inpTags = nmjfe.getInputTags();
		final int sz = inpTags.size();
		final RealVector res = new ArrayRealVector(sz);
		for (int i = 0; i < sz; ++i)
			res.setEntry(i, getConstant(inpTags.get(i)));
		return res;
	}

	private static Map<Object, Double> pointToConstants(NamedMultivariateJacobianFunctionEx nmjfe, RealVector point) {
		Map<Object, Double> res = new HashMap<>();
		final List<? extends Object> inpTags = nmjfe.getInputTags();
		for (int i = 0; i < point.getDimension(); ++i)
			res.put(inpTags.get(i), point.getEntry(i));
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
		assert mCy.getInputTags().equals(this.getOutputTags());

		final List<? extends Object> xExtra = mCy.getInputTags();
		assert Collections.disjoint(getInputTags(), xExtra);
		assert getConstants().keySet().containsAll(xExtra);

		final RealVector yPoint = constantsToPoint(mCy);
		mCy.initializeConstants(pointToConstants(this, point));
		final Pair<RealVector, RealMatrix> pcy = mCy.evaluate(yPoint);

		final int xInDim = mCx.getInputDimension();
		final List<? extends Object> xInTags = mCx.getInputTags();
		final RealVector xPoint = new ArrayRealVector(xInDim);
		for (int i = 0; i < xInDim; ++i) {
			final double val = getValue(xInTags.get(i), point);
			xPoint.setEntry(i, val);
		}
		mCx.initializeConstants(this.getConstants());
		final Pair<RealVector, RealMatrix> pcx = mCx.evaluate(xPoint);

		final RealVector vx = pcx.getFirst(), vy = pcy.getFirst();
		assert vx.getDistance(vy) / (vx.getNorm() + vy.getNorm()) < 1.0e-8 : vx.toString() + " != " + vy.toString();
		if (NAIVE) {
			final RealMatrix cyi = MatrixUtils.inverse(pcy.getSecond());
			final RealMatrix cx = pcx.getSecond();
			final RealMatrix rm = cyi.multiply(cx);
			assert yPoint.getDimension() == getOutputDimension();
			assert rm.getColumnDimension() == getInputDimension();
			assert rm.getRowDimension() == getOutputDimension();
			return Pair.create(yPoint, rm);
		} else {
			assert false;
		}
		return null;
	}

	public String toHTML(Mode mode) {
		Report r = new Report("Implicit Measurement Model");
		r.addHeader("Implicit Measurement Model - " + this.toString());
		r.addHTML(
				"<p>Solve C<sub>y</sub> U<sub>y</sub> C<sub>y</sub><sup>T</sup> = C<sub>x</sub> U<sub>x</sub> C<sub>x</sub><sup>T</sup> ");
		r.addHTML("for U<sub>y</sub> where <b>h</b>(<b>X</b>,<b>Y</b>)=<b>0</b>, ");
		r.addHTML("C<sub>y</sub> is the Jacobian of <b>h</b> with respect to <b>Y</b>, and ");
		r.addHTML("C<sub>x</sub> is the Jacobian of <b>h</b> with respect to <b>X</b>.</p>");
		r.addHTML(super.toHTML(mode));
		r.addSubHeader("C<sub>y</sub> Model");
		r.add(mCy);
		r.addSubHeader("C<sub>x</sub> Model");
		r.add(mCx);
		return r.toHTML(mode);
	}
}
