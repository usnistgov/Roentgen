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
 * <p>
 * To conform to the implementation of explicit measurement model using the
 * {@link NamedMultivariateJacobianFunction}, this class requires the use of
 * constants as implemented in {@link NamedMultivariateJacobianFunctionEx} to
 * pass in the values associated with the output variables (the y values).
 * </p>
 *
 *
 *
 *
 * @author Nicholas
 *
 */
abstract public class ImplicitMeasurementModel //
		extends NamedMultivariateJacobianFunction {

	private final NamedMultivariateJacobianFunction mCy;
	private final NamedMultivariateJacobianFunction mCx;

	private static final boolean NAIVE = true;

	private static Map<Object, Double> buildTempConstants(final List<? extends Object> inps) {
		final Map<Object, Double> res = new HashMap<>();
		for (final Object tag : inps)
			res.put(tag, Double.valueOf(Double.NaN));
		return res;
	}

	/**
	 * <li>C<sub>y</sub> = &delta;h/&delta;C<sub>y</sub>
	 * <li>C<sub>x</sub> = &delta;h/&delta;C<sub>x</sub>
	 *
	 * @param cy
	 * @param cx
	 */

	public ImplicitMeasurementModel( //
			final NamedMultivariateJacobianFunction cy, //
			final NamedMultivariateJacobianFunction cx) {
		super(cx.getInputTags(), cy.getInputTags());
		initializeConstants(buildTempConstants(cy.getInputTags()));
		mCy = cy;
		mCy.initializeConstants(buildTempConstants(cx.getInputTags()));
		mCx = cx;
		mCx.initializeConstants(buildTempConstants(cy.getInputTags()));
	}

	private RealVector constantsToPoint(final NamedMultivariateJacobianFunction nmjfe) {
		final List<? extends Object> inpTags = nmjfe.getInputTags();
		final int sz = inpTags.size();
		final RealVector res = new ArrayRealVector(sz);
		for (int i = 0; i < sz; ++i)
			res.setEntry(i, getConstant(inpTags.get(i)));
		return res;
	}

	private static Map<Object, Double> pointToConstants(final NamedMultivariateJacobianFunction nmjfe,
			final RealVector point) {
		final Map<Object, Double> res = new HashMap<>();
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
			// 1. Form the Cholesky factor Rx of Vx, i.e., the upper triangular matrix such
			// that RxT.Rx = Vx.
			// final CholeskyDecomposition chyd = new
			// CholeskyDecomposition(uv.getCovariances());
			// final RealMatrix Rx = chyd.getL();
			// 2. Factor Jx as the product Jx = Qx.Ux, where Qx is an orthogonal matrix and
			// Ux is upper triangular.
			// final QRDecomposition qrd = new QRDecomposition(jXe.getSecond());
			// final RealMatrix Qx = qrd.getQ(), Ux = qrd.getR();
			// 3. Factor Jy as the product Jy = Ly.Uy, where Ly is lower triangular and Uy
			// is upper triangular.
			// final LUDecomposition lud = new LUDecomposition(jYe.getSecond());
			// ????What about pivoting?????
			// final RealMatrix Ly = lud.getL(), Uy = lud.getU();
			// 4. Solve the matrix equation UyT.M1 = I for M1.
			// final RealMatrix M1 = MatrixUtils.inverse(Uy.transpose());
			// 5. Solve LyT.M2 = M1 for M2,
			// final RealMatrix M2 = MatrixUtils.inverse(Ly.transpose()).multiply(M1);
			// 6. Form M3 = QxT.M2.
			// final RealMatrix M3 = Qx.transpose().multiply(M2);
			// 7. Form M4 = UxT.M3.
			// final RealMatrix M4 = Ux.transpose().multiply(M3);
			// 8. Form M = Rx.M4.
			// final RealMatrix M = Rx.multiply(M4);
			// 9. Orthogonally triangularize M to give the upper triangular matrix R.
			// final RealMatrix R = new QRDecomposition(M).getR();
			// 10. Form Vy = RT.R.
			// return Pair.create(yPoint, (R.transpose()).multiply(R));
		}
		return null;

	}

	@Override
	public String toHTML(final Mode mode) {
		final Report r = new Report("Implicit Measurement Model");
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
