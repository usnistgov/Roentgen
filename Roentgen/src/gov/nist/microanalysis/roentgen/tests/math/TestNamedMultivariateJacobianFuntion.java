package gov.nist.microanalysis.roentgen.tests.math;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import junit.framework.TestCase;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestNamedMultivariateJacobianFuntion extends TestCase {

	public void test1() {
		final List<Object> funcs = Arrays.asList("f1", "f2", "f3");
		final List<Object> vars = Arrays.asList("x0", "x1", "x2", "x3", "x4");
		final LabeledMultivariateJacobianFunction nmvf = new LabeledMultivariateJacobianFunction(vars, funcs) {

			public double f1(final RealVector point) {
				final double x0 = point.getEntry(0);
				final double x1 = point.getEntry(1);
				return x0 + 2.0 * x1 * x1;
			}

			public double f2(final RealVector point) {
				final double x2 = point.getEntry(2);
				final double x4 = point.getEntry(4);
				return x2 + 1.3 * x4 * x4;
			}

			public double f3(final RealVector point) {
				final double x0 = point.getEntry(0);
				final double x3 = point.getEntry(3);
				final double x1 = point.getEntry(1);
				return x0 * x1 + 2.0 * x3 * x3;
			}

			@Override
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				final RealVector vals = MatrixUtils.createRealVector(new double[] { f1(point), f2(point), f3(point) });
				final double[] x = point.toArray();
				final RealMatrix jacob = new Array2DRowRealMatrix(new double[][] { { 1.0, 4.0 * x[1], 0.0, 0.0, 0.0 },
						{ 0.0, 0.0, 1.0, 0.0, 2.6 * x[4] }, { x[1], x[0], 0.0, 4.0 * x[3], 0.0 } });
				return Pair.create(vals, jacob);
			}
		};
		final RealVector point = MatrixUtils.createRealVector(new double[] { 2.0, 3.0, 1.5, 2.3, -3.4 });
		final Pair<RealVector, RealMatrix> res = nmvf.value(point);
		final RealVector val = res.getFirst();
		assertEquals(20.0, val.getEntry(nmvf.outputIndex("f1")), 1.0e-6);
		assertEquals(16.528, val.getEntry(nmvf.outputIndex("f2")), 1.0e-6);
		assertEquals(16.58, val.getEntry(2), 1.0e-6);

		final RealMatrix jacob = res.getSecond();
		assertEquals(1.0, jacob.getEntry(0, 0), 1.e-6);
		assertEquals(12.0, jacob.getEntry(nmvf.outputIndex("f1"), nmvf.inputIndex("x1")), 1.e-6);
		assertEquals(0.0, jacob.getEntry(nmvf.outputIndex("f1"), nmvf.inputIndex("x3")), 1.e-6);

		assertEquals(1.0, jacob.getEntry(1, 2), 1.e-6);
		assertEquals(-8.84, jacob.getEntry(nmvf.outputIndex("f2"), nmvf.inputIndex("x4")), 1.e-6);
		assertEquals(0.0, jacob.getEntry(nmvf.outputIndex("f2"), nmvf.inputIndex("x3")), 1.e-6);

		assertEquals(9.2, jacob.getEntry(2, 3), 1.e-6);
		assertEquals(2.0, jacob.getEntry(nmvf.outputIndex("f3"), nmvf.inputIndex("x1")), 1.e-6);
		assertEquals(0.0, jacob.getEntry(nmvf.outputIndex("f3"), nmvf.inputIndex("x4")), 1.e-6);
	}

}
