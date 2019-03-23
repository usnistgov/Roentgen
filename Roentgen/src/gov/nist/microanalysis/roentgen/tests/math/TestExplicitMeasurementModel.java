package gov.nist.microanalysis.roentgen.tests.math;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import junit.framework.TestCase;

/**
 * <p>
 * Tests {@link ExplicitMeasurementModel} directly and then {@link UncertainValuesCalculator}.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestExplicitMeasurementModel extends TestCase {

	public void test1() throws ArgumentException, IOException {
		final List<String> funcs = Arrays.asList("f1", "f2", "f3");
		final List<String> vars = Arrays.asList("x0", "x1", "x2", "x3", "x4");
		final ExplicitMeasurementModel<String, String> nmvf = new ExplicitMeasurementModel<String, String>(vars,
				funcs) {

			public double f1(
					final RealVector point
			) {
				final double x0 = point.getEntry(0);
				final double x1 = point.getEntry(1);
				return x0 + 2.0 * x1 * x1;
			}

			public double f2(
					final RealVector point
			) {
				final double x2 = point.getEntry(2);
				final double x4 = point.getEntry(4);
				return x2 + 1.3 * x4 * x4;
			}

			public double f3(
					final RealVector point
			) {
				final double x0 = point.getEntry(0);
				final double x3 = point.getEntry(3);
				final double x1 = point.getEntry(1);
				return x0 * x1 + 2.0 * x3 * x3;
			}

			@Override
			public Pair<RealVector, RealMatrix> value(
					final RealVector point
			) {
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

		// 2.0, 3.0, 1.5, 2.3, -3.4
		final RealMatrix covar = MatrixUtils.createRealMatrix(new double[][] {
				{ 0.01 * 0.01, -0.5 * 0.01 * 0.03, 0.3 * 0.01 * 0.15, 0.5 * 0.01 * 0.11, 0.4 * 0.01 * 0.02 },
				{ -0.5 * 0.01 * 0.03, 0.03 * 0.03, -0.2 * 0.03 * 0.15, 0.1 * 0.03 * 0.11, 0.2 * 0.03 * 0.02 },
				{ 0.3 * 0.01 * 0.15, -0.2 * 0.03 * 0.15, 0.15 * 0.15, -0.2 * 0.15 * 0.11, -0.1 * 0.15 * 0.02 },
				{ 0.5 * 0.01 * 0.11, 0.1 * 0.03 * 0.11, -0.2 * 0.15 * 0.11, 0.11 * 0.11, 0.9 * 0.11 * 0.02 },
				{ 0.4 * 0.01 * 0.02, 0.2 * 0.03 * 0.02, -0.1 * 0.15 * 0.02, 0.9 * 0.11 * 0.02, 0.02 * 0.02 } });

		final UncertainValues<String> inputs = new UncertainValues<>(vars, point, covar);

		final UncertainValuesCalculator<String> uvca = UncertainValuesBase.propagateAnalytical(nmvf, inputs);
		final UncertainValuesCalculator<String> uvcf = UncertainValuesBase.propagateFiniteDifference(nmvf, inputs,
				0.001);
		final UncertainValuesCalculator<String> uvcm = UncertainValuesBase.propagateMonteCarlo(nmvf, inputs, 100000);

		final UncertainValuesCalculator<String>.MonteCarlo mc = (UncertainValuesCalculator<String>.MonteCarlo) uvcm
				.getCalculator();
		mc.setParallel(true);

		assertTrue(uvca.equals(uvcf, 0.0001));
		assertTrue(uvca.equals(uvcm, 0.02));

		final List<String> newOrder = Arrays.asList("x2", "x0", "x1", "x4", "x3");
		final UncertainValues<String> reorder = UncertainValues.asUncertainValues(inputs.reorder(newOrder));

		final UncertainValuesCalculator<String> uvcar = UncertainValuesBase.propagateAnalytical(nmvf, reorder);
		final UncertainValuesCalculator<String> uvcfr = UncertainValuesBase.propagateFiniteDifference(nmvf, reorder,
				0.001);
		final UncertainValuesCalculator<String> uvcmr = UncertainValuesBase.propagateMonteCarlo(nmvf, reorder, 100000);

		final UncertainValuesCalculator<String>.MonteCarlo mcr = (UncertainValuesCalculator<String>.MonteCarlo) uvcm
				.getCalculator();
		mcr.setParallel(false);

		final Report r = new Report("Explicit 1");
		r.addHeader("Explicit 1");
		r.addSubHeader("Analytical");
		r.add(uvca, Mode.VERBOSE);
		r.addSubHeader("Finite Difference");
		r.add(uvcf, Mode.VERBOSE);
		r.addSubHeader("Monte Carlo");
		r.add(uvcm, Mode.VERBOSE);

		r.addSubHeader("Analytical - Reorder");
		r.add(uvcar, Mode.VERBOSE);
		r.addSubHeader("Finite Difference - Reorder");
		r.add(uvcfr, Mode.VERBOSE);
		r.addSubHeader("Monte Carlo - Reorder");
		r.add(uvcmr, Mode.VERBOSE);
		r.inBrowser(Mode.NORMAL);

		assertTrue(uvca.equals(uvcar, 0.0001));
		assertTrue(uvca.equals(uvcfr, 0.0001));
		assertTrue(uvca.equals(uvcmr, 0.02));

	}

}
