package gov.nist.microanalysis.roentgen.tests.math;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import junit.framework.TestCase;

/**
 * @author Nicholas
 */
public class TestCompositeMeasurementModel extends TestCase {

	private static final class Step2 extends ExplicitMeasurementModel<String, String> {

		private Step2(
				final List<String> inputLabels, //
				final List<String> outputLabels
		) throws ArgumentException {
			super(inputLabels, outputLabels);
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final double x0 = point.getEntry(0);
			final double x1 = point.getEntry(1);
			final double x2 = point.getEntry(2);
			final double f1_0 = point.getEntry(3);
			final double f1_1 = point.getEntry(4);
			final double f1_2 = point.getEntry(5);
			final RealVector fv = new ArrayRealVector(2);
			final RealMatrix jac = MatrixUtils.createRealMatrix(2, 6);
			fv.setEntry(0, Math.log(7.0 * f1_1 + f1_2) + 3.0 * x1 + 11.0 * x2);
			jac.setEntry(0, 1, 3.0);
			jac.setEntry(0, 2, 11.0);
			jac.setEntry(0, 4, (1.0 / (7.0 * f1_1 + f1_2)) * 7.0);
			jac.setEntry(0, 5, (1.0 / (7.0 * f1_1 + f1_2)) * 1.0);

			fv.setEntry(1, Math.log(f1_2 + 3.0 * f1_0) + 7.0 * x0 + 2.0 * x2);
			jac.setEntry(1, 0, 7.0);
			jac.setEntry(1, 2, 2.0);
			jac.setEntry(1, 3, (1.0 / (f1_2 + 3.0 * f1_0)) * 3.0);
			jac.setEntry(1, 5, (1.0 / (f1_2 + 3.0 * f1_0)) * 1.0);

			return Pair.create(fv, jac);
		}
	}

	private static final class Step1 extends ExplicitMeasurementModel<String, String> {

		private Step1(
				final List<String> inputLabels, //
				final List<String> outputLabels
		) throws ArgumentException {
			super(inputLabels, outputLabels);
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final RealVector fv = new ArrayRealVector(getOutputDimension());
			final RealMatrix jac = MatrixUtils.createRealMatrix(getOutputDimension(), getOutputDimension());
			final int ix0 = inputIndex("X0");
			final int ix1 = inputIndex("X1");
			final int ix2 = inputIndex("X2");
			final int iF1_0 = outputIndex("F1_0");
			fv.setEntry(iF1_0, 13.0 + 5.0 * Math.pow(point.getEntry(ix0), 2.7)
					+ 7.0 * Math.pow(point.getEntry(ix1), 1.3) + 11.0 * Math.pow(point.getEntry(ix2), -3.2));
			jac.setEntry(iF1_0, ix0, 5.0 * 2.7 * Math.pow(point.getEntry(ix0), 1.7));
			jac.setEntry(iF1_0, ix1, 7.0 * 1.3 * Math.pow(point.getEntry(ix1), 0.3));
			jac.setEntry(iF1_0, ix2, 11.0 * (-3.2) * Math.pow(point.getEntry(ix2), -4.2));

			final int iF1_1 = outputIndex("F1_1");
			fv.setEntry(iF1_1, 2.0 + 4.0 * Math.pow(point.getEntry(ix0), 2.0) + 6.0 * Math.pow(point.getEntry(ix1), 2.0)
					+ 8.0 * Math.pow(point.getEntry(ix2), 2.0));
			jac.setEntry(iF1_1, ix0, 2.0 * 4.0 * point.getEntry(ix0));
			jac.setEntry(iF1_1, ix1, 2.0 * 6.0 * point.getEntry(ix1));
			jac.setEntry(iF1_1, ix2, 2.0 * 8.0 * point.getEntry(ix2));

			final int iF1_2 = outputIndex("F1_2");
			fv.setEntry(iF1_2, 1.0 + 0.5 * Math.pow(point.getEntry(ix0), 2.0)
					+ 0.25 * Math.pow(point.getEntry(ix1), 4.0) + 0.125 * Math.pow(point.getEntry(ix2), 8.0));
			jac.setEntry(iF1_2, ix0, point.getEntry(ix0));
			jac.setEntry(iF1_2, ix1, Math.pow(point.getEntry(ix1), 3.0));
			jac.setEntry(iF1_2, ix2, Math.pow(point.getEntry(ix2), 7.0));
			return Pair.create(fv, jac);
		}
	}

	private final List<String> mInputs = Arrays.asList("X0", "X1", "X2");
	private final List<String> mOut1 = Arrays.asList("F1_0", "F1_1", "F1_2");
	private final List<String> mInputs2 = Arrays.asList("X0", "X1", "X2", "F1_0", "F1_1", "F1_2");
	private final List<String> mOut2 = Arrays.asList("F2_0", "F2_1");

	private final RealVector mValues0 = new ArrayRealVector(new double[] { 2.3, 3.5, 1.3 });
	private final RealMatrix mCov0 = new Array2DRowRealMatrix(
			new double[][] { { 0.1, 0.03, -0.013 }, { 0.03, 0.7, 0.002 }, { -0.013, 0.002, 0.3 } });

	private final UncertainValues<String> mInput0 = new UncertainValues<>(mInputs, mValues0, mCov0);

	private final RealVector mValues1 = new ArrayRealVector(new double[] { 2.3, 3.5, 1.3 });
	private final RealMatrix mCov1 = new Array2DRowRealMatrix(
			new double[][] { { 0.01, 0.003, -0.0013 }, { 0.003, 0.07, 0.0002 }, { -0.0013, 0.0002, 0.03 } });

	private final UncertainValues<String> mInput1 = new UncertainValues<>(mInputs, mValues1, mCov1);

	public void test1() throws IOException, ArgumentException {
		final Step1 step1 = new Step1(mInputs, mOut1);

		final UncertainValuesBase<String> uv = UncertainValuesBase.propagateAnalytical(step1, mInput0);
		final Report rep = new Report("Step 1");
		rep.addHeader("Inputs");
		rep.add(mInput0);
		// rep.addHeader("Jacobian");
		// rep.add(new LabeledMultivariateJacobian(mStep1, mValues0));
		rep.addHeader("Outputs 1");
		rep.add(uv);
		rep.inBrowser(Mode.NORMAL);
		// Calculated with Mathematica notebook MultiStepNMJF test1.nb
		final double[] rv = new double[] { 100.812, 110.18, 42.1803 };
		final double[][] rc = new double[][] { { 533.887, 483.696, 455.926 }, { 483.696, 1438.36, 1330.66 },
				{ 455.926, 1330.66, 1305.74 } };
		final UncertainValues<String> res = new UncertainValues<>(mOut1, new ArrayRealVector(rv),
				MatrixUtils.createRealMatrix(rc));
		assertTrue(uv.extract(mOut1).equals(res, 0.01));
	}

	public void test2() throws IOException, ArgumentException {
		final Step1 step1 = new Step1(mInputs, mOut1);
		final Step2 step2 = new Step2(mInputs2, mOut2);

		final List<ExplicitMeasurementModel<String, String>> steps = Arrays.asList(step1, step2);
		final CompositeMeasurementModel<String> msnmjf = //
				new CompositeMeasurementModel<String>("Test1", steps);
		final UncertainValues<String> uv = UncertainValues
				.asUncertainValues(UncertainValuesBase.propagateAnalytical(msnmjf, mInput1));
		final UncertainValuesBase<String> mc = UncertainValuesBase.propagateMonteCarlo(msnmjf, mInput1, 100000);
		final UncertainValuesCalculator<String> delta = UncertainValuesBase.propagateFiniteDifference(msnmjf, mInput1,
				0.001);
		final Report rep = new Report("Step 1 and 2");
		rep.addHeader("Inputs");
		rep.add(mInput1);
		// rep.addHeader("Jacobian");
		// rep.add(LabeledMultivariateJacobian.compute(msnmjf, mValues0));
		rep.addHeader("Outputs 1 and 2");
		rep.add(uv, Mode.VERBOSE);
		rep.addHeader("MC Outputs 1 and 2");
		rep.add(mc, Mode.VERBOSE);
		rep.addHeader("Delta Outputs 1 and 2");
		rep.add(delta, Mode.VERBOSE);
		rep.inBrowser(Mode.NORMAL);

		assertTrue(delta.extract(uv.getLabels()).equals(uv, 0.001));
		assertTrue(mc.extract(uv.getLabels()).equals(uv, 0.5));
	}

	public void test3() throws IOException, ArgumentException {
		final Step1 step1 = new Step1(mInputs, mOut1);
		final Step2 step2 = new Step2(mInputs2, mOut2);

		final List<ExplicitMeasurementModel<String, String>> steps = Arrays.asList(step1, step2);

		final List<String> outputs = Arrays.asList("F2_0", "F1_0", "F1_2");

		final CompositeMeasurementModel<String> trimmed = //
				new CompositeMeasurementModel<String>("Trimmed", steps, outputs);

		final CompositeMeasurementModel<String> full = //
				new CompositeMeasurementModel<String>("Full", steps);

		final UncertainValuesBase<String> trimRes = UncertainValuesBase.propagateAnalytical(trimmed, mInput1);
		final UncertainValuesBase<String> fullRes = UncertainValuesBase.propagateAnalytical(full, mInput1);

		final Report rep = new Report("Step 1 and 2");
		rep.addHeader("Trimmed");
		rep.add(trimRes, Mode.VERBOSE);
		rep.addHeader("Full");
		rep.add(fullRes, Mode.VERBOSE);
		rep.inBrowser(Mode.NORMAL);

		for (final String rLbl : trimRes.getLabels()) {
			testValues(trimRes.getEntry(rLbl), fullRes.getEntry(rLbl), 0.001);
			testValues(trimRes.getUncertainty(rLbl), fullRes.getUncertainty(rLbl), 0.001);
			for (final String cLbl : trimRes.getLabels())
				testValues(trimRes.getCorrelationCoefficient(rLbl, cLbl), //
						fullRes.getCorrelationCoefficient(rLbl, cLbl), 0.01);
		}
	}

	public void testValues(
			final double v1, final double v2, final double eps
			) {
		if (Double.isFinite((v1 - v2) / v2))
			assertEquals(0.0, Math.abs((v1 - v2) / v1), eps);
		else
			assertEquals(v1, v2, eps);
	}

}
