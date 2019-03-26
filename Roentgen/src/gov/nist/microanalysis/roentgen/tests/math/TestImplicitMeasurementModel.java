package gov.nist.microanalysis.roentgen.tests.math;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import gov.nist.juncertainty.ImplicitMeasurementModel;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * Not ready!!! Need test inputs....
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class TestImplicitMeasurementModel {

	public static class Pressure102Example //
			extends ImplicitMeasurementModel<String, String> {

		private static final String PRESSURE = "p";
		private static final String GL = "gl";
		private static final String RHOA = "&rho;<sub>a</sub>";
		private static final String ALPHA = "&alpha;";
		private static final String RHOW = "&rho;<sub>w</sub>";
		private static final String LAMBDA = "&lambda;";
		private static final String A0 = "A0";
		private static final String MW = "m<sub>w</sub>";
		private static final String DTHETA = "&delta;&Theta;";

		static private String idx(
				final String base, final int i
		) {
			return base + "[" + i + "]";
		}

		/**
		 * Build the input (X) labels (used by Pressure102Example constructor)
		 *
		 * @param len The number of pressure measurements
		 * @return List<String> A list of labels
		 */
		static private List<String> xLabels(
				final int len
		) {
			final List<String> res = new ArrayList<>();
			for (int i = 0; i < len; ++i) {
				res.add(idx(DTHETA, i));
				res.add(idx(MW, i));
			}
			res.addAll(Arrays.asList(A0, LAMBDA, ALPHA, RHOA, RHOW, GL));
			return res;

		}

		/**
		 * Build the output (Y) labels (used by Pressure102Example constructor)
		 *
		 * @param len The number of pressure measurements
		 * @return List<String> A list of labels
		 */
		static private List<String> yLabels(
				final int len
		) {
			final List<String> res = new ArrayList<>();
			for (int i = 0; i < len; ++i)
				res.add(idx(PRESSURE, i));
			return res;
		}

		/**
		 * The number of pressure measurements
		 */
		private final int mLength;

		public Pressure102Example(
				final int len
		) throws ArgumentException {
			super(xLabels(len), yLabels(len));
			mLength = len;
		}

		@Override
		public RealMatrix computeCx(
				final RealVector point
		) {
			// Initialize input variables
			final double rhoa = getArg(RHOA, point);
			final double rhow = getArg(RHOW, point);
			final double gl = getArg(GL, point);
			final double a0 = getArg(A0, point);
			final double lambda = getArg(LAMBDA, point);
			;
			final double alpha = getArg(ALPHA, point);
			final double[] mw = new double[mLength];
			final double[] p = new double[mLength];
			final double[] dtheta = new double[mLength];
			for (int j = 0; j < mLength; ++j) {
				mw[j] = getArg(idx(MW, j), point);
				p[j] = getOutputValue(idx(PRESSURE, j));
				dtheta[j] = getArg(idx(DTHETA, j), point);
			}
			// Calculate output values
			final RealMatrix jac = buildEmptyCx();
			for (int i = 0; i < mLength; ++i) {
				final int hIdx = i;
				
				for (final int j = 0; j < mLength; ++i) {
					setCy(hIdx, idx(PRESSURE, j), jac, //
							a0 * (1.0 + 2.0 * lambda * p[i]) * (1.0 + alpha * dtheta[i]));
					setCx(hIdx, idx(DTHETA, j), jac, //
							a0 * p[i] * alpha * (1.0 + p[i] * lambda));
					setCx(hIdx, idx(MW, j), jac, //
							(gl * (rhoa - rhow)) / rhow);
				}
				setCx(hIdx, A0, jac, p[i] * (1.0 + alpha * dtheta[i]) * (1.0 + p[i] * lambda));
				setCx(hIdx, ALPHA, jac, a0 * p[i] * dtheta[i] * (1.0 + p[i] * lambda));
				setCx(hIdx, LAMBDA, jac, a0 * Math.pow(p[i], 2.0) * (1.0 + alpha * dtheta[i]));
				setCx(hIdx, RHOA, jac, (gl * mw[i]) / rhow);
				setCx(hIdx, RHOW, jac, -((gl * mw[i] * rhoa) / Math.pow(rhow, 2.0)));
				setCx(hIdx, GL, jac, (mw[i] * (rhoa - rhow)) / rhow);
			}
			return jac;
		}
		
		// val.setEntry(hIdx, a0 * p[i] * (1.0 + lambda * p[i]) * (1.0 + alpha * dtheta[i]) - mw[i] * (1.0 - rhoa / rhow) * gl);

		@Override
		public RealMatrix computeCy(
				final RealVector point
		) {
			// Initialize input variables
			final double rhoa = point.getEntry(inputIndex(RHOA));
			final double rhow = point.getEntry(inputIndex(RHOW));
			final double gl = point.getEntry(inputIndex(GL));
			final double a0 = point.getEntry(inputIndex(A0));
			final double lambda = point.getEntry(inputIndex(LAMBDA));
			final double alpha = point.getEntry(inputIndex(ALPHA));
			final double[] mw = new double[mLength];
			final double[] p = new double[mLength];
			final double[] dtheta = new double[mLength];
			for (int j = 0; j < mLength; ++j) {
				mw[j] = point.getEntry(inputIndex(idx(MW, j)));
				p[j] = point.getEntry(inputIndex(idx(PRESSURE, j)));
				dtheta[j] = point.getEntry(inputIndex(idx(DTHETA, j)));
			}
			// Calculate output values
			final RealVector val = new ArrayRealVector(getOutputDimension());
			final RealMatrix jac = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			for (int i = 0; i < mLength; ++i) {
				final int hIdx = i;
				val.setEntry(hIdx, a0 * p[i] * (1.0 + lambda * p[i]) * (1.0 + alpha * dtheta[i])
						- mw[i] * (1.0 - rhoa / rhow) * gl);
				for (final int j = 0; j < mLength; ++i) {
					jac.setEntry(hIdx, inputIndex(idx(PRESSURE, j)), //
							a0 * (1.0 + 2.0 * lambda * p[i]) * (1.0 + alpha * dtheta[i]));
					jac.setEntry(hIdx, inputIndex(idx(DTHETA, j)), //
							a0 * p[i] * alpha * (1.0 + p[i] * lambda));
					jac.setEntry(hIdx, inputIndex(idx(MW, j)), //
							(gl * (rhoa - rhow)) / rhow);
				}
				jac.setEntry(hIdx, inputIndex(A0), p[i] * (1.0 + alpha * dtheta[i]) * (1.0 + p[i] * lambda));
				jac.setEntry(hIdx, inputIndex(ALPHA), a0 * p[i] * dtheta[i] * (1.0 + p[i] * lambda));
				jac.setEntry(hIdx, inputIndex(LAMBDA), a0 * Math.pow(p[i], 2.0) * (1.0 + alpha * dtheta[i]));
				jac.setEntry(hIdx, inputIndex(RHOA), (gl * mw[i]) / rhow);
				jac.setEntry(hIdx, inputIndex(RHOW), -((gl * mw[i] * rhoa) / Math.pow(rhow, 2.0)));
				jac.setEntry(hIdx, inputIndex(GL), (mw[i] * (rhoa - rhow)) / rhow);
			}
			return jac;
		}

	}

	@Test
	public void Example102() throws ArgumentException {

		final Pressure102Example imm = new Pressure102Example(3);

		final UncertainValues<String> inputs = new UncertainValues<>(imm.getInputLabels());

		final UncertainValuesCalculator<String> uvc = new UncertainValuesCalculator<>(imm, inputs);
		final UncertainValues<String> jres = UncertainValues.<String>asUncertainValues(uvc);

		uvc.setCalculator(uvc.new FiniteDifference(0.001));
		final UncertainValues<String> fdres = UncertainValues.<String>asUncertainValues(uvc);

		assertTrue("Finite difference does not equal Jacobian", jres.equals(fdres, 0.00001));

	}

}
