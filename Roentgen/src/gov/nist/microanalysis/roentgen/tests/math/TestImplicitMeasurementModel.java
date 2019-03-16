package gov.nist.microanalysis.roentgen.tests.math;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel.HLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;

/**
 * Not ready!!! Need test inputs....
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class TestImplicitMeasurementModel {

	public static class Pressure102Example //
			extends ImplicitMeasurementModel.HModel implements ILabeledMultivariateFunction {

		private static final String PRESSURE = "p";
		private static final String GL = "gl";
		private static final String RHOA = "&rho;<sub>a</sub>";
		private static final String ALPHA = "&alpha;";
		private static final String RHOW = "&rho;<sub>w</sub>";
		private static final String LAMBDA = "&lambda;";
		private static final String A0 = "A0";
		private static final String MW = "m<sub>w</sub>";
		private static final String DTHETA = "&delta;&Theta;";

		static private String idx(String base, int i) {
			return base + "[" + i + "]";
		}

		/**
		 * Build the input (X) labels (used by Pressure102Example constructor)
		 * 
		 * @param len The number of pressure measurements
		 * @return List<String> A list of labels
		 */
		static private List<String> xLabels(int len) {
			List<String> res = new ArrayList<>();
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
		static private List<String> yLabels(int len) {
			List<String> res = new ArrayList<>();
			for (int i = 0; i < len; ++i)
				res.add(idx(PRESSURE, i));
			return res;
		}

		/**
		 * The number of pressure measurements
		 */
		private final int mLength;

		public Pressure102Example(int len) {
			super(xLabels(len), yLabels(len));
			mLength = len;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
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
				final int hIdx = outputIndex(new HLabel(idx(PRESSURE, i)));
				val.setEntry(hIdx, a0 * p[i] * (1.0 + lambda * p[i]) * (1.0 + alpha * dtheta[i])
						- mw[i] * (1.0 - rhoa / rhow) * gl);
				for (int j = 0; j < mLength; ++i) {
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
			return Pair.create(val, jac);
		}

		@Override
		public RealVector optimized(RealVector point) {
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
			for (int i = 0; i < mLength; ++i) {
				final int hIdx = outputIndex(new HLabel(PRESSURE + "[" + i + "]"));
				val.setEntry(hIdx, a0 * p[i] * (1.0 + lambda * p[i]) * (1.0 + alpha * dtheta[i])
						- mw[i] * (1.0 - rhoa / rhow) * gl);
			}
			return val;
		}
	}


	@Test
	public void Example102() throws ArgumentException {

		Pressure102Example hmodel = new Pressure102Example(3);
		ImplicitMeasurementModel imm = new ImplicitMeasurementModel(hmodel, hmodel.getOutputLabels());

		UncertainValues inputs = new UncertainValues(imm.getInputLabels());

		UncertainValuesCalculator uvc = new UncertainValuesCalculator(imm, inputs);
		UncertainValues jres = UncertainValues.force(uvc);
		
		RealVector dinp = inputs.getValues().mapMultiply(0.001);
		uvc.setCalculator(new UncertainValuesCalculator.FiniteDifference(dinp));
		UncertainValues fdres = UncertainValues.force(uvc);
		
		assertTrue("Finite difference does not equal Jacobian", jres.equals(fdres, 0.00001));
		
	}

}
