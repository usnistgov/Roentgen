package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;

/**
 * @author Nicholas WM Ritchie
 *
 */
public class SecondaryFluorescenceModel //
		extends LabeledMultivariateJacobianFunction {

	public static class SecondaryFluorescenceLabel //
			extends BaseLabel<MatrixCorrectionDatum, AtomicShell, Object> {

		public SecondaryFluorescenceLabel(final MatrixCorrectionDatum mcd, final AtomicShell shell) {
			super("Fs", mcd, shell);
		}

		public MatrixCorrectionDatum getDatum() {
			return getObject1();
		}

		public AtomicShell getShell() {
			return getObject2();
		}
	}

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public SecondaryFluorescenceModel(final List<? extends Object> inputLabels,
			final List<? extends Object> outputLabels) {
		super(inputLabels, outputLabels);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(new double[] { 1.0 });
		final double dFs = 1.0e-5;
		final RealMatrix rm = MatrixUtils.createRealMatrix(new double[][] { { dFs * dFs } });
		return Pair.create(rv, rm);
	}

}
