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
 * @author nicho
 *
 */
public class SecondaryFluorescenceModel extends LabeledMultivariateJacobianFunction {

	public static class SecondaryFluorescenceLabel extends BaseLabel<MatrixCorrectionDatum, AtomicShell, Object>{

		public SecondaryFluorescenceLabel(MatrixCorrectionDatum mcd, AtomicShell shell) {
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
	public SecondaryFluorescenceModel(List<? extends Object> inputLabels, List<? extends Object> outputLabels) {
		super(inputLabels, outputLabels);
		// TODO Auto-generated constructor stub
	}

	/* (non-Javadoc)
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		RealVector rv = new ArrayRealVector(new double[] { 1.0 });
		RealMatrix rm = MatrixUtils.createRealMatrix(new double[][] { { 0.02*0.02 } });
		return Pair.create(rv, rm);
	}

}
