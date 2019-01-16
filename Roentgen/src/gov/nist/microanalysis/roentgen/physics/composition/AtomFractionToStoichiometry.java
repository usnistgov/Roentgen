/**
 * 
 */
package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;

/**
 * @author nicho
 *
 */
public class AtomFractionToStoichiometry extends LabeledMultivariateJacobianFunction {

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public AtomFractionToStoichiometry(List<? extends Object> inputLabels, List<? extends Object> outputLabels) {
		super(inputLabels, outputLabels);
		// TODO Auto-generated constructor stub
	}

	/* (non-Javadoc)
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		// TODO Auto-generated method stub
		return null;
	}

}
