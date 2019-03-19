/**
 * 
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.List;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

/**
 * @author nicho
 *
 */
public class RetainInputs<L> extends LabeledMultivariateJacobianFunction<L, L> {

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public RetainInputs(List<L> inputLabels) {
		super(inputLabels, inputLabels);
	}

	/* (non-Javadoc)
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		return Pair.create(point, MatrixUtils.createRealIdentityMatrix(point.getDimension()));
	}

}
