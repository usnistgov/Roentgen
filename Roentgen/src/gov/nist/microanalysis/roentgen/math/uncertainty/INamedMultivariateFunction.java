package gov.nist.microanalysis.roentgen.math.uncertainty;

import org.apache.commons.math3.linear.RealVector;

/**
 * Meant as a complement to NamedMultivariateJacobianFunction. When implemented
 * should compute identically the same value as first return value from
 * {@link INamedMultivariateFunction}.value(...). This acts as an optimization
 * for those cases in which the Jacobian is not required as would be the case
 * when evaluating the function to implement an approximate derivative or
 * a Monte Carlo evaluation.
 * 
 * @author nicholas
 *
 */
public interface INamedMultivariateFunction {

	/**
	 * Evaluate a multivariate function at the point specified. Return the output
	 * values as a RealVector.
	 * 
	 * 
	 * @param point
	 * @return {@link RealVector}
	 */
	public RealVector optimized(RealVector point);

}
