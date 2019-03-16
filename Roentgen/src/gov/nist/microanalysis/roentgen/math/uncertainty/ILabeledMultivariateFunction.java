package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.List;

import org.apache.commons.math3.linear.RealVector;

/**
 * Meant as a complement to LabeledMultivariateJacobianFunction. When
 * implemented should compute identically the same value as first return value
 * from {@link ILabeledMultivariateFunction}.value(...). This acts as an
 * optimization for those cases in which the Jacobian is not required as would
 * be the case when evaluating the function to implement an approximate
 * derivative or a Monte Carlo evaluation.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public interface ILabeledMultivariateFunction {

	/**
	 * Evaluate a multivariate function at the point specified. Return the output
	 * values as a RealVector.
	 *
	 *
	 * @param point
	 * @return {@link RealVector}
	 */
	public RealVector optimized(RealVector point);	
	
	/**
	 * A list containing the labels identifying the input values to optimized(...)
	 * 
	 * @return List&lt;? extends Object&gt;
	 */
	public List<? extends Object> getInputLabels();
	
	/**
	 * A list containing the labels identifying the output values from optimized(...)
	 * 
	 * @return List&lt;? extends Object&gt;
	 */
	public List<? extends Object> getOutputLabels();
}
