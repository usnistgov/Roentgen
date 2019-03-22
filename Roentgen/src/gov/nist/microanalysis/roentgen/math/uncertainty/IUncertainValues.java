/**
 * 
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * @author nicho
 *
 */
public interface IUncertainValues<H> {
	
	/**
	 * Returns a {@link RealVector} containing the L-values.
	 *
	 * @return {@link RealVector}
	 */
	public RealVector getValues();
	
	
	/**
	 * Returns an L x L matrix containing the covariances.
	 * 
	 * @return RealMatrix
	 */
	public RealMatrix getCovariances();

	
	/**
	 * Returns an ordered list of the Object labels associated with the values and
	 * covariances in this object.
	 *
	 * @return List&lt;Object&gt;
	 */
	public List<H> getLabels();
	
	public boolean hasLabel(H lbl);
}
