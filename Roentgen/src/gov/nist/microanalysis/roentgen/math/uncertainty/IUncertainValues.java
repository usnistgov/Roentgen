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
public interface IUncertainValues {
	
	/**
	 * Returns a {@link RealVector} containing the function values associated with this
	 * uncertainty calculation.
	 *
	 * @return {@link RealVector}
	 */
	public RealVector getValues();
	
	
	/**
	 * Returns the matrix containing the covariances associate with this uncertainty calculation.
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
	public List<? extends Object> getLabels();
	
	public boolean hasEntry(Object lbl);
}
