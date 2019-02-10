package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Set;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * An interface defining a minimalist set of properties
 * associated with a material.
 * 
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public interface IMaterial extends IToHTML {
	
	public Set<Element> getElementSet();
	
	public String getHTMLName();

}
