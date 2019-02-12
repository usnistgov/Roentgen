package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * The material class represents the most basic information about a material -
 * its name and the set of constituent elements. Implements the same interface
 * as Composition for those situations in which the mass fractions are not
 * known.
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class Material //
		implements IMaterial, IToHTML {

	private final Set<Element> mElements;
	private final String mHTMLName;

	@Override
	public int hashCode() {
		return Objects.hash(mElements, mHTMLName);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Material other = (Material) obj;
		return mElements.equals(other.mElements) && mHTMLName.equals(other.mHTMLName);
	}

	/**
	 * 
	 */
	public Material(String htmlName, Collection<Element> elms) {
		mHTMLName = htmlName.startsWith("<html>") ? htmlName : "<html>" + htmlName;
		mElements = Collections.unmodifiableSet(new HashSet<>(elms));
	}
	
	public Material(String htmlName, Element ... elms) {
		this(htmlName, Arrays.asList(elms));
	}
	
	

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * gov.nist.microanalysis.roentgen.physics.composition.IMaterial#getElementSet()
	 */
	@Override
	public Set<Element> getElementSet() {
		return mElements;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * gov.nist.microanalysis.roentgen.physics.composition.IMaterial#getHTMLName()
	 */
	@Override
	public String getHTMLName() {
		return mHTMLName;
	}

	@Override
	public String toHTML(Mode mode) {
		switch (mode) {
		case TERSE:
			return mHTMLName.substring(6);
		case NORMAL:
		case VERBOSE:
		default:
			StringBuffer sb = new StringBuffer();
			for (Element elm : getElementSet()) {
				if (sb.length() > 0)
					sb.append(",");
				sb.append(elm.getAbbrev());
			}
			return mHTMLName.substring(6) + "[" + sb.toString() + "]";
		}
	}

	@Override
	public Material asMaterial() {
		return this;
	}

}
