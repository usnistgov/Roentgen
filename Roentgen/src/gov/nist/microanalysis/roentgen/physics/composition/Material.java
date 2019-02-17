package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * The Material class represents the most basic information about a material -
 * its name and the set of constituent elements.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class Material //
		implements IToHTML {

	private final SortedSet<Element> mElements;
	private final String mHTMLName;
	private final int mHashCode;

	@Override
	public int hashCode() {
		return mHashCode;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Material other = (Material) obj;
		return mHTMLName.equals(other.mHTMLName) && mElements.equals(other.mElements);
	}

	public static Set<Material> convert(final Collection<Composition> comps) {
		final Set<Material> res = new HashSet<>();
		for (final Composition comp : comps)
			res.add(comp.getMaterial());
		return res;
	}

	/**
	 *
	 */
	public Material(final String htmlName, final Collection<Element> elms) {
		mHTMLName = htmlName;
		mElements = Collections.unmodifiableSortedSet(new TreeSet<>(elms));
		mHashCode = Objects.hash(mElements, mHTMLName);
	}

	public Material(final String htmlName, final Element... elms) {
		this(htmlName, Arrays.asList(elms));
	}

	public SortedSet<Element> getElementSet() {
		return mElements;
	}

	public String getHTMLName() {
		return mHTMLName;
	}

	@Override
	public String toString() {
		return HTML.stripTags(mHTMLName);
	}

	public boolean contains(final Element elm) {
		return mElements.contains(elm);
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE:
			return mHTMLName;
		case NORMAL:
		case VERBOSE:
		default:
			final StringBuffer sb = new StringBuffer();
			for (final Element elm : getElementSet()) {
				if (sb.length() > 0)
					sb.append(",");
				sb.append(elm.getAbbrev());
			}
			return mHTMLName + "[" + sb.toString() + "]";
		}
	}
}
