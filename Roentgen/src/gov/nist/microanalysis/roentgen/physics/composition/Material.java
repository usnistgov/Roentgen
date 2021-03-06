package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.UUID;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.juncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;

/**
 * The Material class represents the most basic information about a material -
 * its name and the set of constituent elements.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class Material //
		implements IToHTML {

	private final String mHTMLName;
	private final SortedSet<Element> mElements;
	private final Optional<Conductivity> mConductivity;
	private final UUID mUniquizer;

	private final Map<Element, Number> mAtomicWeights;

	private final int mHashCode;

	@Override
	public int hashCode() {
		return mHashCode;
	}

	@Override
	public boolean equals(
			final Object obj
	) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Material other = (Material) obj;
		// Omitting atomic weights
		return mHTMLName.equals(other.mHTMLName) //
				&& mElements.equals(other.mElements) //
				&& mConductivity.equals(other.mConductivity) //
				&& mUniquizer.equals(other.mUniquizer) //
				&& mAtomicWeights.equals(other.mAtomicWeights);
	}

	public static Set<Material> convert(
			final Collection<Composition> comps
	) {
		final Set<Material> res = new HashSet<>();
		for (final Composition comp : comps)
			res.add(comp.getMaterial());
		return res;
	}

	public Material(
			//
			final String htmlName, //
			final Collection<Element> elms, //
			final Conductivity conduct, //
			final Map<Element, ? extends Number> atomicWeights //
	) {
		mHTMLName = htmlName;
		mElements = Collections.unmodifiableSortedSet(new TreeSet<>(elms));
		mConductivity = Optional.ofNullable(conduct);
		mAtomicWeights = new HashMap<>(atomicWeights);
		mUniquizer = UUID.randomUUID();
		mHashCode = Objects.hash(mHTMLName, mElements, mConductivity, mAtomicWeights, mUniquizer);
	}

	public Material(
			final String htmlName, final Collection<Element> elms, final Conductivity conduct
	) {
		this(htmlName, elms, conduct, Collections.emptyMap());
	}

	public Material(
			final String htmlName, final Collection<Element> elms
	) {
		this(htmlName, elms, null);
	}

	public SortedSet<Element> getElementSet() {
		return mElements;
	}

	public String getHTMLName() {
		return mHTMLName.toString();
	}

	public Number getAtomicWeight(
			final Element elm
	) {
		return mAtomicWeights.getOrDefault(elm, elm.getAtomicWeight());
	}

	public UncertainValues<AtomicWeight> getAtomicWeights() {
		final Map<AtomicWeight, Number> unkMap = new HashMap<>();
		for (final Element elm : getElementSet())
			unkMap.put(MaterialLabel.buildAtomicWeightTag(this, elm), getAtomicWeight(elm));
		return new UncertainValues<AtomicWeight>(unkMap);
	}

	public void setAtomicWeight(
			final Element elm, final Number value
	) {
		mAtomicWeights.put(elm, value);
	}

	@Override
	public String toString() {
		return HTML.stripTags(mHTMLName.toString());
	}

	public boolean contains(
			final Element elm
	) {
		return mElements.contains(elm);
	}

	@Override
	public String toHTML(
			final Mode mode
	) {
		switch (mode) {
		default:
		case TERSE:
			return mHTMLName.toString();
		case NORMAL: {
			final StringBuffer sb = new StringBuffer();
			for (final Element elm : getElementSet()) {
				if (sb.length() > 0)
					sb.append(",");
				sb.append(elm.getAbbrev());
			}
			return mHTMLName + "[" + sb.toString() + "]";
		}
		case VERBOSE: {
			final Table t = new Table();
			t.addRow(Table.th("Item"), Table.th("Value"));
			final StringBuffer sb = new StringBuffer();
			for (final Element elm : getElementSet()) {
				if (sb.length() > 0)
					sb.append(",");
				sb.append(elm.getAbbrev());
			}
			t.addRow(Table.td("Name"), Table.td(mHTMLName));
			t.addRow(Table.td("Elements"), Table.td(sb.toString()));
			if (mConductivity.isPresent())
				t.addRow(Table.td("Conductivity"), Table.td(mConductivity.toString()));
			return t.toHTML(Mode.NORMAL);
		}
		}
	}
}
