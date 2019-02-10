package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.IMaterial;

/**
 * Labels to identify Composition-related properties like {@link MassFraction},
 * {@link AtomFraction}, {@link NormalizedMassFraction} and
 * {@link Stoichiometry}.
 * 
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class CompositionalLabel //
		extends BaseLabel<Element, String, Object> {

	private CompositionalLabel(final String prefix, final String html, final Element elm) {
		super(prefix, elm, html);
	}

	public Element getElement() {
		return getObject1();
	}

	public String getHTML() {
		return getObject2();
	}

	public static class MassFraction //
			extends CompositionalLabel {
		private MassFraction(final String html, final Element elm) {
			super("C", html, elm);
		}
	}

	public static MassFraction buildMassFractionTag(final IMaterial comp, final Element elm) {
		return buildMassFractionTag(comp.getHTMLName(), elm);
	}

	public static MassFraction buildMassFractionTag(final String html, final Element elm) {
		return new MassFraction(html, elm);
	}

	public static List<MassFraction> buildMassFractionTags(final String html, final Collection<Element> elms) {
		final List<MassFraction> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(new MassFraction(html, elm));
		return res;
	}

	public static class NormalizedMassFraction extends CompositionalLabel {
		NormalizedMassFraction(final String html, final Element elm) {
			super("N", html, elm);
		}
	}

	public static NormalizedMassFraction buildNormalizedMassFractionTag(String html, Element elm) {
		return new NormalizedMassFraction(html, elm);
	}

	public static List<NormalizedMassFraction> buildNormMassFractionTags(final String html,
			final Collection<Element> elms) {
		final List<NormalizedMassFraction> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(new NormalizedMassFraction(html, elm));
		return res;
	}

	public static class AtomType extends CompositionalLabel {
		private AtomType(final String name, final String html, final Element elm) {
			super(name, html, elm);
		}
	}

	public static class AtomicWeight extends CompositionalLabel {
		private AtomicWeight(final String html, final Element elm) {
			super("A", html, elm);
		}
	}

	public static AtomicWeight buildAtomicWeightTag(final String html, final Element elm) {
		return new AtomicWeight(html, elm);
	}

	public static List<AtomicWeight> buildAtomicWeightTags(final String html, final Collection<Element> elms) {
		final List<AtomicWeight> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(new AtomicWeight(html, elm));
		return res;
	}

	public static class AtomFraction extends AtomType {
		private AtomFraction(final String html, final Element elm) {
			super("f<sub>Atom</sub>", html, elm);
		}
	}

	public static AtomFraction buildAtomFractionTag(final String html, final Element elm) {
		return new AtomFraction(html, elm);
	}

	public static class Stoichiometry extends AtomType {
		Stoichiometry(final String html, final Element elm) {
			super("S", html, elm);
		}
	}

	public static List<Stoichiometry> buildStoichiometryTags(final String html, final Collection<Element> elms) {
		final List<Stoichiometry> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(new Stoichiometry(html, elm));
		return res;
	}

	public static List<AtomFraction> buildAtomFractionTags(final String html, final Collection<Element> elms) {
		final List<AtomFraction> res = new ArrayList<>();
		for (final Element elm : elms)
			res.add(new AtomFraction(html, elm));
		return res;
	}

	public static MaterialMassFraction buildMaterialFractionTag(IMaterial mat) {
		return new MaterialMassFraction(mat);
	}

	/**
	 * Build a tag to associate with the atomic weight of the specified element in
	 * the specified material.
	 *
	 * @param comp
	 * @param elm
	 * @return {@link AtomicWeight}
	 */
	public static AtomicWeight buildAtomicWeightTag(final IMaterial comp, final Element elm) {
		return new AtomicWeight(comp.getHTMLName(), elm);
	}

	public static Stoichiometry buildStoichiometryTag(final IMaterial comp, final Element elm) {
		return new Stoichiometry(comp.getHTMLName(), elm);
	}

	public static Stoichiometry buildStoichiometryTag(final String html, final Element elm) {
		return new Stoichiometry(html, elm);
	}

	public static AtomFraction buildAtomFractionTag(final IMaterial comp, final Element elm) {
		return new AtomFraction(comp.getHTMLName(), elm);
	}

	/**
	 * Returns a list of {@link MassFraction} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;MassFractionTag&gt;
	 */
	public static List<MassFraction> massFractionTags(IMaterial mat) {
		final List<MassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new MassFraction(mat.getHTMLName(), elm));
		return res;
	}

	/**
	 * Returns a list of {@link AtomWeightTag} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;AtomWeightTag&gt;
	 */
	public static List<AtomicWeight> atomWeightTags(IMaterial mat) {
		final List<AtomicWeight> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomicWeight(mat.getHTMLName(), elm));
		return res;
	}

	/**
	 * Returns a list of {@link Stoichiometry} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;StoichiometryTag&gt;
	 */
	public static List<Stoichiometry> stoichiometryTags(IMaterial mat) {
		final List<Stoichiometry> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new Stoichiometry(mat.getHTMLName(), elm));
		return res;
	}

	/**
	 * Returns a list of {@link AtomFraction} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;AtomFractionTag&gt;
	 */
	public List<AtomFraction> atomFractionTags(IMaterial mat) {
		final List<AtomFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomFraction(mat.getHTMLName(), elm));
		return res;
	}

	/**
	 * Returns a list of {@link NormalizedMassFraction} objects one for each element
	 * in the material.
	 *
	 * @return List&lt;NormalizedMassFractionTag&gt;
	 */
	public List<NormalizedMassFraction> normalizedMassFractionTags(IMaterial mat) {
		final List<NormalizedMassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new NormalizedMassFraction(mat.getHTMLName(), elm));
		return res;
	}

}
