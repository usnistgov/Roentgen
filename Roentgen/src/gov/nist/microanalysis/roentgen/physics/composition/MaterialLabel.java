package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * Labels to identify Composition-related properties like {@link MassFraction},
 * {@link AtomFraction}, {@link NormalizedMassFraction} and
 * {@link Stoichiometry}.
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class MaterialLabel //
		extends EPMALabel.BaseLabel<Material, Element, Object> {

	public static class CompositionalStatisticTag extends MaterialLabel {

		protected CompositionalStatisticTag(final String name, final Material mat) {
			super(name, mat);
		}
	}

	public static class AnalyticalTotalTag extends CompositionalStatisticTag {

		private AnalyticalTotalTag(final Material mat) {
			super("&Sigma;C<sub>z</sub>", mat);
		}
	}

	public static class AtomFraction extends AtomType {
		private AtomFraction(final Material mat, final Element elm) {
			super("f<sub>A</sub>", mat, elm);
		}
	}

	public static class AtomicWeight extends MaterialLabel {
		private AtomicWeight(final Material mat, final Element elm) {
			super("A<sub>Z</sub>", mat, elm);
		}
	}

	public static class AtomType extends MaterialLabel {
		private AtomType(final String name, final Material mat, final Element elm) {
			super(name, mat, elm);
		}
	}

	public static class MassFraction //
			extends MaterialLabel {

		private MassFraction(final Material mat, final Element elm) {
			super("C", mat, elm);
		}
	}

	public static class NormalizedMassFraction //
			extends MaterialLabel {
		private NormalizedMassFraction(final Material mat, final Element elm) {
			super("N", mat, elm);
		}
	}

	public static final class MaterialMassFraction //
			extends MaterialLabel {

		public MaterialMassFraction(final Material mat) {
			this(mat, false);
		}

		public MaterialMassFraction(final Material mat, final boolean normalized) {
			super(normalized ? "f<sub>n</sub>" : "f", mat);
		}
	}

	public static class MeanATag extends CompositionalStatisticTag {

		public MeanATag(final Material mat) {
			super("A&#773;", mat);
		}
	}

	public static class MeanZTag extends CompositionalStatisticTag {

		public MeanZTag(final Material mat) {
			super("Z&#773;", mat);
		}
	}

	public static class Stoichiometry extends AtomType {
		Stoichiometry(final Material mat, final Element elm) {
			super("N<sub>atom</sub>", mat, elm);
		}
	}

	/**
	 * Returns a list of {@link AtomWeightTag} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;AtomWeightTag&gt;
	 */
	public static List<AtomicWeight> atomWeightTags(final Material mat) {
		final List<AtomicWeight> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomicWeight(mat, elm));
		return res;
	}

	public static AnalyticalTotalTag buildAnalyticalTotalTag(final Material mat) {
		return new AnalyticalTotalTag(mat);
	}

	public static AtomFraction buildAtomFractionTag(final Material mat, final Element elm) {
		return new AtomFraction(mat, elm);
	}

	public static List<AtomFraction> buildAtomFractionTags(final Material mat) {
		final List<AtomFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomFraction(mat, elm));
		return res;
	}

	/**
	 * Build a tag to associate with the atomic weight of the specified element in
	 * the specified material.
	 *
	 * @param mat
	 * @param elm
	 * @return {@link AtomicWeight}
	 */
	public static AtomicWeight buildAtomicWeightTag(final Material mat, final Element elm) {
		return new AtomicWeight(mat, elm);
	}

	public static List<AtomicWeight> buildAtomicWeightTags(final Material mat) {
		final List<AtomicWeight> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomicWeight(mat, elm));
		return res;
	}

	public static MassFraction buildMassFractionTag(final Material mat, final Element elm) {
		return new MassFraction(mat, elm);
	}

	public static List<MassFraction> buildMassFractionTags(final Material mat) {
		final List<MassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new MassFraction(mat, elm));
		return res;
	}

	public static MaterialMassFraction buildMaterialFractionTag(final Material mat) {
		return new MaterialMassFraction(mat, false);
	}

	public static MeanZTag buildMeanAtomicNumberTag(final Material mat) {
		return new MeanZTag(mat);
	}

	public static MeanATag buildMeanAtomicWeighTag(final Material mat) {
		return new MeanATag(mat);
	}

	public static NormalizedMassFraction buildNormalizedMassFractionTag(final Material mat, final Element elm) {
		return new NormalizedMassFraction(mat, elm);
	}

	public static List<NormalizedMassFraction> buildNormalizedMassFractionTags(final Material mat) {
		final List<NormalizedMassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new NormalizedMassFraction(mat, elm));
		return res;
	}

	public static MaterialMassFraction buildNormalizedMaterialFractionTag(final Material mat) {
		return new MaterialMassFraction(mat, true);
	}

	public static Stoichiometry buildStoichiometryTag(final Material mat, final Element elm) {
		return new Stoichiometry(mat, elm);
	}

	public static List<Stoichiometry> buildStoichiometryTags(final Material mat) {
		final List<Stoichiometry> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new Stoichiometry(mat, elm));
		return res;
	}

	public static <H extends EPMALabel> List<H> convert(List<? extends H> list) {
		List<H> res = new ArrayList<>();
		res.addAll(list);
		return res;
	}

	/**
	 * Returns a list of {@link MassFraction} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;MassFractionTag&gt;
	 */
	public static List<MassFraction> massFractionTags(final Material mat) {
		final List<MassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new MassFraction(mat, elm));
		return res;
	}

	/**
	 * Returns a list of {@link Stoichiometry} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;StoichiometryTag&gt;
	 */
	public static List<Stoichiometry> stoichiometryTags(final Material mat) {
		final List<Stoichiometry> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new Stoichiometry(mat, elm));
		return res;
	}

	private MaterialLabel(final String prefix, final Material material) {
		super(prefix, material);
	}

	private MaterialLabel(final String prefix, final Material material, final Element elm) {
		super(prefix, material, elm);
	}

	/**
	 * Returns a list of {@link AtomFraction} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;AtomFractionTag&gt;
	 */
	public List<AtomFraction> atomFractionTags(final Material mat) {
		final List<AtomFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomFraction(mat, elm));
		return res;
	}

	public Element getElement() {
		return getObject2();
	}

	public Material getMaterial() {
		return getObject1();
	}

}
