package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
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
		extends BaseLabel<Element, Material, Object> {

	private MaterialLabel(final String prefix, final Material material, final Element elm) {
		super(prefix, elm, material);
	}

	public Element getElement() {
		return getObject1();
	}

	public Material getMaterial() {
		return getObject2();
	}

	public static class MassFraction //
			extends MaterialLabel {
		private MassFraction(final Material mat, final Element elm) {
			super("C", mat, elm);
		}
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

	public static class AtomType extends MaterialLabel {
		private AtomType(final String name, final Material mat, final Element elm) {
			super(name, mat, elm);
		}
	}

	public static class AtomicWeight extends MaterialLabel {
		private AtomicWeight(final Material mat, final Element elm) {
			super("A<sub>Z</sub>", mat, elm);
		}
	}

	public static List<AtomicWeight> buildAtomicWeightTags(final Material mat) {
		final List<AtomicWeight> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomicWeight(mat, elm));
		return res;
	}

	public static class AtomFraction extends AtomType {
		private AtomFraction(final Material mat, final Element elm) {
			super("f<sub>atom</sub>", mat, elm);
		}
	}

	public static class Stoichiometry extends AtomType {
		Stoichiometry(final Material mat, final Element elm) {
			super("N<sub>atom</sub>", mat, elm);
		}
	}

	public static List<Stoichiometry> buildStoichiometryTags(final Material mat) {
		final List<Stoichiometry> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new Stoichiometry(mat, elm));
		return res;
	}

	public static List<AtomFraction> buildAtomFractionTags(final Material mat) {
		final List<AtomFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			res.add(new AtomFraction(mat, elm));
		return res;
	}

	public static MaterialMassFraction buildMaterialFractionTag(final Material mat) {
		return new MaterialMassFraction(mat);
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

	public static Stoichiometry buildStoichiometryTag(final Material mat, final Element elm) {
		return new Stoichiometry(mat, elm);
	}

	public static AtomFraction buildAtomFractionTag(final Material mat, final Element elm) {
		return new AtomFraction(mat, elm);
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

	public static class MeanZTag extends BaseLabel<Material, Object, Object> {

		public MeanZTag(final Material mat) {
			super("MeanZ", mat);
		}

		public Material getMaterial() {
			return getObject1();
		}
	}

	public static Object buildMeanAtomicNumberTag(final Material mat) {
		return new MeanZTag(mat);
	}

	public static class MeanATag extends BaseLabel<Material, Object, Object> {

		public MeanATag(final Material mat) {
			super("MeanA", mat);
		}

		public Material getMaterial() {
			return getObject1();
		}
	}

	public static Object buildMeanAtomicWeighTag(final Material mat) {
		return new MeanATag(mat);
	}

	public static class AnalyticalTotalTag extends BaseLabel<Material, Object, Object> {

		private AnalyticalTotalTag(final Material mat) {
			super("Total", mat);
		}

		public Material getMaterial() {
			return getObject1();
		}
	}

	public static Object buildAnalyticalTotalTag(final Material mat) {
		return new AnalyticalTotalTag(mat);
	}

	public static final class MaterialMassFraction //
			extends BaseLabel<Material, Object, Object> {

		public MaterialMassFraction(final Material mat) {
			super("f", mat);
		}

		public Material getMaterial() {
			return getObject1();
		}
	}

}
