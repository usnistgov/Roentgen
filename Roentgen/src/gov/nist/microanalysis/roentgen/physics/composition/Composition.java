package gov.nist.microanalysis.roentgen.physics.composition;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.Normalize;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * The Composition class pulls together all the common representations of
 * compositional data with associated uncertainties into an UncertainValues
 * object. Data is input in the available form (mass fraction, atom fraction,
 * stoichiometry or mixture) and the other representations are computed as
 * possible. The object holds all the available representations (typically mass
 * fraction, atom fraction plus any input representation.) The original input
 * representation is maintained.
 * </p>
 * <p>
 * Compositions are labeled by a String (HTML) name. Compositions with the same
 * name are assumed to be interchangable within calculations.
 * </p>
 *
 * @author Nicholas W.M. Ritchie
 *
 */
public class Composition //
		extends UncertainValues {

	public enum Representation {
		/**
		 * Al2O3 etc
		 */
		Stoichiometry,
		/**
		 * 0.5 atom fraction Qu and 0.5 atom fraction Wz (normalized)
		 */
		AtomFraction,
		/**
		 * 0.29 Qu by mass, 0.68 Wz by mass (not normalized)
		 */
		MassFraction,
		/**
		 * 0.30 Qu by mass, 0.70 Wz by mass (normalized)
		 */
		NormalizedMassFraction,
		/**
		 * 0.58 Al2O3, 0.41 K411 by weight(not normalized)
		 */
		Mixture
	};

	private final Material mMaterial;
	private final Representation mPrimary;
	private Optional<Number> mDensity;
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
		final Composition other = (Composition) obj;
		return super.equals(other) && //
				Objects.equals(mMaterial, other.mMaterial) //
				&& Objects.equals(mPrimary, other.mPrimary) //
				&& mDensity.equals(other.mDensity);
	}

	static public BasicNumberFormat MASS_FRACTION_FORMAT = new BasicNumberFormat("0.0000");
	static public BasicNumberFormat MASS_FRACTION_FORMAT_LONG = new BasicNumberFormat("0.00000");
	static public BasicNumberFormat MEAN_Z_FORMAT = new BasicNumberFormat("0.00");
	static public BasicNumberFormat ATOMIC_FRACTION_FORMAT = new BasicNumberFormat("0.0###");
	static public BasicNumberFormat ATOMIC_FRACTION_FORMAT_LONG = new BasicNumberFormat("0.00E0");

	private static String htmlHelper(final String str) {
		final StringBuffer res = new StringBuffer();
		final int len = str.length();
		boolean inElm = false, inNum = false;
		for (int i = 0; i < len; ++i) {
			final char c = str.charAt(i);
			if (Character.isDigit(c)) {
				inNum = true;
				if (inElm) {
					res.append("<sub>");
					inElm = false;
				}
				res.append(c);
			} else {
				inElm = true;
				if (inNum) {
					res.append("</sub>");
					inNum = false;
				}
				res.append(c);
			}
		}
		if (inNum)
			res.append("</sub>");
		return res.toString();
	}

	private static HashMap<Element, Integer> parseHelper(final String formula) throws ParseException {
		int cx = 0;
		final HashMap<Element, Integer> res = new HashMap<>();
		int next = 0;
		for (int i = 0; i < formula.length(); i = next) {
			next = i + 1;
			int mult = 1;
			final char c0 = formula.charAt(i);
			if (Character.isWhitespace(c0))
				break;
			if (c0 == '(') {
				++cx;
				for (int j = next; j < formula.length(); ++j) {
					if (formula.charAt(j) == '(')
						++cx;
					else if (formula.charAt(j) == ')')
						--cx;
					if (cx == 0) {
						next = j + 1;
						for (int k = next; k < formula.length(); ++k)
							if (Character.isDigit(formula.charAt(k))) {
								mult = Integer.parseInt(formula.substring(next, k + 1));
								next = k + 1;
								break;
							}
						// Parse inside of (...)
						final HashMap<Element, Integer> tmp = parseHelper(formula.substring(i + 1, j).trim());
						for (final Map.Entry<Element, Integer> me : tmp.entrySet()) {
							final int prev = res.containsKey(me.getKey()) ? res.get(me.getKey()) : 0;
							res.put(me.getKey(), Integer.valueOf(prev + (me.getValue() * mult)));
						}
						break;
					}
				}
			} else // Look for element abbreviaton (1 or 2 characters)
			if ((c0 >= 'A') && (c0 <= 'Z')) {
				if (next < formula.length()) {
					// Peak ahead
					final char c1 = formula.charAt(next);
					if ((c1 >= 'a') && (c1 <= 'z'))
						next++;
				}
				final Element elm = Element.parse(formula.substring(i, next));
				if (elm == null)
					throw new ParseException(formula, i);
				// Look for multiplier
				int mult2 = 0;
				if ((next < formula.length()) && Character.isDigit(formula.charAt(next))) {
					for (int k = next; k < formula.length(); ++k)
						if (Character.isDigit(formula.charAt(k))) {
							mult2 = (10 * mult2) + Integer.parseInt(formula.substring(next, k + 1));
							next = k + 1;
						} else
							break;
				} else
					mult2 = 1;
				if (mult2 > 99)
					throw new ParseException(formula, next);
				final int prev = res.containsKey(elm) ? res.get(elm) : 0;
				res.put(elm, Integer.valueOf(prev + mult2));
			} else
				throw new ParseException(formula, i);
		}
		return res;
	}

	/**
	 * @param rep
	 * @param mat
	 * @param uvs
	 * @param density In g/cm<sup>3</sup>
	 */
	private Composition(final Representation rep, final Material mat, final UncertainValues uvs, final Number density) {
		super(uvs.getLabels(), uvs.getValues(), uvs.getCovariances());
		mMaterial = mat;
		mPrimary = rep;
		mDensity = Optional.ofNullable(density);
		// Don't include density in hash
		mHashCode = super.hashCode() ^ Objects.hash(mMaterial, mPrimary);
	}

	private Composition(final Representation rep, final Material mat, final UncertainValues uvs) {
		this(rep, mat, uvs, null);
	}

	private <H> List<H> extractLabels(final Class<H> cls) {
		final List<H> res = new ArrayList<>();
		for (final Object lbl : getLabels())
			if (cls.isInstance(lbl))
				res.add(cls.cast(lbl));
		return res;
	}

	/**
	 * Returns a list of {@link MaterialMassFraction} objects one for each element
	 * in the material.
	 *
	 * @return List&lt;MaterialMassFractionTag&gt;
	 */
	public List<MaterialMassFraction> materialFractionTags() {
		return extractLabels(MaterialMassFraction.class);
	}

	/**
	 * Returns a list of {@link MassFraction} objects one for each element in the
	 * material.
	 *
	 * @return List&lt;MassFractionTag&gt;
	 */
	public List<MaterialLabel.MassFraction> massFractionTags() {
		return MaterialLabel.massFractionTags(this.getMaterial());
	}

	public UncertainValues toMassFraction() {
		return UncertainValues.extract(massFractionTags(), this);
	}

	public UncertainValues getAtomicWeights() {
		return UncertainValues.extract(MaterialLabel.atomWeightTags(this.getMaterial()), this);
	}

	public Optional<Number> getDensity() {
		return mDensity;
	}

	/**
	 * @param n In g/cm<sup>3</sup>
	 */
	public void setDensity(final Number n) {
		mDensity = Optional.of(n);
	}

	public void clearDensity() {
		mDensity = Optional.empty();
	}

	public static class StoichiometryToComposition //
			extends SerialLabeledMultivariateJacobianFunction {
		public StoichiometryToComposition(final Material mat) throws ArgumentException {
			super("StoC", buildSteps(mat));
		}

		private static List<LabeledMultivariateJacobianFunction> buildSteps(final Material mat) {
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			funcs.add(new StoichiometryToAtomFraction(mat));
			funcs.add(new AtomFractionToMassFraction(mat));
			funcs.add(new Normalize(MaterialLabel.buildMassFractionTags(mat)));
			return funcs;
		}

	}

	private static Map<MaterialLabel.AtomicWeight, Number> buildAtomicWeights( //
			final Material mat, final Map<Element, Number> atomicWeights) {
		final Map<MaterialLabel.AtomicWeight, Number> man = new HashMap<>();
		for (final Element elm : mat.getElementSet()) {
			Number A = Double.valueOf(elm.getAtomicWeight());
			if (atomicWeights.containsKey(elm))
				A = atomicWeights.get(elm);
			man.put(MaterialLabel.buildAtomicWeightTag(mat, elm), A);
		}
		return man;
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic, //
			final Map<Element, Number> atomicWeights, //
			final Number density //
	) throws ArgumentException {
		final Material mat = new Material(html, stoic.keySet());
		final Map<MaterialLabel.Stoichiometry, Number> msn = new HashMap<>();
		for (final Map.Entry<Element, Integer> me : stoic.entrySet())
			msn.put(new MaterialLabel.Stoichiometry(mat, me.getKey()), me.getValue());
		final UncertainValues sfrac = new UncertainValues(msn);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(sfrac, wgts);
		final StoichiometryToComposition slmjf = new StoichiometryToComposition(mat);
		final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
		final UncertainValues res = UncertainValues.combine(input, uvs);
		return new Composition(Representation.Stoichiometry, mat, res, density);
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic, //
			final Number density //
	) throws ArgumentException {
		return stoichiometry(html, stoic, Collections.emptyMap(), density);
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic) throws ArgumentException {
		return stoichiometry(html, stoic, Collections.emptyMap(), null);
	}

	public static Composition massFraction(final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		assert elms.size() == vals.getDimension();
		assert elms.size() == vars.getDimension();
		final Material mat = new Material(html, elms);
		final Map<MassFraction, Number> meu = new HashMap<>();
		for (int i = 0; i < elms.size(); ++i)
			meu.put(MaterialLabel.buildMassFractionTag(mat, elms.get(i)),
					new UncertainValue(vals.getEntry(i), Math.sqrt(vars.getEntry(i))));
		final UncertainValues mfracs = new UncertainValues(meu);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(mfracs, wgts);
		try {
			final MassFractionToComposition slmjf = new MassFractionToComposition(mat);
			final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
			final UncertainValues res = UncertainValues.combine(input, uvs);
			return new Composition(Representation.MassFraction, mat, res);
		} catch (final ArgumentException e) {
			e.printStackTrace();
		}
		return null;
	}

	public static Composition massFraction(final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars) throws ArgumentException {
		return massFraction(html, elms, vals, vars, Collections.emptyMap());
	}

	public static class MassFractionToComposition //
			extends SerialLabeledMultivariateJacobianFunction {

		final static List<LabeledMultivariateJacobianFunction> buildSteps(final Material mat) {
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			funcs.add(new MassFractionToAtomFraction(mat));
			funcs.add(new Normalize(MaterialLabel.buildMassFractionTags(mat)));
			return funcs;
		}

		public MassFractionToComposition(final Material mat) throws ArgumentException {
			super("MFtoC", buildSteps(mat));
		}

	}

	public static Composition massFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Map<Element, Number> atomicWeights, //
			final Number density //
	) throws ArgumentException {
		assert elms.size() == vals.getDimension();
		assert elms.size() == cov.getRowDimension();
		assert elms.size() == cov.getColumnDimension();
		final Material mat = new Material(html, elms);
		final List<MaterialLabel.MassFraction> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildMassFractionTag(mat, elm));
		final UncertainValues mfracs = new UncertainValues(labels, vals, cov);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(mfracs, wgts);
		final MassFractionToComposition slmjf = new MassFractionToComposition(mat);
		final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
		final UncertainValues res = UncertainValues.combine(input, uvs);
		return new Composition(Representation.MassFraction, mat, res, density);
	}

	public static Composition massFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Number density //
	) throws ArgumentException {
		return massFraction(html, elms, vals, cov, Collections.emptyMap(), density);
	}

	public static Composition massFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov) throws ArgumentException {
		return massFraction(html, elms, vals, cov, Collections.emptyMap(), null);
	}

	/**
	 * Parses basic strings in common chemical formula notation. It handles simple
	 * formula like H2SO4 and more complex ones like Ca5(PO4)3F. Capitalization is
	 * critical to differentiate PO4 from Po4 and other niceties.
	 *
	 * @param str
	 * @return
	 * @throws ParseException
	 * @throws ArgumentException
	 */
	public static Composition parse( //
			final String str, //
			final Map<Element, Number> atomicWeights //
	) throws ParseException, ArgumentException {
		return Composition.stoichiometry(htmlHelper(str), parseHelper(str), atomicWeights, null);
	}

	/**
	 * Parses basic strings in common chemical formula notation. It handles simple
	 * formula like H2SO4 and more complex ones like Ca5(PO4)3F. Capitalization is
	 * critical to differentiate PO4 from Po4 and other niceties.
	 *
	 * @param str
	 * @return
	 * @throws ParseException
	 * @throws ArgumentException
	 */
	public static Composition parse( //
			final String str) throws ParseException, ArgumentException {
		return Composition.stoichiometry(htmlHelper(str), parseHelper(str), Collections.emptyMap(), null);
	}

	public static Composition atomFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		assert elms.size() == vals.getDimension();
		assert elms.size() == vars.getDimension();
		assert elms.size() == vals.getDimension();
		final Material mat = new Material(html, elms);
		final List<MaterialLabel.AtomFraction> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildAtomFractionTag(mat, elm));
		final UncertainValues sfrac = new UncertainValues(labels, vals, vals);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(sfrac, wgts);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
		final UncertainValues res = UncertainValues.combine(input, uvs);
		return new Composition(Representation.AtomFraction, mat, res);
	}

	public static Composition atomFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars //
	) throws ArgumentException {
		return atomFraction(html, elms, vals, vars, Collections.emptyMap());
	}

	public static Composition atomFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Map<Element, Number> atomicWeights) throws ArgumentException {
		assert elms.size() == vals.getDimension();
		assert elms.size() == cov.getRowDimension();
		assert elms.size() == cov.getColumnDimension();
		final Material mat = new Material(html, elms);
		final List<MaterialLabel.AtomFraction> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildAtomFractionTag(mat, elm));
		final UncertainValues afrac = new UncertainValues(labels, vals, vals);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(afrac, wgts);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
		final UncertainValues res = UncertainValues.combine(input, uvs);
		return new Composition(Representation.AtomFraction, mat, res);
	}

	public static Composition atomFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov //
	) throws ArgumentException {
		return atomFraction(html, elms, vals, cov, Collections.emptyMap());
	}

	/**
	 * Performs the necessary calculations to convert compositional data in atom
	 * fraction into mass fraction and normalized mass fraction.
	 *
	 * @author nicholas
	 *
	 */
	public static class AtomFractionToComposition //
			extends SerialLabeledMultivariateJacobianFunction {

		public AtomFractionToComposition(final Material mat) throws ArgumentException {
			super("AFtoC", buildSteps(mat));
		}

		private static List<LabeledMultivariateJacobianFunction> buildSteps(final Material mat) {
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			funcs.add(new AtomFractionToMassFraction(mat));
			funcs.add(new Normalize(MaterialLabel.buildMassFractionTags(mat)));
			return funcs;
		}

	}

	public static Composition atomFraction(//
			final String html, //
			final Map<Element, ? extends Number> men, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		final Material mat = new Material(html, men.keySet());
		final List<MaterialLabel.AtomFraction> labels = MaterialLabel.buildAtomFractionTags(mat);
		final Map<MaterialLabel.AtomFraction, Number> data = new HashMap<>();
		for (final MaterialLabel.AtomFraction lbl : labels)
			data.put(lbl, men.get(lbl.getElement()));
		final UncertainValues afrac = new UncertainValues(data);
		final UncertainValues wgts = new UncertainValues(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues input = UncertainValues.combine(afrac, wgts);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		final UncertainValues uvs = UncertainValues.propagate(slmjf, input);
		final UncertainValues res = UncertainValues.combine(input, uvs);
		return new Composition(Representation.AtomFraction, mat, res);
	}

	public static Composition atomFraction(//
			final String html, //
			final Map<Element, ? extends Number> men //
	) throws ArgumentException {
		return atomFraction(html, men, Collections.emptyMap());
	}

	/**
	 * A sequence of {@link LabeledMultivariateJacobianFunction}s for converting a
	 * mixture of {@link Composition} objects to mass fractions, atom fractions and
	 * normalized mass fractions.
	 *
	 * @author nicholas
	 *
	 */
	public static class MixtureToComposition //
			extends SerialLabeledMultivariateJacobianFunction {

		private final Material mMaterial;

		public static List<LabeledMultivariateJacobianFunction> buildSteps(//
				final String html, //
				final Set<Material> mats, //
				boolean normalize //
		) {
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			if(normalize) {
				List<MaterialMassFraction> mmfts = new ArrayList<>();
				for(Material mat : mats)
					mmfts.add(MaterialLabel.buildMaterialFractionTag(mat));
				funcs.add(new Normalize(mmfts));
			}
			final MixtureToMassFractions m2mf = new MixtureToMassFractions(html, mats, normalize);
			funcs.add(m2mf);
			final Material newMat = m2mf.getNewMaterial();
			funcs.add(new MassFractionToAtomFraction(newMat));
			funcs.add(new Normalize(MaterialLabel.buildMassFractionTags(newMat)));
			return funcs;
		}

		public MixtureToComposition(final String htmlName, final Set<Material> mats, boolean normalize)
				throws ArgumentException {
			super(htmlName, buildSteps(htmlName, mats, normalize));
			mMaterial = new Material(htmlName, extractElements(mats));
		}

		public Material getMaterial() {
			return mMaterial;
		}
	}

	@SafeVarargs
	private static Map<Composition, Number> convert(final Pair<Composition, Number>... comps) {
		final Map<Composition, Number> res = new HashMap<>();
		for (final Pair<Composition, Number> pr : comps)
			res.put(pr.getFirst(), pr.getSecond());
		return res;
	}

	private static Set<Element> extractElements(final Set<Material> mats) {
		final Set<Element> elms = new HashSet<>();
		for (final Material mat : mats)
			elms.addAll(mat.getElementSet());
		return elms;
	}

	/**
	 * Build a new MassFraction representation of a mixture of the specified
	 * compositions.
	 *
	 * @param htmlName A name for the new material
	 * @param comps    A Pair&lt;Composition, Number&gt; with the composition and
	 *                 fractional quantity.
	 * @return Composition in MassFraction Representation.
	 * @throws ArgumentException
	 */
	@SafeVarargs
	public static Composition combine( //
			final String htmlName, //
			final boolean normalize, //
			final Pair<Composition, Number>... comps) throws ArgumentException {
		return combine(htmlName, convert(comps), normalize);
	}

	/**
	 * Combines fractional amounts of the specified Composition into a Composition
	 * object.
	 * 
	 * There are two uses for this function. The first involves a situation in
	 *
	 * @param htmlName
	 * @param comps
	 * @param normalize
	 * @return {@link Composition}
	 * @throws ArgumentException
	 */
	public static Composition combine( //
			final String htmlName, //
			final Map<Composition, Number> comps, //
			final boolean normalize) throws ArgumentException {
		final List<UncertainValues> input = new ArrayList<>();
		final Set<Material> mats = new HashSet<>();
		final Map<MaterialMassFraction, Number> fracs = new HashMap<>();
		for (final Map.Entry<Composition, Number> pr : comps.entrySet()) {
			input.add(pr.getKey().toMassFraction());
			input.add(pr.getKey().getAtomicWeights());
			fracs.put(new MaterialMassFraction(pr.getKey().getMaterial()), pr.getValue());
			mats.add(pr.getKey().getMaterial());
		}
		input.add(new UncertainValues(fracs));
		final MixtureToComposition combine = new MixtureToComposition(htmlName, mats, normalize);
		final UncertainValues inp = UncertainValues.combine(input);
		final UncertainValues uvs = UncertainValues.propagate(combine, inp);
		final UncertainValues res = UncertainValues.combine(inp, uvs);
		return new Composition(Representation.Mixture, combine.getMaterial(), res);
	}

	/**
	 * @param elm    Element
	 * @param purity [0.0,1.0]
	 * @return Composition
	 * @throws ArgumentException
	 */
	public static Composition pureElement(final Element elm, final double purity) //
			throws ArgumentException {
		final Map<Element, Number> men = new TreeMap<>();
		men.put(elm, new UncertainValue(purity, "d" + elm.toString(), 1.0 - purity));
		return massFraction("Pure " + elm.getAbbrev(),
				Collections.singletonMap(elm, new UncertainValue(purity, 1.0 - purity)), Collections.emptyMap());
	}

	/**
	 * @param elm Element
	 * @return An Composition in MassFraction representation.
	 * @throws ArgumentException
	 */
	public static Composition pureElement(final Element elm) {
		final Map<Object, Number> tags = new HashMap<>();
		final String html = "Pure " + elm.getAbbrev();
		final Material mat = new Material(html, Collections.singleton(elm));
		tags.put(MaterialLabel.buildStoichiometryTag(mat, elm), 1.0);
		MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
		tags.put(mft, 1.0);
		tags.put(Normalize.buildNormalized(mft), 1.0);
		tags.put(MaterialLabel.buildAtomicWeightTag(mat, elm), elm.getAtomicWeight());
		final UncertainValues uvs = new UncertainValues(tags);
		return new Composition(Representation.Stoichiometry, mat, uvs);
	}

	public static Composition massFraction( //
			final String html, //
			final Map<Element, ? extends Number> men, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		final List<Element> elms = new ArrayList<>(men.keySet());
		final RealVector vec = new ArrayRealVector(elms.size());
		final RealVector vars = new ArrayRealVector(elms.size());
		for (int i = 0; i < elms.size(); ++i) {
			final UncertainValue n = new UncertainValue(men.get(elms.get(i)));
			vec.setEntry(i, n.doubleValue());
			vars.setEntry(i, n.variance());
		}
		return Composition.massFraction(html, elms, vec, vars, atomicWeights);
	}

	public static Composition massFraction( //
			final String html, //
			final Map<Element, ? extends Number> men //
	) throws ArgumentException {
		return massFraction(html, men, Collections.emptyMap());
	}

	public static Composition massFraction( //
			final Material mat, //
			final Map<Element, ? extends Number> men //
	) throws ArgumentException {
		return massFraction(mat.getHTMLName(), men, Collections.emptyMap());
	}

	public static Composition massFraction( //
			final Material mat, //
			final UncertainValues uvs //
	) throws ArgumentException {
		final List<Element> elms = new ArrayList<>(mat.getElementSet());
		final List<MassFraction> mfl = MaterialLabel.buildMassFractionTags(mat);
		final RealVector vals = uvs.getValues(mfl);
		final RealMatrix covs = uvs.getCovariances(mfl);
		return massFraction(mat.getHTMLName(), elms, vals, covs);
	}

	public List<Element> getElementList() {
		return new ArrayList<>(mMaterial.getElementSet());
	}

	public Set<Element> getElementSet() {
		return mMaterial.getElementSet();
	}

	public boolean contains(final Element elm) {
		return mMaterial.contains(elm);
	}

	public Representation getNativeRepresentation() {
		return mPrimary;
	}

	public UncertainValue getAtomFraction(final Element elm) {
		final Object tag = MaterialLabel.buildAtomFractionTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValue getMassFraction(final Element elm) {
		final Object tag = MaterialLabel.buildMassFractionTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValue getNomalizedMassFraction(final Element elm) {
		final Object tag = Normalize.buildNormalized(MaterialLabel.buildMassFractionTag(mMaterial, elm));
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValues asStoichiometry() {
		return UncertainValues.extract(MaterialLabel.buildStoichiometryTags(mMaterial), this);
	}

	public UncertainValue getStoichiometry(final Element elm) {
		final Object tag = new MaterialLabel.Stoichiometry(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public static class TotalTag extends BaseLabel<Composition, Object, Object> {

		private TotalTag(final Composition comp) {
			super("Total", comp);
		}

		public Composition getComposition() {
			return getObject1();
		}
	}

	public static Object buildTotalTag(final Composition comp) {
		return new TotalTag(comp);
	}

	public UncertainValue getAnalyticalTotal() {
		final List<MaterialLabel.MassFraction> tags = MaterialLabel.massFractionTags(this.getMaterial());
		final TotalTag tag = new TotalTag(this);
		return UncertainValues.propagateOrdered( //
				LabeledMultivariateJacobianFunction.sum(tags, tag), //
				UncertainValues.extract(tags, this)).getValue(tag);
	}

	public static class MeanZTag extends BaseLabel<Composition, Object, Object> {

		public MeanZTag(final Composition comp) {
			super("MeanZ", comp);
		}

		public Composition getComposition() {
			return getObject1();
		}
	}

	public static Object buildMeanZTag(final Composition comp) {
		return new MeanZTag(comp);
	}

	public UncertainValue getMeanAtomicNumber() {
		final List<? extends Object> tags = MaterialLabel.massFractionTags(this.getMaterial());
		final RealVector z = new ArrayRealVector(tags.size());
		for (int i = 0; i < z.getDimension(); ++i)
			z.setEntry(i, ((MaterialLabel) tags.get(i)).getElement().getAtomicNumber());
		final Object tag = buildMeanZTag(this);
		return UncertainValues.propagateOrdered( //
				LabeledMultivariateJacobianFunction.linear(tags, z, tag), //
				UncertainValues.extract(tags, this)).getValue(tag);
	}

	/**
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.html.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		switch (mPrimary) {
		case MassFraction:
		case NormalizedMassFraction:
			return toHTMLasMassFraction(mode);
		case AtomFraction:
			return toHTMLasAtomFraction(mode);
		case Stoichiometry:
			return toHTMLasStoichiometry(mode);
		case Mixture:
			return toHTMLasMixture(mode);
		default:
			assert false;
			return toHTMLasMassFraction(mode);
		}
	}

	public String toHTML(final Representation rep, final Mode mode) {
		switch (rep) {
		case AtomFraction:
			return toHTMLasAtomFraction(mode);
		case Mixture:
			return toHTMLasMixture(mode);
		case MassFraction:
		case NormalizedMassFraction:
			return toHTMLasMassFraction(mode);
		case Stoichiometry:
			return toHTMLasStoichiometry(mode);
		default:
			return toHTML(mode);
		}
	}

	public Material getMaterial() {
		return mMaterial;
	}

	public List<? extends Object> getLabels(final Representation rep) {
		switch (rep) {
		case AtomFraction:
			return MaterialLabel.buildAtomFractionTags(mMaterial);
		case MassFraction:
			return MaterialLabel.buildMassFractionTags(mMaterial);
		case Mixture: {
			final List<Object> res = new ArrayList<>();
			for (final Object lbl : getLabels())
				if (lbl instanceof MaterialMassFraction)
					res.add(lbl);
			return res;
		}
		case NormalizedMassFraction:
			return Normalize.buildNormalized(MaterialLabel.buildMassFractionTags(mMaterial));
		case Stoichiometry:
			return MaterialLabel.buildStoichiometryTags(mMaterial);
		default:
			assert false;
			return null;
		}

	}

	public boolean hasRepresentation(final Representation rep) {
		final List<? extends Object> lbls = getLabels(rep);
		for (final Object lbl : lbls)
			if (indexOf(lbl) == -1)
				return false;
		return true;
	}

	private String toHTMLasMixture(final Mode mode) {
		assert hasRepresentation(Representation.Mixture);
		switch (mode) {
		default:
		case NORMAL:
			final List<MaterialMassFraction> tags = materialFractionTags();
			final Table t = new Table();
			t.addRow(Table.th("Material"), Table.th("Mass Fraction"), Table.th("Uncertainty"));
			for (final MaterialMassFraction mmft : tags) {
				final UncertainValue uv = getValue(mmft);
				t.addRow( //
						Table.td(mmft.getMaterial()), //
						Table.td(MASS_FRACTION_FORMAT.format(uv.doubleValue())), //
						Table.td(MASS_FRACTION_FORMAT.format(uv.uncertainty())) //
				);
			}
			return HTML.header(mMaterial.getHTMLName()) + t.toHTML(Mode.NORMAL);
		case TERSE:
			return mMaterial.getHTMLName();
		case VERBOSE:
			final StringBuffer sb = new StringBuffer();
			sb.append(toHTMLasMixture(Mode.NORMAL));
			sb.append(toHTMLasMassFraction(Mode.NORMAL));
			return sb.toString();
		}
	}

	private String toHTMLasStoichiometry(final Mode mode) {
		assert hasRepresentation(Representation.Stoichiometry);
		switch (mode) {
		case TERSE:
		case NORMAL:
		default:
			return mMaterial.getHTMLName();
		case VERBOSE: {
			final Table t = new Table();
			t.addRow(Table.th("Element"), Table.thc("Z"), Table.thc("Atoms"), Table.thc("Mass<br/>Fraction"));
			final BasicNumberFormat nf = MASS_FRACTION_FORMAT;
			for (final Object tag : MaterialLabel.stoichiometryTags(this.getMaterial())) {
				final MaterialLabel elementTag = (MaterialLabel) tag;
				final Element elm = elementTag.getElement();
				final UncertainValue uv = getValue(elementTag);
				t.addRow(Table.td(elm.getAbbrev()), //
						Table.tdc(Integer.toString(elm.getAtomicNumber())), //
						Table.tdc(uv.doubleValue()), //
						Table.tdc(nf.format(getMassFraction(elm))));
			}
			if (mDensity.isPresent()) {
				final List<Item> drow = new ArrayList<>();
				drow.add(Table.td("Density")); // Element
				final UncertainValue uv = new UncertainValue(mDensity.get());
				drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", 3));
				t.addRow(drow);
			}
			return HTML.subHeader(mMaterial.getHTMLName()) + t.toHTML(Mode.VERBOSE);
		}
		}
	}

	private String toHTMLasAtomFraction(final Mode mode) {
		assert hasRepresentation(Representation.AtomFraction);
		assert hasRepresentation(Representation.MassFraction);
		switch (mode) {
		case TERSE:
			return mMaterial.getHTMLName();
		case NORMAL: {
			final Table t = new Table();
			final List<Item> hdr = new ArrayList<>();
			hdr.add(Table.th("Element"));
			hdr.add(Table.thc("Z"));
			hdr.add(Table.thc("Mass<br/>Fraction"));
			hdr.add(Table.thc("Uncertainty<br/>(1 &sigma;)"));
			hdr.add(Table.thc("Atomic<br/>Fraction"));
			hdr.add(Table.thc("Uncertainty<br/>(1 &sigma;)"));
			t.addRow(hdr);
			for (final Element elm : getElementSet()) {
				final List<Item> row = new ArrayList<>();
				row.add(Table.td(elm.getAbbrev()));
				row.add(Table.tdc(Integer.toString(elm.getAtomicNumber())));
				final UncertainValue mfVal = getMassFraction(elm);
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(mfVal.doubleValue())));
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(mfVal.uncertainty())));
				final UncertainValue afVal = getAtomFraction(elm);
				row.add(Table.tdc(ATOMIC_FRACTION_FORMAT.formatHTML(afVal.doubleValue())));
				row.add(Table.tdc(ATOMIC_FRACTION_FORMAT.formatHTML(afVal.uncertainty())));
				t.addRow(row);
			}
			if (mDensity.isPresent()) {
				final List<Item> drow = new ArrayList<>();
				drow.add(Table.td("Density")); // Element
				final UncertainValue uv = new UncertainValue(mDensity.get());
				drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", 5));
				t.addRow(drow);
			}
			return t.toHTML(Mode.NORMAL);
		}
		default:
			return toHTML(Mode.VERBOSE, ATOMIC_FRACTION_FORMAT_LONG);
		}
	}

	private String toHTMLasMassFraction(final Mode mode) {
		assert hasRepresentation(Representation.MassFraction);
		switch (mode) {
		case TERSE:
			return mMaterial.getHTMLName();
		case NORMAL: {
			final Table t = new Table();
			try {
				t.addRow(Table.th(mMaterial.getHTMLName(), 8));
				t.addRow(Table.th("Element"), //
						Table.thc("Z"), //
						Table.thc("Mass<br/>Fraction"), //
						Table.thc("Uncertainty<br/>(1 &sigma;)"), //
						Table.thc("Normalized<br/>Mass Fraction"), //
						Table.thc("Uncertainty<br/>(1 &sigma;)"), //
						Table.thc("Atomic<br/>Fraction"), //
						Table.thc("Uncertainty<br/>(1 &sigma;)"));//
				for (final Element elm : getElementSet()) {
					final List<Item> row = new ArrayList<>();
					row.add(Table.td(elm.getAbbrev()));
					row.add(Table.tdc(Integer.toString(elm.getAtomicNumber())));
					final UncertainValue mfVal = getMassFraction(elm);
					row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(mfVal.doubleValue())));
					row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(mfVal.uncertainty())));
					final UncertainValue nmfVal = getNomalizedMassFraction(elm);
					row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(nmfVal.doubleValue())));
					row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(nmfVal.uncertainty())));
					final UncertainValue afVal = getAtomFraction(elm);
					row.add(Table.tdc(ATOMIC_FRACTION_FORMAT.formatHTML(afVal.doubleValue())));
					row.add(Table.tdc(ATOMIC_FRACTION_FORMAT.formatHTML(afVal.uncertainty())));
					t.addRow(row);
				}
				final List<Item> row = new ArrayList<>();
				row.add(Table.td("Total")); // Element
				final UncertainValue meanZ = getMeanAtomicNumber();
				row.add(Table.tdc(meanZ.toHTML(Mode.TERSE))); // Z
				final UncertainValue sum = getAnalyticalTotal();
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(sum.doubleValue()))); // MassFrac
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(sum.uncertainty()))); // Unc
				final UncertainValue nsum = UncertainValue.ONE;
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(nsum.doubleValue()))); // MassFrac
				row.add(Table.tdc(MASS_FRACTION_FORMAT.formatHTML(nsum.uncertainty()))); // Unc
				row.add(Table.td()); // Atomic
				row.add(Table.td()); // Unc
				t.addRow(row);
				if (mDensity.isPresent()) {
					final List<Item> drow = new ArrayList<>();
					drow.add(Table.td("Density")); // Element
					final UncertainValue uv = new UncertainValue(mDensity.get());
					drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", 7));
					t.addRow(drow);
				}
			} catch (final Exception e) {
				e.printStackTrace();
			}
			return t.toHTML(Mode.NORMAL);
		}
		default:
			return toHTML(Mode.VERBOSE, MASS_FRACTION_FORMAT);
		}
	}

	@Override
	public String toString() {
		return "C[" + mMaterial.toString() + "]";
	}

}
