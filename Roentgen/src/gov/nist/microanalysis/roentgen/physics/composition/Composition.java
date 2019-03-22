package gov.nist.microanalysis.roentgen.physics.composition;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.CompositeLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.math.uncertainty.models.Normalize;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MaterialMassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.Stoichiometry;
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
		extends UncertainValuesCalculator<MaterialLabel> {

	/**
	 * Performs the necessary calculations to convert compositional data in atom
	 * fraction into mass fraction and normalized mass fraction.
	 *
	 * @author nicholas
	 *
	 */
	public static class AtomFractionToComposition //
			extends CompositeLabeledMultivariateJacobianFunction<MaterialLabel> {

		private static List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> buildSteps(
				final Material mat) {
			final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> funcs = new ArrayList<>();
			funcs.add(new AtomFractionToMassFraction(mat));
			funcs.add(new MassFractionToNormalizedMassFraction(mat));
			funcs.add(new CompositionStatistics(mat));
			return funcs;
		}

		public AtomFractionToComposition(final Material mat) throws ArgumentException {
			super("AFtoC", buildSteps(mat));
		}

	};

	public static class ElementByDifferenceToComposition //
			extends CompositeLabeledMultivariateJacobianFunction<MaterialLabel> {

		final static List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> buildSteps(
				final Material mat, final Element elm) {
			assert mat.getElementSet().contains(elm);
			final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> funcs = new ArrayList<>();
			funcs.add(new ElementByDifference(mat, elm));
			funcs.add(new MassFractionToAtomFraction(mat));
			funcs.add(new MassFractionToNormalizedMassFraction(mat));
			funcs.add(new CompositionStatistics(mat));
			return funcs;
		}

		public ElementByDifferenceToComposition(final Material mat, final Element elm) throws ArgumentException {
			super("EbDtoC", buildSteps(mat, elm));
		}

	}

	public static class MassFractionToComposition //
			extends CompositeLabeledMultivariateJacobianFunction<MaterialLabel> {

		final static List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> buildSteps(
				final Material mat,
				final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElmRule) {
			final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> funcs = new ArrayList<>();
			if (otherElmRule != null)
				funcs.add(otherElmRule);
			funcs.add(new MassFractionToAtomFraction(mat));
			funcs.add(new MassFractionToNormalizedMassFraction(mat));
			funcs.add(new CompositionStatistics(mat));
			return funcs;
		}

		public MassFractionToComposition(final Material mat) throws ArgumentException {
			super("MFtoC", buildSteps(mat, null));
		}

		public MassFractionToComposition(final Material mat,
				final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElmRule)
				throws ArgumentException {
			super("MFtoC", buildSteps(mat, otherElmRule));
		}

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
			extends CompositeLabeledMultivariateJacobianFunction<MaterialLabel> {

		public static List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> buildSteps(//
				final Material newMat, //
				final Set<Material> mats, //
				final boolean normalize //
		) {
			final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> funcs = new ArrayList<>();
			if (normalize) {
				final List<MaterialMassFraction> mmfts = new ArrayList<>();
				final List<MaterialMassFraction> nmfts = new ArrayList<>();
				for (final Material mat : mats) {
					mmfts.add(MaterialLabel.buildMaterialFractionTag(mat));
					nmfts.add(MaterialLabel.buildNormalizedMaterialFractionTag(mat));
				}
				funcs.add(new Normalize<MaterialMassFraction, MaterialMassFraction>(mmfts, nmfts));
			}
			funcs.add(new MixtureToMassFractions(newMat, mats, normalize));
			funcs.add(new MassFractionToAtomFraction(newMat));
			funcs.add(new MassFractionToNormalizedMassFraction(newMat));
			funcs.add(new CompositionStatistics(newMat));
			return funcs;
		}

		public MixtureToComposition(final Material newMat, final Set<Material> mats, final boolean normalize)
				throws ArgumentException {
			super("MtoC", buildSteps(newMat, mats, normalize));
		}
	}

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
	}

	public static class StoichiometryToComposition //
			extends CompositeLabeledMultivariateJacobianFunction<MaterialLabel> {

		private static List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> buildSteps(
				final Material mat) {
			final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> funcs = new ArrayList<>();
			funcs.add(new StoichiometryToAtomFraction(mat));
			funcs.add(new AtomFractionToMassFraction(mat));
			funcs.add(new MassFractionToNormalizedMassFraction(mat));
			funcs.add(new CompositionStatistics(mat));
			return funcs;
		}

		public StoichiometryToComposition(final Material mat) throws ArgumentException {
			super("StoC", buildSteps(mat));
		}

	}

	static public BasicNumberFormat MASS_FRACTION_FORMAT = new BasicNumberFormat("0.0000");

	static public BasicNumberFormat MASS_FRACTION_FORMAT_LONG = new BasicNumberFormat("0.00000");
	static public BasicNumberFormat MEAN_Z_FORMAT = new BasicNumberFormat("0.00");
	static public BasicNumberFormat ATOMIC_FRACTION_FORMAT = new BasicNumberFormat("0.0###");
	static public BasicNumberFormat ATOMIC_FRACTION_FORMAT_LONG = new BasicNumberFormat("0.00E0");

	public static Composition atomFraction( //
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov //
	) throws ArgumentException {
		return atomFraction(html, elms, vals, cov, Collections.emptyMap());
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
		final List<MaterialLabel> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildAtomFractionTag(mat, elm));
		final UncertainValues<MaterialLabel> afrac = new UncertainValues<MaterialLabel>(labels, vals, vals);
		final UncertainValues<MaterialLabel> wgts = new UncertainValues<MaterialLabel>(
				buildAtomicWeights(mat, atomicWeights));
		final UncertainValuesBase<MaterialLabel> input = UncertainValuesBase.combine(Arrays.asList(afrac, wgts), false);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		return new Composition(Representation.AtomFraction, mat, slmjf, input);
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
			final RealVector vars, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		assert elms.size() == vals.getDimension();
		assert elms.size() == vars.getDimension();
		assert elms.size() == vals.getDimension();
		final Material mat = new Material(html, elms);
		final List<MaterialLabel> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildAtomFractionTag(mat, elm));
		final UncertainValues<MaterialLabel> sfrac = new UncertainValues<MaterialLabel>(labels, vals, vals);
		final UncertainValues<MaterialLabel> wgts = new UncertainValues<MaterialLabel>(
				buildAtomicWeights(mat, atomicWeights));
		final UncertainValuesBase<MaterialLabel> input = UncertainValuesBase.combine(Arrays.asList(sfrac, wgts), false);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		return new Composition(Representation.AtomFraction, mat, slmjf, input);
	}

	public static Composition atomFraction(//
			final String html, //
			final Map<Element, ? extends Number> men //
	) throws ArgumentException {
		return atomFraction(html, men, Collections.emptyMap());
	}

	public static Composition atomFraction(//
			final String html, //
			final Map<Element, ? extends Number> men, //
			final Map<Element, Number> atomicWeights //
	) throws ArgumentException {
		final Material mat = new Material(html, men.keySet());
		final List<AtomFraction> labels = MaterialLabel.buildAtomFractionTags(mat);
		final Map<MaterialLabel, Number> data = new HashMap<>();
		for (final MaterialLabel lbl : labels)
			data.put(lbl, men.get(lbl.getElement()));
		final UncertainValues<MaterialLabel> afrac = new UncertainValues<MaterialLabel>(data);
		final UncertainValues<MaterialLabel> wgts = new UncertainValues<MaterialLabel>(
				buildAtomicWeights(mat, atomicWeights));
		final UncertainValuesBase<MaterialLabel> input = UncertainValuesBase.combine(Arrays.asList(afrac, wgts), false);
		final AtomFractionToComposition slmjf = new AtomFractionToComposition(mat);
		return new Composition(Representation.AtomFraction, mat, slmjf, input);
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
		final List<UncertainValuesBase<MaterialLabel>> input = new ArrayList<>();
		final Set<Material> mats = new HashSet<>();
		final Map<MaterialLabel, Number> fracs = new HashMap<>();
		final Set<Element> elms = new HashSet<>();
		final List<LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel>> steps = new ArrayList<>();
		for (final Map.Entry<Composition, Number> pr : comps.entrySet()) {
			final Composition comp = pr.getKey();
			input.add(comp.getInputs());
			steps.add(comp.getFunction());
			fracs.put(new MaterialMassFraction(comp.getMaterial()), pr.getValue());
			mats.add(comp.getMaterial());
			elms.addAll(comp.getElementSet());
		}
		input.add(new UncertainValues<MaterialLabel>(fracs));
		final Material newMat = new Material(htmlName, elms);
		steps.add(new MixtureToComposition(newMat, mats, normalize));
		final CompositeLabeledMultivariateJacobianFunction<MaterialLabel> combine = //
				new CompositeLabeledMultivariateJacobianFunction<MaterialLabel>("Mixture", steps);
		final UncertainValues<MaterialLabel> inp = UncertainValues
				.<MaterialLabel>asUncertainValues(UncertainValuesBase.combine(input, true));
		return new Composition(Representation.Mixture, newMat, combine, inp);
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
	public static Composition combineDelta( //
			final String htmlName, //
			final Map<Composition, Number> comps, //
			final boolean normalize) throws ArgumentException {
		final Composition res = combine(htmlName, comps, normalize);
		final FiniteDifference delta = res.new FiniteDifference(res.getInputValues().mapMultiply(0.001));
		res.setCalculator(delta);
		return res;
	}

	public static Composition elementByDifference(//
			final Material mat, //
			final Element diffElm, //
			final Map<Element, Number> men, //
			final Map<Object, Number> weights) throws ArgumentException {
		final Map<MaterialLabel, Number> tmp = new HashMap<>();
		for (final Element elm : men.keySet())
			if (!elm.equals(diffElm))
				tmp.put(MaterialLabel.buildMassFractionTag(mat, elm), men.get(elm));
		for (final Element elm : mat.getElementSet())
			tmp.put(MaterialLabel.buildAtomicWeightTag(mat, elm), mat.getAtomicWeight(elm));
		final UncertainValues<MaterialLabel> inp = new UncertainValues<MaterialLabel>(tmp);
		final LabeledMultivariateJacobianFunction<MaterialLabel, MaterialLabel> func = new ElementByDifferenceToComposition(
				mat, diffElm);
		return new Composition(Representation.MassFraction, mat, func, inp);
	}

	public static Composition massFraction(final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov) throws ArgumentException {
		return massFraction(mat, elms, vals, cov, Collections.emptyMap(), null, null);
	}

	/**
	 * All massFraction implementations eventually pass through this one.
	 *
	 *
	 * @param mat           The material
	 * @param elms          An ordered list of elements in the material
	 * @param vals          Mass fraction values in the same order as elms
	 * @param cov           Covariance matrix in the same row/col order as elms
	 * @param atomicWeights A Map&lt;Element, Number&gt; with optional atomic weight
	 *                      values. Otherwise book values are used. (May be
	 *                      Collections.emptyMap())
	 * @param density       Optional density in g/cm<sup>3</sup>
	 * @param otherElmRules A list of {@link LabeledMultivariateJacobianFunction}
	 * @return
	 * @throws ArgumentException
	 */
	public static Composition massFraction(//
			final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Map<Element, Number> atomicWeights, //
			final Number density, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElmRules //
	) throws ArgumentException {
		assert mat.getElementSet().containsAll(elms);
		assert allElements(elms, otherElmRules).containsAll(mat.getElementSet());
		assert elms.size() == vals.getDimension();
		assert elms.size() == cov.getRowDimension();
		assert elms.size() == cov.getColumnDimension();
		final List<MaterialLabel.MassFraction> labels = new ArrayList<>();
		for (final Element elm : elms)
			labels.add(MaterialLabel.buildMassFractionTag(mat, elm));
		final UncertainValues<MaterialLabel> mfracs = new UncertainValues<MaterialLabel>(labels, vals, cov);
		final UncertainValues<MaterialLabel> wgts = new UncertainValues<MaterialLabel>(
				buildAtomicWeights(mat, atomicWeights));
		final UncertainValues<MaterialLabel> input = UncertainValues.<MaterialLabel>asUncertainValues(//
				UncertainValuesBase.combine(Arrays.asList(mfracs, wgts), false));
		final MassFractionToComposition slmjf = new MassFractionToComposition(mat, otherElmRules);
		return new Composition(Representation.MassFraction, mat, slmjf, input, density);
	}

	public static Composition massFraction(final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Number density //
	) throws ArgumentException {
		return massFraction(mat, elms, vals, cov, Collections.emptyMap(), density, null);
	}

	public static Composition massFraction(//
			final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars //
	) throws ArgumentException {
		return massFraction(mat, elms, vals, vars, Collections.emptyMap(), null);
	}

	public static Composition massFraction(//
			final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars, //
			final Map<Element, Number> atomicWeights) throws ArgumentException {
		return massFraction(mat, elms, vals, vars, atomicWeights, null);
	}

	public static Composition massFraction(//
			final Material mat, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars, //
			final Map<Element, Number> atomicWeights, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElmRules //
	) throws ArgumentException {
		return massFraction(mat, elms, vals, MatrixUtils.createRealDiagonalMatrix(vars.toArray()), atomicWeights, null,
				otherElmRules);
	}

	public static Composition massFraction( //
			final Material mat, //
			final Map<Element, ? extends Number> men, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElmRule //
	) throws ArgumentException {
		final List<Element> elms = new ArrayList<>(men.keySet());
		final RealVector vec = new ArrayRealVector(elms.size());
		final RealMatrix cov = MatrixUtils.createRealMatrix(elms.size(), elms.size());
		for (int i = 0; i < elms.size(); ++i) {
			final UncertainValue n = new UncertainValue(men.get(elms.get(i)));
			vec.setEntry(i, n.doubleValue());
			cov.setEntry(i, i, n.variance());
		}
		return massFraction(mat, elms, vec, cov, Collections.emptyMap(), null, otherElmRule);
	}

	public static Composition massFraction(//
			final Material mat, //
			final UncertainValuesBase<EPMALabel> uvs //
	) throws ArgumentException {
		final List<Element> elms = new ArrayList<>(mat.getElementSet());
		final List<EPMALabel> inputs = new ArrayList<>();
		inputs.addAll(MaterialLabel.buildMassFractionTags(mat));
		inputs.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		final RealVector vals = uvs.extractValues(inputs);
		final RealMatrix covs = uvs.extractCovariances(inputs);
		return massFraction(mat, elms, vals, covs, null);
	}

	public static Composition massFraction(//
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov) throws ArgumentException {
		return massFraction(new Material(html, elms), elms, vals, cov, Collections.emptyMap(), null, null);
	}

	public static Composition massFraction(//
			final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealMatrix cov, //
			final Number density //
	) throws ArgumentException {
		return massFraction(new Material(html, elms), elms, vals, cov, Collections.emptyMap(), density, null);
	}

	public static Composition massFraction(final String html, //
			final List<Element> elms, //
			final RealVector vals, //
			final RealVector vars) throws ArgumentException {
		return massFraction(new Material(html, elms), elms, vals, vars, Collections.emptyMap());
	}

	public static Composition massFraction(final String html, //
			final Map<Element, ? extends Number> men //
	) throws ArgumentException {
		return massFraction(new Material(html, men.keySet()), men, null);
	}

	public static Composition massFraction( //
			final String html, //
			final Map<Element, ? extends Number> men, //
			final Map<Element, Number> atomicWeights) throws ArgumentException {
		final List<Element> elms = new ArrayList<>(men.keySet());
		final RealVector vec = new ArrayRealVector(elms.size());
		final RealVector vars = new ArrayRealVector(elms.size());
		for (int i = 0; i < elms.size(); ++i) {
			final UncertainValue n = new UncertainValue(men.get(elms.get(i)));
			vec.setEntry(i, n.doubleValue());
			vars.setEntry(i, n.variance());
		}
		return Composition.massFraction(new Material(html, elms), elms, vec, vars, atomicWeights);
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
	 * @param elm Element
	 * @return An Composition in MassFraction representation.
	 * @throws ArgumentException
	 */
	public static Composition pureElement(final Element elm) {
		try {
			return atomFraction("Pure " + elm.getAbbrev(), Collections.singletonMap(elm, 1.0));
		} catch (final ArgumentException e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * @param elm    Element
	 * @param purity [0.0,1.0]
	 * @return Composition
	 * @throws ArgumentException
	 */
	public static Composition pureElement(final Element elm, final double purity) {
		final Map<Element, Number> men = new TreeMap<>();
		men.put(elm, new UncertainValue(purity,  1.0 - purity));
		try {
			return massFraction("Pure " + elm.getAbbrev(),
					Collections.singletonMap(elm, new UncertainValue(purity, 1.0 - purity)), Collections.emptyMap());
		} catch (final ArgumentException e) {
			e.printStackTrace();
		}
		return null;
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic) throws ArgumentException {
		return stoichiometry(html, stoic, Collections.emptyMap(), null);
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic, //
			final Map<Element, Number> atomicWeights, //
			final Number density //
	) throws ArgumentException {
		final Material mat = new Material(html, stoic.keySet());
		final Map<MaterialLabel, Number> msn = new HashMap<>();
		for (final Map.Entry<Element, Integer> me : stoic.entrySet())
			msn.put(new MaterialLabel.Stoichiometry(mat, me.getKey()), me.getValue());
		msn.putAll(buildAtomicWeights(mat, atomicWeights));
		final UncertainValues<MaterialLabel> input = new UncertainValues<MaterialLabel>(msn);
		final StoichiometryToComposition slmjf = new StoichiometryToComposition(mat);
		return new Composition(Representation.Stoichiometry, mat, slmjf, input, density);
	}

	public static Composition stoichiometry( //
			final String html, //
			final Map<Element, Integer> stoic, //
			final Number density //
	) throws ArgumentException {
		return stoichiometry(html, stoic, Collections.emptyMap(), density);
	}

	static private Set<Element> allElements(final Collection<Element> elms,
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> otherElm) {
		final Set<Element> allElms = new HashSet<>();
		allElms.addAll(elms);
		if (otherElm != null)
			for (final Object lbl : otherElm.getOutputLabels())
				if (lbl instanceof MassFraction)
					allElms.add(((MassFraction) lbl).getElement());
		return allElms;
	}

	private static Map<MaterialLabel, Number> buildAtomicWeights( //
			final Material mat, final Map<Element, Number> atomicWeights) {
		final Map<MaterialLabel, Number> man = new HashMap<>();
		for (final Element elm : mat.getElementSet())
			man.put(MaterialLabel.buildAtomicWeightTag(mat, elm), mat.getAtomicWeight(elm));
		return man;
	}

	@SafeVarargs
	private static Map<Composition, Number> convert(final Pair<Composition, Number>... comps) {
		final Map<Composition, Number> res = new HashMap<>();
		for (final Pair<Composition, Number> pr : comps)
			res.put(pr.getFirst(), pr.getSecond());
		return res;
	}

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

	private final Material mMaterial;

	private final Representation mPrimary;

	private Optional<Number> mDensity;

	private final int mHashCode;

	private Composition(final Representation rep, final Material mat,
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> func,
			final UncertainValuesBase<MaterialLabel> inputs) throws ArgumentException {
		this(rep, mat, func, inputs, null);
	}

	private Composition(//
			final Representation rep, //
			final Material mat, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> func, //
			final UncertainValuesBase<MaterialLabel> inputs, //
			final Number density //
	) throws ArgumentException {
		super(func, inputs);
		mMaterial = mat;
		mPrimary = rep;
		mDensity = Optional.ofNullable(density);
		// Don't include density in hash
		mHashCode = super.hashCode() ^ Objects.hash(mMaterial, mPrimary);
	}

	public UncertainValues<Stoichiometry> asStoichiometry() {
		try {
			return UncertainValues.<Stoichiometry>extract(MaterialLabel.buildStoichiometryTags(mMaterial), this);
		} catch (final ArgumentException e) {
			// Should never happen!!!
			e.printStackTrace();
		}
		return null;
	}

	public void clearDensity() {
		mDensity = Optional.empty();
	}

	public boolean contains(final Element elm) {
		return mMaterial.contains(elm);
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
		return Objects.equals(mMaterial, other.mMaterial) //
				&& Objects.equals(mPrimary, other.mPrimary) //
				&& mDensity.equals(other.mDensity) //
				&& super.equals(other); // ;
	}

	public UncertainValue getAnalyticalTotal() {
		return getValue(MaterialLabel.buildAnalyticalTotalTag(getMaterial()));
	}

	public UncertainValue getAtomFraction(final Element elm) {
		final MaterialLabel tag = MaterialLabel.buildAtomFractionTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValue getAtomicWeight(final Element elm) {
		final MaterialLabel tag = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValues<AtomicWeight> getAtomicWeights() {
		try {
			return UncertainValues.extract(MaterialLabel.atomWeightTags(this.getMaterial()), this);
		} catch (final ArgumentException e) {
			// Should never happen!!!
			e.printStackTrace();
		}
		return null;
	}

	public Optional<Number> getDensity() {
		return mDensity;
	}

	public List<Element> getElementList() {
		return new ArrayList<>(mMaterial.getElementSet());
	}

	public Set<Element> getElementSet() {
		return mMaterial.getElementSet();
	}

	public List<? extends MaterialLabel> getLabels(final Representation rep) {
		switch (rep) {
		case AtomFraction:
			return MaterialLabel.buildAtomFractionTags(mMaterial);
		case MassFraction:
			return MaterialLabel.buildMassFractionTags(mMaterial);
		case Mixture: {
			final List<MaterialLabel> res = new ArrayList<>();
			for (final MaterialLabel lbl : getLabels())
				if (lbl instanceof MaterialMassFraction)
					res.add(lbl);
			return res;
		}
		case NormalizedMassFraction:
			return MaterialLabel.buildNormalizedMassFractionTags(mMaterial);
		case Stoichiometry:
			return MaterialLabel.buildStoichiometryTags(mMaterial);
		default:
			assert false;
			return null;
		}

	}

	public UncertainValue getMassFraction(final Element elm) {
		final MaterialLabel tag = MaterialLabel.buildMassFractionTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public Material getMaterial() {
		return mMaterial;
	}

	public UncertainValue getMeanAtomicNumber() {
		return getValue(MaterialLabel.buildMeanAtomicNumberTag(getMaterial()));
	}

	public UncertainValue getMeanAtomicWeight() {
		return getValue(MaterialLabel.buildMeanAtomicWeighTag(getMaterial()));
	}

	public Representation getNativeRepresentation() {
		return mPrimary;
	}

	public UncertainValue getNomalizedMassFraction(final Element elm) {
		final MaterialLabel tag = MaterialLabel.buildNormalizedMassFractionTag(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	public UncertainValue getStoichiometry(final Element elm) {
		final MaterialLabel tag = new MaterialLabel.Stoichiometry(mMaterial, elm);
		return indexOf(tag) == -1 ? UncertainValue.ZERO : getValue(tag);
	}

	@Override
	public int hashCode() {
		return mHashCode;
	}

	public boolean hasRepresentation(final Representation rep) {
		final List<? extends MaterialLabel> lbls = getLabels(rep);
		for (final MaterialLabel lbl : lbls)
			if (indexOf(lbl) == -1)
				return false;
		return true;
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
	 * @param n In g/cm<sup>3</sup>
	 */
	public void setDensity(final Number n) {
		mDensity = Optional.of(n);
	}

	@Override
	public String toHTML(final Mode mode) {
		return toHTML(mPrimary, mode);
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

	public UncertainValues<MassFraction> toMassFraction() {
		try {
			return UncertainValues.extract(massFractionTags(), this);
		} catch (final ArgumentException e) {
			// Should never happen
			e.printStackTrace();
		}
		return null;
	}

	public String toMassFractionString() {
		final StringBuffer sb = new StringBuffer();
		for (final Element elm : getElementList()) {
			if (sb.length() == 0)
				sb.append("[(");
			else
				sb.append("),(C[");
			sb.append(elm.getAbbrev());
			sb.append("]=");
			sb.append(MASS_FRACTION_FORMAT.format(getMassFraction(elm).doubleValue()));
		}
		sb.append(")]");
		return sb.toString();
	}

	@Override
	public String toString() {
		return "C[" + mMaterial.toString() + "]";
	}

	private <H> List<H> extractLabels(final Class<H> cls) {
		final List<H> res = new ArrayList<>();
		for (final Object lbl : getLabels())
			if (cls.isInstance(lbl))
				res.add(cls.cast(lbl));
		return res;
	}

	private String toHTMLasAtomFraction(final Mode mode) {
		assert hasRepresentation(Representation.AtomFraction);
		assert hasRepresentation(Representation.MassFraction);
		switch (mode) {
		case TERSE:
			return mMaterial.getHTMLName();
		case NORMAL: {
			final Table t = new Table();
			final int ROW_LEN = 6;
			final List<Item> hdr = new ArrayList<>();
			hdr.add(Table.th("Element"));
			hdr.add(Table.thc("Z"));
			hdr.add(Table.thc("Mass<br/>Fraction"));
			hdr.add(Table.thc("Uncertainty<br/>(1 &sigma;)"));
			hdr.add(Table.thc("Atomic<br/>Fraction"));
			hdr.add(Table.thc("Uncertainty<br/>(1 &sigma;)"));
			assert hdr.size() == ROW_LEN;
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
				assert row.size() == ROW_LEN;
				t.addRow(row);
			}
			if (mDensity.isPresent()) {
				final List<Item> drow = new ArrayList<>();
				drow.add(Table.td("Density")); // Element
				final UncertainValue uv = new UncertainValue(mDensity.get());
				drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", ROW_LEN - 1));
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
				final int ROW_LEN = 8;
				t.addRow(Table.th(mMaterial.getHTMLName(), ROW_LEN));
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
					assert row.size() == ROW_LEN;
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
				assert row.size() == ROW_LEN;
				t.addRow(row);
				if (mDensity.isPresent()) {
					final List<Item> drow = new ArrayList<>();
					drow.add(Table.td("Density")); // Element
					final UncertainValue uv = new UncertainValue(mDensity.get());
					drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", ROW_LEN - 1));
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

	private String toHTMLasMixture(final Mode mode) {
		assert hasRepresentation(Representation.Mixture);
		switch (mode) {
		default:
		case NORMAL:
			final Table t = new Table();
			t.addRow(Table.th("Material"), Table.th("Mass Fraction"), Table.th("Uncertainty"));
			for (final MaterialMassFraction mmft : materialFractionTags()) {
				final UncertainValue uv = getValue(mmft);
				t.addRow( //
						Table.td(mmft.toHTML(Mode.TERSE)), //
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
			final int ROW_LEN = 4;
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
				drow.add(Table.td(uv.toHTML(Mode.TERSE) + " g/cm<sup>3</sup>", ROW_LEN - 1));
				t.addRow(drow);
			}
			return HTML.subHeader(mMaterial.getHTMLName()) + t.toHTML(Mode.VERBOSE);
		}
		}
	}
}
