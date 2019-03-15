package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.DataStore.UniqueString;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * A matrix correction model is a model that computes the values associated with
 * {@link MatrixCorrectionLabel} and {@link KRatioLabel} for a set of
 * {@link ElementXRaySet} objects associated with {@link MatrixCorrectionDatum}s
 * associated with Standards relative to a {@link MatrixCorrectionDatum}
 * associated with an unknown.
 *
 * A matrix correction model may also calculate values for {@link KRatioLabel}
 * and other related tags.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
abstract public class MatrixCorrectionModel2 //
		extends SerialLabeledMultivariateJacobianFunction {

	private static class ElementLabel extends BaseLabel<Element, Object, Object> {

		private ElementLabel(final String name, final Element obj) {
			super(name, obj);
		}
	}

	public static class MatrixCorrectionDatumLabel extends BaseLabel<MatrixCorrectionDatum, Object, Object> {
		public MatrixCorrectionDatumLabel(final String name, final MatrixCorrectionDatum mcd) {
			super(name, mcd);
		}
	}

	public static class RoughnessLabel extends MatrixCorrectionDatumLabel {
		public RoughnessLabel(final MatrixCorrectionDatum mcd) {
			super("dz", mcd);
		}
	}

	public static class MaterialBasedLabel extends BaseLabel<Material, Object, Object> {
		public MaterialBasedLabel(final String name, final Material mcd) {
			super(name, mcd);
		}
	}

	public static class MatrixCorrectionDatumTag2<H> extends BaseLabel<MatrixCorrectionDatum, H, Object> {

		MatrixCorrectionDatumTag2(final String name, final MatrixCorrectionDatum mcd, final H obj2) {
			super(name, mcd, obj2);
		}
	}

	public static class ChiLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		ChiLabel(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr) {
			super("&chi;", mcd, cxr);
		}

	}

	public static class FofChiLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		private FofChiLabel(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr) {
			super("F(&chi;)", mcd, cxr);
		}
	}

	public static class FofChiReducedLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		private FofChiReducedLabel(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr) {
			super("F<sub>r</sub>(&chi;)", mcd, cxr);
		}

	}

	public static class Phi0Label extends MatrixCorrectionDatumTag2<AtomicShell> {

		Phi0Label(final MatrixCorrectionDatum mcd, final AtomicShell shell) {
			super("&phi;<sub>0</sub>", mcd, shell);
		}

	}

	public static class ZAFLabel extends BaseLabel<MatrixCorrectionDatum, MatrixCorrectionDatum, CharacteristicXRay> {

		private ZAFLabel(final String name, final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std,
				final CharacteristicXRay cxr) {
			super(name, unk, std, cxr);
		}
	}

	public static class XRayWeightLabel extends BaseLabel<CharacteristicXRay, Object, Object> {

		public XRayWeightLabel(final CharacteristicXRay cxr) {
			super("w", cxr);
		}
	}

	public static class IonizationExponentLabel extends BaseLabel<AtomicShell, Object, Object> {

		public IonizationExponentLabel(final AtomicShell sh) {
			super("m", sh);
		}
	}

	/**
	 * A list of types of variables that can be included in the uncertainty
	 * calculation.
	 *
	 * @author Nicholas W. M. Ritchie
	 *
	 */
	public enum VariateX {

		MeanIonizationPotential("Mean ionization potential"), //
		MassAbsorptionCofficient("Mass absorption coefficent"), //
		StandardComposition("Standard composition"), //
		UnknownComposition("Unknown composition"), //
		BeamEnergy("Beam energy"), //
		TakeOffAngle("Take-off angle"), //
		WeightsOfLines("Weights-of-Lines"), //
		IonizationExponent("Ionization exponent"), //
		SurfaceRoughness("Surface roughness"), //
		SecondaryFluorescence("Secondary fluorescence"), //
		Coating("Coating"); //

		private final String mName;

		private VariateX(final String name) {
			mName = name;
		}

		@Override
		public String toString() {
			return mName;
		}

	}

	protected final Set<KRatioLabel> mKRatios;
	// The types of variables to compute Jacobian elements.

	private final Material mUnknownMaterial;

	public MatrixCorrectionModel2(//
			final String name, //
			final Set<KRatioLabel> kratios, //
			final List<LabeledMultivariateJacobianFunction> steps, //
			List<? extends Object> outputLabels //
	) throws ArgumentException {
		super(name, steps, outputLabels);
		mKRatios = kratios;
		mUnknownMaterial = mKRatios.iterator().next().getUnknown().getMaterial();
		// Validate the inputs...
		final List<? extends Object> outputTags = getOutputLabels();
		final List<? extends Object> inputTags = getInputLabels();

		for (final KRatioLabel krl : mKRatios) {
			if (!krl.getUnknown().getMaterial().equals(mUnknownMaterial))
				throw new ArgumentException("The k-ratios can only be associated with one material.");
			final MatrixCorrectionLabel mct = new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(),
					krl.getXRaySet());
			if (!outputTags.contains(mct))
				throw new ArgumentException(toString() + " does not calculate the required output " + mct.toString());
			final Composition stdComp = krl.getStandard().getComposition();

			for (final Object inpLabel : stdComp.getInputLabels()) {
				if (!inputTags.contains(inpLabel))
					throw new ArgumentException(toString() + " must take " + inpLabel.toString() + " as an argument.");
			}
			/*
			 * final Material unkMat = krl.getUnknown().getMaterial(); for (final Element
			 * elm : unkMat.getElementSet()) { final MaterialLabel.MassFraction mft =
			 * MaterialLabel.buildMassFractionTag(unkMat, elm); if
			 * (!inputTags.contains(mft)) throw new ArgumentException(toString() +
			 * " must take " + mft.toString() + " as an argument."); }
			 */ }
	}

	public MatrixCorrectionModel2(String string, Set<KRatioLabel> kratios,
			List<LabeledMultivariateJacobianFunction> buildSteps) throws ArgumentException {
		this(string, kratios,buildSteps, SerialLabeledMultivariateJacobianFunction.allOutputs(buildSteps, false));
	}

	public Set<KRatioLabel> getKRatios() {
		return Collections.unmodifiableSet(mKRatios);
	}

	public Material getUnknownMaterial() {
		return mUnknownMaterial;
	}

	public Set<KRatioLabel> getKRatios(final Element elm) {
		final Set<KRatioLabel> res = new HashSet<>();
		for (final KRatioLabel krl : mKRatios)
			if (krl.getElement().equals(elm))
				res.add(krl);
		return Collections.unmodifiableSet(res);
	}

	public Set<Element> getElementSet() {
		final Set<Element> res = new HashSet<>();
		for (final KRatioLabel krl : mKRatios)
			res.add(krl.getElement());
		return Collections.unmodifiableSet(res);
	}

	/**
	 * Only StandardComposition and UnknownComposition.
	 *
	 * @return
	 */

	static public Object aLabel(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std,
			final CharacteristicXRay cxr) {
		return new MatrixCorrectionModel2.ZAFLabel("A", unk, std, cxr);
	}

	public static MatrixCorrectionModel2.RoughnessLabel roughnessLabel(final MatrixCorrectionDatum mcd) {
		return new MatrixCorrectionModel2.RoughnessLabel(mcd);
	}

	static public Object chiLabel(final MatrixCorrectionDatum comp, final CharacteristicXRay other) {
		return new MatrixCorrectionModel2.ChiLabel(comp, other);
	}

	static public Object coatingMassThickness(final Layer layer) {
		return new BaseLabel<UniqueString, Object, Object>("Thickness", layer.getName());
	}

	static public Object FofChiLabel(final MatrixCorrectionDatum comp, final CharacteristicXRay other) {
		return new MatrixCorrectionModel2.FofChiLabel(comp, other);
	}

	static public Object FofChiReducedLabel(final MatrixCorrectionDatum comp, final CharacteristicXRay other) {
		return new MatrixCorrectionModel2.FofChiReducedLabel(comp, other);
	}

	static public Object phi0Label(final MatrixCorrectionDatum comp, final AtomicShell other) {
		return new MatrixCorrectionModel2.Phi0Label(comp, other);
	}

	static public Object FxFLabel(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr) {
		return characterisiticLabel("F(&chi;)/F", mcd, cxr);
	}

	static public Object atomicNumberLabel(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr) {
		return characterisiticLabel("Z", mcd, cxr);
	}

	static public Object beamEnergyLabel(final MatrixCorrectionDatum datum) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumLabel("E<sub>0</sub>", datum);
	}

	static public Object takeOffAngleLabel(final MatrixCorrectionDatum datum) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumLabel("TOA", datum);
	}

	static public Object meanIonizationLabel(final Element elm) {
		return new MatrixCorrectionModel2.ElementLabel("J", elm);
	}

	static public Object matMacLabel(final Material mat, final CharacteristicXRay cxr) {
		return new MaterialMACFunction.MaterialMAC(mat, cxr);
	}

	static public Object shellLabel(final String name, final MatrixCorrectionDatum mcd, final AtomicShell other) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumTag2<>(name, mcd, other);
	}

	static public Object characterisiticLabel(final String name, final MatrixCorrectionDatum mcd,
			final CharacteristicXRay other) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumTag2<CharacteristicXRay>(name, mcd, other);
	}

	static public Object zafLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final CharacteristicXRay cxr) {
		return new MatrixCorrectionLabel(unk, std, cxr);
	}

	static public Object zafLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final ElementXRaySet exrs) {
		return new MatrixCorrectionLabel(unk, std, exrs);
	}

	static public Object zafLabel(final KRatioLabel krl) {
		return new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet());
	}

	static public Object zLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final CharacteristicXRay cxr) {
		return new MatrixCorrectionModel2.ZAFLabel("Z", unk, std, cxr);
	}

	abstract public UncertainValuesBase buildInput(Map<MassFraction, ? extends Number> estUnknown)
			throws ArgumentException;

}
