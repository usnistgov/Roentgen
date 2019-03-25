package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.EPMALabel.MaterialMAC;
import gov.nist.microanalysis.roentgen.DataStore.UniqueString;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.Layer;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;

/**
 * A matrix correction model is a model that computes the values associated with
 * {@link ZAFMultiLineLabel} and {@link KRatioLabel} for a set of
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
		extends CompositeMeasurementModel<EPMALabel> {

	public static class ChiLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		ChiLabel(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
		) {
			super("&chi;", mcd, cxr);
		}

	}

	public static class CoatingThickness extends EPMALabel.BaseLabel<UniqueString, Object, Object> {

		private CoatingThickness(
				final Layer layer
		) {
			super("t<sub>Coating</sub>", layer.getName());
		}
	}

	public static class FofChiLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		private FofChiLabel(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
		) {
			super("F(&chi;)", mcd, cxr);
		}
	}

	public static class FofChiReducedLabel extends MatrixCorrectionDatumTag2<CharacteristicXRay> {

		private FofChiReducedLabel(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
		) {
			super("F<sub>r</sub>(&chi;)", mcd, cxr);
		}

	}

	public static class IonizationExponentLabel extends EPMALabel.BaseLabel<AtomicShell, Object, Object> {

		public IonizationExponentLabel(
				final AtomicShell sh
		) {
			super("m", sh);
		}
	}

	public static class MaterialBasedLabel extends EPMALabel.BaseLabel<Material, Object, Object> {
		public MaterialBasedLabel(
				final String name, final Material mcd
		) {
			super(name, mcd);
		}
	}

	public static class MatrixCorrectionDatumLabel extends EPMALabel.BaseLabel<MatrixCorrectionDatum, Object, Object> {
		public MatrixCorrectionDatumLabel(
				final String name, final MatrixCorrectionDatum mcd
		) {
			super(name, mcd);
		}
	}

	public static class MatrixCorrectionDatumTag2<H> extends EPMALabel.BaseLabel<MatrixCorrectionDatum, H, Object> {

		MatrixCorrectionDatumTag2(
				final String name, final MatrixCorrectionDatum mcd, final H obj2
		) {
			super(name, mcd, obj2);
		}
	}

	public static class Phi0Label extends MatrixCorrectionDatumTag2<AtomicShell> {

		Phi0Label(
				final MatrixCorrectionDatum mcd, final AtomicShell shell
		) {
			super("&phi;<sub>0</sub>", mcd, shell);
		}

	}

	public static class RoughnessLabel extends MatrixCorrectionDatumLabel {
		public RoughnessLabel(
				final MatrixCorrectionDatum mcd
		) {
			super("dz", mcd);
		}
	}

	public static class XRayWeightLabel extends EPMALabel.BaseLabel<CharacteristicXRay, Object, Object> {

		public XRayWeightLabel(
				final CharacteristicXRay cxr
		) {
			super("w", cxr);
		}
	}

	public static class ZAFLabel
			extends EPMALabel.BaseLabel<MatrixCorrectionDatum, MatrixCorrectionDatum, CharacteristicXRay> {

		private ZAFLabel(
				final String name, final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std,
				final CharacteristicXRay cxr
		) {
			super(name, unk, std, cxr);
		}
	}

	public static class ZAFMultiLineLabel
			extends EPMALabel.BaseLabel<UnknownMatrixCorrectionDatum, StandardMatrixCorrectionDatum, ElementXRaySet> {

		private ZAFMultiLineLabel(
				final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final CharacteristicXRay cxr
		) {
			super("ZAF", unk, std, new ElementXRaySet(cxr));
		}

		private ZAFMultiLineLabel(
				final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final ElementXRaySet exrs
		) {
			super("ZAF", unk, std, exrs);
		}

		public ElementXRaySet getElementXRaySet() {
			return getObject3();
		}

		public StandardMatrixCorrectionDatum getStandard() {
			return getObject2();
		}

		public UnknownMatrixCorrectionDatum getUnknown() {
			return getObject1();
		}
	}

	private static class ElementLabel extends EPMALabel.BaseLabel<Element, Object, Object> {

		private ElementLabel(
				final String name, final Element obj
		) {
			super(name, obj);
		}
	}

	/**
	 * Only StandardComposition and UnknownComposition.
	 *
	 * @return
	 */

	static public ZAFLabel aLabel(
			final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final CharacteristicXRay cxr
	) {
		return new MatrixCorrectionModel2.ZAFLabel("A", unk, std, cxr);
	}

	static public MatrixCorrectionDatumTag2<?> atomicNumberLabel(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
	) {
		return characterisiticLabel("Z", mcd, cxr);
	}

	static public MatrixCorrectionDatumLabel beamEnergyLabel(
			final MatrixCorrectionDatum datum
	) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumLabel("E<sub>0</sub>", datum);
	}

	static public MatrixCorrectionDatumTag2<?> characterisiticLabel(
			final String name, final MatrixCorrectionDatum mcd, final CharacteristicXRay other
	) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumTag2<CharacteristicXRay>(name, mcd, other);
	}

	static public ChiLabel chiLabel(
			final MatrixCorrectionDatum comp, final CharacteristicXRay other
	) {
		return new MatrixCorrectionModel2.ChiLabel(comp, other);
	}

	static public CoatingThickness coatingMassThickness(
			final Layer layer
	) {
		return new CoatingThickness(layer);
	}

	static public FofChiLabel FofChiLabel(
			final MatrixCorrectionDatum comp, final CharacteristicXRay other
	) {
		return new MatrixCorrectionModel2.FofChiLabel(comp, other);
	}

	static public FofChiReducedLabel FofChiReducedLabel(
			final MatrixCorrectionDatum comp, final CharacteristicXRay other
	) {
		return new MatrixCorrectionModel2.FofChiReducedLabel(comp, other);
	}

	static public MatrixCorrectionDatumTag2<?> FxFLabel(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
	) {
		return characterisiticLabel("F(&chi;)/F", mcd, cxr);
	}

	static public MaterialMAC matMacLabel(
			final Material mat, final CharacteristicXRay cxr
	) {
		return new MaterialMAC(mat, cxr);
	}

	static public ElementLabel meanIonizationLabel(
			final Element elm
	) {
		return new MatrixCorrectionModel2.ElementLabel("J", elm);
	}

	static public Phi0Label phi0Label(
			final MatrixCorrectionDatum comp, final AtomicShell other
	) {
		return new MatrixCorrectionModel2.Phi0Label(comp, other);
	}

	public static MatrixCorrectionModel2.RoughnessLabel roughnessLabel(
			final MatrixCorrectionDatum mcd
	) {
		return new MatrixCorrectionModel2.RoughnessLabel(mcd);
	}

	static public MatrixCorrectionDatumTag2<?> shellLabel(
			final String name, final MatrixCorrectionDatum mcd, final AtomicShell other
	) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumTag2<>(name, mcd, other);
	}

	static public MatrixCorrectionDatumLabel takeOffAngleLabel(
			final MatrixCorrectionDatum datum
	) {
		return new MatrixCorrectionModel2.MatrixCorrectionDatumLabel("TOA", datum);
	}

	static public ZAFMultiLineLabel zafLabel(
			final KRatioLabel krl
	) {
		return new ZAFMultiLineLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet());
	}

	static public ZAFMultiLineLabel zafLabel(
			final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final CharacteristicXRay cxr
	) {
		return new ZAFMultiLineLabel(unk, std, cxr);
	}

	static public ZAFMultiLineLabel zafLabel(
			final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std, final ElementXRaySet exrs
	) {
		return new ZAFMultiLineLabel(unk, std, exrs);
	}

	static public ZAFLabel zLabel(
			final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final CharacteristicXRay cxr
	) {
		return new MatrixCorrectionModel2.ZAFLabel("Z", unk, std, cxr);
	}

	protected final Set<KRatioLabel> mKRatios;
	// The types of variables to compute Jacobian elements.

	private final Material mUnknownMaterial;

	public MatrixCorrectionModel2(
			final String string, final Set<KRatioLabel> kratios, //
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps //
	) throws ArgumentException {
		this(string, kratios, buildSteps, CompositeMeasurementModel.allOutputs(buildSteps));
	}

	public MatrixCorrectionModel2(
			final String name, //
			final Set<KRatioLabel> kratios, //
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> list, //
			final List<? extends EPMALabel> outputLabels //
	) throws ArgumentException {
		super(name, list, outputLabels);
		mKRatios = kratios;
		mUnknownMaterial = mKRatios.iterator().next().getUnknown().getMaterial();
		// Validate the inputs...
		final List<? extends Object> outputTags = getOutputLabels();
		final List<? extends Object> inputTags = getInputLabels();

		for (final KRatioLabel krl : mKRatios) {
			if (!krl.getUnknown().getMaterial().equals(mUnknownMaterial))
				throw new ArgumentException("The k-ratios must all be associated with the same material.");
			final ZAFMultiLineLabel mct = zafLabel(krl);
			if (!outputTags.contains(mct))
				throw new ArgumentException(toString() + " does not calculate the required output " + mct.toString());
			final Composition stdComp = krl.getStandard().getComposition();

			for (final Object inpLabel : stdComp.getInputLabels()) {
				if (!inputTags.contains(inpLabel))
					throw new ArgumentException(toString() + " must take " + inpLabel.toString() + " as an argument.");
			}
		}
	}

	abstract public UncertainValuesBase<EPMALabel> buildInput(
			Material unknownMat //
	) throws ArgumentException;

	public Set<Element> getElementSet() {
		final Set<Element> res = new HashSet<>();
		for (final KRatioLabel krl : mKRatios)
			res.add(krl.getElement());
		return Collections.unmodifiableSet(res);
	}

	public Set<KRatioLabel> getKRatios() {
		return Collections.unmodifiableSet(mKRatios);
	}

	public Set<KRatioLabel> getKRatios(
			final Element elm
	) {
		final Set<KRatioLabel> res = new HashSet<>();
		for (final KRatioLabel krl : mKRatios)
			if (krl.getElement().equals(elm))
				res.add(krl);
		return Collections.unmodifiableSet(res);
	}

	public Material getUnknownMaterial() {
		return mUnknownMaterial;
	}

}
