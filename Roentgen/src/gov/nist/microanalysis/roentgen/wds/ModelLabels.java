package gov.nist.microanalysis.roentgen.wds;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author Nicholas W. M. Ritchie
 * @param <H>
 *
 */
public class ModelLabels<H, J> extends BaseLabel<MatrixCorrectionDatum, H, J> {

	private ModelLabels(final String name, final MatrixCorrectionDatum obj1, final H obj2) {
		super(name, obj1, obj2);
	}

	private ModelLabels(final String name, final MatrixCorrectionDatum obj1, final H obj2, final J obj3) {
		super(name, obj1, obj2, obj3);
	}

	public static Object buildNormCharacteristicIntensity(final MatrixCorrectionDatum mcd, final ElementXRaySet exrs) {
		return new ModelLabels<ElementXRaySet, Object>("I<sub>norm,char</sub>", mcd, exrs);
	}

	public static Object buildNormCharacteristicIntensity(final MatrixCorrectionDatum mcd,
			final CharacteristicXRay cxr) {
		return buildNormCharacteristicIntensity(mcd, new ElementXRaySet(cxr));
	}

	public static Object buildSpectrometerPosition(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr,
			final int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("R", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildRawIntensity(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr,
			final int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("I", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildLiveTime(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("LT", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildProbeCurrent(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr,
			final int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("PC", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildNormalizedIntensity(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr,
			final int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("I<sub>norm</sub>", mcd, cxr, Integer.valueOf(index));
	}

}
