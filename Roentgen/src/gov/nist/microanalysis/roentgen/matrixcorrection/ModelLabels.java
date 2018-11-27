/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author nicholas
 * @param <H>
 *
 */
public class ModelLabels<H, J> extends BaseLabel<MatrixCorrectionDatum, H, J> {

	private ModelLabels(String name, MatrixCorrectionDatum obj1, H obj2) {
		super(name, obj1, obj2);
	}

	private ModelLabels(String name, MatrixCorrectionDatum obj1, H obj2, J obj3) {
		super(name, obj1, obj2, obj3);
	}

	public static Object buildNormCharacteristicIntensity(MatrixCorrectionDatum mcd, ElementXRaySet exrs) {
		return new ModelLabels<ElementXRaySet, Object>("I<sub>norm,char</sub>", mcd, exrs);
	}

	public static Object buildNormCharacteristicIntensity(MatrixCorrectionDatum mcd, CharacteristicXRay cxr) {
		return buildNormCharacteristicIntensity(mcd, new ElementXRaySet(cxr));
	}

	public static Object buildSpectrometerPosition(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("R", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildRawIntensity(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("I", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildLiveTime(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("LT", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildProbeCurrent(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("PC", mcd, cxr, Integer.valueOf(index));
	}

	public static Object buildNormalizedIntensity(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		return new ModelLabels<CharacteristicXRay, Integer>("I<sub>norm</sub>", mcd, cxr, Integer.valueOf(index));
	}

}
