package gov.nist.microanalysis.roentgen.wds;

import gov.nist.microanalysis.roentgen.EPMALabel.BaseLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;

/**
 * @author Nicholas W. M. Ritchie
 * @param <H>
 *
 */
public class ModelLabels<H, J> extends BaseLabel<MatrixCorrectionDatum, H, J> {

	private ModelLabels(
			final String name, final MatrixCorrectionDatum obj1, final H obj2
	) {
		super(name, obj1, obj2);
	}

	private ModelLabels(
			final String name, final MatrixCorrectionDatum obj1, final H obj2, final J obj3
	) {
		super(name, obj1, obj2, obj3);
	}

	public static class NormalizedCharacteristicIntensity extends ModelLabels<CharacteristicXRay, Object> {
		private NormalizedCharacteristicIntensity(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
		) {
			super("I<sub>norm,char</sub>", mcd, cxr);
		}
	}

	public static NormalizedCharacteristicIntensity buildNormCharacteristicIntensity(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
	) {
		return new NormalizedCharacteristicIntensity(mcd, cxr);
	}

	public static class SpectrometerPosition extends ModelLabels<CharacteristicXRay, Integer> {
		private SpectrometerPosition(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			super("L", mcd, cxr, Integer.valueOf(index));
		}
	}

	public static SpectrometerPosition buildSpectrometerPosition(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		return new SpectrometerPosition(mcd, cxr, Integer.valueOf(index));
	}

	public static class RawIntensity extends ModelLabels<CharacteristicXRay, Integer> {
		private RawIntensity(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			super("I", mcd, cxr, index);
		}

	}

	public static RawIntensity buildRawIntensity(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		return new RawIntensity(mcd, cxr, Integer.valueOf(index));
	}

	public static class LiveTime extends ModelLabels<CharacteristicXRay, Integer> {
		private LiveTime(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			super("LT", mcd, cxr, Integer.valueOf(index));
		}
	}

	public static LiveTime buildLiveTime(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		return new LiveTime(mcd, cxr, Integer.valueOf(index));
	}

	public static class ProbeCurrent extends ModelLabels<CharacteristicXRay, Integer> {

		private ProbeCurrent(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			super("PC", mcd, cxr, Integer.valueOf(index));
		}

	}

	public static ProbeCurrent buildProbeCurrent(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		return new ProbeCurrent(mcd, cxr, Integer.valueOf(index));
	}

	public static class NormalizedIntensity extends ModelLabels<CharacteristicXRay, Integer> {

		private NormalizedIntensity(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			super("I<sub>norm</sub>", mcd, cxr, Integer.valueOf(index));
		}

	}

	public static NormalizedIntensity buildNormalizedIntensity(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		return new NormalizedIntensity(mcd, cxr, Integer.valueOf(index));
	}

}
