package gov.nist.microanalysis.roentgen.wds;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.CompositeMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ParallelMeasurementModelBuilder;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;

/**
 * Calculates the normalized (counts per (nanoamp second)) and continuum
 * corrected intensity from a three point (on-peak, low background and high
 * background) WDS measurement.
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class TwoPointContinuumModel //
		extends CompositeMeasurementModel<ModelLabels<?, ?>> {

	public static class ComputePeakContinuum
			extends ExplicitMeasurementModel<ModelLabels<?, ?>, ModelLabels<?, ?>> {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;

		public ComputePeakContinuum(
				final MatrixCorrectionDatum mcd, //
				final CharacteristicXRay cxr
		) throws ArgumentException {
			super(buildPCInputs(mcd, cxr),
					Collections.singletonList(ModelLabels.buildNormCharacteristicIntensity(mcd, cxr)));
			mMcd = mcd;
			mCxr = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final int[] iNormI = { //
					inputIndex(ModelLabels.buildNormalizedIntensity(mMcd, mCxr, LOW_BACK)), //
					inputIndex(ModelLabels.buildNormalizedIntensity(mMcd, mCxr, ON_PEAK)), //
					inputIndex(ModelLabels.buildNormalizedIntensity(mMcd, mCxr, HIGH_BACK)) //
			};
			final int[] rI = { //
					inputIndex(ModelLabels.buildSpectrometerPosition(mMcd, mCxr, LOW_BACK)), //
					inputIndex(ModelLabels.buildSpectrometerPosition(mMcd, mCxr, ON_PEAK)), //
					inputIndex(ModelLabels.buildSpectrometerPosition(mMcd, mCxr, HIGH_BACK)) //
			};
			final int LOW = 0, PEAK = 1, HIGH = 2;
			final double[] iNorm = { //
					point.getEntry(iNormI[LOW]), //
					point.getEntry(iNormI[PEAK]), //
					point.getEntry(iNormI[HIGH]) };
			final double[] r = { //
					point.getEntry(rI[LOW]), //
					point.getEntry(rI[PEAK]), //
					point.getEntry(rI[HIGH]) //
			};

			final RealVector rv = new ArrayRealVector(1);
			final double iPeak = iNorm[PEAK]
					- (iNorm[LOW] + ((iNorm[HIGH] - iNorm[LOW]) / (r[HIGH] - r[LOW])) * (r[PEAK] - r[LOW]));
			rv.setEntry(0, iPeak);
			final RealMatrix rm = MatrixUtils.createRealMatrix(1, 6);
			// WRT normalized intensity
			rm.setEntry(0, iNormI[LOW], (r[PEAK] - r[HIGH]) / (r[HIGH] - r[LOW])); // OK
			rm.setEntry(0, iNormI[PEAK], 1.0); // OK
			rm.setEntry(0, iNormI[HIGH], (r[PEAK] - r[LOW]) / (r[HIGH] - r[LOW]));
			// WRT spectrometer position
			rm.setEntry(0, rI[LOW], (iNorm[HIGH] - iNorm[LOW]) * (r[HIGH] - r[PEAK]) / Math.pow(r[HIGH] - r[LOW], 2.0)); // OK
			rm.setEntry(0, rI[PEAK], (iNorm[LOW] - iNorm[HIGH]) / (r[HIGH] - r[LOW])); // OK
			rm.setEntry(0, rI[HIGH], (iNorm[LOW] - iNorm[HIGH]) * (r[LOW] - r[PEAK]) / Math.pow(r[HIGH] - r[LOW], 2.0)); // OK
			return Pair.create(rv, rm);
		}

	}

	private static class NormalizeIntensity
			extends ExplicitMeasurementModel<ModelLabels<?, ?>, ModelLabels<?, ?>> {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;
		private final int mIndex;

		private NormalizeIntensity(
				final MatrixCorrectionDatum mcd, //
				final CharacteristicXRay cxr, //
				final int index
		) throws ArgumentException {
			super(buildNIInputs(mcd, cxr, index),
					Collections.singletonList(ModelLabels.buildNormalizedIntensity(mcd, cxr, index)));
			mMcd = mcd;
			mCxr = cxr;
			mIndex = index;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final int iI = inputIndex(ModelLabels.buildRawIntensity(mMcd, mCxr, mIndex));
			final int ltI = inputIndex(ModelLabels.buildLiveTime(mMcd, mCxr, mIndex));
			final int pcI = inputIndex(ModelLabels.buildProbeCurrent(mMcd, mCxr, mIndex));

			final double i = point.getEntry(iI);
			final double lt = point.getEntry(ltI);
			final double pc = point.getEntry(pcI);

			final RealVector rv = new ArrayRealVector(1);
			rv.setEntry(0, i / (lt * pc));

			final RealMatrix rm = MatrixUtils.createRealMatrix(1, 3);
			rm.setEntry(0, iI, 1.0 / (lt * pc));
			rm.setEntry(0, ltI, -i / (lt * lt * pc));
			rm.setEntry(0, pcI, -1 / (lt * pc * pc));

			return Pair.create(rv, rm);
		}
	}

	public static final int HIGH_BACK = 1;

	public static final int ON_PEAK = 0;

	public static final int LOW_BACK = -1;

	static public List<? extends ExplicitMeasurementModel<? extends ModelLabels<?, ?>, ? extends ModelLabels<?, ?>>> buildModel(
			//
			final MatrixCorrectionDatum mcd, //
			final CharacteristicXRay cxr //
	) throws ArgumentException {
		final List<ExplicitMeasurementModel<? extends ModelLabels<?, ?>, ? extends ModelLabels<?, ?>>> res = new ArrayList<>();
		final ParallelMeasurementModelBuilder<ModelLabels<?, ?>, ModelLabels<?, ?>> builder = //
				new ParallelMeasurementModelBuilder<ModelLabels<?, ?>, ModelLabels<?, ?>>("Normalizer");
		builder.add(new NormalizeIntensity(mcd, cxr, LOW_BACK));
		builder.add(new NormalizeIntensity(mcd, cxr, ON_PEAK));
		builder.add(new NormalizeIntensity(mcd, cxr, HIGH_BACK));
		res.add(builder.build());
		res.add(new ComputePeakContinuum(mcd, cxr));
		return res;
	}

	static private List<ModelLabels<?, ?>> buildNIInputs(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
	) {
		final List<ModelLabels<?, ?>> res = new ArrayList<>();
		res.add(ModelLabels.buildRawIntensity(mcd, cxr, index));
		res.add(ModelLabels.buildLiveTime(mcd, cxr, index));
		res.add(ModelLabels.buildProbeCurrent(mcd, cxr, index));
		return res;
	}

	static private List<ModelLabels<?, ?>> buildPCInputs(
			final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
	) {
		final List<ModelLabels<?, ?>> res = new ArrayList<>();
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, LOW_BACK));
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, ON_PEAK));
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, HIGH_BACK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, LOW_BACK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, ON_PEAK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, HIGH_BACK));
		return res;
	}

	/**
	 * @param mcd Details the material and conditions
	 * @param cxr The characteristic x-ray
	 * @throws ArgumentException
	 */
	public TwoPointContinuumModel(
			//
			final MatrixCorrectionDatum mcd, //
			final CharacteristicXRay cxr //
	) throws ArgumentException {
		super("TwoPoint", buildModel(mcd, cxr));
	}
}
