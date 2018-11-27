package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunctionBuilder;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;

/**
 * Calculates the normalized (counts per (nanoamp second)) and continuum
 * corrected intensity from a three point (on-peak, low background and high
 * background) WDS measurement.
 * 
 * 
 * @author nicholas
 *
 */
public class TwoPointContinuumModel extends SerialLabeledMultivariateJacobianFunction {

	public static final int HIGH_BACK = 1;
	public static final int ON_PEAK = 0;
	public static final int LOW_BACK = -1;

	static private List<? extends Object> buildNIInputs(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
		List<Object> res = new ArrayList<>();
		res.add(ModelLabels.buildRawIntensity(mcd, cxr, index));
		res.add(ModelLabels.buildLiveTime(mcd, cxr, index));
		res.add(ModelLabels.buildProbeCurrent(mcd, cxr, index));
		return res;
	}

	static private List<? extends Object> buildPCInputs(MatrixCorrectionDatum mcd, CharacteristicXRay cxr) {
		List<Object> res = new ArrayList<>();
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, LOW_BACK));
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, ON_PEAK));
		res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, HIGH_BACK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, LOW_BACK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, ON_PEAK));
		res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, HIGH_BACK));
		return res;
	}

	private static class NormalizeIntensity extends LabeledMultivariateJacobianFunction {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;
		private final int mIndex;

		private NormalizeIntensity(MatrixCorrectionDatum mcd, CharacteristicXRay cxr, int index) {
			super(buildNIInputs(mcd, cxr, index),
					Collections.singletonList(ModelLabels.buildNormalizedIntensity(mcd, cxr, index)));
			mMcd = mcd;
			mCxr = cxr;
			mIndex = index;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
			int iI = inputIndex(ModelLabels.buildRawIntensity(mMcd, mCxr, mIndex));
			int ltI = inputIndex(ModelLabels.buildLiveTime(mMcd, mCxr, mIndex));
			int pcI = inputIndex(ModelLabels.buildProbeCurrent(mMcd, mCxr, mIndex));

			double i = point.getEntry(iI);
			double lt = point.getEntry(ltI);
			double pc = point.getEntry(pcI);

			RealVector rv = new ArrayRealVector(1);
			rv.setEntry(0, i / (lt * pc));

			RealMatrix rm = MatrixUtils.createRealMatrix(1, 3);
			rm.setEntry(0, iI, 1.0 / (lt * pc));
			rm.setEntry(0, ltI, -i / (lt * lt * pc));
			rm.setEntry(0, pcI, -1 / (lt * pc * pc));

			return Pair.create(rv, rm);
		}
	}

	public static class ComputePeakContinuum extends LabeledMultivariateJacobianFunction {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;

		public ComputePeakContinuum(MatrixCorrectionDatum mcd, CharacteristicXRay cxr) {
			super(buildPCInputs(mcd, cxr),
					Collections.singletonList(ModelLabels.buildNormCharacteristicIntensity(mcd, cxr)));
			mMcd = mcd;
			mCxr = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
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

			RealVector rv = new ArrayRealVector(1);
			double iPeak = iNorm[PEAK] - (iNorm[LOW] + ((iNorm[HIGH]-iNorm[LOW])/(r[HIGH]-r[LOW]))*(r[PEAK]-r[LOW]));
			rv.setEntry(0, iPeak);
			RealMatrix rm = MatrixUtils.createRealMatrix(1, 6);
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

	static public List<LabeledMultivariateJacobianFunction> buildModel( //
			MatrixCorrectionDatum mcd, //
			CharacteristicXRay cxr //
	) throws ArgumentException {
		final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
		LabeledMultivariateJacobianFunctionBuilder builder = new LabeledMultivariateJacobianFunctionBuilder(
				"Normalizer");
		builder.add(new NormalizeIntensity(mcd, cxr, LOW_BACK));
		builder.add(new NormalizeIntensity(mcd, cxr, ON_PEAK));
		builder.add(new NormalizeIntensity(mcd, cxr, HIGH_BACK));
		res.add(builder.build());
		res.add(new ComputePeakContinuum(mcd, cxr));
		return res;
	}

	/**
	 * @param mcd Details the material and conditions
	 * @param cxr The characteristic x-ray
	 * @throws ArgumentException
	 */
	public TwoPointContinuumModel(//
			MatrixCorrectionDatum mcd, //
			CharacteristicXRay cxr //
	) throws ArgumentException {
		super("TwoPoint", buildModel(mcd, cxr));
	}
}
