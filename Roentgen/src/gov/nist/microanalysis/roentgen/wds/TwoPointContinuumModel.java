package gov.nist.microanalysis.roentgen.wds;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ParallelMeasurementModelBuilder;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.LiveTime;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.NormalizedCharacteristicIntensity;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.NormalizedIntensity;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.ProbeCurrent;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.RawIntensity;
import gov.nist.microanalysis.roentgen.wds.ModelLabels.SpectrometerPosition;

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
		extends CompositeMeasurementModel<EPMALabel> {

	public static class ComputePeakContinuum
			extends ExplicitMeasurementModel<EPMALabel, NormalizedCharacteristicIntensity> {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;

		static private List<EPMALabel> buildInputs(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, LOW_BACK));
			res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, ON_PEAK));
			res.add(ModelLabels.buildNormalizedIntensity(mcd, cxr, HIGH_BACK));
			res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, LOW_BACK));
			res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, ON_PEAK));
			res.add(ModelLabels.buildSpectrometerPosition(mcd, cxr, HIGH_BACK));
			return res;
		}

		public ComputePeakContinuum(
				final MatrixCorrectionDatum mcd, //
				final CharacteristicXRay cxr
		) throws ArgumentException {
			super(buildInputs(mcd, cxr),
					Collections.singletonList(ModelLabels.buildNormCharacteristicIntensity(mcd, cxr)));
			mMcd = mcd;
			mCxr = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final int LOW = 0, PEAK = 1, HIGH = 2;
			final NormalizedIntensity[] normI = { //
					ModelLabels.buildNormalizedIntensity(mMcd, mCxr, LOW_BACK), //
					ModelLabels.buildNormalizedIntensity(mMcd, mCxr, ON_PEAK), //
					ModelLabels.buildNormalizedIntensity(mMcd, mCxr, HIGH_BACK) //
			};
			final SpectrometerPosition[] rI = { //
					ModelLabels.buildSpectrometerPosition(mMcd, mCxr, LOW_BACK), //
					ModelLabels.buildSpectrometerPosition(mMcd, mCxr, ON_PEAK), //
					ModelLabels.buildSpectrometerPosition(mMcd, mCxr, HIGH_BACK) //
			};
			final double[] iNorm = { //
					getArg(normI[LOW], point), //
					getArg(normI[PEAK], point), //
					getArg(normI[HIGH], point) //
			};
			final double[] r = { //
					getArg(rI[LOW], point), //
					getArg(rI[PEAK], point), //
					getArg(rI[HIGH], point) //
			};

			final RealVector rv = buildResult();
			final RealMatrix rm = buildJacobian();
			final NormalizedCharacteristicIntensity nci = ModelLabels.buildNormCharacteristicIntensity(mMcd, mCxr);
			final double iPeak = iNorm[PEAK]
					- (iNorm[LOW] + (iNorm[HIGH] - iNorm[LOW]) * (r[PEAK] - r[LOW]) / (r[HIGH] - r[LOW]));
			setResult(nci, rv, iPeak);
			// WRT normalized intensity
			setJacobian(normI[LOW], nci, rm, (r[PEAK] - r[HIGH]) / (r[HIGH] - r[LOW])); // OK
			setJacobian(normI[PEAK], nci, rm, 1.0); // OK
			setJacobian(normI[HIGH], nci, rm, (r[PEAK] - r[LOW]) / (r[HIGH] - r[LOW])); // OK
			// WRT spectrometer position
			setJacobian(rI[LOW], nci, rm,
					(iNorm[HIGH] - iNorm[LOW]) * (r[HIGH] - r[PEAK]) / Math.pow(r[HIGH] - r[LOW], 2.0)); // OK
			setJacobian(rI[PEAK], nci, rm, (iNorm[LOW] - iNorm[HIGH]) / (r[HIGH] - r[LOW])); // OK
			setJacobian(rI[HIGH], nci, rm,
					(iNorm[LOW] - iNorm[HIGH]) * (r[LOW] - r[PEAK]) / Math.pow(r[HIGH] - r[LOW], 2.0)); // OK
			return Pair.create(rv, rm);
		}

	}

	private static class NormalizeIntensity //
			extends ExplicitMeasurementModel<EPMALabel, NormalizedIntensity> {

		private final MatrixCorrectionDatum mMcd;
		private final CharacteristicXRay mCxr;
		private final int mIndex;

		static private List<EPMALabel> buildInputs(
				final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr, final int index
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(ModelLabels.buildRawIntensity(mcd, cxr, index));
			res.add(ModelLabels.buildLiveTime(mcd, cxr, index));
			res.add(ModelLabels.buildProbeCurrent(mcd, cxr, index));
			return res;
		}

		private NormalizeIntensity(
				final MatrixCorrectionDatum mcd, //
				final CharacteristicXRay cxr, //
				final int index
		) throws ArgumentException {
			super(buildInputs(mcd, cxr, index),
					Collections.singletonList(ModelLabels.buildNormalizedIntensity(mcd, cxr, index)));
			mMcd = mcd;
			mCxr = cxr;
			mIndex = index;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final RawIntensity iL = ModelLabels.buildRawIntensity(mMcd, mCxr, mIndex);
			final LiveTime ltL = ModelLabels.buildLiveTime(mMcd, mCxr, mIndex);
			final ProbeCurrent pcL = ModelLabels.buildProbeCurrent(mMcd, mCxr, mIndex);

			final double i = getArg(iL, point);
			final double lt = getArg(ltL, point);
			final double pc = getArg(pcL, point);

			final RealVector rv = buildResult();
			final RealMatrix rm = buildJacobian();
			final NormalizedIntensity ni = ModelLabels.buildNormalizedIntensity(mMcd, mCxr, mIndex);
			final double val = i / (lt * pc);
			setResult(ni, rv, val);

			setJacobian(iL, ni, rm, val / i);
			setJacobian(ltL, ni, rm, -val / lt);
			setJacobian(pcL, ni, rm, -val / pc);

			return Pair.create(rv, rm);
		}
	}

	public static class ComputeKRatio //
			extends ExplicitMeasurementModel<NormalizedCharacteristicIntensity, KRatioLabel> {

		private final UnknownMatrixCorrectionDatum mUnknown;
		private final StandardMatrixCorrectionDatum mStandard;
		private final CharacteristicXRay mXRay;

		private static List<NormalizedCharacteristicIntensity> buildInputs(
				final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final CharacteristicXRay cxr
		) {
			final List<NormalizedCharacteristicIntensity> res = new ArrayList<>();
			res.add(ModelLabels.buildNormCharacteristicIntensity(unk, cxr));
			res.add(ModelLabels.buildNormCharacteristicIntensity(std, cxr));
			return res;
		}

		public ComputeKRatio(
				final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final CharacteristicXRay cxr
		) throws ArgumentException {
			super(buildInputs(unk, std, cxr),
					Collections.singletonList(new KRatioLabel(unk, std, cxr, Method.Measured)));
			mUnknown = unk;
			mStandard = std;
			mXRay = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final NormalizedCharacteristicIntensity nciu = //
					ModelLabels.buildNormCharacteristicIntensity(mUnknown, mXRay);
			final NormalizedCharacteristicIntensity ncis = //
					ModelLabels.buildNormCharacteristicIntensity(mStandard, mXRay);

			final KRatioLabel krl = new KRatioLabel(mUnknown, mStandard, mXRay, Method.Measured);
			final double iu = getArg(nciu, point);
			final double is = getArg(ncis, point);

			final RealVector val = buildResult();
			final RealMatrix jac = buildJacobian();

			setResult(krl, val, iu / is);
			setJacobian(nciu, krl, jac, 1 / is);
			setJacobian(ncis, krl, jac, -iu / (is * is));

			return Pair.create(val, jac);
		}

	}

	public static final int HIGH_BACK = 1;

	public static final int ON_PEAK = 0;

	public static final int LOW_BACK = -1;

	static public List<? extends ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildModel(
			final UnknownMatrixCorrectionDatum unk, //
			final StandardMatrixCorrectionDatum std, //
			final CharacteristicXRay cxr //
	) throws ArgumentException {
		final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
		{
			final ParallelMeasurementModelBuilder<EPMALabel, NormalizedIntensity> builder = //
					new ParallelMeasurementModelBuilder<EPMALabel, NormalizedIntensity>("Normalizer");
			builder.add(new NormalizeIntensity(unk, cxr, LOW_BACK));
			builder.add(new NormalizeIntensity(unk, cxr, ON_PEAK));
			builder.add(new NormalizeIntensity(unk, cxr, HIGH_BACK));
			builder.add(new NormalizeIntensity(std, cxr, LOW_BACK));
			builder.add(new NormalizeIntensity(std, cxr, ON_PEAK));
			builder.add(new NormalizeIntensity(std, cxr, HIGH_BACK));
			res.add(builder.build());
		}
		{
			final ParallelMeasurementModelBuilder<EPMALabel, NormalizedCharacteristicIntensity> builder = //
					new ParallelMeasurementModelBuilder<EPMALabel, NormalizedCharacteristicIntensity>(
							"NormCharacteristicIntensity");
			builder.add(new ComputePeakContinuum(unk, cxr));
			builder.add(new ComputePeakContinuum(std, cxr));
			res.add(builder.build());
		}
		res.add(new ComputeKRatio(unk, std, cxr));
		return res;
	}

	/**
	 * @param mcd Details the material and conditions
	 * @param cxr The characteristic x-ray
	 * @throws ArgumentException
	 */
	public TwoPointContinuumModel(
			final UnknownMatrixCorrectionDatum unk, //
			final StandardMatrixCorrectionDatum std, //
			final CharacteristicXRay cxr //
	) throws ArgumentException {
		super("TwoPoint", buildModel(unk, std, cxr));
	}
}
