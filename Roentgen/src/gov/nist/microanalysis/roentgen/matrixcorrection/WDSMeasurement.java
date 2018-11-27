package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * This class implements an abstract class that associates WDS data for on-peak
 * measurements with background measurements. It allows but does not require
 * background measurements to be shared between on-peak measurements. It is
 * designed to serve as the basis for either classic peak + low + high
 * measurements or shared multipoint backgrounds. It also can be used
 * 
 * simple two-point background corrected intensity measurement. The data points
 * are the pairs of measured intensity and live-times at specific spectrometer
 * positions. One position represents the on-peak position and other two
 * represent the background positions.
 * 
 * @author Nicholas
 *
 */
public class WDSMeasurement extends LabeledMultivariateJacobianFunction {

	enum Model {
		TwoPoint
	};

	protected static class WDSIntensityTag extends BaseLabel<Double,Integer,Object> {

		protected WDSIntensityTag(String name, double lPos, int rep) {
			super(name, Double.valueOf(lPos), Integer.valueOf(rep));
		}

		public double getL() {
			return ((Double) getObject1()).doubleValue();
		}

		public int getRepetition() {
			return ((Integer) getObject2()).intValue();
		}
	}

	protected static class WDSLiveTimeTag extends BaseLabel<Double,Integer,Object> {

		protected WDSLiveTimeTag(String name, double lPos, int rep) {
			super(name, Double.valueOf(lPos), Integer.valueOf(rep));
		}

		public double getL() {
			return ((Double) getObject1()).doubleValue();
		}

		public int getRepetition() {
			return ((Integer) getObject2()).intValue();
		}
	}

	public static class OnPeakIntensityTag extends WDSIntensityTag {

		private final CharacteristicXRay mXRay;
		private final List<BackgroundPacket> mBackgrounds;

		public OnPeakIntensityTag(CharacteristicXRay cxr, double pos, int rep, List<BackgroundPacket> backgrounds) {
			super("I<sub>" + cxr.toHTML(Mode.TERSE) + "</sub>", pos, rep);
			mBackgrounds = new ArrayList<>(backgrounds);
			mXRay = cxr;
		}
	}

	private final static BasicNumberFormat sFormat = new BasicNumberFormat("0.00");

	public static class BackgroundWDSIntensityTag extends WDSIntensityTag {

		public BackgroundWDSIntensityTag(double lPos, int rep) {
			super("I<sub>" + sFormat.format(lPos) + "," + rep + "</sub>", Double.valueOf(lPos), Integer.valueOf(rep));
		}
	}

	public static class OnPeakLiveTimeTag extends WDSLiveTimeTag {

		private final CharacteristicXRay mCharacteristic;

		public OnPeakLiveTimeTag(CharacteristicXRay cxr, double lPos, int rep) {
			super("LT<sub>" + cxr.toHTML(Mode.TERSE) + "</sub>", lPos, rep);
			mCharacteristic = cxr;
		}

		public CharacteristicXRay getXRay() {
			return mCharacteristic;
		}
	}

	public static class BackgroundLiveTimeTag extends WDSLiveTimeTag {

		public BackgroundLiveTimeTag(double lPos, int rep) {
			super("LT<sub>" + sFormat.format(lPos) + "," + rep + "</sub>", Double.valueOf(lPos), Integer.valueOf(rep));
		}
	}

	/**
	 * @author Nicholas
	 *
	 */
	public static class OnPeakPacket {

		private final OnPeakIntensityTag mOnPeakI;
		private final OnPeakLiveTimeTag mOnPeakLT;
		private final ProbeCurrentLabel mProbeCurrent;
		private final Model mModel;
		private final CharacteristicXRay mXRay;

		OnPeakPacket(Model model, CharacteristicXRay cxr, double lPos, int rep, ProbeCurrentLabel pcTag,
				List<BackgroundPacket> backs) {
			assert model == Model.TwoPoint;
			assert backs.size() == 2;
			mModel = model;
			mXRay = cxr;
			mOnPeakI = new OnPeakIntensityTag(cxr, lPos, rep, backs);
			mOnPeakLT = new OnPeakLiveTimeTag(cxr, lPos, rep);
			mProbeCurrent = pcTag;
		}
	}

	public static class BackgroundPacket {

		@Override
		public int hashCode() {
			return Objects.hash(mBackI, mBackLT);
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			BackgroundPacket other = (BackgroundPacket) obj;
			return Objects.equals(mBackI, other.mBackI) && //
					Objects.equals(mBackLT, other.mBackLT);
		}

		private final BackgroundWDSIntensityTag mBackI;
		private final BackgroundLiveTimeTag mBackLT;

		public BackgroundPacket(double lPos, int rep) {
			mBackI = new BackgroundWDSIntensityTag(lPos, rep);
			mBackLT = new BackgroundLiveTimeTag(lPos, rep);
		}
	}

	/**
	 * This tag represents a single background corrected and live-time scaled
	 * mesurement of x-ray intensity.
	 * 
	 * 
	 * @author Nicholas
	 *
	 */
	public static class CorrectedIntensity extends BaseLabel<CharacteristicXRay,Integer,Object> {

		public CorrectedIntensity(CharacteristicXRay cxr, int rep) {
			super("I<sub>BC," + cxr.toHTML(Mode.TERSE) + "</sub>", cxr, Integer.valueOf(rep));
		}
	}

	protected final ArrayList<OnPeakPacket> mMeasurements;

	private static List<BaseLabel<?,?,?>> inputTags(List<OnPeakPacket> dataTags) {
		HashSet<BackgroundPacket> backs = new HashSet<>();
		ArrayList<BaseLabel<?,?,?>> res = new ArrayList<>();
		for (OnPeakPacket pack : dataTags) {
			res.add(pack.mOnPeakI);
			res.add(pack.mOnPeakLT);
			for (BackgroundPacket back : pack.mOnPeakI.mBackgrounds)
				if (!backs.contains(back)) {
					res.add(back.mBackI);
					res.add(back.mBackLT);
					backs.add(back);
				}
		}
		return res;
	}

	private static List<BaseLabel<?,?,?>> outputTags(List<OnPeakPacket> dataTags) {
		ArrayList<BaseLabel<?,?,?>> res = new ArrayList<>();
		for (OnPeakPacket pack : dataTags)
			res.add(new CorrectedIntensity(pack.mOnPeakI.mXRay, pack.mOnPeakI.getRepetition()));
		return res;
	}

	public WDSMeasurement(ArrayList<OnPeakPacket> dataTags) {
		super(inputTags(dataTags), outputTags(dataTags));
		mMeasurements = new ArrayList<>(dataTags);
	}
	
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		final RealVector rvres = new ArrayRealVector(getOutputDimension());
		final RealMatrix rmres = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		for (OnPeakPacket opp : mMeasurements) {
			final OnPeakIntensityTag opi = opp.mOnPeakI;
			final OnPeakLiveTimeTag opl = opp.mOnPeakLT;
			final ProbeCurrentLabel oppc = opp.mProbeCurrent;
			final List<BackgroundPacket> backs = opi.mBackgrounds;
			final Model model = opp.mModel;
			if (model == Model.TwoPoint) {
				/*
				 * Assumptions: 1) Probe current same for on-peak and background measurements 2)
				 * Spectrometer position without error
				 */
				// Two background points
				final BackgroundPacket bp0 = backs.get(0);
				final BackgroundPacket bp1 = backs.get(1);
				// Input values
				final int ii_opi = inputIndex(opi);
				final int ii_oplt = inputIndex(opl);
				final int ii_oppc = inputIndex(oppc);
				final int ii_b0i = inputIndex(bp0.mBackI);
				final int ii_b0lt = inputIndex(bp0.mBackLT);
				final int ii_b1i = inputIndex(bp1.mBackI);
				final int ii_b1lt = inputIndex(bp1.mBackLT);
				// Measured values (on peak)
				final double vopi = point.getEntry(ii_opi);
				final double voplt = point.getEntry(ii_oplt);
				final double voppc = point.getEntry(ii_oppc);
				final double vopL = opi.getL();
				// Measured values (background 0)
				final double vb0i = point.getEntry(ii_b0i);
				final double vb0lt = point.getEntry(ii_b0lt);
				final double vb0L = bp0.mBackI.getL();
				// Measured values (background 1)
				final double vb1i = point.getEntry(ii_b1i);
				final double vb1lt = point.getEntry(ii_b1lt);
				final double vb1L = bp1.mBackI.getL();
				// Dose normalized intensities
				final double vinp = vopi / (voplt * voppc);
				final double vinb0 = vb0i / (vb0lt * voppc);
				final double vinb1 = vb1i / (vb1lt * voppc);

				final double kk = (vopL - vb0L) / (vb1L - vb0L);
				final double bci = vinp - kk * (vinb1 - vinb0);
				CorrectedIntensity ci = new CorrectedIntensity(opp.mXRay, opi.getRepetition());
				final int row = outputIndex(ci);
				// Background corrected intensity value
				rvres.setEntry(row, bci);
				// Partial derivatives
				rmres.setEntry(row, ii_opi, 1.0 / (voplt * voppc));
				rmres.setEntry(row, ii_oplt, -vopi / (voplt * voplt * voppc));
				rmres.setEntry(row, ii_oppc, -vopi / (voplt * voppc * voppc));
				rmres.setEntry(row, ii_b0i, kk / (vb0lt * voppc));
				rmres.setEntry(row, ii_b0lt, kk * (vb0i / voppc));
				rmres.setEntry(row, ii_b1i, -kk / (vb1lt * voppc));
				rmres.setEntry(row, ii_b1lt, -kk * (vb1i / voppc));
			} else
				assert false;

		}
		return null;
	}
}
