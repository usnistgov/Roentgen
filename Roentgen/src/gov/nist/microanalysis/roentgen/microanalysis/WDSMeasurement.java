package gov.nist.microanalysis.roentgen.microanalysis;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
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
public class WDSMeasurement extends NamedMultivariateJacobianFunction {

	enum Model {
		TwoPoint
	};

	protected static class WDSIntensity extends BaseTag {

		protected WDSIntensity(String name, double lPos, int rep) {
			super(name, Double.valueOf(lPos), Integer.valueOf(rep));
		}

		public double getL() {
			return ((Double) getObject1()).doubleValue();
		}

		public int getRepeatition() {
			return ((Integer) getObject2()).intValue();
		}
	}

	protected static class WDSLiveTime extends BaseTag {

		protected WDSLiveTime(String name, double lPos, int rep) {
			super(name, Double.valueOf(lPos), Integer.valueOf(rep));
		}

		public double getL() {
			return ((Double) getObject1()).doubleValue();
		}

		public int getRepeatition() {
			return ((Integer) getObject2()).intValue();
		}
	}

	public static class OnPeakIntensity extends WDSIntensity {

		private final CharacteristicXRay mXRay;
		private final List<BackgroundPacket> mBackgrounds;

		public OnPeakIntensity(CharacteristicXRay cxr, double pos, int rep, List<BackgroundPacket> backgrounds) {
			super("I<sub>" + cxr.toHTML(Mode.TERSE) + "</sub>", pos, rep);
			mBackgrounds = new ArrayList<>(backgrounds);
			mXRay = cxr;
		}
	}

	private final static BasicNumberFormat sFormat = new BasicNumberFormat("0.00");

	public static class BackgroundWDSIntensity extends WDSIntensity {

		public BackgroundWDSIntensity(double lPos, int rep) {
			super("I<sub>" + sFormat.format(lPos) + "," + rep + "</sub>", Double.valueOf(lPos), Integer.valueOf(rep));
		}
	}

	public static class OnPeakLiveTime extends WDSLiveTime {

		private final CharacteristicXRay mCharacteristic;

		public OnPeakLiveTime(CharacteristicXRay cxr, double lPos, int rep) {
			super("LT<sub>" + cxr.toHTML(Mode.TERSE) + "</sub>", lPos, rep);
			mCharacteristic = cxr;
		}

		public CharacteristicXRay getXRay() {
			return mCharacteristic;
		}
	}

	public static class BackgroundLiveTime extends WDSLiveTime {

		public BackgroundLiveTime(double lPos, int rep) {
			super("LT<sub>" + sFormat.format(lPos) + "," + rep + "</sub>", Double.valueOf(lPos), Integer.valueOf(rep));
		}
	}

	/**
	 * @author Nicholas
	 *
	 */
	public static class OnPeakPacket {

		private final OnPeakIntensity mOnPeakI;
		private final OnPeakLiveTime mOnPeakLT;
		private final Model mModel;
		private final CharacteristicXRay mXRay;

		OnPeakPacket(Model model, CharacteristicXRay cxr, double lPos, int rep, List<BackgroundPacket> backs) {
			assert model == Model.TwoPoint;
			assert backs.size() == 2;
			mModel = model;
			mXRay = cxr;
			mOnPeakI = new OnPeakIntensity(cxr, lPos, rep, backs);
			mOnPeakLT = new OnPeakLiveTime(cxr, lPos, rep);
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

		private final BackgroundWDSIntensity mBackI;
		private final BackgroundLiveTime mBackLT;

		public BackgroundPacket(double lPos, int rep) {
			mBackI = new BackgroundWDSIntensity(lPos, rep);
			mBackLT = new BackgroundLiveTime(lPos, rep);
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
	public static class CorrectedIntensity extends BaseTag {

		public CorrectedIntensity(CharacteristicXRay cxr, int rep) {
			super("I<sub>BC," + cxr.toHTML(Mode.TERSE) + "</sub>", cxr, Integer.valueOf(rep));
		}
	}

	protected final ArrayList<OnPeakPacket> mMeasurements;

	private static List<BaseTag> inputTags(List<OnPeakPacket> dataTags) {
		HashSet<BackgroundPacket> backs = new HashSet<>();
		ArrayList<BaseTag> res = new ArrayList<>();
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

	private static List<BaseTag> outputTags(List<OnPeakPacket> dataTags) {
		ArrayList<BaseTag> res = new ArrayList<>();
		for (OnPeakPacket pack : dataTags)
			res.add(new CorrectedIntensity(pack.mOnPeakI.mXRay, pack.mOnPeakI.getRepeatition()));
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
		for(OnPeakPacket opp : mMeasurements) {
			final OnPeakIntensity opi = opp.mOnPeakI;
			final OnPeakLiveTime opl = opp.mOnPeakLT;
			final List<BackgroundPacket> backs= opi.mBackgrounds;
			final Model model = opp.mModel;
			if(model==Model.TwoPoint) {
				int ii_opi = inputIndex(opi);
				int ii_oplt = inputIndex(opl);
				final double vopi = point.getEntry(ii_opi);
				final double voplt = point.getEntry(ii_oplt);
				final double vopL = opi.getL();
				final BackgroundPacket bp0 = backs.get(0);
				final BackgroundPacket bp1 = backs.get(1);
				int ii_b0i = inputIndex(bp0.mBackI);
				int ii_b0lt = inputIndex(bp0.mBackLT);
				final double vb0L = bp0.mBackI.getL();
				final double vb0i = point.getEntry(ii_b0i);
				final double vb0lt = point.getEntry(ii_b0lt);
				final double vb1L = bp1.mBackI.getL();
				int ii_b1i = inputIndex(bp1.mBackI);
				int ii_b1lt = inputIndex(bp1.mBackLT);
				final double vb1i = point.getEntry(ii_b1i);
				final double vb1lt = point.getEntry(ii_b1lt);
				final double backEst = (vb0i/vb0lt) + (((vb1i/vb1lt)-(vb0i/vb0lt))/(vb1L-vb0L))*(vopL-vb0L);
				CorrectedIntensity ci = new CorrectedIntensity(opp.mXRay, opi.getRepeatition());
				final double i = (vopi/voplt) - backEst;
				final int row = outputIndex(ci);
				rvres.setEntry(row, i);
				rmres.setEntry(row, ii_opi, Double.NaN);
				rmres.setEntry(row, ii_oplt, Double.NaN);
				rmres.setEntry(row, ii_b0i, Double.NaN);
				rmres.setEntry(row, ii_b0lt, Double.NaN);
				rmres.setEntry(row, ii_b1i, Double.NaN);
				rmres.setEntry(row, ii_b1lt, Double.NaN);
			} else
				assert false;
			
			
			
			
			
		}
		return null;
	}
}
