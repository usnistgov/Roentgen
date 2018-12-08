package gov.nist.microanalysis.roentgen.spectrum;

import java.util.Collections;
import java.util.Date;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class EDSSpectrum {

	/**
	 * (Required) User friendly name
	 */
	private final String mName;
	/**
	 * (Required) Spectrum channel data
	 */
	private final RealVector mData;
	// Approximate gain and offset

	/**
	 * (Required) eV/channel gain
	 */
	private final double mGain;
	/**
	 * (Required) eV zero offset
	 */
	private final double mOffset;

	/**
	 * (Optional) Calibrated energy scale
	 */
	private EnergyCalibration mEnergy;
	/**
	 * (Optional) Calibrated lineshape
	 */
	private LineshapeCalibration mLineshape;

	/**
	 * (Optional) Date acquired
	 */
	private Date mAcquired;

	/**
	 * (Optional) Live time in seconds
	 */
	private double mLiveTime = Double.NaN;

	/**
	 * (Optional) Real time in seconds
	 */
	private double mRealTime = Double.NaN;

	/**
	 * (Optional) Beam energy in eV
	 */
	private double mBeamEnergy = Double.NaN;

	/**
	 * (Optional) Probe current in nA
	 */
	private double mProbeCurrent = Double.NaN;

	/**
	 * (Optional) Stage position
	 */
	private Map<String, Double> mPosition;

	/**
	 * (Optional) Set independently known composition of the material from which
	 * this spectrum was collected
	 */
	private Composition mComposition;

	/**
	 * Creates an EDSSpectrum object with the specified
	 *
	 * @param name         User friendly
	 * @param eVperCh      Channe width in eV
	 * @param zeroOffset   Zero offset in eV
	 * @param data         Channel data
	 * @param unmodifiable Whether the channel data should be modifiable
	 */
	public EDSSpectrum(final String name, final double eVperCh, final double zeroOffset, final double[] data,
			final boolean unmodifiable) {
		if (unmodifiable)
			mData = RealVector.unmodifiableRealVector(new ArrayRealVector(data));
		else
			mData = new ArrayRealVector(data);
		mName = name;
		mGain = eVperCh;
		mOffset = zeroOffset;
		mEnergy = EnergyCalibration.Linear(mOffset, mGain, data.length);

	}

	/**
	 * The energy associated with the low energy side of the channel.
	 *
	 * @param ch
	 * @return double in eV
	 */
	public double energy(final int ch) {
		return ch * mGain + mOffset;
	}

	public int channel(final double e) {
		return (int) Math.floor((e - mOffset) / mGain);
	}

	public RealVector getData() {
		return mData;
	}

	@Override
	public String toString() {
		return mName;
	}

	public int size() {
		return mData.getDimension();
	}

	/**
	 * (Optional) Get the {@link EnergyCalibration}.
	 */
	public EnergyCalibration getEnergyCalibration() {
		return mEnergy;
	}

	/**
	 * (Optional) Set the {@link EnergyCalibration}.
	 *
	 * @param energy {@link EnergyCalibration}
	 */
	public void setEnergyCalibration(final EnergyCalibration energy) {
		mEnergy = energy;
	}

	/**
	 * (Optional) Get the {@link LineshapeCalibration}
	 *
	 * @return {@link LineshapeCalibration}
	 */
	public LineshapeCalibration getLineshapeCalibration() {
		return mLineshape;
	}

	/**
	 * (Optional) Set the {@link LineshapeCalibration}
	 *
	 * @param lineshape
	 */
	public void setLineshapeCalibration(final LineshapeCalibration lineshape) {
		mLineshape = lineshape;
	}

	/**
	 * (Optional) Time the spectrum was acquired.
	 *
	 * @return
	 */
	public Date getAcquired() {
		return mAcquired;
	}

	/**
	 * (Optional) Time the spectrum was acquired.
	 *
	 * @param acquired
	 */
	public void setAcquired(final Date acquired) {
		mAcquired = acquired;
	}

	/**
	 * (Optional) Live time in seconds
	 *
	 * @return double seconds
	 */
	public double getLiveTime() {
		return mLiveTime;
	}

	/**
	 * (Optional) Live time in seconds
	 *
	 * @param liveTime seconds
	 */
	public void setLiveTime(final double liveTime) {
		mLiveTime = liveTime;
	}

	/**
	 * (Optional) Real time in seconds
	 *
	 * @return seconds
	 */
	public double getRealTime() {
		return mRealTime;
	}

	/**
	 * (Optional) Real time in seconds
	 *
	 * @param realTime seconds
	 */
	public void setRealTime(final double realTime) {
		mRealTime = realTime;
	}

	/**
	 * (Optional) Beam energy in eV
	 *
	 * @return double in eV
	 */
	public double getBeamEnergy() {
		return mBeamEnergy;
	}

	/**
	 * (Optional) Set the configured beam energy in eV
	 *
	 * @param beamEnergy eV
	 */
	public void setBeamEnergy(final double beamEnergy) {
		mBeamEnergy = beamEnergy;
	}

	/**
	 * (Optional) Probe current (Faraday current) in nA
	 *
	 * @return double in nA
	 */
	public double getProbeCurrent() {
		return mProbeCurrent;
	}

	/**
	 * (Optional) Probe current (Faraday current) in nA
	 *
	 * @param probeCurrent in nA
	 */
	public void setProbeCurrent(final double probeCurrent) {
		mProbeCurrent = probeCurrent;
	}

	public Map<String, Double> getPosition() {
		return mPosition;
	}

	static final public String X_AXIS = "X";
	static final public String Y_AXIS = "Y";
	static final public String Z_AXIS = "Z";
	static final public String R_AXIS = "R";
	static final public String X_TILT = "X Tilt";
	static final public String Y_TILT = "Y Tilt";

	public void setPosition(final Map<String, Double> position) {
		mPosition = Collections.unmodifiableMap(position);
	}

	/**
	 * Gets the current value assigned to composition
	 *
	 * @return Returns the composition.
	 */
	public Composition getComposition() {
		return mComposition;
	}

	/**
	 * Sets the value assigned to composition.
	 *
	 * @param composition The value to which to set composition.
	 */
	public void setComposition(final Composition composition) {
		mComposition = composition;
	}
}
