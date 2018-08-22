package gov.nist.microanalysis.roentgen.matrixcorrection;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.google.common.base.Objects;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class MeasurementDatum implements IToHTML {

	private static int mNextIndex = 0;

	private final int mIndex;
	private final Sample mSample;
	private Vector3D mSurfaceOrientation;
	private final double mNominalLiveTime;

	public final static Vector3D BEAM_NORMAL_TO_SURFACE = Vector3D.MINUS_K;

	private MeasurementDatum(Sample sample, Vector3D orient, double liveTime) {
		mIndex = (++mNextIndex);
		mSample = sample;
		mSurfaceOrientation = orient;
		mNominalLiveTime = liveTime;
	}

	@Override
	public int hashCode() {
		return Objects.hashCode(mSample, mSurfaceOrientation, mNominalLiveTime, mIndex);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final MeasurementDatum other = (MeasurementDatum) obj;
		return Objects.equal(mSample, other.mSample) && //
				Objects.equal(mSurfaceOrientation, other.mSurfaceOrientation) && //
				(mNominalLiveTime == other.mNominalLiveTime) && //
				(mIndex == other.mIndex);
	}

	/**
	 * Creates a {@link MeasurementDatum} object that is suitable for use as an
	 * unknown.
	 *
	 * @param liveTime
	 * @param probeCurrent
	 * @param beamEnergy
	 * @return MeasurementDatum
	 */
	static public MeasurementDatum unknown(final Sample samp, final double liveTime) {
		final MeasurementDatum res = new MeasurementDatum(samp, BEAM_NORMAL_TO_SURFACE, liveTime);
		assert res.isSuitableAsUnknown();
		return res;
	}

	/**
	 * Creates a {@link MeasurementDatum} object that is suitable for use as a
	 * standard.
	 *
	 * @param liveTime
	 * @param probeCurrent
	 * @param beamEnergy
	 * @param comp
	 * @return MeasurementDatum
	 */
	static public MeasurementDatum standard(Sample samp, Element elm, double liveTime) {
		final MeasurementDatum res = new MeasurementDatum(samp, BEAM_NORMAL_TO_SURFACE, liveTime);
		assert res.isSuitableAsStandard(elm);
		return res;
	}

	/**
	 * Modifies a {@link MeasurementDatum} object to specify an orientation of the
	 * sample relative to the beam.
	 *
	 * @param orient
	 * @return {@link MeasurementDatum}
	 */
	public MeasurementDatum withOrientation(final Vector3D normal) {
		mSurfaceOrientation = normal.normalize();
		return this;
	}

	/**
	 * Modifies a {@link MeasurementDatum} object to specify that the sample is
	 * oriented normal to the beam.
	 *
	 * @return {@link MeasurementDatum}
	 */
	public MeasurementDatum normalToBeam() {
		mSurfaceOrientation = BEAM_NORMAL_TO_SURFACE;
		return this;
	}

	/**
	 * Does the object contain all the data necessary to be used as an unknown?
	 *
	 * @return boolean
	 */
	public boolean isSuitableAsUnknown() {
		assert (mNominalLiveTime > 0.0);
		return true;
	}

	/**
	 * Does the object contain all the data necessary to be used as a standard?
	 *
	 * @return
	 */
	public boolean isSuitableAsStandard(Element elm) {
		return isSuitableAsUnknown() && mSample.isSuitableAsStandard(elm);
	}

	/**
	 * @return An {@link UncertainValue} in seconds
	 */
	public double getLiveTime() {
		return mNominalLiveTime;
	}

	/**
	 * Returns the orientation of the sample relative to the beam which is assumed
	 * to propagate down towards larger Z.
	 *
	 * @return Vector3D (normalized to unity length)
	 */
	public Vector3D getSurfaceNormal() {
		return mSurfaceOrientation;
	}

	/**
	 * The sample is oriented such that the electron beam is normal to the surface.
	 *
	 * @return boolean
	 */
	public boolean isNormal() {
		return mSurfaceOrientation.equals(BEAM_NORMAL_TO_SURFACE);
	}

	@Override
	public String toHTML(final Mode mode) {
		Table table = new Table();
		BasicNumberFormat bnf = new BasicNumberFormat("0.00");
		table.addRow(Table.td("Sample"), Table.td(mSample));
		table.addRow(Table.td("Orientation"), Table.td(HTML.toHTML(mSurfaceOrientation, mode)));
		table.addRow(Table.td("Live Time"), Table.td(bnf.formatHTML(mNominalLiveTime)));
		return table.toHTML(mode);
	}
}
