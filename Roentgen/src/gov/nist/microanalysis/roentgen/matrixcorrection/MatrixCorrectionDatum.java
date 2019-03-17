package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Objects;
import java.util.Optional;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Material;

/**
 * A package describing a set of arguments to the matrix correction algorithm.
 *
 * {@link MatrixCorrectionDatum} come in three different flavors:
 * <ol>
 * <li>Standards - with a defined Composition and isStandard()==true</li>
 * <li>Unknown estimates - with a defined Composition and
 * isStandard()==false</li>
 * <li>True unknowns - with only an Element set</li>
 * </ol>
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public abstract class MatrixCorrectionDatum //
		implements IToHTML {

	// Either mComposition or mElements is defined (never both!)
	protected final UncertainValue mBeamEnergy;
	protected final UncertainValue mTakeOffAngle;
	protected final Optional<Double> mRoughness;
	protected final Optional<Layer> mCoating;

	@Override
	public int hashCode() {
		return Objects.hash(mBeamEnergy, mTakeOffAngle, mRoughness, mCoating);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final MatrixCorrectionDatum other = (MatrixCorrectionDatum) obj;
		return Objects.equals(mBeamEnergy, other.mBeamEnergy) && //
				Objects.equals(mTakeOffAngle, other.mTakeOffAngle) && //
				Objects.equals(mRoughness, other.mRoughness) && //
				Objects.equals(mCoating, other.mCoating);
	}

	/**
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 */
	protected MatrixCorrectionDatum(//
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle //
	) {
		this(beamEnergy, takeOffAngle, Double.NaN, null);
	}

	/**
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 */
	protected MatrixCorrectionDatum(final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness //
	) {
		this(beamEnergy, takeOffAngle, roughness, null);
	}

	/**
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 * @param coating      A Layer object
	 */
	protected MatrixCorrectionDatum(final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness, final Layer coating//
	) {
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mRoughness = Double.isNaN(roughness) || (roughness == 0.0) ? Optional.empty()
				: Optional.of(Double.valueOf(roughness));
		mCoating = Optional.ofNullable(coating);
	}

	/**
	 * Helper to simplify computing roughness.
	 *
	 * @param dimension In nanometers
	 * @param density   g/cm<sup>3</sup>
	 * @return double cm (g/cm<sup>3</sup>)
	 */
	static public double roughness(final double dimension, final double density) {
		return 1.0e-7 * dimension * density;
	}

	public UncertainValue getBeamEnergy() {
		return mBeamEnergy;
	}

	public UncertainValue getTakeOffAngle() {
		return mTakeOffAngle;
	}

	/**
	 * Returns the max(roughness, 1.0 nm g/cm<sup>3</cm>)
	 *
	 * @return double
	 */
	public double getRoughness() {
		return mRoughness.orElse(roughness(1.0, 1.0));
	}

	public boolean hasRoughness() {
		return mRoughness.isPresent();
	}

	public boolean hasCoating() {
		return mCoating.isPresent();
	}

	public Layer getCoating() {
		return mCoating.orElse(null);
	}

	public abstract Material getMaterial();
}
