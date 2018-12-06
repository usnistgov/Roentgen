package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

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
 * @author nicholas
 *
 */
public abstract class MatrixCorrectionDatum //
		implements IToHTML {

	// Either mComposition or mElements is defined (never both!)
	protected final UncertainValue mBeamEnergy;
	protected final UncertainValue mTakeOffAngle;
	protected final Optional<Double> mRoughness;

	@Override
	public int hashCode() {
		return Objects.hash(mBeamEnergy, mTakeOffAngle, mRoughness);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MatrixCorrectionDatum other = (MatrixCorrectionDatum) obj;
		return Objects.equals(mBeamEnergy, other.mBeamEnergy) && //
				Objects.equals(mTakeOffAngle, other.mTakeOffAngle) && //
				Objects.equals(mRoughness, other.mRoughness);
	}

	/**
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 */
	protected MatrixCorrectionDatum(//
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle //
	) {
		this(beamEnergy, takeOffAngle, Double.NaN);
	}

	/**
	 * @param comp
	 * @param isStandard
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 */
	protected MatrixCorrectionDatum(
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle, //
			double roughness //
	) {
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mRoughness = Double.isNaN(roughness) ? Optional.empty() : Optional.of(Double.valueOf(roughness));
	}

	/**
	 * Helper to simplify computing roughness.
	 * 
	 * @param dimension In nanometers
	 * @param density   g/cm<sup>3</sup>
	 * @return double cm (g/cm<sup>3</sup>)
	 */
	static public double roughness(double dimension, double density) {
		return 1.0e-7 * dimension * density;
	}

	public UncertainValue getBeamEnergy() {
		return mBeamEnergy;
	}

	public UncertainValue getTakeOffAngle() {
		return mTakeOffAngle;
	}

	public double getRoughness() {
		return mRoughness.orElse(0.0);
	}

	public boolean hasRoughness() {
		return mRoughness.isPresent();
	}
	
	/**
	 * Returns the Composition of the standard or the estimated Composition of
	 * the unknown.
	 * 
	 * @return {@link Composition}
	 */
	abstract Composition getComposition();
}
