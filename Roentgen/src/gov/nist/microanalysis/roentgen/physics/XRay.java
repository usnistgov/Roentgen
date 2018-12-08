package gov.nist.microanalysis.roentgen.physics;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import com.duckandcover.html.IToHTML;
import com.google.common.base.Objects;
import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * Represents a generic x-ray. The base class for {@link CharacteristicXray}.
 * Comparison (sorting order) of x-rays is only determined by the energy (not
 * direction, source or intensity)
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 253 $
 */
public class XRay implements IToHTML, Comparable<XRay> {

	public static final BasicNumberFormat ENERGY_FORMAT = new BasicNumberFormat("0,000.0");

	/**
	 * The photon energy in eV
	 */
	private final double mEnergy;

	/**
	 * Constructs a XRay of the specified energy (in eV) and the specified
	 * intensity.
	 *
	 * @param energy in eV
	 */
	public XRay(final double energy) {
		Preconditions.checkArgument(energy > 0, "XRay energy must be greater than zero. " + energy);
		mEnergy = energy;
	}

	/**
	 * Used to compute the angular dependence of the emission profile. This is
	 * critical for Bremsstrahlung for which emission is focused in the direction of
	 * the electron trajectory.
	 *
	 * @param photonDirection
	 * @return double
	 */
	public double angularDependence(final Vector3D photonDirection) {
		return 1.0;
	}

	/**
	 * The x-ray energy in eV
	 *
	 * @return double in eV
	 */
	public double getEnergy() {
		return mEnergy;
	}

	/**
	 * Returns the x-ray wavelength in picometer.
	 *
	 * @return double (in pm.)
	 */
	public double getWavelength() {
		return 1.0e12 * (PhysicalConstants.PlanckConstant * PhysicalConstants.SpeedOfLight)
				/ (getEnergy() * PhysicalConstants.ElectronCharge);
	}

	@Override
	public String toString() {
		final StringBuffer sb = new StringBuffer();
		sb.append("x-ray(");
		sb.append(mEnergy);
		sb.append(" eV)");
		return sb.toString();
	}

	/**
	 * Hash code generated from energy (but not intensity.)
	 *
	 * @return int
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return super.hashCode() ^ Objects.hashCode(mEnergy);
	}

	/**
	 * Equals checks energy only. (Does not check intensity.)
	 *
	 * @param obj
	 * @return boolean
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final XRay other = (XRay) obj;
		return (Double.doubleToLongBits(mEnergy) == Double.doubleToLongBits(other.mEnergy));
	}

	/**
	 * X-rays are sorted by energy only!!!! This is regardless of their superclass.
	 *
	 * @param o
	 * @return int
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(final XRay o) {
		return Double.compare(mEnergy, o.mEnergy);
	}

	/**
	 * @param mode
	 * @return
	 * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
	 */
	@Override
	public String toHTML(final Mode mode) {
		return ENERGY_FORMAT.formatHTML(mEnergy) + " eV X-ray";
	}
}
