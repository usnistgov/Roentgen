package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Objects;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;

public class Layer //
		implements IToHTML {

	private final Composition mComposition;
	private final UncertainValue mThickness; // In cm
	private final UncertainValue mDensity; // In g/cm^3

	@Override
	public int hashCode() {
		return Objects.hash(mComposition, mDensity, mThickness);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Layer other = (Layer) obj;
		return Objects.equals(mComposition, other.mComposition) //
				&& Objects.equals(mDensity, other.mDensity) //
				&& Objects.equals(mThickness, other.mThickness);
	}

	/**
	 * @param comp
	 * @param thickness_nm In nanometers
	 * @param density      In g/cm<sup>3</sup>
	 */
	public Layer(final Composition comp, final UncertainValue thickness_nm, final UncertainValue density) {
		assert comp != null;
		assert (thickness_nm != null) && (thickness_nm.doubleValue() > 0.0);
		assert (density != null) && (density.doubleValue() > 0.0);
		mComposition = comp;
		mThickness = thickness_nm.multiply(1.0e-7); // 1 nm = 1.0e-7 cm
		mDensity = density;
	}

	/**
	 * Makes a Layer consisting of thickness_nm (in nanometers) of carbon with
	 * density of 1.5 &pm; 0.1 g/cm<sup>3</sup>.
	 * 
	 * @param thickness_nm
	 * @return Layer object
	 */
	public static Layer carbonCoating(UncertainValue thickness_nm) {
		return new Layer(Composition.pureElement(Element.Carbon), thickness_nm, new UncertainValue(1.5, 0.1));
	}

	/**
	 * Makes a Layer consisting of thickness_nm (in nanometers) of gold with density
	 * of 19.3 &pm; 0.1 g/cm<sup>3</sup>.
	 * 
	 * @param thickness_nm
	 * @return Layer object
	 */
	public static Layer goldCoating(UncertainValue thickness_nm) {
		return new Layer(Composition.pureElement(Element.Gold), thickness_nm, new UncertainValue(19.3, 0.1));
	}

	public Composition getComposition() {
		return mComposition;
	}
	
	public Material getMaterial() {
		return mComposition.getMaterial();
	}

	/**
	 * Thickness in meters
	 * 
	 * @return in meters
	 */
	public UncertainValue getThickness() {
		return mThickness;
	}

	/**
	 * @return in g/cm<sup>3</sup>
	 */
	public UncertainValue getDensity() {
		return mDensity;
	}

	public UncertainValue getMassThickness() {
		return new UncertainValue(mThickness.doubleValue() * mDensity.doubleValue(), //
				(mThickness.doubleValue() * mDensity.uncertainty()
						+ mDensity.doubleValue() * mThickness.uncertainty()));
	}

	@Override
	public String toHTML(final Mode mode) {
		return mComposition.toHTML(mode.demote()) + "[" + mThickness.toHTML(mode) + "cm," + mDensity.toHTML(mode)
				+ " g/cm<sup>3</sup>]";
	}

	@Override
	public String toString() {
		return mComposition.toString() + "[" + mThickness.toString() + "cm," + mDensity.toString() + " g/cm^3]";
	}

}
