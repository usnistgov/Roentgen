package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Objects;

import com.duckandcover.html.IToHTML;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.DataStore.UniqueString;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

public class Layer //
		implements IToHTML {

	private final UniqueString mName;
	private final Composition mComposition;
	private final UncertainValue mThickness; // In cm
	private final UncertainValue mDensity; // In g/cm^3

	@Override
	public int hashCode() {
		return Objects.hash(mComposition, mDensity, mThickness, mName);
	}

	@Override
	public boolean equals(
			final Object obj
	) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Layer other = (Layer) obj;
		return Objects.equals(mDensity, other.mDensity) //
				&& Objects.equals(mThickness, other.mThickness) //
				&& Objects.equals(mName, other.mName) //
				&& Objects.equals(mComposition, other.mComposition); //
	}

	/**
	 * @param comp
	 * @param thickness_nm In nanometers
	 * @param density      In g/cm<sup>3</sup>
	 */
	public Layer(
			final Composition comp, final UncertainValue thickness_nm, final UncertainValue density
	) {
		assert comp != null;
		assert (thickness_nm != null) && (thickness_nm.doubleValue() > 0.0);
		assert (density != null) && (density.doubleValue() > 0.0);
		mComposition = comp;
		mThickness = thickness_nm.multiply(1.0e-7); // 1 nm = 1.0e-7 cm
		mDensity = density;
		final HalfUpFormat nf = new HalfUpFormat("0.0");
		final String name = mComposition.toString() + "[" + nf.format(1.0e7 * mThickness.doubleValue()) + " nm, "
				+ nf.format(mDensity.doubleValue()) + " g/cm^3]";
		mName = new UniqueString(name);
	}

	/**
	 * Makes a Layer consisting of thickness_nm (in nanometers) of carbon with
	 * density of 1.5 &pm; 0.1 g/cm<sup>3</sup>.
	 *
	 * @param thickness_nm Thickness in nanometers
	 * @return Layer object
	 * @throws ArgumentException
	 */
	public static Layer carbonCoating(
			final UncertainValue thickness_nm
	) throws ArgumentException {
		return new Layer(Composition.pureElement(Element.Carbon), thickness_nm, new UncertainValue(1.5, 0.1));
	}

	/**
	 * Makes a Layer consisting of thickness_nm (in nanometers) of gold with density
	 * of 19.3 &pm; 0.1 g/cm<sup>3</sup>.
	 *
	 * @param thickness_nm
	 * @return Layer object
	 * @throws ArgumentException
	 */
	public static Layer goldCoating(
			final UncertainValue thickness_nm
	) throws ArgumentException {
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
	 * @return in centimeters
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

	public UniqueString getName() {
		return mName;
	}

	/**
	 * @return in g/cm<sup>2</sup>
	 */
	public UncertainValue getMassThickness() {
		return new UncertainValue(mThickness.doubleValue() * mDensity.doubleValue(), //
				(mThickness.doubleValue() * mDensity.uncertainty()
						+ mDensity.doubleValue() * mThickness.uncertainty()));
	}

	@Override
	public String toHTML(
			final Mode mode
	) {

		final HalfUpFormat nf = new HalfUpFormat("0.0");
		return mComposition.toHTML(mode.demote()) + "[" + nf.format(1.0e7 * mThickness.doubleValue()) + " nm, "
				+ nf.format(mDensity.doubleValue()) + " g/cm^3]";
	}

	@Override
	public String toString() {
		return mName.toString();
	}

}
