package gov.nist.microanalysis.roentgen.physics.composition;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;

public class Layer implements IToHTML {

	private final Composition mComposition;
	private final UncertainValue mThickness;
	private final UncertainValue mDensity;

	/**
	 * @param comp
	 * @param thickness In micrometers
	 * @param density   In g/cm<sup>3</sup>
	 */
	public Layer(final Composition comp, final UncertainValue thickness, final UncertainValue density) {
		assert comp != null;
		assert (thickness != null) && (thickness.doubleValue() > 0.0);
		assert (density != null) && (density.doubleValue() > 0.0);
		mComposition = comp;
		mThickness = thickness;
		mDensity = density;
	}

	public Composition getComposition() {
		return mComposition;
	}

	/**
	 * @return in micrometers
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

	@Override
	public String toHTML(final Mode mode) {
		return mComposition.toHTML(mode.demote()) + "[" + mThickness.toHTML(mode) + "&mu;m," + mDensity.toHTML(mode)
				+ " g/cm<sup>3</sup>]";
	}
}
