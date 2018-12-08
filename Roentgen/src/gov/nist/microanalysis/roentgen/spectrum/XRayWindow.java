package gov.nist.microanalysis.roentgen.spectrum;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>
 * Implements an x-ray window.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class XRayWindow {

	private final PolynomialSplineFunction mTransparency;

	/**
	 * Constructs an XRayWindow from a set of x,y data points that define the
	 * general outline of the curve. The intermediate values are interpolated using
	 * a spline interpolation algorithm.
	 */
	private XRayWindow(final double[] xvals, final double[] yvals) {
		final SplineInterpolator si = new SplineInterpolator();
		mTransparency = si.interpolate(xvals, yvals);
	}

	private XRayWindow(final PolynomialSplineFunction psf) {
		mTransparency = psf;
	}

	public static XRayWindow create(final Composition[] mf, final double[] massThickness) {

		return null;

	}

	public double getTransmission(final double energy) {
		return mTransparency.value(energy);
	}

	public double[] getKnots() {
		return mTransparency.getKnots();
	}

}
