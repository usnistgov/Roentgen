package gov.nist.microanalysis.roentgen.spectrum;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class EfficencyCalibration {

	private final PolynomialSplineFunction mCalibration;

	/**
	 * Constructs a EfficencyCalibration
	 */
	public EfficencyCalibration(final double[] xval, final double[] yval) {
		mCalibration = (new SplineInterpolator()).interpolate(xval, yval);
	}

	public static EfficencyCalibration siliconDetector(final double thickness, final XRayWindow window) {
		final Element si = Element.Silicon;
		final ElementalMAC mac = new ElementalMAC();

		return null;

	}

}
