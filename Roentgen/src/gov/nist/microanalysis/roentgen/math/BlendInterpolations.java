package gov.nist.microanalysis.roentgen.math;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * <p>
 * Takes a number of {@link PolynomialSplineFunction} interpolated functions and
 * creates a {@link PolynomialSplineFunction} that represents the linear
 * combination. Care is taken to ensure that the resulting interpolation uses
 * essentially the same knots (to within a tolerance) to ensure the best
 * possible fit.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class BlendInterpolations {

	private final Map<PolynomialSplineFunction, Double> mFunctions;

	/**
	 * Constructs a BlendInterpolations
	 */
	public BlendInterpolations() {
		mFunctions = new HashMap<>();
	}

	public void add(final PolynomialSplineFunction psf, final double scale) {
		mFunctions.put(psf, Double.valueOf(scale));
	}

	public PolynomialSplineFunction compute(final double tol) {
		final Set<Double> xvals = new TreeSet<>();
		for (final Map.Entry<PolynomialSplineFunction, Double> me : mFunctions.entrySet()) {
			final PolynomialSplineFunction psf = me.getKey();
			for (final double knot : psf.getKnots())
				xvals.add(knot);
		}
		double[] xv = new double[xvals.size()];
		{
			int i = 0;
			final double prev = -Double.MAX_VALUE;
			for (final double dv : xvals)
				if (dv - prev > tol) {
					xv[i] = dv;
					++i;
				}
			xv = Arrays.copyOf(xv, i);
		}
		final double[] yv = new double[xv.length];
		for (int i = 0; i < xv.length; ++i)
			for (final Map.Entry<PolynomialSplineFunction, Double> me : mFunctions.entrySet())
				yv[i] = me.getValue().doubleValue() * me.getKey().value(xv[i]);
		final SplineInterpolator si = new SplineInterpolator();
		return si.interpolate(xv, yv);
	}
}
