package gov.nist.microanalysis.roentgen.physics;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * <p>
 * Calculates the mass absorption coefficient by performing log-log
 * interpolation on the Chantler-2005 database.
 * </p>
 * <p>
 * Chantler, C.T., Olsen, K., Dragoset, R.A., Chang, J., Kishore, A.R.,
 * Kotochigova, S.A., and Zucker, D.S. (2005), <i>X-Ray Form Factor, Attenuation
 * and Scattering Tables (version 2.1).</i> [Online] Available:
 * http://physics.nist.gov/ffast [Wednesday, 12-Oct-2016 08:02:04 EDT]. National
 * Institute of Standards and Technology, Gaithersburg, MD.
 * </p>
 * <p>
 * Originally published as Chantler, C.T., J. Phys. Chem. Ref. Data
 * <b>29</b>(4), 597-1048 (2000); and Chantler, C.T., J. Phys. Chem. Ref. Data
 * <b>24</b>, 71-643 (1995).
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 280 $
 */
public class ElementalMAC {

	private double mLimit = 1.0;

	// Only is created once!
	private static Implementation mImplementation = new Implementation();

	private static class Implementation {

		private final Map<Element, PolynomialSplineFunction> mData;

		private Implementation() {
			mData = loadData();
		}

		private final Map<Element, PolynomialSplineFunction> loadData() {
			final InputStream is = ElementalMAC.class.getResourceAsStream("Chantler2005MAC.csv");
			final BufferedReader br = new BufferedReader(new InputStreamReader(is));
			final Map<Element, PolynomialSplineFunction> res = new HashMap<>();
			try {
				String keVline = br.readLine();
				String macLine = keVline != null ? br.readLine() : null;
				Element elm = Element.Hydrogen;
				final LinearInterpolator uvi = new LinearInterpolator();
				while ((keVline != null) && (macLine != null)) {
					final String[] eVItems = keVline.split(",");
					final String[] macItems = macLine.split(",");
					assert eVItems.length == macItems.length;
					final double[] eVres = new double[eVItems.length];
					final double[] macRes = new double[eVItems.length];
					for (int i = 0; i < eVItems.length; ++i) {
						eVres[i] = Math.log(1000.0 * Double.parseDouble(eVItems[i].trim()));
						macRes[i] = Math.log(Double.parseDouble(macItems[i].trim()));
					}
					// The tables are linearly interpolated in log-log space.
					final PolynomialSplineFunction psf = uvi.interpolate(eVres, macRes);
					res.put(elm, psf);
					keVline = br.readLine();
					macLine = keVline != null ? br.readLine() : null;
					elm = elm.next();
				}
			} catch (final Exception e) {
				final String msg = "Unable to read Chantler2005 MAC data file.";
				Globals.getLogger().error(msg, e);
				throw new Error(msg, e);
			}
			return res;
		}

	}

	private double limit(final double frac) {
		return frac < mLimit ? frac : mLimit;
	}

	public static class ElementMAC extends BaseLabel<Element, XRay, Object> {

		public ElementMAC(final Element elm, final XRay xr) {
			super("[&mu;/&rho;]", elm, xr);
		}

		public Element getElement() {
			return getObject1();
		}

		public XRay getXRay() {
			return getObject2();
		}

	}

	public static final double MIN_E = 0.0;
	public static final double MAX_E = 432.9451e3;

	/**
	 * Constructs a MassAbsorptionCoefficient2
	 */
	public ElementalMAC() {
		mLimit = 2.0;
	}

	public ElementalMAC(final double limit) {
		mLimit = Math.max(0.0, limit);
	}

	/**
	 * @param el The Element
	 * @param eV The energy in eV
	 * @return double
	 */
	private double fractionalUncertainty(final Element el, final double eV) {
		double err = 0.0;
		if (eV < 200.0)
			err = 1.50; // 100-200%
		else if (eV < 500)
			err = 0.5 + (((1.0 - 0.5) * (eV - 200)) / (500.0 - 200.0)); // 50-100%
		else if (eV < 1000)
			err = 0.05 + (((0.20 - 0.05) * (eV - 500)) / (1000.0 - 500.0));
		for (final AtomicShell sh : AtomicShell.forElement(el)) {
			final double ee = sh.getEdgeEnergy();
			final double delta = (eV - ee) / eV;
			if ((Math.abs(delta) < 0.001) || (Math.abs(eV - ee) < 5.0)) // Near edges (within 0.1%)
				err = Math.max(err, 0.5);
			else
				switch (sh.getShell()) {
				case K:
					if (Math.abs(delta) < 0.1)
						err = Math.max(err, 0.15); // 10-20%
					else if ((eV > ee) && (eV < (1.1 * ee)))
						err = Math.max(err, 0.03); // 3%
					else
						err = Math.max(err, 0.01); // 1%
					break;
				case L1:
				case M1:
				case M2:
				case M3:
				case N1: // Special case for N shells added by
					// NWMR
				case N3:
				case N4:
				case N5:
					if (Math.abs(delta) < 0.15)
						err = Math.max(err, 0.225); // 15-30%
					else if ((delta > 0) && (delta < 0.4))
						err = Math.max(err, 0.04); // 4%
					else
						err = Math.max(err, 0.01); // 1%
					break;
				case L2:
				case L3:
				case M4:
				case M5:
				case N6: // Special case for N shells added
					// by NWMR
				case N7:
					if (Math.abs(delta) < 0.15)
						err = Math.max(err, 0.3); // 20-40%
					else if ((delta > 0.0) && (delta < 0.40))
						err = Math.max(err, 0.04); // 4%
					else
						err = Math.max(err, 0.01); // 1%
					break;
				default:
					break;
				}
		}
		return limit(err);
	}

	/**
	 * @param el The Element
	 * @param eV The energy in eV
	 * @return double
	 */
	public String fractionalUncertaintySource(final Element el, final double eV) {
		String result = "None";
		final HalfUpFormat df = new HalfUpFormat("0.0%");
		double err = 0.0;
		if (eV < 200.0) {
			result = "eV < 200.0"; // 100-200%
			err = 1.5;
		} else if (eV < 500) {
			result = "eV < 500";
			err = 0.5 + (((1.0 - 0.5) * (eV - 200)) / (500.0 - 200.0)); // 50-100%
		} else if (eV < 1000) {
			result = "eV < 1000";
			err = 0.05 + (((0.20 - 0.05) * (eV - 500)) / (1000.0 - 500.0));
		}
		for (final AtomicShell sh : AtomicShell.forElement(el)) {
			final double ee = sh.getEdgeEnergy();
			final double delta = (eV - ee) / eV;
			if ((Math.abs(delta) < 0.001) || (Math.abs(eV - ee) < 5.0)) { // Near edges (within 0.1%)
				if (err < 0.5) {
					result = "Near " + sh.toString() + " edge (within 0.1%)";
					err = 0.5;
				}
			} else
				switch (sh.getShell()) {
				case K:
					if (Math.abs(delta) < 0.1) {
						if (err < 0.15) {
							result = "Near " + sh.toString() + " edge (delta < 0.1)";
							err = 0.15; // 10-20%
						}
					} else if ((eV > ee) && (eV < (1.1 * ee))) {
						if (err < 0.03) {
							result = "Near " + sh.toString() + " edge ((eV > ee) && (eV < (1.1 * ee))";
							err = 0.03;
						}
					} else {
						if (err < 0.01) {
							result = "Above " + sh.toString();
							err = 0.01;
						}
					}
					break;
				case L1:
				case M1:
				case M2:
				case M3:
				case N1: // Special case for N shells added by
					// NWMR
				case N3:
				case N4:
				case N5:
					if (Math.abs(delta) < 0.15) {
						if (err < 0.225) {
							result = "Math.abs(delta) < 0.15 " + sh.toString();
							err = 0.225;
						}
					} else if ((delta > 0) && (delta < 0.4)) {
						if (err < 0.04) {
							result = "(delta > 0) && (delta < 0.4) " + sh.toString();
							err = 0.04;
						}
					} else {
						if (err < 0.01) {
							result = "Above " + sh.toString();
							err = 0.01;
						}
					}
					break;
				case L2:
				case L3:
				case M4:
				case M5:
				case N6: // Special case for N shells added
					// by NWMR
				case N7:
					if (Math.abs(delta) < 0.15) {
						if (err < 0.3) {
							result = "Math.abs(delta) < 0.15 " + sh.toString();
							err = 0.3;
						}
					} else if ((delta > 0.0) && (delta < 0.40)) {
						if (err < 0.04) {
							result = "(delta > 0) && (delta < 0.4) " + sh.toString();
							err = 0.04;
						}
					} else {
						if (err < 0.01) {
							result = "Above " + sh.toString();
							err = 0.01;
						}
					}
					break;
				default:
					break;
				}
		}
		return result + " -> " + df.format(limit(err));
	}

	/**
	 * Computes the mass absorption coefficient for the specified element at the
	 * specified x-ray energy.
	 *
	 * @param el     Element
	 * @param energy in eV
	 * @return UncertainValue In cm^2/g
	 */
	public UncertainValue compute(final Element el, final double energy) {
		assert mImplementation.mData.get(el) != null;
		assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
		assert energy < MAX_E : "Energy too high.";
		final double res = Math.exp(mImplementation.mData.get(el).value(Math.log(energy)));
		return new UncertainValue(res, new ElementMAC(el, new XRay(energy)), res * fractionalUncertainty(el, energy));
	}

	/**
	 * Returns a PolynomialSplineFunction that implements a cubic spline
	 * interpolation function over the mass absorption coefficent values for the
	 * specified element. The values returned by this function are likely to agree
	 * well but not perfectly with those computed using the other compute functions
	 * in this class which use a log-log interpolations.
	 *
	 * @param elm
	 * @return UnivariateFunction
	 */
	public PolynomialSplineFunction getMACFunction(final Element elm) {
		final PolynomialSplineFunction uvf = mImplementation.mData.get(elm);
		final double[] x = uvf.getKnots();
		final double[] y = new double[x.length];
		for (int i = 0; i < x.length; ++i)
			y[i] = uvf.value(x[i]);
		return (new SplineInterpolator()).interpolate(x, y);
	}

	/**
	 * Returns an array containing the energies at which values of the MACs were
	 * provided that served as the basis for the interpolation.
	 *
	 * @param elm
	 * @return double[]
	 */
	public double[] getKnots(final Element elm) {
		final double[] res = mImplementation.mData.get(elm).getKnots();
		for (int i = 0; i < res.length; ++i)
			res[i] = Math.exp(res[i]);
		return res;
	}

	/**
	 * @param el     Element
	 * @param energy in eV
	 * @return true if a MAC is available for this element
	 */
	public boolean isAvailable(final Element el, final XRay xr) {
		assert xr != null;
		final double energy = xr.getEnergy();
		assert el != null;
		if (!mImplementation.mData.containsKey(el))
			return false;
		return (energy > MIN_E) && (energy <= MAX_E);
	}

	/**
	 * Computes the mass absorption coefficient for the specified element at the
	 * specified x-ray energy.
	 *
	 * @param el         Element
	 * @param XRayPhoton
	 * @return UncertainValue In cm^2/g
	 */
	public UncertainValue compute(final Element el, final XRay xr) {
		assert mImplementation.mData.get(el) != null;
		final double energy = xr.getEnergy();
		assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
		assert energy < MAX_E : "Energy too high.";
		final double res = Math.exp(mImplementation.mData.get(el).value(Math.log(energy)));
		return new UncertainValue(res, new ElementMAC(el, xr), res * fractionalUncertainty(el, energy));
	}

	/**
	 * Returns the elemental MAC associated with the specified material and xray
	 * energy in eV.
	 *
	 * @param element
	 * @param energy  (eV)
	 * @return MAC in cm^2/g
	 */
	public double computeQ(final Element el, final double energy) {
		assert mImplementation.mData.get(el) != null;
		assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
		assert energy < MAX_E : "Energy too high.";
		return Math.exp(mImplementation.mData.get(el).value(Math.log(energy)));
	}

	/**
	 * Computes the MAC associated with the specified material and xray energy in
	 * eV.
	 *
	 * @param material
	 * @param energy   (eV)
	 * @return MAC in cm^2/g
	 */
	public double computeQ(final Composition comp, final double energy) {
		Preconditions.checkArgument(energy > MIN_E, "Can't compute a MAC for a negative or zero x-ray energy");
		Preconditions.checkArgument(energy < MAX_E, "Energy too high.");
		final Set<Element> elms = comp.getElementSet();
		double res = 0.0;
		for (final Element elm : elms)
			res += computeQ(elm, energy) * comp.getValue(elm).doubleValue();
		return res;
	}

	/**
	 * Computes the MAC associated with the specified material and xray.
	 *
	 * @param material
	 * @param xr
	 * @return MAC in cm^2/g
	 */
	public double computeQ(final Composition material, final XRay xr) {
		return computeQ(material, xr.getEnergy());
	}

	/**
	 * Set a maximum limit for the fractional uncertainty. Default = 0.9 or 90%.
	 *
	 * @param limit
	 */
	public void setLimit(final double limit) {
		mLimit = Math.max(0.0, limit);
	}

	/**
	 * Returns the largest fractional uncertainty that the function
	 * fractionalUncertainty(...) will return.
	 *
	 *
	 * @return limit
	 */
	public double getLimit() {
		return mLimit;
	}

}
