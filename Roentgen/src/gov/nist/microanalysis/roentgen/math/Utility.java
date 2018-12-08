package gov.nist.microanalysis.roentgen.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.Function;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.MathArrays;

import com.google.common.base.Preconditions;

/**
 * <p>
 * A collection of numeric utility functions.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 285 $
 */
public class Utility {

	public static final double bound(final double a, final double min, final double max) {
		Preconditions.checkArgument(min <= max);
		return a < min ? min : (a > max ? max : a);
	}

	/**
	 * Returns min if val is less than min, max-1 if val is greater than or equal to
	 * max, otherwise val.
	 *
	 * @param val
	 * @param min Inclusive
	 * @param max Exclusive
	 * @return val bounded to between [min, max)
	 */
	public static int bound(final int val, final int min, final int max) {
		Preconditions.checkArgument(min <= max);
		return val < min ? min : (val >= max ? max - 1 : val);

	}

	/**
	 * Ensures that all numbers in the array vals are between min and max by
	 * increasing smaller values to min and decreasing larger values to max. All
	 * vals that are NaN are replaced with defVal.
	 *
	 * @param vals
	 * @param min    Lower limit
	 * @param max    Upper limit
	 * @param defVal The value which NaN values will be replaced with
	 * @return double[] - Freshly allocated
	 */
	public static double[] bound(final double[] vals, final double min, final double max, final double defVal) {
		Preconditions.checkArgument(min <= max);
		final double[] res = new double[vals.length];
		for (int i = 0; i < res.length; ++i) {
			double tmp = vals[i];
			if (Double.isNaN(tmp))
				tmp = defVal;
			else if (tmp < min)
				tmp = min;
			else if (tmp > max)
				tmp = max;
			res[i] = tmp;
		}
		return res;
	}

	/**
	 * Returns min if val is less than min, max-1 if val is greater than or equal to
	 * max, otherwise val.
	 *
	 * @param val
	 * @param min Inclusive
	 * @param max Exclusive
	 * @return val bounded to between [min, max)
	 */
	public static long bound(final long val, final long min, final long max) {
		Preconditions.checkArgument(min <= max);
		return val < min ? min : (val >= max ? max - 1 : val);

	}

	/**
	 * A numerically stable method for solving the quadratic equation
	 * a&middot;x<sup>2</sup> + b&middot;x + c = 0.
	 *
	 * @param a
	 * @param b
	 * @param c
	 * @return double[2] containing the solutions (smaller, larger) or null if there
	 *         is no real solution.
	 */
	static final public double[] quadraticSolver(final double a, final double b, final double c) {
		final double r = (b * b) - (4.0 * a * c);
		if (r < 0.0)
			return null;
		final double q = -0.5 * (b + (Math.signum(b) * Math.sqrt(r)));
		final double[] res = new double[] { q / a, c / q };
		if (res[0] > res[1]) {
			final double tmp = res[0];
			res[0] = res[1];
			res[1] = tmp;
		}
		return res;
	}

	static final public double[] multiply(final double a, final double[] v) {
		return MathArrays.scale(a, v);
	}

	static final public double[] divide(final double[] v, final double den) {
		return MathArrays.scale(1.0 / den, v);
	}

	static final public double dot(final double[] a, final double[] b) {
		return MathArrays.linearCombination(a, b);
	}

	static final public double sum(final double... a) {
		double sum = 0.0;
		for (final double v : a)
			sum += v;
		return sum;
	}

	static final public double sum(final Double... a) {
		double sum = 0.0;
		for (final double v : a)
			sum += v;
		return sum;
	}

	static final public int sum(final int... a) {
		int sum = 0;
		for (final int v : a)
			sum += v;
		return sum;
	}

	static final public double mean(final double... a) {
		return sum(a) / a.length;
	}

	static final public double max(final double... a) {
		double max = a[0];
		for (final double v : a)
			max = Math.max(max, v);
		return max;
	}

	static final public int max(final int... a) {
		int max = a[0];
		for (final int v : a)
			max = Math.max(max, v);
		return max;
	}

	static final public double max(final Collection<Double> a) {
		double max = a.iterator().next();
		for (final double v : a)
			max = Math.max(max, v);
		return max;
	}

	static final public double min(final double... a) {
		double min = a[0];
		for (final double v : a)
			min = Math.min(min, v);
		return min;
	}

	static final public int min(final int... a) {
		int min = a[0];
		for (final int v : a)
			min = Math.min(min, v);
		return min;
	}

	static final public double min(final Collection<Double> a) {
		double min = a.iterator().next();
		for (final double v : a)
			min = Math.min(min, v);
		return min;
	}

	static final public double[] asArray(final Collection<Double> cd) {
		final double[] res = new double[cd.size()];
		int i = 0;
		for (final double d : cd)
			res[i++] = d;
		return res;
	}

	static final public double sumOfSquares(final double... a) {
		double sum = 0.0;
		for (final double v : a)
			sum += v * v;
		return sum;
	}

	static final public double variance(final double... a) {
		return sumOfSquares(a) / a.length;
	}

	static final public double sampleStdDev(final double... a) {
		return Math.sqrt(sumOfSquares(a) / (a.length - 1));
	}

	static final public double stdDev(final double... a) {
		return Math.sqrt(sumOfSquares(a) / a.length);
	}

	/**
	 * Performs an element-by-element divide
	 *
	 * @param num
	 * @param den
	 * @return double[] res
	 */
	static final public double[] ebeDivide(final double[] num, final double[] den) {
		return MathArrays.ebeDivide(num, den);
	}

	/**
	 * expRand - Selects a random value from an exponential distibution. The mean
	 * value returned is 1.0.
	 *
	 * @return double - Returns a random variable in the range [0,infinity)
	 */
	static final public double expRand() {
		return -Math.log(1.0 - ThreadLocalRandom.current().nextDouble());
	}

	/**
	 * expRand - Selects a random value from an exponential distibution. The mean
	 * value returned is 1.0.
	 *
	 * @param rand A uniformly distributed random variant between [0 (Inclusive) ,1
	 *             (Exclusive))
	 * @return double - Returns a an exponential distributed random variable in the
	 *         range [0,infinity)
	 */
	static final public double expRand(final double rand) {
		assert rand >= 0.0;
		assert rand < 1.0;
		return -Math.log(1.0 - rand);
	}

	/**
	 * Computes a random 3-vector uniform in solid angle using the algorithm of
	 * Robert Knop in Commun. ACM, ACM, 1970, 13, 326. I have tested this routine
	 * and it does seem to produce directions that are truly distributed evenly
	 * about the unit sphere.
	 *
	 * @return Vector3D
	 */
	static final public Vector3D randomDirection() {
		double x, y, s;
		do {
			x = 2.0 * (Math.random() - 0.5);
			y = 2.0 * (Math.random() - 0.5);
			s = (x * x) + (y * y);
		} while (s > 1.0);
		final double z = (2.0 * s) - 1.0;
		s = Math.sqrt((1 - (z * z)) / s);
		return new Vector3D(s * x, s * y, z);
	}

	/**
	 * Scatters the initial direction by a polar scatter angle theta and an
	 * azimuthal scatter angle phi. See Salvat 2011 1.131 and 1.132
	 *
	 * @param initial Initial direction (normalized)
	 * @param theta   0 to &pi; Polar angle
	 * @param phi     0 to 2&pi; Azimuthal angle
	 * @return Vector3D normalized
	 */
	static final public Vector3D scatter(final Vector3D initial, final double theta, final double phi) {
		assert theta >= 0.0;
		assert theta <= Math.PI;
		assert phi >= 0.0;
		assert phi <= (2.0 * Math.PI);
		assert Math.abs(initial.getNorm() - 1.0) < 1.0e-8;
		final double z0 = initial.getZ();
		final double ct = Math.cos(theta), st = Math.sin(theta);
		final double cp = Math.cos(phi), sp = Math.sin(phi);
		if (Math.abs(1.0 - z0) > 1.0e-8) {
			// Not along Z-axis
			final double x0 = initial.getX(), y0 = initial.getY();
			final double sw = Math.sqrt(1.0 - (z0 * z0));
			final double up = (x0 * ct) + ((st * ((x0 * z0 * cp) - (y0 * sp))) / sw);
			final double vp = (y0 * ct) + ((st * ((y0 * z0 * cp) + (x0 * sp))) / sw);
			final double wp = (z0 * ct) - (sw * st * cp);
			final double n = Math.sqrt((up * up) + (vp * vp) + (wp * wp));
			final Vector3D newDir = new Vector3D(up / n, vp / n, wp / n);
			assert (Vector3D.angle(newDir, initial) - theta) < 1.0e-6 : "theta = " + theta + " ~ "
					+ Vector3D.angle(newDir, initial) + "   phi = " + phi + " == " + initial;
			return newDir;
		} else {
			// Along Z-axis
			final double signum = z0 > 0.0 ? 1.0 : -1.0;
			return new Vector3D(signum * st * cp, signum * st * sp, signum * ct);
		}
	}

	/**
	 * Scatters the initial direction by a polar scatter angle theta and an
	 * azimuthal scatter angle phi. See Mathematica notebook scatter.nb in
	 * resources.
	 *
	 * @param initial Initial direction (normalized)
	 * @param phi     0 to &pi; Polar angle
	 * @param eta     0 to 2&pi; Azimuthal angle
	 * @return Vector3D normalized
	 */
	static final public Vector3D scatter2(final Vector3D initial, final double theta, final double phi) {
		assert theta >= 0.0;
		assert theta <= Math.PI;
		assert phi >= 0.0;
		assert phi <= (2.0 * Math.PI);
		assert Math.abs(initial.getNorm() - 1.0) < 1.0e-8;
		final double ct = Math.cos(theta), st = Math.sin(theta);
		final double cp = Math.cos(phi), sp = Math.sin(phi);
		final double z0 = initial.getZ();
		if (Math.abs(1 - z0) > 1.0e-8) {
			final double x0 = initial.getX(), y0 = initial.getY();
			final double den = Math.sqrt((x0 * x0) + (y0 * y0));
			final double x = (x0 * ct) + ((st * ((y0 * sp) + (x0 * z0 * cp))) / den);
			final double y = (y0 * ct) + ((st * ((y0 * z0 * cp) - (x0 * sp))) / den);
			final double z = (z0 * ct) - (den * cp * st);
			final double norm = Math.sqrt((x * x) + (y * y) + (z * z));
			assert Math.abs(norm - 1.0) < 1.0e-6;
			final Vector3D newDir = new Vector3D(x / norm, y / norm, z / norm);
			assert (Vector3D.angle(newDir, initial) - theta) < 1.0e-6 : "theta = " + theta + " ~ "
					+ Vector3D.angle(newDir, initial) + "   phi = " + phi + " == " + initial;
			return newDir;

		} else {
			final double signum = z0 > 0.0 ? 1.0 : -1.0;
			return new Vector3D(signum * st * cp, signum * st * sp, signum * ct);
		}
	}

	/**
	 * Computes a equally spaced array of double values between minVal and maxVal.
	 *
	 * @param minVal The initial value res[0]
	 * @param maxVal The final value - res[res.length-1]
	 * @param len    Array length
	 * @return double[len]
	 */
	public static double[] linearArray(final double minVal, final double maxVal, final int len) {
		final double[] res = new double[len];
		for (int i = 0; i < len; ++i)
			res[i] = minVal + (((maxVal - minVal) * i) / (len - 1));
		return res;
	}

	/**
	 * <p>
	 * Computes an array of double of length 'len' in which the i-th value is
	 * computed using the function 'function' whose argument is 'i'. Use with lambda
	 * functions as in:
	 * </p>
	 * <p>
	 * <code>Utility.functionalArray(100, i -> Math.log(i))
	 * </code>
	 * </p>
	 *
	 * @param len      Number of items in the array
	 * @param function Function&lt;Integer,Double&gt;
	 * @return double[len]
	 */
	public static double[] functionalArray(final int len, final Function<Integer, Double> function) {
		final double[] res = new double[len];
		for (int i = 0; i < len; ++i)
			res[i] = function.apply(i);
		return res;
	}

	/**
	 * Creates an array of length 'len' that is computed by applying 'function' to
	 * an equally spaced set of values in the range minVal to maxVal.
	 *
	 * @param minVal   double
	 * @param maxVal   double
	 * @param len      int
	 * @param function Function&lt;Double,Double&gt;
	 * @return double[]
	 */
	public static double[] functionalArray(final double minVal, final double maxVal, final int len,
			final Function<Double, Double> function) {
		return apply(linearArray(minVal, maxVal, len), function);
	}

	/**
	 * Takes an array da of type double and creates a new array by applying
	 * 'function' to the elements 'd' in 'da' in index order.
	 *
	 * @param da       array of double
	 * @param function Function&lt;Double,Double&gt;
	 * @return double[]
	 */
	public static double[] apply(final double[] da, final Function<Double, Double> function) {
		final double[] res = new double[da.length];
		for (int i = 0; i < da.length; ++i)
			res[i] = function.apply(da[i]);
		return res;
	}

	/**
	 * Takes an array ta of type T objects and creates a new array by applying
	 * function to the elements 't' in 'ta' in index order.
	 *
	 * @param ta       An array of T
	 * @param function A function from T to Double/double
	 * @return double[] An array of double of length ta.length.
	 */
	public static <T> double[] apply(final T[] ta, final Function<T, Double> function) {
		final double[] res = new double[ta.length];
		for (int i = 0; i < res.length; ++i)
			res[i] = function.apply(ta[i]).doubleValue();
		return res;
	}

	/**
	 * Takes an array ta of type T objects and creates a new S[] based on 'function'
	 * applied to 't' in 'ta'.
	 *
	 * @param ta
	 * @param function
	 * @return S[]
	 */
	@SuppressWarnings("unchecked")
	public static <S, T> S[] apply2(final T[] ta, final Function<T, S> function) {
		final ArrayList<S> res = new ArrayList<>();
		for (final T t : ta)
			res.add(function.apply(t));
		return (S[]) res.toArray();
	}

	/**
	 * Computes a matrix of numbers
	 *
	 * @param xVals An array of x values
	 * @param yVals An array of y values
	 * @param func  Function taking a Double[2] { x, y} as an argument returning a
	 *              Double
	 * @return double[xVals.length][yVals.length] A matrix computed from func( { x,
	 *         y})
	 */
	public static double[][] fillMesh(final double[] xVals, final double[] yVals,
			final Function<Double[], Double> func) {
		final double[][] mesh = new double[xVals.length][yVals.length];
		for (int i = 0; i < xVals.length; ++i)
			for (int j = 0; j < yVals.length; ++j) {
				final Double[] args = new Double[] { xVals[i], yVals[j] };
				mesh[i][j] = func.apply(args);
			}
		return mesh;
	}

	/**
	 * Convert an array of double values into a tab-separated table. The numbers
	 * always use "." as the decimal spacer regardless of locale.
	 *
	 * @param vals
	 * @param spacer - Nominally "\t" or ", "
	 * @return String
	 */
	public static String tabulate(final double[][] vals, final String spacer) {
		final StringBuffer sb = new StringBuffer();
		for (final double[] val : vals) {
			for (int j = 0; j < val.length; ++j) {
				if (j > 0)
					sb.append(spacer);
				sb.append(Double.toString(val[j]));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	/**
	 * Convert an array of double values into a tab-separated table. The numbers
	 * always use "." as the decimal spacer regardless of locale.
	 *
	 * @param xVals  - Row labels
	 * @param yVals  - Col labels
	 * @param vals
	 * @param spacer - Nominally "\t" or ", "
	 * @return String
	 */
	public static String tabulate(final double[] xVals, final double[] yVals, final double[][] vals,
			final String spacer) {
		final StringBuffer sb = new StringBuffer();
		assert vals.length == xVals.length;
		// Header of yVals
		sb.append("");
		for (final double yVal : yVals) {
			sb.append(spacer);
			sb.append(Double.toString(yVal));
		}
		sb.append("\n");
		for (int i = 0; i < vals.length; ++i) {
			final double[] val = vals[i];
			assert val.length == yVals.length;
			sb.append(Double.toString(xVals[i]));
			for (final double element : val) {
				sb.append(spacer);
				sb.append(Double.toString(element));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	public enum CONVOLVE_END_MODE {
		/**
		 * Assume -1, -2, and size(), size()+1 etc are replaced with zero.
		 */
		ASSUME_ZERO,
		/**
		 * Assume -1, -2,.. are replaced with v[0] and size(), size()+1 etc are replaced
		 * with v[size()-1].
		 */
		COPY_END_VALUE,
		/**
		 * Assume v[-1]= v[1], v[-2]=v[2] and v[size()]=v[size()-2], etc.
		 */
		REVERSE_ORDER_COPIES
	};

	/**
	 * Convolves the input array with the specified kernel using the specified mode
	 * to determine how the ends are handled.
	 *
	 * @param inp    double[]
	 * @param kernel double[]
	 * @param cm     {@link CONVOLVE_END_MODE}
	 * @return
	 */
	public static double[] convolve(final double[] inp, final double[] kernel, final CONVOLVE_END_MODE cm) {
		return convolve(inp, kernel, cm, 1);
	}

	/**
	 * Convolves the input array with the specified kernel using the specified mode
	 * to determine how the ends are handled.
	 *
	 * @param inp    double[]
	 * @param kernel double[]
	 * @param cm     {@link CONVOLVE_END_MODE}
	 * @param reps   int Number of times to apply convolution kernel
	 * @return
	 */
	public static double[] convolve(final double[] inp, final double[] kernel, final CONVOLVE_END_MODE cm,
			final int reps) {
		final int n = kernel.length / 2;
		final double[] remap = new double[inp.length + (2 * n)];
		System.arraycopy(inp, 0, remap, n, inp.length);
		switch (cm) {
		case ASSUME_ZERO:
			// Don't do anything
			break;
		case COPY_END_VALUE:
			// Replicate the first and last value
			for (int i = 0; i < n; ++i) {
				remap[i] = remap[n];
				remap[remap.length - (i + 1)] = remap[remap.length - (n + 1)];
			}
			break;
		case REVERSE_ORDER_COPIES:
			for (int i = 1; i <= n; ++i) {
				remap[n - i] = remap[n + i];
				remap[remap.length - (n + 1 + i)] = remap[(remap.length - (n + 1)) + i];
			}
			break;
		}
		double[] tmp = remap;
		for (int r = 0; r < reps; r++)
			tmp = MathArrays.convolve(tmp, kernel);
		final int start = n + ((reps * (kernel.length - 1)) / 2);
		return Arrays.copyOfRange(tmp, start, inp.length + start);
	}

	/**
	 * Apply a Gaussian Kernel Smoother to the data.
	 *
	 * @param data
	 * @param width Width of the kernel parameter in channels
	 * @param cm    {@link CONVOLVE_END_MODE}
	 * @return double[] of the same length as the input.
	 */
	public static double[] gaussianKernelSmooth(final double[] data, final double width, final CONVOLVE_END_MODE cm) {
		final int dcw = (int) Math.ceil(3.0 * width);
		final double[] kernel = new double[(2 * dcw) + 1];
		for (int i = -dcw; i <= dcw; ++i)
			kernel[i + dcw] = Math.exp(-0.5 * Math.pow(i / width, 2.0));
		return convolve(data, MathArrays.normalizeArray(kernel, 1.0), cm);
	}

	/*
	 * Coefficients extracted on 31-Oct-2016 from
	 * https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter&section=17
	 */

	/**
	 * Savitsky-Golay convolution coefficients for smoothing - Quadratic/Cubic
	 * degree 5
	 */
	static public final double[] SAVITZKY_GOLAY_QC5 = { -3.0 / 35.0, 12.0 / 35.0, 17.0 / 35.0, 12.0 / 35.0,
			-3.0 / 35.0 };

	/**
	 * Savitsky-Golay convolution coefficients for smoothing - Quadratic/Cubic
	 * degree 7
	 */
	static public final double[] SAVITZKY_GOLAY_QC7 = { -2.0 / 21.0, 3.0 / 21.0, 6.0 / 21.0, 7.0 / 21.0, 6.0 / 21.0,
			3.0 / 21.0, -2.0 / 21.0 };

	/**
	 * Savitsky-Golay convolution coefficients for smoothing - Quadratic/Cubic
	 * degree 9
	 */
	static public final double[] SAVITZKY_GOLAY_QC9 = { -21.0 / 231.0, 14.0 / 231.0, 39.0 / 231.0, 54.0 / 231.0,
			59.0 / 231.0, 54.0 / 231.0, 39.0 / 231.0, 14.0 / 231.0, -21.0 / 231.0 };

	/**
	 * Savitsky-Golay convolution coefficients for smoothing - Quartic/Quintic
	 * degree 7
	 */
	static public final double[] SAVITZKY_GOLAY_QQ7 = { 5.0 / 231.0, -30.0 / 231.0, 75.0 / 231.0, 131.0 / 231.0,
			75.0 / 231.0, -30.0 / 231.0, 5.0 / 231.0 };

	/**
	 * Savitsky-Golay convolution coefficients for smoothing - Quartic/Quintic
	 * degree 9
	 */
	static public final double[] SAVITZKY_GOLAY_QQ9 = { 15.0 / 429.0, -55.0 / 429.0, 30.0 / 429.0, 135.0 / 429.0,
			179.0 / 429.0, 135.0 / 429.0, 30.0 / 429.0, -55.0 / 429.0, 15.0 / 429.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Linear/Quadratic degree 3
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_LQ3 = { -1.0 / 2.0, 0.0 / 2.0, 1.0 / 2.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Linear/Quadratic degree 5
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_LQ5 = { -2.0 / 10.0, -1.0 / 10.0, 0.0 / 10.0, 1.0 / 10.0,
			2.0 / 10.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Linear/Quadratic degree 7
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_LQ7 = { -3.0 / 28.0, -2.0 / 28.0, -1.0 / 28.0, 0.0,
			1.0 / 28.0, 2.0 / 28.0, 3.0 / 28.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Linear/Quadratic degree 9
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_LQ9 = { -4.0 / 60.0, -3.0 / 60.0, -2.0 / 60.0, -1.0 / 60.0,
			0.0, 1.0 / 60.0, 2.0 / 60.0, 3.0 / 60.0, 4.0 / 60.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Quartic/Quintic degree 5
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_QC5 = { 1.0 / 12.0, -8.0 / 12.0, 0.0, 8.0 / 12.0,
			-1.0 / 12.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Quartic/Quintic degree 7
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_QC7 = { 22.0 / 252.0, -67.0 / 252.0, -58.0 / 252.0, 0.0,
			58.0 / 252.0, 67.0 / 252.0, -22.0 / 252.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 1st derivative estimation -
	 * Quartic/Quintic degree 9
	 */
	static public final double[] SAVITSKY_GOLAY_1st_DERIV_QC9 = { 86.0 / 1188.0, -142.0 / 1188.0, -193.0 / 1188.0,
			-126.0 / 1188.0, 0.0, 126.0 / 1188.0, 193.0 / 1188.0, 142.0 / 1188.0, -86.0 / 1188.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quadratic/Cubic degree 5
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QC5 = { 2.0 / 7.0, -1.0 / 7.0, -2.0 / 7.0, -1.0 / 7.0,
			2.0 / 7.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quadratic/Cubic degree 7
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QC7 = { 5.0 / 42.0, 0.0, -3.0 / 42.0, -4.0 / 42.0,
			-3.0 / 42.0, 0.0, 5.0 / 42.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quadratic/Cubic degree 9
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QC9 = { 28.0 / 462.0, 7.0 / 462.0, -8.0 / 462.0,
			-17.0 / 462.0, -20.0 / 462.0, -17.0 / 462.0, -8.0 / 462.0, 7.0 / 462.0, 28.0 / 462.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quartic/Quintic degree 5
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QQ5 = { -1.0 / 12.0, 16.0 / 12.0, -30.0 / 12.0, 16.0 / 12.0,
			-1.0 / 12.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quartic/Quintic degree 7
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QQ7 = { -13.0 / 132.0, 67.0 / 132.0, -19.0 / 132.0,
			-70.0 / 132.0, -19.0 / 132.0, 67.0 / 132.0, -13.0 / 132.0 };

	/**
	 * Savitsky-Golay convolution coefficients for 2st derivative estimation -
	 * Quartic/Quintic degree 9
	 */
	static public final double[] SAVITSKY_GOLAY_2nd_DERIV_QQ9 = { -126.0 / 1716.0, 371.0 / 1716.0, 151.0 / 1716.0,
			-211.0 / 1716.0, -370.0 / 1716.0, -211.0 / 1716.0, 151.0 / 1716.0, 371.0 / 1716.0, -126.0 / 1716.0 };

}
