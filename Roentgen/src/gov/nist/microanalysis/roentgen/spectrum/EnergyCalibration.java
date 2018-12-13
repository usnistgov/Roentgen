package gov.nist.microanalysis.roentgen.spectrum;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * Implements various different models mapping channel number into energy bins.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 232 $
 */
abstract public class EnergyCalibration extends LabeledMultivariateJacobianFunction implements IToHTML {

	private static List<Object> buildParams(final String prefix, final int nParams) {
		final List<Object> res = new ArrayList<>();
		for (int i = 0; i < nParams; ++i)
			res.add(prefix + "[" + Integer.toString(i) + "]");
		return res;
	}

	protected EnergyCalibration(final List<Object> inputParams, final int nCh) {
		super(inputParams, buildParams("E", nCh));
	}

	/**
	 * Compute the energy for the specified channel.
	 *
	 * @param ch
	 * @param params
	 * @return double
	 */
	abstract public double compute(double channel);

	public double averageEnergyForChannel(final int ch) {
		return 0.5 * (compute(ch) + compute(ch + 1));
	}

	public double minEnergyForChannel(final int ch) {
		return compute(ch);
	}

	public double maxEnergyForChannel(final int ch) {
		return compute(ch + 1);
	}

	public double channelWidth(final int ch) {
		return compute(ch + 1) - compute(ch);
	}

	public int channelIndex(final double energy) {
		return (int) Math.floor(invert(energy));
	}

	public RealVector range(final int low, final int high) {
		final RealVector res = new ArrayRealVector(high - low);
		for (int ch = low; ch < high; ++ch)
			res.setEntry(ch, compute(ch));
		return res;
	}

	/**
	 * Compute the channel for the specified energy where the returned channel may
	 * be non-integral.
	 *
	 * @param ch
	 * @param params
	 * @return double
	 */
	abstract public double invert(double energy);

	public static Polynomial Linear(final double zeroOffset, final double gain, final int nCh) {
		return new Polynomial(new double[] { zeroOffset, gain }, nCh);
	}

	public static Polynomial Quadratic(final double zeroOffset, final double gain, final double quad, final int nCh) {
		return new Polynomial(new double[] { zeroOffset, gain, quad }, nCh);
	}

	static public class Polynomial extends EnergyCalibration {

		private static final String[] ORDER_NAMES = new String[] { "Offset", "Linear", "Quadratic", "Cubic", "Quartic",
				"Quintic", "Sextic", "Septic", "Octic" };

		private final LaguerreSolver mSolver;
		private final PolynomialFunction mPolynomial;

		public Polynomial(final double[] params, final int nCh) {
			super(Arrays.asList((Object[]) Arrays.copyOf(ORDER_NAMES, params.length)), nCh);
			assert params.length >= 2 : "EnergyScaleFunction.Polynomial order must be >1.";
			assert params.length < 9 : "EnergyScaleFunction.Polynomial order must be 9.";
			mPolynomial = new PolynomialFunction(params);
			mSolver = new LaguerreSolver();
		}

		@Override
		public double compute(final double ch) {
			return mPolynomial.value(ch);
		}

		@Override
		public double invert(final double energy) {
			final double[] params = mPolynomial.getCoefficients();
			final int dim = params.length;
			if (dim == 2)
				return (energy - params[0]) / params[1];
			else if (params.length == 3) {
				// A la Press et al final
				final double coeff0 = params[0] - energy;
				final double q = params[1]
						+ (Math.signum(params[1]) * Math.sqrt((params[1] * params[1]) - (4.0 * coeff0 * params[2])));
				return (-2.0 * coeff0) / q;
			} else {
				// Since we know the calibration is essentially linear with
				// higher order corrections, we'll estimate the answer assuming
				// linearity.
				final double[] coeffs = params.clone();
				coeffs[0] -= energy;
				final double estimate = -coeffs[0] / coeffs[1];
				final PolynomialFunction pf = new PolynomialFunction(coeffs);
				return mSolver.solve(100, pf, estimate - 100.0, estimate + 100.0);
			}
		}

		@Override
		public int hashCode() {
			return Objects.hash(mPolynomial) ^ super.hashCode();
		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			if (!super.equals(obj))
				return false;
			final Polynomial other = (Polynomial) obj;
			if (!mPolynomial.equals(other.mPolynomial))
				return false;
			return true;
		}

		/**
		 * @param point
		 * @return
		 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
		 */
		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector vals = new ArrayRealVector(point.getDimension());
			final RealMatrix cov = new Array2DRowRealMatrix(point.getDimension(), mPolynomial.degree());
			for (int i = 0; i < point.getDimension(); ++i) {
				final double ch = point.getEntry(i);
				vals.setEntry(i, compute(ch));
				double x = 1.0;
				for (int j = 0; j < mPolynomial.degree(); ++j) {
					cov.setEntry(i, j, x);
					x *= ch;
				}
			}
			return Pair.create(vals, cov);
		}

		@Override
		public String toHTML(final Mode mode) {
			final BasicNumberFormat bnf = new BasicNumberFormat();
			final StringBuffer sb = new StringBuffer();
			final double[] coeff = mPolynomial.getCoefficients();
			for (int i = coeff.length - 1; i >= 0; i--) {
				sb.append(bnf.formatHTML(coeff[i]));
				if (i > 0) {
					sb.append("&middot;" + HTML.i("ch"));
					if (i > 1)
						sb.append(HTML.sup(Integer.toString(i)));
					sb.append(" + ");
				}
			}
			sb.append(" eV");
			if (mode != Mode.TERSE)
				sb.append(" where " + HTML.i("ch") + " is the channel index");
			return sb.toString();
		}
	}

	static public class QuadraticInSqrt extends EnergyCalibration implements IToHTML {

		private static final double MIN_B = 0.001 / Math.sqrt(2000);

		private final double[] mParameters;

		private static List<Object> inputTags() {
			return Arrays.asList(new Object[] { "Offset", "Sqrt", "Linear" });
		}

		public QuadraticInSqrt(final double offset, final double sqrt, final double linear, final int nChannels) {
			super(inputTags(), nChannels);
			mParameters = new double[] { offset, sqrt, linear };
		}

		@Override
		public double compute(final double ch) {
			assert ch >= 0.0;
			final double a = mParameters[2];
			final double b = mParameters[1];
			final double c = mParameters[0];
			return (a * ch) + (b * Math.sqrt(ch)) + c;
		}

		@Override
		public double invert(final double energy) {
			assert energy >= 0;
			double ch;
			final double a = mParameters[2];
			final double b = mParameters[1];
			final double c = mParameters[0];
			final double de = energy - c;
			if (Math.abs(b) > MIN_B)
				// (b^2 - 2 a c + 2 a e0 - Sqrt[b^4 - 4 a b^2 c + 4 a b^2
				// e0])/(2
				// a^2)
				ch = (((b * b) + (2.0 * a * de)) - Math.sqrt(b * b * ((b * b) + (4.0 * a * de)))) / (2.0 * a * a);
			else
				ch = de / a;
			assert Math.abs(compute(ch) - energy) < 0.1;
			return ch;

		}

		/**
		 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
		 */
		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector vals = new ArrayRealVector(point.getDimension());
			final RealMatrix cov = new Array2DRowRealMatrix(point.getDimension(), 3);
			for (int i = 0; i < point.getDimension(); ++i) {
				final double ch = point.getEntry(i);
				vals.setEntry(i, compute(ch));
				cov.setEntry(i, 0, 1);
				cov.setEntry(i, 1, Math.sqrt(ch));
				cov.setEntry(i, 2, ch);
			}
			return Pair.create(vals, cov);
		}

		@Override
		public String toHTML(final Mode terse) {
			final Number a = mParameters[1];
			final Number b = mParameters[0];
			final Number c = mParameters[2];
			final BasicNumberFormat bnf = new BasicNumberFormat();
			return bnf.formatHTML(a) + "&middot;i + " + bnf.formatHTML(b) + "&middot;i" + HTML.sup("1/2") + " "
					+ bnf.formatHTML(c);
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = super.hashCode();
			result = prime * result + Arrays.hashCode(mParameters);
			return result;
		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (!super.equals(obj))
				return false;
			if (getClass() != obj.getClass())
				return false;
			if (!super.equals(obj))
				return false;
			final QuadraticInSqrt other = (QuadraticInSqrt) obj;
			if (!Arrays.equals(mParameters, other.mParameters))
				return false;
			return true;
		}
	}
}