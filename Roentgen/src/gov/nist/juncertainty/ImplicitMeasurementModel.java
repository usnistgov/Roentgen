package gov.nist.juncertainty;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>
 * Implements a mechanism to handle implicit measurement model steps as a step
 * within an {@link ExplicitMeasurementModel}.
 * </p>
 * <p>
 * Implicit measurement models are measurement models defined by the equation
 * f(X,Y)=0 rather than the more common Y = f(X). Implicit measurement models
 * are more complex because they can not be solved explicitly for Y and it is
 * necessary to perform a non-linear optimization to estimate Y.
 * </p>
 * <p>
 * The ImplicitMeasurementModel class breaks evaluating the model into the
 * computation of three parts - the Jacobian matrix of J<sub>x</sub> and
 * J<sub>Y</sub> and the value of h(X,Y). J<sub>x</sub> equals the Jacobian of
 * the model with respect to the input variables and J<sub>Y</sub> equals the
 * Jacobian of the model with respect to the output variables. Finally h(X,Y) is
 * the actual value of the model which should after non-linear optimization be a
 * zero array.
 * </p>
 * 
 * <p>
 * See JCGM 102 for more details.
 * </p>
 *
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
abstract public class ImplicitMeasurementModel<X, Y extends X> //
		extends ExplicitMeasurementModel<X, Y> {

	public static class HLabel {

		private final int mIndex;

		public HLabel(
				final int index
		) {
			mIndex = index;
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
			final HLabel other = (HLabel) obj;
			return mIndex == other.mIndex;
		}

		@Override
		public int hashCode() {
			return Objects.hash(mIndex);
		}

		@Override
		public String toString() {
			return "H[" + mIndex + "]";
		}
	}

	/**
	 * A {@link List} of the specified size
	 * 
	 * @param size The number of HLabels to build
	 * @return List&lt;HLabel&gt;
	 */
	static public List<HLabel> buildHLabels(
			final int size
	) {
		final ArrayList<HLabel> res = new ArrayList<>();
		for (int i = 0; i < size; ++i)
			res.add(new HLabel(i));
		return res;
	}

	/**
	 * @param inputLabels A {@link List} of input variables
	 * @param outputLabels A {@link List} of output variables
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public ImplicitMeasurementModel(
			final List<X> inputLabels, //
			final List<Y> outputLabels
	) throws ArgumentException {
		super(inputLabels, outputLabels);
	}

	/**
	 * Implement this by computing the Jacobian of the model with respect to the
	 * input variable X.
	 *
	 * @param point A {@link List} of input variables
	 * @return RealMatrix A matrix of dimension getOutputDimension() &times;
	 *         getInputDimension().
	 */
	abstract public RealMatrix computeJx(
			final RealVector point
	);

	/**
	 * Implement this by computing the Jacobian of the model with respect to the
	 * output variable Y.
	 *
	 * @param point c
	 * @return RealMatrix A matrix of dimension getOutputDimension() &times;
	 *         getOutputDimension().
	 */
	abstract public RealMatrix computeJy(
			final RealVector point
	);

	/**
	 * Computes the h(X,Y) for x provided in point and y provided via the
	 * alternative input mechanism. h<sub>i</sub>(X,Y) should be very close to zero.
	 *
	 * @param x The point X at which to evaluate h(X,Y)
	 * @return RealVector The result h(x,Y)
	 */
	abstract public RealVector computeH(
			final RealVector x
	);

	/**
	 * Computes the Jacobian as Cy<sup>-1</sup>&middot;Cx.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealMatrix cy = computeJy(point);
		final RealMatrix cx = computeJx(point);
		// Uses a niave algorithm. See appendix B of JGGM 102 for a better way.
		// Unfortunately we can't implement the more sophisticated model here because
		// we need to calculate the Jacobian to fit within the ExplicitMeasurementModel
		// framework.
		return Pair.create(computeValue(point.toArray()), MatrixUtils.inverse(cy).multiply(cx));
	}

	/**
	 * Implements the optimized computeValue(...) function.
	 * 
	 * @see gov.nist.juncertainty.ExplicitMeasurementModel#computeValue(double[])
	 */
	@Override
	public RealVector computeValue(
			final double[] point
			) {
		final RealVector vals = buildResult();
		for (int i = 0; i < vals.getDimension(); ++i)
			vals.setEntry(i, getArg(getOutputLabel(i), point));
		return vals;
	}

	/**
	 * Builds a {@link RealMatrix} of the correct dimensions to hold the Jacobian
	 * with respect to the input variables. Useful to implement computeJx(...).
	 * 
	 * @return {@link RealMatrix}
	 */
	protected RealMatrix buildEmptyJx() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
	}

	/**
	 * Builds a {@link RealMatrix} of the correct dimensions to hold the Jacobian
	 * with respect to the output variables. Useful to implement computeJy(...).
	 * 
	 * @return {@link RealMatrix}
	 */
	protected RealMatrix buildEmptyJy() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getOutputDimension());
	}

	/**
	 * Computes a {@link RealVector} of the correct dimension to hold the results of
	 * h(X,Y).
	 * 
	 * @return {@link RealVector}
	 */
	protected RealVector buildEmptyH() {
		return new ArrayRealVector(getOutputDimension());
	}

	/**
	 * Set the J<sub>X</sub> Jacobian element at row outputLabel and column
	 * inputLabel.
	 * 
	 * @param outputLabel
	 * @param inputLabel
	 * @param jx          The model Jacobian with respect to the input variables
	 * @param value       The value to assign jx.setEntry(output,input)
	 */
	protected void setJx(
			final int outputLabel, //
			final X inputLabel, //
			final RealMatrix jx, //
			final double value
			) {
		jx.setEntry(outputLabel, inputIndex(inputLabel), value);
	}

	/**
	 * Set the J<sub>Y</sub> Jacobian element at row outputLabel and column
	 * inputLabel.
	 * 
	 * @param hIndex
	 * @param outputLabel
	 * @param jy          The model Jacobian with respect to the output variables
	 * @param value       The value to assign jy.setEntry(hIndex,output)
	 */
	protected void setJy(
			final int hIndex, //
			final Y outputLabel, //
			final RealMatrix jy, //
			final double value
			) {
		jy.setEntry(hIndex, outputIndex(outputLabel), value);
	}

	/**
	 * Sets the value associated with the hIndex in the h vector.
	 * 
	 * @param hIndex
	 * @param res
	 * @param value
	 */
	protected void setH(
			final int hIndex, //
			final RealVector res, //
			final double value
			) {
		res.setEntry(hIndex, value);
	}

}
