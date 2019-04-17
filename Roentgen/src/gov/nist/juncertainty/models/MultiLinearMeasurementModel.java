package gov.nist.juncertainty.models;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>
 * A extension of {@link ExplicitMeasurementModel} for the special
 * case of multiple functions whose outputs are linear combinations of the input
 * vector elements.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public abstract class MultiLinearMeasurementModel<G, H> //
		extends ExplicitMeasurementModel<G, H> //
		implements ILabeledMultivariateFunction<G, H> {

	private final SimplyLazy<RealMatrix> mJacobian = new SimplyLazy<RealMatrix>() {

		@Override
		protected RealMatrix initialize() {
			return buildLinearTransform(getOutputDimension(), getInputDimension());
		}

	};

	/**
	 * Constructs a MultiLinearJacobianFunction
	 *
	 * @param inputPrefix
	 * @param outputPrefix
	 * @throws ArgumentException
	 */
	public MultiLinearMeasurementModel(
			final List<G> inputLabels, //
			final List<H> outputLabels
	) throws ArgumentException {
		super(inputLabels, outputLabels);
	}

	/**
	 * <p>
	 * Implement this function to build a matrix that implements the linear
	 * transform.
	 * </p>
	 * <p>
	 * Each row in the matrix represents a single linear function of the input
	 * point. The width of the matrix is the same as the length of the input point.
	 * </p>
	 *
	 * @param nRows
	 * @param nCols
	 * @return RealMatrix
	 */
	public abstract RealMatrix buildLinearTransform(
			int nRows, int nCols
	);

	/**
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealMatrix jac = mJacobian.get();
		return Pair.create(jac.operate(point), jac);
	}

	@Override
	public RealVector optimized(
			final RealVector point
	) {
		return mJacobian.get().operate(point);
	}
}
