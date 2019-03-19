package gov.nist.microanalysis.roentgen.math.uncertainty.models;

import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;

/**
 * Normalize takes the input values and normalizes them to a sum of unity.
 * Negative values are truncated to zero. If the resulting sum of the truncated
 * values is zero, then the inputs are returned and a unity Jacobian matrix.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class Normalize<H, K> //
		extends LabeledMultivariateJacobianFunction<H, K> //
		implements ILabeledMultivariateFunction<H, K> {

	/**
	 * @param inputLabels
	 */
	public Normalize(final List<H> inputLabels, final List<K> outputLabels) {
		super(inputLabels, outputLabels);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction
	 * #optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(final RealVector point) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.getDimension();
		assert getInputDimension() == dim;
		double norm = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			norm += Math.max(0.0, point.getEntry(i));
		return point.mapDivide(norm);
	}

	@Override
	public String toString() {
		return "Normalize" + getInputLabels();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.getDimension();
		assert getInputDimension() == dim;
		double norm = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			norm += Math.max(0.0, point.getEntry(i));
		if (norm == 0.0) {
			return Pair.create(point.copy(), MatrixUtils.createRealIdentityMatrix(dim));
		} else {
			final RealMatrix j = MatrixUtils.createRealMatrix(dim, dim);
			final RealVector resV = new ArrayRealVector(dim);
			for (int r = 0; r < dim; ++r) {
				final double a = point.getEntry(r);
				final double n = a / norm;
				if (a > 0.0) {
					final double dnadx = -(n * n) / a;
					final double dnada = (n * (1.0 - n)) / a;
					for (int c = 0; c < dim; ++c)
						j.setEntry(r, c, r == c ? dnada : dnadx);
				}
				resV.setEntry(r, a / norm);
			}
			return Pair.create(resV, j);
		}
	}

}
