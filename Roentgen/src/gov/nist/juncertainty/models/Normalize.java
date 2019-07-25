package gov.nist.juncertainty.models;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * Normalize takes the input values and normalizes them to a sum of unity.
 * Doesn't perform special checks for zero or negative entries.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class Normalize<H, K> //
		extends ExplicitMeasurementModel<H, K> {

	/**
	 * @param inputLabels
	 * @throws ArgumentException
	 */
	public Normalize(
			final List<H> inputLabels, //
			final List<K> outputLabels
	) throws ArgumentException {
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
	public RealVector computeValue(
			final double[] point
	) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.length;
		assert getInputDimension() == dim;
		double norm = 0.0;
		final RealVector res = buildResult();
		for (int i = 0; i < dim; ++i) {
			res.setEntry(i, point[i]);
			norm += point[i];
		}
		return res.mapDivideToSelf(norm);
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
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.getDimension();
		assert getInputDimension() == dim;
		double norm = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			norm += point.getEntry(i);
		final RealMatrix j = buildJacobian();
		final RealVector resV = buildResult();
		for (int r = 0; r < dim; ++r) {
			final double a = point.getEntry(r);
			final double na = a / norm;
			final double dnadx = -(na * na) / a;
			final double dnada = (na * (1.0 - na)) / a;
			for (int c = 0; c < dim; ++c)
				j.setEntry(r, c, r == c ? dnada : dnadx);
			resV.setEntry(r, na);
		}
		return Pair.create(resV, j);

	}

}
