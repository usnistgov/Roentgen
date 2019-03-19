package gov.nist.microanalysis.roentgen.math.uncertainty.models;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;

/**
 * COmputes the weighted sum of the input label values times a set of
 * coefficients.
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class WeightedSum<G,H> extends LabeledMultivariateJacobianFunction<G,H> //
		implements ILabeledMultivariateFunction<G,H> {

	private final RealVector mCoeffs;

	public WeightedSum( //
			final List<G> inLabels, //
			final RealVector coeffs, //
			final H outLabel //
	) {
		super(inLabels, Collections.singletonList(outLabel));
		assert inLabels.size() == coeffs.getDimension();
		mCoeffs = coeffs;
	}

	public static <G,H> WeightedSum<G,H> buildSum(List<G> inLabels, H outLabel) {
		final ArrayRealVector coeffs = new ArrayRealVector(inLabels.size());
		coeffs.set(1.0);
		return new WeightedSum<G,H>(inLabels, coeffs, outLabel);
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(1);
		final RealMatrix rm = MatrixUtils.createRealMatrix(1, point.getDimension());
		rm.setRow(0, mCoeffs.toArray());
		rv.setEntry(0, point.dotProduct(mCoeffs));
		return Pair.create(rv, rm);
	}

	@Override
	public RealVector optimized(RealVector point) {
		RealVector rv = new ArrayRealVector(1);
		rv.setEntry(0, point.dotProduct(mCoeffs));
		return rv;

	}
}