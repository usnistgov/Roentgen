package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

/**
 * Normalize takes the input values and normalizes them to a sum of unity.
 * Negative values are truncated to zero. If the resulting sum of the truncated
 * values is zero, then the inputs are returned and a unity Jacobian matrix.
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class Normalize //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	public static class Normalized extends BaseLabel<Object, Object, Object> {

		private Normalized(Object obj) {
			super("N", obj);
		}
	}

	public final static List<Object> buildNormalized(List<? extends Object> inputs) {
		List<Object> res = new ArrayList<>();
		for (Object lbl : inputs)
			res.add(buildNormalized(lbl));
		return res;
	}

	public static Normalized buildNormalized(Object lbl) {
		if (lbl instanceof BaseLabel)
			return new Normalized((BaseLabel<?, ?, ?>) lbl);
		else
			return new Normalized(lbl);
	}

	/**
	 * @param inputLabels
	 */
	public Normalize(List<? extends Object> inputLabels) {
		super(inputLabels, buildNormalized(inputLabels));
	}

	/**
	 * <p>
	 * If obj is an instance of Normalized then unwrap returns the object inside of
	 * Normalized (getObject1()). Otherwise unwrap just returns obj.
	 * </p>
	 * <p>
	 * Useful for functions that can take either normalized or unnormalized
	 * equivalents of tags.
	 * </p>
	 * 
	 * @param obj
	 * @return
	 */
	public static Object unwrap(Object obj) {
		if (obj instanceof Normalized)
			return ((Normalized) obj).getObject1();
		else
			return obj;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.getDimension();
		assert getInputDimension() == dim;
		double norm = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			norm += Math.max(0.0, point.getEntry(i));
		if (norm == 0.0) {
			return Pair.create(point.copy(), MatrixUtils.createRealIdentityMatrix(dim));
		} else {
			RealMatrix j = MatrixUtils.createRealMatrix(dim, dim);
			RealVector resV = new ArrayRealVector(dim);
			for (int r = 0; r < dim; ++r) {
				final double a = point.getEntry(r);
				if (a > 0.0) {
					final double dnadx = -a / (norm * norm);
					final double dnada = 1.0 / norm + dnadx;
					for (int c = 0; c < dim; ++c)
						j.setEntry(r, c, r == c ? dnada : dnadx);
				}
				resV.setEntry(r, a / norm);
			}
			return Pair.create(resV, j);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction
	 * #optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(RealVector point) {
		assert getInputDimension() == getOutputDimension();
		final int dim = point.getDimension();
		assert getInputDimension() == dim;
		double norm = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			norm += Math.max(0.0, point.getEntry(i));
		return point.mapDivide(norm);
	}
}
