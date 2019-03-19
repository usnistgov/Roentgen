package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

/**
 * <p>
 * A extension of {@link LabeledMultivariateJacobianFunction} for the special
 * case of multiple functions whose outputs are linear combinations of the input
 * vector elements.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public abstract class MultiLinearJacobianFunction //
		extends LabeledMultivariateJacobianFunction<String, String> {

	private final SimplyLazy<RealMatrix> mJacobian = new SimplyLazy<RealMatrix>() {

		@Override
		protected RealMatrix initialize() {
			return buildLinearTransform(getInputDimension());
		}

	};

	private static List<String> buildLabels(final String prefix, final int nCh) {
		final List<String> labels = new ArrayList<>();
		for (int i = 0; i < nCh; ++i)
			labels.add(prefix + "[" + Integer.toString(i) + "]");
		return labels;
	}

	/**
	 * Constructs a MultiLinearJacobianFunction
	 *
	 * @param inputPrefix
	 * @param outputPrefix
	 */
	public MultiLinearJacobianFunction(final int nInput, final String inputPrefix, final String outputPrefix) {
		super(buildLabels(inputPrefix, nInput), buildLabels(outputPrefix, nInput));
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
	 * @param nCh
	 * @return RealMatrix
	 */
	public abstract RealMatrix buildLinearTransform(int nCh);

	/**
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealMatrix jac = mJacobian.get();
		return Pair.create(jac.operate(point), jac);
	}
}
