/**
 * 
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * <p>Extracts a subset of the rows and columns associated with the provided
 * base {@link NamedMultivariateJacobianFunction}.  Use with care as this function
 * will evaluate the base {@link NamedMultivariateJacobianFunction} each time
 * the value(...) method is called.<p>
 * 
 * <p>If you want to avoid the cost of repeated evaluations, it is possible to pass a 
 * {@link NamedMultivariateJacobian} (a {@link NamedMultivariateJacobianFunction}
 * evaluated at a specific point) in place of the {@link NamedMultivariateJacobianFunction}
 * so long as the point of evaluation doesn't change.</p>
 * 
 * 
 * @author nicho
 *
 */
public class SubsetNamedMultivariateJacobianFunction extends NamedMultivariateJacobianFunction {

	private final NamedMultivariateJacobianFunction mBase;

	public SubsetNamedMultivariateJacobianFunction( //
			NamedMultivariateJacobianFunction nmvjf, //
			List<? extends Object> inputTags, //
			List<? extends Object> outputTags //
	) throws ArgumentException {
		super(inputTags, outputTags);
		mBase = nmvjf;
		{
			List<? extends Object> inputs = getInputTags();
			for (int i = 0; i < inputs.size(); ++i)
				if (mBase.inputIndex(inputs.get(i)) == -1)
					throw new ArgumentException("The input argument " + inputs.get(i) + " is missing.");
		}
		final int[] bOutIdx = new int[getOutputDimension()];
		{
			List<? extends Object> outputs = getInputTags();
			for (int i = 0; i < bOutIdx.length; ++i)
				if (mBase.outputIndex(outputs.get(i)) == -1)
					throw new ArgumentException("The output argument " + outputs.get(i) + " is missing.");
		}

	}

	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {
		final int[] bInIdx = new int[getInputDimension()];
		{
			List<? extends Object> inputs = getInputTags();
			for (int i = 0; i < bInIdx.length; ++i)
				bInIdx[i] = mBase.inputIndex(inputs.get(i));
		}
		final int[] bOutIdx = new int[getOutputDimension()];
		{
			List<? extends Object> outputs = getInputTags();
			for (int i = 0; i < bOutIdx.length; ++i)
				bOutIdx[i] = mBase.outputIndex(outputs.get(i));
		}
		final double[] vals = new double[bOutIdx.length];
		final double[][] cov = new double[bOutIdx.length][bInIdx.length];
		final Pair<RealVector, RealMatrix> tmp = mBase.value(point);
		final RealVector bVals = tmp.getFirst();
		final RealMatrix bCov = tmp.getSecond();
		for (int r = 0; r < bOutIdx.length; ++r) {
			vals[r] = bVals.getEntry(bOutIdx[r]);
			for (int c = 0; c < bInIdx.length; ++c)
				cov[r][c] = bCov.getEntry(bOutIdx[r], bInIdx[c]);

		}
		return Pair.create( //
				new ArrayRealVector(vals), //
				new Array2DRowRealMatrix(cov) //
		);
	}

}
