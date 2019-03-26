/**
 *
 */
package gov.nist.juncertainty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
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

	static public List<HLabel> buildHLabels(
			final int size
	) {
		final ArrayList<HLabel> res = new ArrayList<>();
		for (int i = 0; i < size; ++i)
			res.add(new HLabel(i));
		return res;
	}

	private final Map<Y, Double> mOutputValues;

	/**
	 * @param inputLabels
	 * @param outputLabels
	 * @throws ArgumentException
	 */
	public ImplicitMeasurementModel(
			final List<X> inputLabels, //
			final List<Y> outputLabels
	) throws ArgumentException {
		super(inputLabels, outputLabels);
		mOutputValues = new HashMap<>();
	}
	
	abstract public RealMatrix computeCx(final RealVector point);

	abstract public RealMatrix computeCy(final RealVector point);

	public void setOutputValues(Map<Y, Double> outputValues) {
		mOutputValues.putAll(outputValues);
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
		final RealMatrix cy = computeCy(point);
		final RealMatrix cx = computeCx(point);

		final RealVector vals = buildResult();
		for (int i = 0; i < vals.getDimension(); ++i)
			vals.setEntry(i, getArg(getOutputLabel(i), point));

		return Pair.create(vals, MatrixUtils.inverse(cy).multiply(cx));
	}
	
	protected RealMatrix buildEmptyCx() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
	}
	
	protected RealMatrix buildEmptyCy() {
		return MatrixUtils.createRealMatrix(getOutputDimension(), getOutputDimension());
	}
	
	protected double getOutputValue(
			final Y label
	) {
		return mOutputValues.get(label);
	}

	protected void setCx(int hIndex, X xLabel, RealMatrix cx, double value) {
		cx.setEntry(hIndex, inputIndex(xLabel), value);
	}

	protected void setCy(int hIndex, Y yLabel, RealMatrix cy, double value) {
		cy.setEntry(hIndex, outputIndex(yLabel), value);
	}

}
