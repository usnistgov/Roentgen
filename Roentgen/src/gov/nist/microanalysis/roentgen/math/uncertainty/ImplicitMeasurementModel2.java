/**
 *
 */
package gov.nist.microanalysis.roentgen.math.uncertainty;

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
 * @author Nicholas W. M. Ritchie
 *
 */
public class ImplicitMeasurementModel2<X, Y extends X> //
		extends ExplicitMeasurementModel<X, Y> {

	public static class HLabel {

		private final int mIndex;

		public HLabel(
				final int index
		) {
			mIndex = index;
		}

		@Override
		public String toString() {
			return "H[" + mIndex + "]";
		}

		@Override
		public int hashCode() {
			return Objects.hash(mIndex);
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
	}

	private final ExplicitMeasurementModel<Y, HLabel> mCy;
	private final ExplicitMeasurementModel<X, HLabel> mCx;

	private static <X, Y> List<X> combine(
			final List<? extends X> inLabels, final List<? extends X> outLabels
	) {
		final List<X> res = new ArrayList<>();
		res.addAll(inLabels);
		res.addAll(outLabels);
		return res;
	}

	/**
	 * @param inputLabels
	 * @param outputLabels
	 * @throws ArgumentException
	 */
	public ImplicitMeasurementModel2(
			final ExplicitMeasurementModel<Y, HLabel> cy, //
			final ExplicitMeasurementModel<X, HLabel> cx //
	) throws ArgumentException {
		super(combine(cx.getInputLabels(), cy.getInputLabels()), cy.getInputLabels());
		assert cy.getInputDimension() == cy.getOutputDimension();
		assert cx.getOutputDimension() == cy.getOutputDimension();
		mCx = cx;
		mCy = cy;
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

		final RealVector cyPoint = new ArrayRealVector(mCy.getInputDimension());
		for (int i = 0; i < mCy.getInputDimension(); ++i) {
			final Y inLabel = mCy.getInputLabel(i);
			cyPoint.setEntry(i, getArg(inLabel, point));
		}
		final Pair<RealVector, RealMatrix> cyPr = mCy.evaluate(cyPoint);
		final RealMatrix cyJ = cyPr.getSecond();

		final RealVector cxPoint = new ArrayRealVector(mCy.getInputDimension());
		for (int i = 0; i < mCx.getInputDimension(); ++i) {
			final X inLabel = mCx.getInputLabel(i);
			cxPoint.setEntry(i, getArg(inLabel, point));
		}
		final Pair<RealVector, RealMatrix> cxPr = mCy.evaluate(cxPoint);
		final RealMatrix cxJ = cxPr.getSecond();

		return Pair.create(cyPoint, MatrixUtils.inverse(cyJ).multiply(cxJ));
	}

}
