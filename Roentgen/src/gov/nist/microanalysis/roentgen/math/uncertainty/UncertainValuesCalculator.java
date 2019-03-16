package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.utility.FastIndex;

/**
 * A lazy mechanism to perform uncertainty propagation. The
 * {@link UncertainValuesCalculator} class delays calculating the resulting
 * outputs until as late as possible. When a new set of input
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class UncertainValuesCalculator extends UncertainValuesBase {

	public interface ICalculator {

		/**
		 * Give a {@link LabeledMultivariateJacobianFunction} and and {@link UncertainValues}
		 * computes the values and the covariance in some manner or another.
		 * 
		 * @param func
		 * @param point
		 * @return Pair&lt;RealVector, RealMatrix&gt;
		 */
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction func, //
				final RealVector point //
		);

	}

	public static class Jacobian implements ICalculator {

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction func, //
				final RealVector point //
		) {
			return func.evaluate(point);
		}

		@Override
		public String toString() {
			return "Jacobian";
		}
	}

	public static class FiniteDifference implements ICalculator {

		@Override
		public int hashCode() {
			return Objects.hash(mDeltaInputs);
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			FiniteDifference other = (FiniteDifference) obj;
			return Objects.equals(mDeltaInputs, other.mDeltaInputs);
		}

		final RealVector mDeltaInputs;

		public FiniteDifference(final RealVector dinp) {
			mDeltaInputs = dinp;
		}

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction func, //
				final RealVector point //
		) {
			final int inDim = func.getInputDimension(), outDim = func.getOutputDimension();
			assert point.getDimension() == inDim;
			final RealMatrix rm = NullableRealMatrix.build(outDim, inDim);
			for (int c = 0; c < inDim; ++c) {
				final RealVector pt0 = new ArrayRealVector(point), pt1 = new ArrayRealVector(point);
				final double deltaX = Math.abs(mDeltaInputs.getEntry(c));
				pt0.setEntry(c, pt0.getEntry(c) + 0.5 * deltaX);
				pt1.setEntry(c, pt1.getEntry(c) - 0.5 * deltaX);
				final RealVector output0 = func.compute(pt0), output1 = func.compute(pt1);
				for (int r = 0; r < outDim; ++r) {
					final double value = (output0.getEntry(r) - output1.getEntry(r)) / deltaX;
					rm.setEntry(r, c, Double.isNaN(value) ? 0.0 : value);
				}
			}
			return Pair.create(func.compute(point), rm);
		}

		@Override
		public String toString() {
			return "Finite Difference";
		}
	}

	private final LabeledMultivariateJacobianFunction mFunction;
	private UncertainValues mInputs;
	private UncertainValuesBase mRawInputs;
	private ICalculator mCalculator;
	private final boolean mRetainInputs;

	@Override
	public int hashCode() {
		return Objects.hash(mFunction, mInputs, mRawInputs, mCalculator, mRetainInputs);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValuesCalculator other = (UncertainValuesCalculator) obj;
		return Objects.equals(mFunction, other.mFunction) //
				&& Objects.equals(mInputs, other.mInputs) //
				&& Objects.equals(mRawInputs, other.mRawInputs) //
				&& Objects.equals(mCalculator, other.mCalculator) //
				&& (mRetainInputs == other.mRetainInputs);

	}

	private final SimplyLazy<Pair<RealVector, RealMatrix>> mEvaluated = new SimplyLazy<Pair<RealVector, RealMatrix>>() {

		@Override
		protected Pair<RealVector, RealMatrix> initialize() {
			final Pair<RealVector, RealMatrix> tmp = mCalculator.compute(mFunction, mInputs.getValues());
			if (mRetainInputs) {
				// Does not change the column dimension, just adds the inputs to the covariances
				final RealVector rv1 = tmp.getFirst();
				final RealMatrix rm1 = tmp.getSecond();
				final int inpLen = rm1.getColumnDimension();
				assert inpLen == mInputs.getDimension();
				final RealVector rv2 = new ArrayRealVector(inpLen + rv1.getDimension());
				assert inpLen == mInputs.getDimension();
				final RealMatrix rm2 = MatrixUtils.createRealMatrix(rv2.getDimension(), inpLen);
				rv2.setSubVector(0, mInputs.getValues());
				rv2.setSubVector(inpLen, rv1);
				rm2.setSubMatrix(MatrixUtils.createRealIdentityMatrix(inpLen).getData(), 0, 0);
				rm2.setSubMatrix(rm1.getData(), inpLen, 0);
				return Pair.create(rv2, rm2);
			}
			return tmp;

		}

	};

	private final SimplyLazy<UncertainValues> mFullOutputs = new SimplyLazy<UncertainValues>() {

		@Override
		protected UncertainValues initialize() {
			final Pair<RealVector, RealMatrix> eval = mEvaluated.get();
			final RealMatrix jac = eval.getSecond();
			final RealMatrix cov = jac.multiply(mInputs.getCovariances().multiply(jac.transpose()));
			return new UncertainValues(//
					getOutputLabels(), //
					eval.getFirst(), //
					cov);
		}
	};

	transient private SimplyLazy<RealVector> mOutputValues = new SimplyLazy<RealVector>() {
		@Override
		protected RealVector initialize() {
			if (mFullOutputs.initialized())
				return mFullOutputs.get().getValues();
			else {
				final RealVector rv1 = mFunction.compute(mInputs.getValues());
				if (mRetainInputs) {
					final int inpLen = mInputs.getDimension();
					final RealVector rv2 = new ArrayRealVector(inpLen + rv1.getDimension());
					rv2.setSubVector(0, mInputs.getValues());
					rv2.setSubVector(inpLen, rv1);
					return rv2;
				}
				return rv1;
			}
		}
	};

	private final static List<? extends Object> buildOutputs(//
			final LabeledMultivariateJacobianFunction func, //
			final boolean retainInputs //
	) {
		final FastIndex<Object> res = new FastIndex<>();
		if (retainInputs)
			res.addMissing(func.getInputLabels());
		res.addMissing(func.getOutputLabels());
		return res;
	}

	public UncertainValuesCalculator(//
			final LabeledMultivariateJacobianFunction func, //
			final UncertainValuesBase inputs, //
			final boolean retainInputs //
	) throws ArgumentException {
		super(buildOutputs(func, retainInputs));
		mFunction = func;
		if (inputs != null) {
			mRawInputs = inputs;
			mInputs = UncertainValues.force(inputs.reorder(mFunction.getInputLabels()));
		}
		mCalculator = new Jacobian();
		mRetainInputs = retainInputs;
	}

	/**
	 * @throws ArgumentException
	 *
	 */
	public UncertainValuesCalculator(final LabeledMultivariateJacobianFunction func, final UncertainValuesBase inputs) //
			throws ArgumentException {
		this(func, inputs, false);
	}

	public static UncertainValuesCalculator buildDelta(//
			final LabeledMultivariateJacobianFunction func, final UncertainValuesBase inps, final double sigma//
	) throws ArgumentException {
		final UncertainValuesCalculator res = new UncertainValuesCalculator(func, inps, false);
		final RealVector rv = inps.getValues(func.getInputLabels()).mapMultiply(sigma);
		res.setCalculator(new FiniteDifference(rv));
		return res;
	}
	
	public static UncertainValuesCalculator buildDelta(//
			final LabeledMultivariateJacobianFunction func, final UncertainValuesBase inps, final RealVector delta//
	) throws ArgumentException {
		final UncertainValuesCalculator res = new UncertainValuesCalculator(func, inps, false);
		res.setCalculator(new FiniteDifference(delta));
		return res;
	}

	public void setInputs(final UncertainValuesBase inputs) //
			throws ArgumentException {
		synchronized (this) {
			mInputs = UncertainValues.force(inputs.reorder(mFunction.getInputLabels()));
			mEvaluated.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	public void reset() {
		synchronized (this) {
			mEvaluated.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	@Override
	public String toString() {
		return mFunction.toString() + (mInputs != null ? "(" + mInputs.toString() + ")" : "(No inputs)");
	}

	public LabeledMultivariateJacobianFunction getFunction() {
		return mFunction;
	}

	/**
	 * Returns the input {@link UncertainValues} in the correct order for the
	 * associated {@link LabeledMultivariateJacobianFunction}.
	 * 
	 * 
	 * @return {@link UncertainValues}
	 */
	public UncertainValues getInputs() {
		return mInputs;
	}

	/**
	 * Returns the input {@link RealVector} in the correct order for the associated
	 * {@link LabeledMultivariateJacobianFunction}.
	 * 
	 * 
	 * @return {@link RealVector}
	 */
	public RealVector getInputValues() {
		return mInputs.getValues();
	}

	/**
	 * Returns the input {@link UncertainValuesBase} in the form in which it was
	 * passed to the constructor - any order with any extra values.
	 * 
	 * @return {@link UncertainValuesBase}
	 */
	public UncertainValuesBase getRawInputs() {
		return mRawInputs;
	}

	public RealMatrix getJacobian() {
		return mEvaluated.get().getSecond();
	}

	@Override
	public RealVector getValues() {
		return mOutputValues.get();
	}

	@Override
	public RealMatrix getCovariances() {
		return mFullOutputs.get().getCovariances();
	}

	public void setCalculator(final ICalculator calc) {
		synchronized (this) {
			mCalculator = calc;
			mEvaluated.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	public int getInputDimension() {
		return mFunction.getInputDimension();
	}

	public int getOutputDimension() {
		return getDimension();
	}

	public List<? extends Object> getOutputLabels() {
		return getLabels();
	}

	public List<? extends Object> getInputLabels() {
		return mFunction.getInputLabels();
	}

	public double getJacobianEntry(final int r, final int c) {
		return mEvaluated.get().getSecond().getEntry(r, c);
	}

	public double getJacobianEntry(final Object rLbl, final Object cLbl) {
		final int r = mFunction.inputIndex(rLbl);
		assert r >= 0 : //
		"Row " + rLbl + " is missing.";
		final int c = mFunction.outputIndex(cLbl);
		assert c >= 0 : //
		"Column " + cLbl + " is missing.";
		return getJacobianEntry(r, c);
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects. Equivalent to
	 * <code>getOutputValues(uvs, 1.0e-6)</code>
	 *
	 * @param uvs Input {@link UncertainValues}
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues() throws ArgumentException {
		return getOutputValues(1.0e-6);
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects.
	 *
	 * @param uvs Input {@link UncertainValues}
	 * @param tol Tolerance relative to value
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues( //
			final double tol//
	) throws ArgumentException {
		final Pair<RealVector, RealMatrix> pvm = mEvaluated.get();
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<Object, UncertainValue> res = new HashMap<>();
		final List<? extends Object> inLabels = getInputLabels();
		final List<? extends Object> outLabels = getOutputLabels();
		for (int i = 0; i < outLabels.size(); ++i)
			res.put(outLabels.get(i), new UncertainValue(vals.getEntry(i)));
		
		RealMatrix covs = mInputs.getCovariances();
		RealMatrix zeroed = MatrixUtils.createRealMatrix(covs.getRowDimension(), covs.getColumnDimension());
		for (int inIdx = 0; inIdx < inLabels.size(); ++inIdx) {
			final Object inLabel = inLabels.get(inIdx);
			final int idx = mInputs.indexOf(inLabel);
			zeroed.setEntry(idx, idx, covs.getEntry(idx, idx));
			final RealMatrix covLabel = jac.multiply(zeroed).multiply(jac.transpose());
			zeroed.setEntry(idx, idx, 0.0);
			for (int outIdx = 0; outIdx < outLabels.size(); ++outIdx) {
				final UncertainValue val = res.get(outLabels.get(outIdx));
				final double sigma = Math.sqrt(covLabel.getEntry(outIdx, outIdx));
				if (sigma > Math.abs(tol * val.doubleValue()))
					val.assignComponent(inLabel, sigma);
			}
		}
		return res;
	}

	/**
	 * <p>
	 * Get the resulting UncertainValue for each output quantity with sources
	 * grouped and named according to the map of names to a collection of labels.
	 * </p>
	 * </p>
	 * Returns a map from output value to an UncertainValue with discretely broken
	 * out uncertainty components according to the labels map.
	 * </p>
	 *
	 * @param uvs
	 * @param labels Group according to this mapping
	 * @param tol    Minimum size uncertainty to include
	 * @return HashMap&lt;? extends Object, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<? extends Object, UncertainValue> getOutputValues(
			final Map<String, Collection<? extends Object>> labels, //
			final double tol //
	) throws ArgumentException {
		final UncertainValues ordered = mFullOutputs.get();
		final Pair<RealVector, RealMatrix> pvm = mEvaluated.get();
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<Object, UncertainValue> res = new HashMap<>();
		final List<? extends Object> outLabels = getOutputLabels();
		for (int i = 0; i < outLabels.size(); ++i)
			res.put(outLabels.get(i), new UncertainValue(vals.getEntry(i)));
		for (final Map.Entry<String, Collection<? extends Object>> me : labels.entrySet()) {
			final UncertainValues zeroed = UncertainValues.zeroBut(me.getValue(), ordered);
			final RealMatrix covLabel = jac.multiply(zeroed.getCovariances()).multiply(jac.transpose());
			for (int outIdx = 0; outIdx < outLabels.size(); ++outIdx) {
				final UncertainValue val = res.get(outLabels.get(outIdx));
				final double sigma = Math.sqrt(covLabel.getEntry(outIdx, outIdx));
				if (sigma > Math.abs(tol * val.doubleValue()))
					val.assignComponent(me.getKey(), sigma);
			}
		}
		return res;
	}

	/**
	 * Perform the uncertainty calculation using a Monte Carlo evaluation model.
	 * 
	 * @param nEvals
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public UncertainValues propagateMC(int nEvals) throws ArgumentException {
		final MCPropagator mcp = new MCPropagator(this);
		return mcp.compute(nEvals);
	}

}
