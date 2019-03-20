package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.EstimateUncertainValues;
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
public class UncertainValuesCalculator<H, K> extends UncertainValuesBase<K> {

	public interface ICalculator {

		/**
		 * Give a {@link LabeledMultivariateJacobianFunction} and and
		 * {@link UncertainValues} computes the values and the covariance in some manner
		 * or another.
		 *
		 * @param func
		 * @param point
		 * @return Pair&lt;RealVector, RealMatrix&gt;
		 */
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction<?, ?> func, //
				final RealVector point //
		);

	}

	public static class Jacobian implements ICalculator {

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction<?, ?> func, //
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
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			final FiniteDifference other = (FiniteDifference) obj;
			return Objects.equals(mDeltaInputs, other.mDeltaInputs);
		}

		final RealVector mDeltaInputs;

		public FiniteDifference(final RealVector dinp) {
			mDeltaInputs = dinp;
		}

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction<?, ?> func, //
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

	private final LabeledMultivariateJacobianFunction<? extends H, ? extends K> mFunction;
	private UncertainValues<H> mInputs;
	private UncertainValuesBase<? extends H> mRawInputs;
	private ICalculator mCalculator;
	private boolean mRetainInputs;

	@Override
	public int hashCode() {
		return Objects.hash(mFunction, mInputs, mRawInputs, mCalculator);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValuesCalculator<?, ?> other = (UncertainValuesCalculator<?, ?>) obj;
		return Objects.equals(mFunction, other.mFunction) //
				&& Objects.equals(mInputs, other.mInputs) //
				&& Objects.equals(mRawInputs, other.mRawInputs) //
				&& Objects.equals(mCalculator, other.mCalculator);

	}

	private final SimplyLazy<Pair<RealVector, RealMatrix>> mEvaluated = new SimplyLazy<Pair<RealVector, RealMatrix>>() {

		@Override
		protected Pair<RealVector, RealMatrix> initialize() {
			return mCalculator.compute(mFunction, mInputs.getValues());
		}

	};

	private final SimplyLazy<UncertainValues<K>> mFullOutputs = new SimplyLazy<UncertainValues<K>>() {

		@Override
		protected UncertainValues<K> initialize() {
			final Pair<RealVector, RealMatrix> eval = mEvaluated.get();
			RealMatrix jac = eval.getSecond();
			final int extra = mInputs.getDimension();
			assert jac.getColumnDimension() == extra;
			RealVector val = eval.getFirst();
			if (mRetainInputs) {
				final int base = val.getDimension();
				RealVector val2 = new ArrayRealVector(base + extra);
				val2.setSubVector(0, val);
				RealMatrix jac2 = MatrixUtils.createRealMatrix(base + extra, jac.getColumnDimension());
				jac2.setSubMatrix(jac.getData(), 0, 0);
				// The inputs are in the order of mRawInput not mInputs
				for (int i = 0; i < extra; ++i) {
					final int outIdx = indexOf(mInputs.getLabel(i));
					jac2.setEntry(outIdx, i, 1.0);
					val2.setEntry(outIdx, mInputs.getEntry(i));
				}
				jac = jac2;
				val = val2;
			}
			final RealMatrix cov = jac.multiply(mInputs.getCovariances().multiply(jac.transpose()));
			return new UncertainValues<K>(//
					getOutputLabels(), //
					val, //
					cov);
		}
	};

	transient private SimplyLazy<RealVector> mOutputValues = new SimplyLazy<RealVector>() {
		@Override
		protected RealVector initialize() {
			RealVector val;
			if (mFullOutputs.initialized())
				val = mFullOutputs.get().getValues();
			else {
				val = mFunction.compute(mInputs.getValues());
				if (mRetainInputs) {
					RealVector val2 = new ArrayRealVector(val.getDimension() + mInputs.getDimension());
					val2.setSubVector(0, val);
					final int extra = mInputs.getDimension();
					// The inputs are in the order of mRawInput not mInputs
					for (int i = 0; i < extra; ++i) {
						final int outIdx = indexOf(mInputs.getLabel(i));
						val2.setEntry(outIdx, mInputs.getEntry(i));
					}
					val = val2;
				}
			}
			return val;
		}
	};

	private final static <H, K> List<K> buildOutputs1(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func //
	) {
		final FastIndex<K> res = new FastIndex<>();
		res.addAll(func.getOutputLabels());
		return (FastIndex<K>) res;
	}

	private final static <H> List<H> buildOutputs2(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
			final UncertainValuesBase<H> inputs, //
			final boolean retainInputs //
	) {
		final FastIndex<H> res = new FastIndex<>();
		res.addAll(func.getOutputLabels());
		if (retainInputs)
			res.addAll(inputs.getLabels());
		return res;
	}

	/**
	 * It is intentional that both the generic parameters on func are the same type.
	 * This is necessary so that when the inputs and outputs are combined, the
	 * result is compatible with both types.
	 * 
	 * 
	 * @param func
	 * @param inputs
	 * @param retainInputs
	 * @throws ArgumentException
	 */
	@SuppressWarnings("unchecked")
	public UncertainValuesCalculator(//
			final LabeledMultivariateJacobianFunction<? extends K, ? extends K> func, //
			final UncertainValuesBase<K> inputs, //
			final boolean retainInputs //
	) throws ArgumentException {
		super(buildOutputs2(func, inputs, retainInputs));
		mFunction = (LabeledMultivariateJacobianFunction<? extends H, ? extends K>) func;
		if (inputs != null) {
			mRawInputs = (UncertainValuesBase<? extends H>) inputs;
			setInputs((UncertainValuesBase<H>) inputs);
		}
		mCalculator = new Jacobian();
		mRetainInputs = retainInputs;
	}

	public UncertainValuesCalculator(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends K> func, //
			final UncertainValuesBase<H> inputs //
	) throws ArgumentException {
		super(buildOutputs1(func));
		mFunction = func;
		if (inputs != null) {
			mRawInputs = inputs;
			setInputs(inputs);
		}
		mCalculator = new Jacobian();
		mRetainInputs = false;
	}

	public static <J, L> UncertainValuesCalculator<J, L> buildDelta(//
			final LabeledMultivariateJacobianFunction<? extends J, ? extends L> func, //
			final UncertainValuesBase<J> inps, //
			final double sigma, //
			final boolean retainInputs //
	) throws ArgumentException {
		final UncertainValuesCalculator<J, L> res = new UncertainValuesCalculator<J, L>(func, inps);
		final RealVector rv = inps.getValues(func.getInputLabels()).mapMultiply(sigma);
		res.setCalculator(new FiniteDifference(rv));
		return res;
	}

	public static <J, L> UncertainValuesCalculator<J, L> buildDelta(//
			final LabeledMultivariateJacobianFunction<J, L> func, //
			final UncertainValuesBase<J> inps, //
			final RealVector delta //
	) throws ArgumentException {
		final UncertainValuesCalculator<J, L> res = new UncertainValuesCalculator<J, L>(func, inps);
		final RealVector reordered = new ArrayRealVector(delta.getDimension());
		for (int i = 0; i < res.getInputDimension(); ++i)
			reordered.setEntry(i, delta.getEntry(inps.indexOf(res.getInputLabels().get(i))));
		res.setCalculator(new FiniteDifference(reordered));
		return res;
	}

	public void setInputs( //
			final UncertainValuesBase<H> inputs //
	) throws ArgumentException {
		synchronized (this) {
			final List<? extends H> inputLabels = mFunction.getInputLabels();
			final UncertainValuesBase<? extends H> reorder = inputs.reorder(inputLabels);
			mInputs = UncertainValues.force(reorder);
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

	public LabeledMultivariateJacobianFunction<? extends H, ? extends K> getFunction() {
		return mFunction;
	}

	/**
	 * Returns the input {@link UncertainValues} in the correct order for the
	 * associated {@link LabeledMultivariateJacobianFunction}.
	 *
	 *
	 * @return {@link UncertainValues}
	 */
	public UncertainValues<H> getInputs() {
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
	public UncertainValuesBase<? extends H> getRawInputs() {
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

	public List<K> getOutputLabels() {
		return getLabels();
	}

	public List<? extends H> getInputLabels() {
		return mFunction.getInputLabels();
	}

	public double getJacobianEntry(final int r, final int c) {
		return mEvaluated.get().getSecond().getEntry(r, c);
	}

	public double getJacobianEntry(final H rLbl, final H cLbl) {
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
	 * @return HashMap&lt;H, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<K, UncertainValue> getOutputValues() throws ArgumentException {
		return getOutputValues(1.0e-6);
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects.
	 *
	 * @param uvs Input {@link UncertainValues}
	 * @param tol Tolerance relative to value
	 * @return HashMap&lt;H, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<K, UncertainValue> getOutputValues( //
			final double tol //
	) throws ArgumentException {
		final Pair<RealVector, RealMatrix> pvm = mEvaluated.get();
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<K, UncertainValue> res = new HashMap<>();
		final List<? extends H> inLabels = getInputLabels();
		final List<? extends K> outLabels = getOutputLabels();
		for (int i = 0; i < outLabels.size(); ++i)
			res.put(outLabels.get(i), new UncertainValue(vals.getEntry(i)));

		final RealMatrix covs = mInputs.getCovariances();
		final RealMatrix zeroed = MatrixUtils.createRealMatrix(covs.getRowDimension(), covs.getColumnDimension());
		for (int inIdx = 0; inIdx < inLabels.size(); ++inIdx) {
			final H inLabel = inLabels.get(inIdx);
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
	 * @return HashMap&lt;H, UncertainValue&gt;
	 * @throws ArgumentException
	 */
	public HashMap<K, UncertainValue> getOutputValues( //
			final Map<String, Collection<? extends K>> labels, //
			final double tol //
	) throws ArgumentException {
		final UncertainValues<K> ordered = mFullOutputs.get();
		final Pair<RealVector, RealMatrix> pvm = mEvaluated.get();
		final RealVector vals = pvm.getFirst();
		final RealMatrix jac = pvm.getSecond();
		final HashMap<K, UncertainValue> res = new HashMap<>();
		final List<? extends K> outLabels = getOutputLabels();
		for (int i = 0; i < outLabels.size(); ++i)
			res.put(outLabels.get(i), new UncertainValue(vals.getEntry(i)));
		for (final Entry<String, Collection<? extends K>> me : labels.entrySet()) {
			final UncertainValues<K> zeroed = UncertainValues.zeroBut(me.getValue(), ordered);
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
	public EstimateUncertainValues<K> propagateMC(final int nEvals) //
			throws ArgumentException {
		final MCPropagator<H, K> mcp = new MCPropagator<H, K>(this);
		return mcp.compute(nEvals);
	}

}
