package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import org.apache.commons.math3.distribution.MultivariateRealDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.lazy.SimplyLazy;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.Constraint;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;
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
public class UncertainValuesCalculator<H> //
		extends UncertainValuesBase<H> {

	public class Analytical extends JacobianEvaluator<H> {

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
				final RealVector point //
		) {
			return func.evaluate(point);
		}

		@Override
		public String toString() {
			return "Analytical";
		}
	}

	public class FiniteDifference extends JacobianEvaluator<H> {

		final RealVector mDeltaInputs;

		/**
		 * A RealVector containing the delta size used by the finite difference
		 * differentiator. The ordering dinp should
		 * 
		 * @param dinp
		 */
		public FiniteDifference(final RealVector dinp) {
			// dinp is in the order of mRawInputs not mInputs
			final List<? extends H> labels = getInputLabels();
			final RealVector reordered = new ArrayRealVector(labels.size());
			for (int i = 0; i < reordered.getDimension(); ++i)
				reordered.setEntry(i, dinp.getEntry(mRawInputs.indexOf(labels.get(i))));
			mDeltaInputs = reordered;
		}

		@Override
		public Pair<RealVector, RealMatrix> compute( //
				final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
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

	public interface ICalculator<H> {

		public Pair<RealVector, RealMatrix> compute(//
				LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, RealVector point //
		);

		public UncertainValuesBase<H> evaluate(//
				LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
				List<? extends H> outputLabels //
		);

	}

	public class MonteCarlo implements ICalculator<H> {

		private final class MultiEvaluate<J> extends RecursiveAction {

			private static final int MAX_ITERATIONS = 2000;

			private static final long serialVersionUID = 2121417483581902926L;

			private final int mNEvals;

			private MultiEvaluate(final int nEvals //
			) {
				super();
				mNEvals = nEvals;
			}

			@Override
			protected void compute() {
				if (mNEvals <= MAX_ITERATIONS) {
					for (int i = 0; i < mNEvals; ++i)
						evaluate();
				} else
					invokeAll(new MultiEvaluate<J>(mNEvals / 2), new MultiEvaluate<J>(mNEvals - (mNEvals / 2)));
			}
		}

		private MultivariateRealDistribution mDistribution;

		private final int mNEvals;

		private final EstimateUncertainValues<H> mOutputs;

		public MonteCarlo( //
				final int nEvals //
		) {
			this(new SafeMultivariateNormalDistribution(getInputs()), nEvals);
		}

		public MonteCarlo( //
				final MultivariateRealDistribution dist, //
				final int nEvals //
		) {
			mDistribution = dist;
			mNEvals = nEvals;
			mOutputs = new EstimateUncertainValues<H>(getOutputLabels());
		}

		@Override
		public UncertainValuesBase<H> evaluate( //
				final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
				final List<? extends H> outputLabels //
		) {
			assert (mDistribution != null);
			final ForkJoinPool fjp = new ForkJoinPool();
			fjp.invoke(new MultiEvaluate<H>(mNEvals));
			return mOutputs;
		}

		/**
		 * Perform one evaluation of the {@link LabeledMultivariateJacobianFunction}
		 */
		private void evaluate() {
			final RealVector pt = MatrixUtils.createRealVector(mDistribution.sample());
			assert !pt.isNaN();
			{
				final UncertainValues<H> inputs = mInputs.get();
				for (Map.Entry<H, Constraint> con : mConstraints.entrySet()) {
					final int idx = inputs.indexOf(con.getKey());
					if (idx != -1)
						pt.setEntry(idx, con.getValue().limit(pt.getEntry(idx)));
				}
			}
			final int base = mFunction.getOutputDimension();
			final int extra = mFunction.getInputDimension();
			final RealVector val2 = new ArrayRealVector(base + extra);
			val2.setSubVector(0, mFunction.compute(pt));
			// The inputs are in the order of mRawInput not mInputs
			final List<? extends H> inputs = getInputLabels();
			final List<H> outputLabels = getOutputLabels();
			for (int i = 0; i < extra; ++i) {
				final int outIdx = outputLabels.indexOf(inputs.get(i));
				val2.setEntry(outIdx, pt.getEntry(i));
			}
			mOutputs.add(val2);
		}

		@Override
		public Pair<RealVector, RealMatrix> compute(//
				final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
				final RealVector point //
		) {
			// The Monte Carlo method does not compute the Jacobian
			return Pair.create(mFunction.compute(point), (RealMatrix) null);
		}

		public String toString() {
			return "Monte Carlo[evals = " + mNEvals + "]";
		}

	}

	private interface ICompute<H> {

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
				final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
				final RealVector point //
		);

	}

	private abstract class JacobianEvaluator<J> //
			implements ICalculator<J>, ICompute<J> {

		public JacobianEvaluator() {
		}

		@Override
		public UncertainValues<J> evaluate( //
				final LabeledMultivariateJacobianFunction<? extends J, ? extends J> func, //
				final List<? extends J> outputLabels //
		) {
			final Pair<RealVector, RealMatrix> eval = mEvaluation.get();
			final RealMatrix jac = eval.getSecond();
			UncertainValues<H> inputs = mInputs.get();
			final int extra = inputs.getDimension();
			assert jac.getColumnDimension() == extra;
			final RealVector val = eval.getFirst();
			final int base = val.getDimension();
			assert base == func.getOutputDimension();
			final RealVector val2 = new ArrayRealVector(getOutputDimension());
			val2.setSubVector(0, val);
			final RealMatrix jac2 = MatrixUtils.createRealMatrix(getOutputDimension(), jac.getColumnDimension());
			jac2.setSubMatrix(jac.getData(), 0, 0);
			// The inputs are in the order of mRawInput not mInputs
			for (int i = 0; i < extra; ++i) {
				final int outIdx = indexOf(inputs.getLabel(i));
				assert outIdx < jac2.getRowDimension() : //
				"Oops! " + outIdx + " for " + inputs.getLabel(i);
				val2.setEntry(outIdx, inputs.getEntry(i));
				jac2.setEntry(outIdx, i, 1.0);
			}
			final RealMatrix cov2 = jac2.multiply(inputs.getCovariances().multiply(jac2.transpose()));
			return new UncertainValues<J>(//
					outputLabels, //
					val2, //
					cov2);
		}
	}

	private final static <H> List<H> buildOutputs1(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
			final UncertainValuesBase<H> inputs) {
		final FastIndex<H> res = new FastIndex<>();
		res.addAll(func.getOutputLabels());
		res.addMissing(inputs.getLabels());
		return res;
	}

	private final LabeledMultivariateJacobianFunction<? extends H, ? extends H> mFunction;

	private final SimplyLazy<UncertainValues<H>> mInputs = new SimplyLazy<UncertainValues<H>>() {

		@Override
		protected UncertainValues<H> initialize() {
			final List<? extends H> inputLabels = mFunction.getInputLabels();
			final UncertainValuesBase<H> reorder = mRawInputs.reorder(inputLabels);
			final UncertainValues<H> res = UncertainValues.asUncertainValues(reorder).copy();
			for (final Map.Entry<H, Constraint> me : mConstraints.entrySet()) {
				if (res.hasLabel(me.getKey())) {
					final double val = res.getEntry(me.getKey());
					final double con = me.getValue().limit(val);
					if (con != val)
						res.set(me.getKey(), con);
				}
			}
			return res;
		}

	};

	private final Map<H, Constraint> mConstraints;

	private UncertainValuesBase<H> mRawInputs;

	private ICalculator<H> mCalculator;

	private final SimplyLazy<Pair<RealVector, RealMatrix>> mEvaluation = new SimplyLazy<Pair<RealVector, RealMatrix>>() {

		@Override
		protected Pair<RealVector, RealMatrix> initialize() {
			return mCalculator.compute(mFunction, mInputs.get().getValues());
		}
	};

	private final SimplyLazy<UncertainValuesBase<H>> mFullOutputs = new SimplyLazy<UncertainValuesBase<H>>() {

		@Override
		protected UncertainValuesBase<H> initialize() {
			final UncertainValuesBase<H> res = mCalculator.evaluate(mFunction, getOutputLabels());
			return res;
		}
	};

	transient private SimplyLazy<RealVector> mOutputValues = new SimplyLazy<RealVector>() {
		@Override
		protected RealVector initialize() {
			RealVector val;
			if (mFullOutputs.initialized())
				val = mFullOutputs.get().getValues();
			else {
				final UncertainValues<H> inputs = mInputs.get();
				val = mFunction.compute(inputs.getValues());
				final RealVector val2 = new ArrayRealVector(val.getDimension() + inputs.getDimension());
				val2.setSubVector(0, val);
				final int extra = inputs.getDimension();
				// The inputs are in the order of mRawInput not mInputs
				for (int i = 0; i < extra; ++i) {
					final int outIdx = indexOf(inputs.getLabel(i));
					val2.setEntry(outIdx, inputs.getEntry(i));
				}
				return val2;
			}
			return val;
		}
	};

	public UncertainValuesCalculator(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends H> func, //
			final UncertainValuesBase<H> inputs //
	) throws ArgumentException {
		super(buildOutputs1(func, inputs));
		mFunction = func;
		if (inputs != null) {
			mRawInputs = inputs;
			setInputs(inputs);
		}
		mCalculator = new Analytical();
		mConstraints = new HashMap<>();
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UncertainValuesCalculator<?> other = (UncertainValuesCalculator<?>) obj;
		return Objects.equals(mFunction, other.mFunction) //
				&& Objects.equals(mInputs, other.mInputs) //
				&& Objects.equals(mRawInputs, other.mRawInputs) //
				&& Objects.equals(mCalculator, other.mCalculator);

	}

	@Override
	public RealMatrix getCovariances() {
		return mFullOutputs.get().getCovariances();
	}

	public LabeledMultivariateJacobianFunction<? extends H, ? extends H> getFunction() {
		return mFunction;
	}

	public int getInputDimension() {
		return mFunction.getInputDimension();
	}

	public List<? extends H> getInputLabels() {
		return mFunction.getInputLabels();
	}

	/**
	 * Returns the input {@link UncertainValues} in the correct order for the
	 * associated {@link LabeledMultivariateJacobianFunction}. If there are any
	 * Constraints imposed on this {@link UncertainValuesCalculator} then the
	 * returned values will reflect the constraints.
	 *
	 * @return {@link UncertainValues}
	 */
	public UncertainValues<H> getInputs() {
		return mInputs.get();
	}

	/**
	 * Returns the input {@link RealVector} in the correct order for the associated
	 * {@link LabeledMultivariateJacobianFunction}.
	 *
	 *
	 * @return {@link RealVector}
	 */
	public RealVector getInputValues() {
		return mInputs.get().getValues();
	}

	public Optional<Jacobian<H, H>> getJacobian() {
		final RealMatrix jac = mEvaluation.get().getSecond();
		if (jac != null)
			return Optional.of(new Jacobian<H, H>(getInputLabels(), mFunction.getOutputLabels(), jac));
		else
			return Optional.empty();
	}

	public int getOutputDimension() {
		return getDimension();
	}

	public List<H> getOutputLabels() {
		return getLabels();
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
	public HashMap<H, UncertainValueEx<H>> getOutputValues() //
			throws ArgumentException {
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
	public HashMap<H, UncertainValueEx<H>> getOutputValues( //
			final double tol //
	) throws ArgumentException {
		final UncertainValuesBase<H> pvm = mFullOutputs.get();
		final HashMap<H, UncertainValueEx<H>> res = new HashMap<>();
		for (final H outLbl : getOutputLabels()) {
			final Map<H, Double> vars = new HashMap<>();
			for (final H inLbl : getInputLabels()) {
				if (pvm.getCovariance(inLbl, outLbl) > tol)
					vars.put(inLbl, pvm.getCovariance(inLbl, outLbl));
			}
			final int outIdx = pvm.indexOf(outLbl);
			res.put(outLbl, new UncertainValueEx<H>(pvm.getEntry(outIdx), pvm.getUncertainty(outIdx), vars));
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
	public HashMap<H, UncertainValueEx<String>> getOutputValues( //
			final Map<String, Collection<? extends H>> labels, //
			final double tol //
	) throws ArgumentException {
		final UncertainValuesBase<H> pvm = mFullOutputs.get();
		final HashMap<H, UncertainValueEx<String>> res = new HashMap<>();
		for (final H outLbl : getOutputLabels()) {
			final Map<String, Double> vars = new HashMap<>();
			for (final Map.Entry<String, Collection<? extends H>> me : labels.entrySet()) {
				final String varLbl = me.getKey();
				double covSum = 0.0;
				for (final H inLbl : me.getValue())
					covSum += pvm.getCovariance(inLbl, outLbl);
				vars.put(varLbl, covSum);
			}
			final int outIdx = pvm.indexOf(outLbl);
			res.put(outLbl, new UncertainValueEx<String>(pvm.getEntry(outIdx), pvm.getUncertainty(outIdx), vars));
		}
		return res;
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

	public UncertainValuesBase<H> getUncertainValues() {
		return mFullOutputs.get();
	}

	@Override
	public RealVector getValues() {
		return mOutputValues.get();
	}

	/**
	 * Forces the recalculation of all lazy evaluated quantities
	 */
	private void reset() {
		synchronized (this) {
			mInputs.reset();
			mEvaluation.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	public void setCalculator(final ICalculator<H> calc) {
		synchronized (this) {
			mCalculator = calc;
			mEvaluation.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	/**
	 * Sets the input value and covariance matrix for the calculation. The input
	 * {@link UncertainValuesBase} is evaluated, reordered and any
	 * {@link Constraint}s are applied.
	 *
	 *
	 * @param inputs
	 * @throws ArgumentException
	 */
	public void setInputs( //
			final UncertainValuesBase<H> inputs //
	) throws ArgumentException {
		final List<H> missing = inputs.missing(getInputLabels());
		if (missing.size() > 0)
			throw new ArgumentException("The input values " + missing + " are not available in setInputs(...).");
		synchronized (this) {
			mRawInputs = inputs;
			reset();
		}
	}

	/**
	 * Forces the otherwise lazy evaluation of the calculation.
	 *
	 * @return this
	 */
	public UncertainValuesCalculator<H> force() {
		mFullOutputs.get();
		return this;
	}

	public void setConstraint(final H label, final Constraint constraint) {
		mConstraints.put(label, constraint);
		reset();
	}

	public void addConstraints(Map<H, Constraint> constraints) {
		mConstraints.putAll(constraints);
		reset();
	}

	public void clearConstraints() {
		mConstraints.clear();
		reset();
	}

	@Override
	public String toString() {
		return mFunction.toString() + (mInputs != null ? "(" + mInputs.toString() + ")" : "(No inputs)");
	}
}