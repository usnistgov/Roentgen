package gov.nist.juncertainty;

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

import gov.nist.juncertainty.utility.FastIndex;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.Constraint;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;

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

	public interface ICalculator<H> {

		public Pair<RealVector, RealMatrix> compute(
				ExplicitMeasurementModel<? extends H, ? extends H> func, //
				RealVector point //
		);

		public UncertainValuesBase<H> evaluate(
				ExplicitMeasurementModel<? extends H, ? extends H> func, //
				List<? extends H> outputLabels //
		);
	}

	public class Analytical extends JacobianEvaluator<H> {

		@Override
		public Pair<RealVector, RealMatrix> compute(
				final ExplicitMeasurementModel<? extends H, ? extends H> func, //
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
		 * dinp is a RealVector containing the delta size used by the finite difference
		 * differentiator. The ordering dinp should be the same as
		 * {@link UncertainValuesCalculator}.getInputValues().
		 *
		 * @param dinp The RealVector with differences
		 */
		public FiniteDifference(
				final RealVector dinp
		) {
			mDeltaInputs = dinp;
		}

		/**
		 * Builds a FiniteDifference evaluator that uses
		 * getInputValues().mapMultiply(frac). Be careful as an zeros in
		 * getInputValues() may lead to mischief.
		 *
		 * @param frac The fractional difference size
		 */
		public FiniteDifference(
				final double frac
		) {
			mDeltaInputs = getInputValues().mapMultiply(frac);
		}

		@Override
		public Pair<RealVector, RealMatrix> compute(
				final ExplicitMeasurementModel<? extends H, ? extends H> func, //
				final RealVector point //
		) {
			final int inDim = func.getInputDimension(), outDim = func.getOutputDimension();
			assert point.getDimension() == inDim;
			final RealMatrix rm = MatrixUtils.createRealMatrix(outDim, inDim);
			for (int c = 0; c < inDim; ++c) {
				final RealVector pt0 = new ArrayRealVector(point), pt1 = new ArrayRealVector(point);
				final double deltaX = Math.abs(mDeltaInputs.getEntry(c));
				pt0.setEntry(c, pt0.getEntry(c) + 0.5 * deltaX);
				pt1.setEntry(c, pt1.getEntry(c) - 0.5 * deltaX);
				final RealVector output0 = func.compute(pt0), output1 = func.compute(pt1);
				for (int r = 0; r < outDim; ++r) {
					final double value = (output0.getEntry(r) - output1.getEntry(r)) / deltaX;
					rm.setEntry(r, c, value);
				}
			}
			return Pair.create(func.compute(point), rm);
		}

		@Override
		public String toString() {
			return "Finite Difference";
		}
	}

	public class MonteCarlo implements ICalculator<H> {

		private final class MultiEvaluate<J> extends RecursiveAction {

			private static final int MAX_ITERATIONS = 2000;

			private static final long serialVersionUID = 2121417483581902926L;

			private final int mNEvals;

			private MultiEvaluate(
					final int nEvals //
			) {
				super();
				mNEvals = nEvals;
			}

			@Override
			protected void compute() {
				assert mParallel;
				if (mNEvals <= MAX_ITERATIONS) {
					for (int i = 0; i < mNEvals; ++i)
						performOneEvaluation();
				} else
					invokeAll(new MultiEvaluate<J>(mNEvals / 2), new MultiEvaluate<J>(mNEvals - (mNEvals / 2)));
			}
		}

		private MultivariateRealDistribution mDistribution;

		private final int mNEvals;

		private boolean mParallel;

		private final EstimateUncertainValues<H> mOutputs;

		/**
		 * Constructs a MonteCarlo ICalculator object.
		 * 
		 * @param nEvals The number of evaluations to perform
		 */
		public MonteCarlo(
				final int nEvals //
		) {
			this(nEvals, true);
		}

		/**
		 * Builds a Monte Carlo evaluator around the an instance of
		 * {@link SafeMultivariateNormalDistribution}.
		 *
		 * @param nEvals   Number of evaluations to execute
		 * @param parallel Perform in parallel or in series.
		 */
		public MonteCarlo(
				final int nEvals, final boolean parallel //
		) {
			this(new SafeMultivariateNormalDistribution(getInputs()), nEvals, parallel);
		}

		/**
		 * Builds a Monte Carlo evaluator around the specified uncertainty distribution
		 * which is presumably representative of getInputs().
		 *
		 * @param dist     A {@link MultivariateRealDistribution} in the order of
		 *                 {@link UncertainValuesCalculator}.getInputs().
		 * @param nEvals   Number of evaluations to execute
		 * @param parallel Perform in parallel or sequentially.
		 */
		public MonteCarlo(
				final MultivariateRealDistribution dist, //
				final int nEvals, //
				final boolean parallel //
		) {
			mDistribution = dist;
			mNEvals = nEvals;
			mOutputs = new EstimateUncertainValues<H>(getOutputLabels());
			mParallel = parallel;
		}

		@Override
		public UncertainValuesBase<H> evaluate(
				final ExplicitMeasurementModel<? extends H, ? extends H> func, //
				final List<? extends H> outputLabels //
		) {
			assert (mDistribution != null);
			if (mParallel) {
				final ForkJoinPool fjp = new ForkJoinPool();
				fjp.invoke(new MultiEvaluate<H>(mNEvals));
			} else {
				for (int i = 0; i < mNEvals; ++i)
					performOneEvaluation();
			}
			return mOutputs;
		}

		/**
		 * Perform one evaluation of the {@link ExplicitMeasurementModel}
		 */
		private void performOneEvaluation() {
			// pt is in the correct order for mFunction
			final RealVector pt = MatrixUtils.createRealVector(mDistribution.sample());
			assert !pt.isNaN();
			final int base = mFunction.getOutputDimension();
			final RealVector val2 = new ArrayRealVector(getOutputDimension());
			// Constraints are applied to pt in compute(...)
			val2.setSubVector(0, mFunction.compute(pt));
			val2.setSubVector(base, pt);
			mOutputs.add(val2);
		}

		@Override
		public Pair<RealVector, RealMatrix> compute(
				final ExplicitMeasurementModel<? extends H, ? extends H> func, //
				final RealVector point //
		) {
			// The Monte Carlo method does not compute the Jacobian
			return Pair.create(mFunction.compute(point), (RealMatrix) null);
		}

		@Override
		public String toString() {
			return "Monte Carlo[evals = " + mNEvals + "]";
		}

		/**
		 * Are the evaluations performed in a parallel multi-threaded style or in
		 * sequence
		 * 
		 * @return true in parallel, false otherwise
		 */
		public boolean isParallel() {
			return mParallel;
		}

		/**
		 * Set to true to use multi-threaded evaluation.
		 * 
		 * @param parallel true to multi-thread
		 */
		public void setParallel(
				final boolean parallel
		) {
			mParallel = parallel;
		}

	}

	private abstract class JacobianEvaluator<J> //
			implements ICalculator<J> {

		public JacobianEvaluator() {
		}

		@Override
		public UncertainValues<J> evaluate(
				final ExplicitMeasurementModel<? extends J, ? extends J> func, //
				final List<? extends J> outputLabels //
		) {
			final Pair<RealVector, RealMatrix> eval = mEvaluation.get();
			final RealVector vals = eval.getFirst();
			final RealMatrix jac = eval.getSecond();
			final int extra = getInputDimension();
			assert jac.getColumnDimension() == extra;
			final int base = vals.getDimension();
			assert base == func.getOutputDimension();
			assert getOutputDimension() == extra + base;
			final RealVector val2 = new ArrayRealVector(getOutputDimension());
			final RealMatrix jac2 = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			final UncertainValues<? extends H> inputs = mInputs.get();
			assert inputs.getDimension() == extra;
			val2.setSubVector(0, vals);
			val2.setSubVector(vals.getDimension(), inputs.getValues());
			jac2.setSubMatrix(jac.getData(), 0, 0);
			// Tack an identity matrix on the end...
			for (int i = 0; i < extra; ++i)
				jac2.setEntry(i + base, i, 1.0);
			final RealMatrix cov2 = jac2.multiply(inputs.getCovariances().multiply(jac2.transpose()));
			return new UncertainValues<J>(//
					outputLabels, //
					val2, //
					cov2);
		}
	}

	/***
	 * Builds the outputs in the order of func.getOutputLabels() +
	 * func.getInputLabels(). Makes sure there are no duplicates.
	 * 
	 * @param emm The measurement model
	 * @return List&lt;H&gt;
	 */

	private final static <H> List<H> buildOutputs1(
			final ExplicitMeasurementModel<? extends H, ? extends H> emm
	) {
		final FastIndex<H> res = new FastIndex<>();
		res.addAll(emm.getOutputLabels());
		res.addAll(emm.getInputLabels());
		return res;
	}

	private final ExplicitMeasurementModel<? extends H, ? extends H> mFunction;

	private final SimplyLazy<UncertainValues<? extends H>> mInputs = new SimplyLazy<UncertainValues<? extends H>>() {

		@Override
		protected UncertainValues<? extends H> initialize() {
			final List<? extends H> inputLabels = mFunction.getInputLabels();
			final UncertainValuesBase<? extends H> reorder = mRawInputs.reorder(inputLabels);
			if (reorder instanceof UncertainValues<?>)
				return ((UncertainValues<? extends H>) reorder).copy();
			else
				return UncertainValues.asUncertainValues(reorder);
		}

	};

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
			return mCalculator.evaluate(mFunction, getOutputLabels());
		}
	};

	transient private SimplyLazy<RealVector> mOutputValues = new SimplyLazy<RealVector>() {
		@Override
		protected RealVector initialize() {
			RealVector val;
			if (mFullOutputs.initialized())
				return mFullOutputs.get().getValues();
			else {
				final RealVector inputs = mInputs.get().getValues();
				val = mFunction.compute(inputs);
				final RealVector val2 = new ArrayRealVector(getOutputDimension());
				val2.setSubVector(0, val);
				val2.setSubVector(val.getDimension(), inputs);
				return val2;
			}
		}
	};

	/**
	 * Constructs a {@link UncertainValuesCalculator} to propagate the measurement
	 * model over the specified set of input values to estimate the output
	 * uncertainties.
	 * 
	 * @param emm    The measurement model
	 * @param inputs The input values
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public UncertainValuesCalculator(
			final ExplicitMeasurementModel<? extends H, ? extends H> emm, //
			final UncertainValuesBase<H> inputs //
	) throws ArgumentException {
		super(buildOutputs1(emm));
		mFunction = emm;
		if (inputs != null) {
			mRawInputs = inputs;
			setInputs(inputs);
		}
		mCalculator = new Analytical();
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

	/**
	 * Returns the model as presented to the constructor.
	 * 
	 * @return The model
	 */
	public ExplicitMeasurementModel<? extends H, ? extends H> getFunction() {
		return mFunction;
	}

	/**
	 * The number of input variables.
	 * 
	 * @return int
	 */
	public int getInputDimension() {
		return mFunction.getInputDimension();
	}

	/**
	 * An ordered list of input variable labels.
	 * 
	 * @return List&lt;? extends H&gt;
	 */
	public List<? extends H> getInputLabels() {
		return mFunction.getInputLabels();
	}

	/**
	 * Returns the input {@link UncertainValues} in the correct order for the
	 * associated {@link ExplicitMeasurementModel}. If there are any Constraints
	 * imposed on this {@link UncertainValuesCalculator} then the returned values
	 * will reflect the constraints.
	 *
	 * @return {@link UncertainValues}
	 */
	public UncertainValues<? extends H> getInputs() {
		return mInputs.get();
	}

	/**
	 * Returns the input {@link RealVector} in the correct order for the associated
	 * {@link ExplicitMeasurementModel}.
	 *
	 *
	 * @return {@link RealVector}
	 */
	public RealVector getInputValues() {
		return mInputs.get().getValues();
	}

	/**
	 * Returns the Jacobian when the {@link ICalculator} associated with this {@link UncertainValuesCalculator}
	 * support calculating the Jacobian (the {@link MonteCarlo} doesn't).
	 * 
	 * @return Optional&lt;Jacobian&lt;H,H&gt;&gt;
	 */
	public Optional<Jacobian<H, H>> getJacobian() {
		final RealMatrix jac = mEvaluation.get().getSecond();
		if (jac != null)
			return Optional.of(new Jacobian<H, H>(getInputLabels(), mFunction.getOutputLabels(), jac));
		else
			return Optional.empty();
	}

	/**
	 * @return The number of output variables
	 */
	public int getOutputDimension() {
		return getDimension();
	}

	/**
	 * An ordered list of the labels associated with the output variables.
	 * 
	 * @return List&lt;H&gt;
	 */
	public List<H> getOutputLabels() {
		return getLabels();
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects. Equivalent to
	 * <code>getOutputValues(1.0e-6)</code>
	 *
	 * @return HashMap&lt;H, UncertainValue&gt;
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public HashMap<H, UncertainValueEx<H>> getOutputValues() //
			throws ArgumentException {
		return getOutputValues(1.0e-6);
	}

	/**
	 * Returns a Map of the labels associated with output values expressed as
	 * {@link UncertainValue} objects.
	 *
	 * @param tol Tolerance relative to value
	 * @return HashMap&lt;H, UncertainValueEx&lt;H&gt;&gt;
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public HashMap<H, UncertainValueEx<H>> getOutputValues(
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
	 * <p>
	 * Returns a map from output value to an UncertainValue with discretely broken
	 * out uncertainty components according to the labels map.
	 * </p>
	 *
	 * @param labels Group according to this mapping
	 * @param tol    Minimum size uncertainty to include
	 * @return HashMap&lt;H, UncertainValue&gt;
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public HashMap<H, UncertainValueEx<String>> getOutputValues(
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

	/**
	 * Specify an {@link ICalculator} and force a recalculation.
	 * 
	 * @param calc An {@link ICalculator}
	 */
	public void setCalculator(
			final ICalculator<H> calc
	) {
		synchronized (this) {
			mCalculator = calc;
			mEvaluation.reset();
			mFullOutputs.reset();
			mOutputValues.reset();
		}
	}

	/**
	 * The calculator used to propagate the uncertainties.
	 * 
	 * @return An instance of {@link ICalculator}
	 */
	public ICalculator<H> getCalculator() {
		return mCalculator;
	}

	/**
	 * Sets the input value and covariance matrix for the calculation. The input
	 * {@link UncertainValuesBase} is evaluated, reordered and any
	 * {@link Constraint}s are applied.
	 *
	 *
	 * @param inputs The input values as an {@link UncertainValuesBase}
	 * @throws ArgumentException When there is an inconsistency in the function
	 *                           arguments
	 */
	public void setInputs(
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

	@Override
	public String toString() {
		return mFunction.toString() + (mInputs != null ? "(" + mInputs.toString() + ")" : "(No inputs)");
	}
}