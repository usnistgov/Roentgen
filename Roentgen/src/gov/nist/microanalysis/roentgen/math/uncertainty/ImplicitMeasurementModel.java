package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.utility.FastIndex;

/**
 * Solve C<sub>y</sub> U<sub>y</sub> C<sub>y</sub><sup>T</sup> = C<sub>x,z</sub>
 * U<sub>x,z</sub> C<sub>x,z</sub><sup>T</sup> for
 * C<sub>y</sub><sup>-1</sup>C<sub>x,z</sub> and <b>Y</b>.
 * <ul>
 * <li><b>h</b>(<b>X</b>,<b>Y</b>,<b>Z</b>) = <b>0</b>
 * <li>C<sub>y</sub> = &delta;h/&delta;C<sub>y</sub>
 * <li>C<sub>x,z</sub> = &delta;h/&delta;C<sub>x</sub> and
 * &delta;h/&delta;C<sub>z</sub>,
 * </ul>
 *
 * <p>
 * The values associated with the output labels are passed in as constant values
 * to ensure that the input and output labels are disjoint.
 * </p>
 *
 *
 * @author Nicholas
 *
 */
public class ImplicitMeasurementModel<G> //
		extends LabeledMultivariateJacobianFunction<G, G>//
		implements ILabeledMultivariateFunction<G, G> {

	public static abstract class HModel<G> //
			extends LabeledMultivariateJacobianFunction<G, G> {
		
		private final List<G> mHInputLabels;
		private final List<G> mHOutputLabels;
	
		private static <J> List<J> combine(final List<? extends J> inp, final List<? extends J> outp) {
			final List<J> res = new ArrayList<>();
			res.addAll(inp);
			res.addAll(outp);
			return res;
		}

		public HModel(final List<? extends G> inputLabels, final List<? extends G> outputLabels, final List<? extends G> hLabels) {
			super(combine(inputLabels, outputLabels), new ArrayList<G>(hLabels));
			assert hLabels.size() == outputLabels.size();
			mHInputLabels = Collections.unmodifiableList(new FastIndex<>(inputLabels));
			mHOutputLabels = Collections.unmodifiableList(new FastIndex<>(outputLabels));
		}

		/**
		 * @return the mHInputLabels
		 */
		public List<G> getHInputLabels() {
			return mHInputLabels;
		}

		/**
		 * @return the hOutputLabels
		 */
		public List<G> getHOutputLabels() {
			return mHOutputLabels;
		}
	}

	private final HModel<G> mHFunction;
	private RealVector mHValues;
	private final ILabeledMultivariateFunction<G, G> mAlternativeModel;

	/***
	 * A set of constant values for use evaluating the MulitvariateJacobianFunction.
	 */
	private final Map<G, Double> mConstants = new HashMap<>();

	private static <H> List<H> buildInputs(//
			final LabeledMultivariateJacobianFunction<? extends H, ? extends H> h, //
			final List<? extends H> outputLabels //
	) {
		final List<H> res = new ArrayList<>(h.getInputLabels());
		res.removeAll(outputLabels);
		return res;
	}

	public ImplicitMeasurementModel(//
			final HModel<G> h, //
			final List<? extends G> outputLabels, //
			final ILabeledMultivariateFunction<G, G> altModel //
	) {
		super(buildInputs(h, outputLabels), new ArrayList<>(outputLabels));
		mHFunction = h;
		final ArrayList<G> inputLabels = new ArrayList<>(h.getInputLabels());
		inputLabels.removeAll(outputLabels);
		assert altModel.getInputLabels().containsAll(inputLabels);
		assert altModel.getOutputLabels().containsAll(outputLabels);
		assert inputLabels.containsAll(altModel.getInputLabels());
		assert outputLabels.containsAll(altModel.getOutputLabels());
		mAlternativeModel = altModel;
	}

	public ImplicitMeasurementModel(//
			final HModel<G> h, //
			final List<G> outputLabels //
	) {
		this(h, outputLabels, null);
	}

	public boolean hasAlternativeModel() {
		return mAlternativeModel != null;
	}

	public ILabeledMultivariateFunction<G, G> getAlternativeModel() {
		return mAlternativeModel;
	}

	@Override
	public String toString() {
		return "Implicit[" + mHFunction + "]";
	}

	public HModel<G> getHModel() {
		return mHFunction;
	}

	/**
	 * Initializes the constant labels with the specified values. The label for
	 * value <code>vals.getEntry(i)</code> is <code>list.get(i)</code>. Checks first
	 * to see is the label is being used as an input variable and won't define it as
	 * a constant if it is an input item.
	 *
	 * @param list
	 * @param vals
	 */
	public void initializeConstants(final List<G> list, final RealVector vals) {
		for (int i = 0; i < list.size(); ++i)
			if (inputIndex(list.get(i)) == -1)
				mConstants.put(list.get(i), vals.getEntry(i));
	}

	/**
	 * Initializes the constant labels with the associated values.
	 *
	 * @param consts Map&lt;Object,Double&gt; where Object is a label
	 */
	public void initializeConstants(final Map<G, ? extends Number> consts) {
		for (final Entry<G, ? extends Number> me : consts.entrySet())
			if (inputIndex(me.getKey()) == -1)
				mConstants.put(me.getKey(), me.getValue().doubleValue());
	}

	/**
	 * Returns the constant value associated with <code>label</code>.
	 *
	 * @param label
	 * @return double
	 */
	public double getConstant(final G label) {
		assert inputIndex(label) == -1 : "Label " + label + " is a variable.";
		assert mConstants.containsKey(label) : "Label " + label + " is not a constant";
		return mConstants.get(label).doubleValue();
	}

	/**
	 * Check whether <code>label</code> is defined as a constant value.
	 *
	 * @param label A label
	 * @return true if <code>label</code> is initialized as a constant, false
	 *         otherwise.
	 */
	public boolean isConstant(final Object label) {
		return mConstants.containsKey(label);
	}

	public int getConstantDimension() {
		return mConstants.size();
	}

	public void assertHValues(final RealVector h, final double tol) {
		for (int i = 0; i < h.getDimension(); ++i)
			assert h.getEntry(i) < tol : //
			mHFunction.getInputLabel(i) + " = " + h.getEntry(i) + " (> " + tol;
	}

	public RealVector getHValues() {
		return mHValues;
	}

	/**
	 * Implements value(...) for an implicit measurement model by evaluating the
	 * h-function splitting it apart into <b>J<sub>y</sub></b> and
	 * <b>J<sub>x</sub></b> parts from which is computed <b>J</b> =
	 * <b>J<sub>y</sub></b><sup>-1</sup><b>J<sub>x</sub></b>
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 *      value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector hPoint = new ArrayRealVector(mHFunction.getInputDimension());
		for (int r = 0; r < hPoint.getDimension(); ++r) {
			final G hLbl = mHFunction.getInputLabel(r);
			final int ii = inputIndex(hLbl);
			if (ii != -1)
				hPoint.setEntry(r, point.getEntry(ii));
			else {
				assert isConstant(hLbl) : //
				hLbl + " is missing in " + this;
				hPoint.setEntry(r, getConstant(hLbl));
			}
		}
		final Pair<RealVector, RealMatrix> hpr = mHFunction.evaluate(hPoint);
		mHValues = hpr.getFirst();
		// hpr.getFirst() values should all be close to zero...
		final RealMatrix hrm = hpr.getSecond();
		// Partials wrt outputs
		final RealMatrix mcy = extractCy(hrm);
		// Partials wrt inputs
		final RealMatrix mcx = extractCx(hrm);
		final RealMatrix mcyInv = MatrixUtils.inverse(mcy);
		final RealMatrix jac = mcyInv.multiply(mcx);
		// Output values are input values
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (int r = 0; r < rv.getDimension(); ++r)
			rv.setEntry(r, getConstant(getOutputLabel(r)));
		return Pair.create(rv, jac);
	}

	private RealMatrix extractCx(final RealMatrix hrm) {
		final List<? extends G> yLabels = getOutputLabels();
		final List<Object> xLabels = new ArrayList<>();
		xLabels.addAll(getInputLabels());
		xLabels.removeAll(yLabels);
		final int n = xLabels.size(), m = yLabels.size();
		final List<? extends G> hIn = mHFunction.getInputLabels();
		// output x inputs
		final RealMatrix mc = MatrixUtils.createRealMatrix(m, n);
		for (int r = 0; r < m; ++r) {
			for (int c = 0; c < n; ++c) {
				final int cIdx = hIn.indexOf(xLabels.get(c));
				assert cIdx >= 0 : xLabels.get(c);
				mc.setEntry(r, c, hrm.getEntry(r, cIdx));
			}
		}
		return mc;
	}

	private RealMatrix extractCy(final RealMatrix hrm) {
		final List<? extends G> labels = getOutputLabels();
		final int m = labels.size();
		final List<? extends G> hIn = mHFunction.getInputLabels();
		// output x inputs
		final RealMatrix mc = MatrixUtils.createRealMatrix(m, m);
		for (int r = 0; r < m; ++r) {
			for (int c = 0; c < m; ++c) {
				final int cIdx = hIn.indexOf(labels.get(c));
				assert cIdx >= 0 : labels.get(c);
				mc.setEntry(r, c, hrm.getEntry(r, cIdx));
			}
		}
		return mc;
	}

	/**
	 * This is customed to handle the specific case in which a alternative model has
	 * been provided to perform the same calculation as the implicit model but in an
	 * explicit manner. This is useful for the Monte Carlo and Finite Difference
	 * calculations.
	 *
	 *
	 * @see gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction#optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		if (mAlternativeModel != null) {
			final List<? extends G> inLabels = mAlternativeModel.getInputLabels();
			final RealVector altPt = new ArrayRealVector(inLabels.size());
			for (int i = 0; i < inLabels.size(); ++i) {
				final int ii = inputIndex(inLabels.get(i));
				assert ii >= 0 : //
				inLabels.get(i) + " is missing in IMM.optimized(..)";
				altPt.setEntry(i, point.getEntry(ii));
			}
			final RealVector res = mAlternativeModel.optimized(altPt);
			final List<? extends G> outLabels = mAlternativeModel.getOutputLabels();
			for (int i = 0; i < outLabels.size(); ++i)
				rv.setEntry(outputIndex(outLabels.get(i)), res.getEntry(i));
		} else {
			// Trivial but not very useful implemention....
			for (int r = 0; r < rv.getDimension(); ++r)
				rv.setEntry(r, getConstant(getOutputLabel(r)));
		}
		return rv;
	}
}
