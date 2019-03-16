package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

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
public class ImplicitMeasurementModel //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	public static class HLabel extends BaseLabel<Object, Object, Object> {

		public HLabel(final Object label) {
			super("h", label);
		}
	}

	public static abstract class HModel //
			extends LabeledMultivariateJacobianFunction {

		private static List<? extends Object> combine(final List<? extends Object> inp,
				final List<? extends Object> outp) {
			final List<Object> res = new ArrayList<>();
			res.addAll(inp);
			res.addAll(outp);
			return res;
		}

		/**
		 * Builds a special set of labels to represent the function
		 * <b>h</b>(<b>X</b>,<b>Y</b>,<b>Z</b>) = <b>0</b>.
		 *
		 * @param outlabels
		 * @return List&lt;? extends Object&gt;
		 */
		private static List<? extends Object> buildHLabels(//
				final List<? extends Object> outlabels //
		) {
			final List<Object> res = new ArrayList<>();
			for (final Object label : outlabels)
				res.add(new HLabel(label));
			return res;
		}

		public HModel(final List<? extends Object> inputLabels, final List<? extends Object> outputLabels) {
			super(combine(inputLabels, outputLabels), buildHLabels(outputLabels));
		}
	}

	private final HModel mHFunction;
	private RealVector mHValues;
	private final ILabeledMultivariateFunction mAlternativeModel;

	/***
	 * A set of constant values for use evaluating the MulitvariateJacobianFunction.
	 */
	private final Map<Object, Double> mConstants = new HashMap<>();

	private static List<? extends Object> buildInputs(//
			final LabeledMultivariateJacobianFunction h, //
			final List<? extends Object> outputLabels //
	) {
		final List<? extends Object> res = new ArrayList<>(h.getInputLabels());
		res.removeAll(outputLabels);
		return res;
	}

	public ImplicitMeasurementModel(//
			final HModel h, //
			final List<? extends Object> outputLabels, //
			final ILabeledMultivariateFunction altModel //
	) {
		super(buildInputs(h, outputLabels), outputLabels);
		mHFunction = h;
		ArrayList<? extends Object> inputLabels = new ArrayList<>(h.getInputLabels());
		inputLabels.removeAll(outputLabels);
		assert altModel.getInputLabels().containsAll(inputLabels);
		assert altModel.getOutputLabels().containsAll(outputLabels);
		assert inputLabels.containsAll(altModel.getInputLabels());
		assert outputLabels.containsAll(altModel.getOutputLabels());
		mAlternativeModel = altModel;
	}

	public ImplicitMeasurementModel(//
			final HModel h, //
			final List<? extends Object> outputLabels //
	) {
		this(h, outputLabels, null);
	}

	public boolean hasAlternativeModel() {
		return mAlternativeModel != null;
	}

	public ILabeledMultivariateFunction getAlternativeModel() {
		return mAlternativeModel;
	}

	@Override
	public String toString() {
		return "Implicit[" + mHFunction + "]";
	}

	public HModel getHModel() {
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
	public void initializeConstants(final List<? extends Object> list, final RealVector vals) {
		for (int i = 0; i < list.size(); ++i)
			if (inputIndex(list.get(i)) == -1)
				mConstants.put(list.get(i), vals.getEntry(i));
	}

	/**
	 * Initializes the constant labels with the associated values.
	 *
	 * @param consts Map&lt;Object,Double&gt; where Object is a label
	 */
	public void initializeConstants(final Map<Object, ? extends Number> consts) {
		for (final Entry<Object, ? extends Number> me : consts.entrySet())
			if (inputIndex(me.getKey()) == -1)
				mConstants.put(me.getKey(), me.getValue().doubleValue());
	}

	/**
	 * Returns the constant value associated with <code>label</code>.
	 *
	 * @param label
	 * @return double
	 */
	public double getConstant(final Object label) {
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
			final Object hLbl = mHFunction.getInputLabel(r);
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
		final List<? extends Object> yLabels = getOutputLabels();
		final List<Object> xLabels = new ArrayList<>();
		xLabels.addAll(getInputLabels());
		xLabels.removeAll(yLabels);
		final int n = xLabels.size(), m = yLabels.size();
		final List<? extends Object> hIn = mHFunction.getInputLabels();
		final List<? extends Object> hOut = mHFunction.getOutputLabels();
		// output x inputs
		final RealMatrix mc = MatrixUtils.createRealMatrix(m, n);
		for (int r = 0; r < m; ++r) {
			final int rIdx = hOut.indexOf(new HLabel(yLabels.get(r)));
			assert rIdx >= 0 : yLabels.get(r);
			for (int c = 0; c < n; ++c) {
				final int cIdx = hIn.indexOf(xLabels.get(c));
				assert cIdx >= 0 : xLabels.get(c);
				mc.setEntry(r, c, hrm.getEntry(rIdx, cIdx));
			}
		}
		return mc;
	}

	private RealMatrix extractCy(final RealMatrix hrm) {
		final List<? extends Object> labels = getOutputLabels();
		final int m = labels.size();
		final List<? extends Object> hIn = mHFunction.getInputLabels();
		final List<? extends Object> hOut = mHFunction.getOutputLabels();
		// output x inputs
		final RealMatrix mc = MatrixUtils.createRealMatrix(m, m);
		for (int r = 0; r < m; ++r) {
			final int rIdx = hOut.indexOf(new HLabel(labels.get(r)));
			assert rIdx >= 0 : labels.get(r);
			for (int c = 0; c < m; ++c) {
				final int cIdx = hIn.indexOf(labels.get(c));
				assert cIdx >= 0 : labels.get(c);
				mc.setEntry(r, c, hrm.getEntry(rIdx, cIdx));
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
			final List<? extends Object> inLabels = mAlternativeModel.getInputLabels();
			final RealVector altPt = new ArrayRealVector(inLabels.size());
			for (int i = 0; i < inLabels.size(); ++i) {
				final int ii = inputIndex(inLabels.get(i));
				assert ii >= 0 : //
				inLabels.get(i) + " is missing in IMM.optimized(..)";
				altPt.setEntry(i, point.getEntry(ii));
			}
			final RealVector res = mAlternativeModel.optimized(altPt);
			final List<? extends Object> outLabels = mAlternativeModel.getOutputLabels();
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
