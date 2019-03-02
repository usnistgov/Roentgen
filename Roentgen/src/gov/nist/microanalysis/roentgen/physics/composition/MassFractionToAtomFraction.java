package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;

public class MassFractionToAtomFraction //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	private static List<? extends Object> buildInputTags(final Material mat) {
		final List<Object> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildMassFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}

	/**
	 * Constructs a MassFractionToAtomFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public MassFractionToAtomFraction(//
			final Material mat, final Collection<Element> atomicWeightElms) {
		super(buildInputTags(mat), MaterialLabel.buildAtomFractionTags(mat));
	}

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public MassFractionToAtomFraction(final Material mat) {
		this(mat, Collections.emptySet());
	}

	private double getAtomicWeight(final MaterialLabel tag, final RealVector point) {
		return getValue(MaterialLabel.buildAtomicWeightTag(tag.getMaterial(), tag.getElement()), point);
	}

	private double denom(final RealVector point) {
		double res = 0.0;
		for (final Object tag : getInputLabels())
			if (tag instanceof MaterialLabel.MassFraction)
				res += getValue(tag, point) / getAtomicWeight((MaterialLabel) tag, point);
		return res;
	}

	private double delta(final MaterialLabel a1, final MaterialLabel a2) {
		return a1.getElement().equals(a2.getElement()) ? 1.0 : 0.0;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getOutputDimension(), getInputDimension());
		final double den = denom(point);
		for (int i1 = 0; i1 < getOutputDimension(); ++i1) {
			final MaterialLabel.AtomFraction aft1 = (MaterialLabel.AtomFraction) getOutputLabel(i1);
			final double w1 = getAtomicWeight(aft1, point);
			final double c1 = getValue(MaterialLabel.buildMassFractionTag(aft1.getMaterial(), aft1.getElement()),
					point);
			final double a1 = (c1 / w1) / den;
			vals.setEntry(i1, a1);
			for (int i2 = 0; i2 < getInputDimension(); ++i2) {
				final Object il = getInputLabel(i2);
				if (il instanceof MaterialLabel.MassFraction) {
					final MaterialLabel.MassFraction mft2 = (MaterialLabel.MassFraction) getInputLabel(i2);
					final double w2 = getAtomicWeight(mft2, point);
					jac.setEntry(i1, i2, (delta(aft1, mft2) - (c1 / (w2 * den))) / (w1 * den));
				} else {
					assert il instanceof MaterialLabel.AtomicWeight;
					final MaterialLabel.AtomicWeight awt2 = (MaterialLabel.AtomicWeight) getInputLabel(i2);
					final double w2 = point.getEntry(i2);
					final MaterialLabel.MassFraction mft2 = MaterialLabel.buildMassFractionTag(awt2.getMaterial(),
							awt2.getElement());
					final double c2 = getValue(mft2, point);
					assert (w1 - w2) * delta(aft1, mft2) == 0.0 : w1 + "!=" + w2;
					jac.setEntry(i1, i2, (c2 / (den * w2 * w2)) * (c1 / (w1 * den) - delta(aft1, mft2)));
				}
			}
		}
		return Pair.create(vals, jac);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final double den = denom(point);
		for (int i1 = 0; i1 < getOutputDimension(); ++i1) {
			final MaterialLabel.AtomFraction aft_o = (MaterialLabel.AtomFraction) getOutputLabel(i1);
			final double w1 = getAtomicWeight(aft_o, point);
			final double c1 = getValue(MaterialLabel.buildMassFractionTag(aft_o.getMaterial(), aft_o.getElement()),
					point);
			final double a1 = (c1 / w1) / den;
			vals.setEntry(i1, a1);
		}
		return vals;
	}
}