package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Collection;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.AtomFractionTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.AtomTypeTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.StoichiometryTag;

public class StoichiometryToAtomFraction //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public StoichiometryToAtomFraction(final String html, final Collection<Element> elms) {
		super(Composition.buildStoichiometryTags(html, elms), Composition.buildAtomFractionTags(html, elms));
	}

	private double denom(final RealVector point) {
		double res = 0.0;
		for (final Object tag : getInputLabels())
			res += getValue(tag, point);
		return res;
	}

	private double delta(final AtomTypeTag a1, final AtomTypeTag a2) {
		return a1.getElement().equals(a2.getElement()) ? 1.0 : 0.0;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getInputDimension(), getOutputDimension());
		final double den = denom(point);
		for (int outIdx = 0; outIdx < getOutputDimension(); ++outIdx) {
			final AtomFractionTag aft1 = (AtomFractionTag) getOutputLabel(outIdx);
			final StoichiometryTag st1 = Composition.buildStoichiometryTag(aft1.getHTML(), aft1.getElement());
			final double s1 = getValue(st1, point);
			vals.setEntry(outIdx, s1 / den);
			for (int inIdx = 0; inIdx < getInputDimension(); ++inIdx) {
				final StoichiometryTag st2 = (StoichiometryTag) getInputLabel(inIdx);
				jac.setEntry(outIdx, inIdx, (delta(aft1, st2) - s1 / den) / den);
			}
		}
		return Pair.create(vals, jac);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final double den = denom(point);
		for (int outIdx = 0; outIdx < getOutputDimension(); ++outIdx) {
			final AtomFractionTag aft1 = (AtomFractionTag) getOutputLabel(outIdx);
			final StoichiometryTag st1 = Composition.buildStoichiometryTag(aft1.getHTML(), aft1.getElement());
			final double s1 = getValue(st1, point);
			vals.setEntry(outIdx, s1 / den);
		}
		return vals;
	}

}