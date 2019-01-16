package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.AtomTypeTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;

public class AtomFractionToMassFraction extends LabeledMultivariateJacobianFunction {

	/**
	 * Nominally the tabulated atomic weights but can be customized for special
	 * situations.
	 */
	private final Map<Element, Double> mAtomicWeights;

	private static Map<Element, Double> buildAtomicWeights(final Set<Element> elms) {
		final Map<Element, Double> res = new HashMap<>();
		for (final Element elm : elms)
			res.put(elm, elm.getAtomicWeight());
		return res;
	}

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public AtomFractionToMassFraction(final Composition comp, final Map<Element, Double> atomicWeights) {
		super(comp.getLabels(), Composition.buildTags(comp, Representation.MassFraction));
		assert (comp.mRepresentation == Representation.AtomFraction)
				|| (comp.mRepresentation == Representation.Stoichiometry);
		mAtomicWeights = new HashMap<>(atomicWeights);
	}

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public AtomFractionToMassFraction(final Composition comp) {
		this(comp, buildAtomicWeights(comp.getElementSet()));
	}

	private double denom(final RealVector point) {
		double res = 0.0;
		for (final Object tag : getInputLabels())
			res += getValue(tag, point) * mAtomicWeights.get(((AtomTypeTag) tag).getElement());
		return res;
	}

	private AtomTypeTag find(MassFractionTag mft) {
		for (Object inp : getInputLabels())
			if (inp instanceof AtomTypeTag) {
				AtomTypeTag att = (AtomTypeTag) inp;
				if (att.getElement() == mft.getElement())
					return att;
			}
		return null;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getInputDimension(), getOutputDimension());
		final double den = denom(point);
		for (int i = 0; i < getOutputDimension(); ++i) {
			final MassFractionTag mft = (MassFractionTag) getOutputLabels().get(i);
			final double zm = mAtomicWeights.get(mft.getElement());
			for (int j = 0; j < getInputDimension(); ++j) {
				final AtomTypeTag aft = (AtomTypeTag) getInputLabels().get(j);
				final double za = mAtomicWeights.get(aft.getElement());
				double dCdA;
				if (mft.getElement() == aft.getElement()) {
					final double atomFrac = getValue(aft, point);
					final double mf = atomFrac * mAtomicWeights.get(aft.getElement()) / den;
					vals.setEntry(i, mf);
					dCdA = zm / den - atomFrac * Math.pow(zm / den, 2.0);
				} else {
					final double aa = getValue(find(mft), point);
					dCdA = -aa * zm * za / Math.pow(den, 2.0);
				}
				writeJacobian(i, aft, dCdA, jac);
			}
		}
		return Pair.create(vals, jac);
	}
}