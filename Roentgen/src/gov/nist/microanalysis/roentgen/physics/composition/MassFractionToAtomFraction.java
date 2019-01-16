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

public class MassFractionToAtomFraction extends LabeledMultivariateJacobianFunction {

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
	public MassFractionToAtomFraction(final Composition comp, final Map<Element, Double> atomicWeights) {
		super(comp.getLabels(), Composition.buildTags(comp, Representation.AtomFraction));
		assert (comp.mRepresentation == Representation.MassFraction);
		mAtomicWeights = new HashMap<>(atomicWeights);
	}

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public MassFractionToAtomFraction(final Composition comp) {
		this(comp, buildAtomicWeights(comp.getElementSet()));
	}

	private double denom(final RealVector point) {
		double res = 0.0;
		for (final Object tag : getInputLabels())
			res += getValue(tag, point) / mAtomicWeights.get(((MassFractionTag) tag).getElement());
		return res;
	}

	private MassFractionTag find(AtomTypeTag att) {
		for (Object inp : getInputLabels())
			if (inp instanceof MassFractionTag) {
				MassFractionTag mft = (MassFractionTag) inp;
				if (att.getElement() == mft.getElement())
					return mft;
			}
		return null;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getInputDimension(), getOutputDimension());
		final double den = denom(point);
		for (int i = 0; i < getOutputDimension(); ++i) {
			final AtomTypeTag att = (AtomTypeTag) getOutputLabels().get(i);
			final double za = mAtomicWeights.get(att.getElement());
			for (int j = 0; j < getInputDimension(); ++j) {
				final MassFractionTag mft = (MassFractionTag) getInputLabels().get(j);
				final double zm = mAtomicWeights.get(mft.getElement());
				double dCdA;
				if (mft.getElement() == att.getElement()) {
					final double massFrac = getValue(mft, point);
					final double atomFrac = massFrac / (zm * den);
					vals.setEntry(i, atomFrac);
					dCdA = 1.0 / (za * den) - massFrac / Math.pow(za * den, 2.0);
				} else {
					final double massFrac = getValue(find(att), point);
					dCdA = -massFrac / (zm * za * Math.pow(den, 2.0));
				}
				writeJacobian(i, mft, dCdA, jac);
			}
		}
		return Pair.create(vals, jac);
	}
}