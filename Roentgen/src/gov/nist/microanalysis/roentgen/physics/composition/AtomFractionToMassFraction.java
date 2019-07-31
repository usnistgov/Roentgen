package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections4.map.HashedMap;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * Converts compositional data in atom fraction into mass fraction. Requires
 * also the atomic weights for each element in the material.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class AtomFractionToMassFraction //
		extends ExplicitMeasurementModel<MaterialLabel, MassFraction> {

	private static List<MaterialLabel> buildInputTags(final Material mat) {
		final List<MaterialLabel> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildAtomFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}

	public final Material mMaterial;

	/**
	 * Constructs a AtomicFractionToMassFraction instance.
	 *
	 * @param String                    html
	 * @param Collection&lt;Element&gt; The elements present in the material.
	 * @throws ArgumentException
	 */
	public AtomFractionToMassFraction(final Material mat //
	) throws ArgumentException {
		super(buildInputTags(mat), MaterialLabel.buildMassFractionTags(mat));
		mMaterial = mat;
	}

	@Override
	public RealVector computeValue(final double[] point) {
		final RealVector res = buildResult();
		Map<Element, Double> awafs = new HashedMap<>();
		double den = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
			final AtomFraction aft = MaterialLabel.buildAtomFractionTag(mMaterial, elm);
			final double awaf = getArg(awt, point) * getArg(aft, point);
			den += awaf;
			awafs.put(elm, Double.valueOf(awaf));
		}
		for (final MassFraction mft : getOutputLabels()) {
			assert mft.getMaterial() == mMaterial;
			setResult(mft, res, awafs.get(mft.getElement()).doubleValue() / den);
		}
		return res;
	}

	@Override
	public String toString() {
		return "Atom Fraction-to-Mass Fraction[" + mMaterial + "]";
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector res = buildResult();
		final RealMatrix jac = buildJacobian();
		Map<Element, Pair<AtomicWeight, AtomFraction>> tags = new HashedMap<>();
		double den = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
			final AtomFraction aft = MaterialLabel.buildAtomFractionTag(mMaterial, elm);
			tags.put(elm, Pair.create(awt, aft));
			final double aw = getArg(awt, point);
			final double af = getArg(aft, point);
			den += aw * af;
		}
		for (final MassFraction mft1 : getOutputLabels()) {
			assert mft1.getMaterial() == mMaterial;
			Pair<AtomicWeight, AtomFraction> tag1 = tags.get(mft1.getElement());
			final double w1 = getArg(tag1.getFirst(), point);
			final double a1 = getArg(tag1.getSecond(), point);
			final double c1 = a1 * w1 / den;
			setResult(mft1, res, c1);
			for (final MassFraction mft2 : getOutputLabels()) {
				assert mft2.getMaterial() == mMaterial;
				Pair<AtomicWeight, AtomFraction> tag2 = tags.get(mft2.getElement());
				final AtomicWeight awt2 = tag2.getFirst();
				final AtomFraction aft2 = tag2.getSecond();
				final double w2 = getArg(awt2, point);
				final double a2 = getArg(aft2, point);
				if (mft1.getElement().equals(mft2.getElement())) {
					double ff = (1.0 / den) * (1.0 - a1 * w1 / den);
					setJacobian(aft2, mft1, jac, w1 * ff); // ok
					setJacobian(awt2, mft1, jac, a1 * ff); // ok
				} else {
					double ff = -a1 * w1 / (den * den);
					setJacobian(aft2, mft1, jac, w2 * ff); // ok
					setJacobian(awt2, mft1, jac, a2 * ff); // ok
				}
			}
		}
		return Pair.create(res, jac);
	}

}