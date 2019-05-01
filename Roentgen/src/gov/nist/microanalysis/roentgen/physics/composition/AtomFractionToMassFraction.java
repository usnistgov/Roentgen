package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

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

	private static List<MaterialLabel> buildInputTags(
			final Material mat
	) {
		final List<MaterialLabel> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildAtomFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}

	public final Material mMaterial;

	/**
	 * Constructs a AtomicFractionToMassFraction instance.
	 *
	 * @param String html
	 * @param        Collection&lt;Element&gt; The elements present in the material.
	 * @throws ArgumentException
	 */
	public AtomFractionToMassFraction(
			//
			final Material mat //
	) throws ArgumentException {
		super(buildInputTags(mat), MaterialLabel.buildMassFractionTags(mat));
		mMaterial = mat;
	}

	@Override
	public RealVector computeValue(
			final double[] point
	) {
		final RealVector res = buildResult();
		double den = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			final double aw = getArg(MaterialLabel.buildAtomicWeightTag(mMaterial, elm), point);
			final double af = getArg(MaterialLabel.buildAtomFractionTag(mMaterial, elm), point);
			den += aw * af;
		}
		for (final MassFraction mft : getOutputLabels()) {
			final Material mat = mft.getMaterial();
			final Element elm = mft.getElement();
			final double w1 = getArg(MaterialLabel.buildAtomicWeightTag(mat, elm), point);
			final double a1 = getArg(MaterialLabel.buildAtomFractionTag(mat, elm), point);
			setResult(mft, res, a1 * w1 / den);
		}
		return res;
	}

	@Override
	public String toString() {
		return "Atom Fraction-to-Mass Fraction[" + mMaterial + "]";
	}

	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealVector res = buildResult();
		final RealMatrix jac = buildJacobian();
		double den = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			final double aw = getArg(MaterialLabel.buildAtomicWeightTag(mMaterial, elm), point);
			final double af = getArg(MaterialLabel.buildAtomFractionTag(mMaterial, elm), point);
			den += aw * af;
		}
		for (final MassFraction mft1 : getOutputLabels()) {
			final Material mat1 = mft1.getMaterial();
			final Element elm1 = mft1.getElement();
			final AtomicWeight awt1 = MaterialLabel.buildAtomicWeightTag(mat1, elm1);
			final AtomFraction aft1 = MaterialLabel.buildAtomFractionTag(mat1, elm1);
			final double w1 = getArg(awt1, point);
			final double a1 = getArg(aft1, point);
			final double c1 = a1 * w1 / den;
			setResult(mft1, res, c1);
			for (final MassFraction mft2 : getOutputLabels()) {
				final Material mat2 = mft2.getMaterial();
				final Element elm2 = mft2.getElement();
				final AtomicWeight awt2 = MaterialLabel.buildAtomicWeightTag(mat2, elm2);
				final AtomFraction aft2 = MaterialLabel.buildAtomFractionTag(mat2, elm2);
				final double w2 = getArg(awt2, point);
				final double a2 = getArg(aft2, point);
				if (elm1.equals(elm2)) {
					setJacobian(aft2, mft1, jac, (w1 / den) * (1.0 - a1 * w2 / den));
					setJacobian(awt2, mft1, jac, (a1 / den) * (1.0 - w1 * a2 / den));
				} else {
					setJacobian(aft2, mft1, jac, (w1 / den) * (-a1 * w2 / den));
					setJacobian(awt2, mft1, jac, (a1 / den) * (-w1 * a2 / den));
				}
			}
		}
		return Pair.create(res, jac);
	}

}