package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.models.Normalize;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MaterialMassFraction;

/**
 * Converts a mixture of Composition objects into the mass fraction of the
 * constituent elements.
 *
 * @author nicholas
 *
 */
final class MixtureToMassFractions //
		extends LabeledMultivariateJacobianFunction //
		implements ILabeledMultivariateFunction {

	private final List<Object> mInputs;
	private final Material mNewMaterial;

	static private List<? extends Object> buildOutputs(//
			final Material newMat, //
			final Set<Material> mats //
	) {
		final Set<Element> elms = new HashSet<>();
		final List<Object> res = new ArrayList<>();
		for (final Material mat : mats)
			elms.addAll(mat.getElementSet());
		res.addAll(MaterialLabel.buildMassFractionTags(newMat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(newMat));
		return res;
	}

	static private List<Object> buildInputs(//
			final Set<Material> mats, //
			final boolean normalize //
	) {
		final List<Object> res = new ArrayList<>();
		for (final Material mat : mats) {
			res.addAll(MaterialLabel.massFractionTags(mat));
			res.addAll(MaterialLabel.atomWeightTags(mat));
			if (normalize)
				res.add(Normalize.buildNormalized(MaterialLabel.buildMaterialFractionTag(mat)));
			else
				res.add(MaterialLabel.buildMaterialFractionTag(mat));
		}
		return res;
	}

	/**
	 * Converts a mixture of Composition objects into the mass fraction of the
	 * constituent elements.
	 *
	 * @param newMat    A Material to define
	 * @param mats      Set&lt;Material&gt;
	 * @param normalize
	 */
	public MixtureToMassFractions(final Material newMat, final Set<Material> mats, final boolean normalize) {
		super(buildInputs(mats, normalize), buildOutputs(newMat, mats));
		mInputs = new ArrayList<>();
		mNewMaterial = newMat;
		for (final Material mat : mats) {
			final MaterialMassFraction mmft = MaterialLabel.buildMaterialFractionTag(mat);
			mInputs.add(normalize ? Normalize.buildNormalized(mmft) : mmft);
		}
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MaterialLabel.MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final MaterialLabel.AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (int mmfti = 0; mmfti < mInputs.size(); ++mmfti) {
				final Object mmft = mInputs.get(mmfti);
				final Material mat = ((MaterialMassFraction) Normalize.unwrap(mmft)).getMaterial();
				final int mfti = inputIndex(MaterialLabel.buildMassFractionTag(mat, elm));
				// The material may or may not have the element...
				final double mmf = point.getEntry(mmfti);
				if (mfti != -1) {
					final int awti = inputIndex(MaterialLabel.buildAtomicWeightTag(mat, elm));
					assert awti != -1;
					final double mf = point.getEntry(mfti);
					final double aw = point.getEntry(awti);
					tmpCz += mf * mmf;
					tmpAz += mf * mmf / aw;
				}
			}
			final double cZ = tmpCz;
			final double aZ = tmpCz / tmpAz;
			for (int mmfti = 0; mmfti < mInputs.size(); ++mmfti) {
				final Object mmft = mInputs.get(mmfti);
				final Material mat = ((MaterialMassFraction) Normalize.unwrap(mmft)).getMaterial();
				final MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
				final int mfti = inputIndex(mft);
				// The material may or may not have the element...
				final double mI = point.getEntry(mmfti);
				if (mfti != -1) {
					rm.setEntry(iCz, mfti, mI);
					final int awti = inputIndex(MaterialLabel.buildAtomicWeightTag(mat, elm));
					assert awti != -1;
					final double cIZ = point.getEntry(mfti);
					rm.setEntry(iCz, mmfti, cIZ);
					rm.setEntry(iCz, mfti, mI);

					final double aIZ = point.getEntry(awti);
					final double kk = aZ * (1.0 - aZ / aIZ) / cZ;
					rm.setEntry(iAz, mfti, mI * kk);
					rm.setEntry(iAz, awti, Math.pow(aZ / aIZ, 2.0) * (cIZ / cZ) * mI);
					rm.setEntry(iAz, mmfti, cIZ * kk);
				}
			}
			rv.setEntry(iCz, cZ);
			rv.setEntry(iAz, aZ);

		}
		return Pair.create(rv, rm);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MaterialLabel.MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final MaterialLabel.AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (int mmfti = 0; mmfti < mInputs.size(); ++mmfti) {
				final Object mmft = mInputs.get(mmfti);
				final Material mat = ((MaterialMassFraction) Normalize.unwrap(mmft)).getMaterial();
				final int mfti = inputIndex(MaterialLabel.buildMassFractionTag(mat, elm));
				// The material may or may not have the element...
				final double mmf = point.getEntry(mmfti);
				if (mfti != -1) {
					final int awti = inputIndex(MaterialLabel.buildAtomicWeightTag(mat, elm));
					assert awti != -1;
					final double mf = point.getEntry(mfti);
					final double aw = point.getEntry(awti);
					tmpCz += mf * mmf;
					tmpAz += mf * mmf / aw;
				}
			}
			final double cZ = tmpCz;
			final double aZ = tmpCz / tmpAz;
			rv.setEntry(iCz, cZ);
			rv.setEntry(iAz, aZ);

		}
		return rv;
	}

	@Override
	public String toString() {
		return "Mixture-to-" + mNewMaterial.toString();
	}
}