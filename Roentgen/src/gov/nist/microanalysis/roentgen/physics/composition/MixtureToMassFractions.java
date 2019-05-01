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

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
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
		extends ExplicitMeasurementModel<MaterialLabel, MaterialLabel> {

	static private List<MaterialLabel> buildInputs(
			//
			final Set<Material> mats, //
			final boolean normalize //
	) {
		final List<MaterialLabel> res = new ArrayList<>();
		for (final Material mat : mats) {
			res.addAll(MaterialLabel.massFractionTags(mat));
			res.addAll(MaterialLabel.atomWeightTags(mat));
			if (normalize)
				res.add(MaterialLabel.buildNormalizedMaterialFractionTag(mat));
			else
				res.add(MaterialLabel.buildMaterialFractionTag(mat));
		}
		return res;
	}

	static private List<MaterialLabel> buildOutputs(
			//
			final Material newMat, //
			final Set<Material> mats //
	) {
		final Set<Element> elms = new HashSet<>();
		final List<MaterialLabel> res = new ArrayList<>();
		for (final Material mat : mats)
			elms.addAll(mat.getElementSet());
		res.addAll(MaterialLabel.buildMassFractionTags(newMat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(newMat));
		return res;
	}

	private final List<MaterialMassFraction> mInputs;

	private final Material mNewMaterial;

	/**
	 * Converts a mixture of Composition objects into the mass fraction of the
	 * constituent elements.
	 *
	 * @param newMat    A Material to define
	 * @param mats      Set&lt;Material&gt;
	 * @param normalize
	 * @throws ArgumentException
	 */
	public MixtureToMassFractions(
			final Material newMat, //
			final Set<Material> mats, //
			final boolean normalize
	) throws ArgumentException {
		super(buildInputs(mats, normalize), buildOutputs(newMat, mats));
		mInputs = new ArrayList<>();
		mNewMaterial = newMat;
		for (final Material mat : mats) {
			if (normalize)
				mInputs.add(MaterialLabel.buildNormalizedMaterialFractionTag(mat));
			else
				mInputs.add(MaterialLabel.buildMaterialFractionTag(mat));
		}
	}

	@Override
	public RealVector computeValue(
			final double[] point
	) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (int i = 0; i < mInputs.size(); ++i) {
				final MaterialMassFraction mmft = mInputs.get(i);
				final Material mat = mmft.getMaterial();
				final MassFraction mfti = MaterialLabel.buildMassFractionTag(mat, elm);
				// The material may or may not have the element...
				final double mmf = getArg(mmft, point);
				if (inputIndex(mfti) != -1) {
					final AtomicWeight awti = MaterialLabel.buildAtomicWeightTag(mat, elm);
					final double mf = getArg(mfti, point);
					final double aw = getArg(awti, point);
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

	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (int i = 0; i < mInputs.size(); ++i) {
				final MaterialMassFraction mmft = mInputs.get(i);
				final Material mat = mmft.getMaterial();
				final int mfi = inputIndex(MaterialLabel.buildMassFractionTag(mat, elm));
				// The material may or may not have the element...
				final int mmfti = inputIndex(mmft);
				final double mmf = point.getEntry(mmfti);
				if (mfi != -1) {
					final int awi = inputIndex(MaterialLabel.buildAtomicWeightTag(mat, elm));
					assert awi != -1;
					final double mf = point.getEntry(mfi);
					final double aw = point.getEntry(awi);
					tmpCz += mf * mmf;
					tmpAz += mf * mmf / aw;
				}
			}
			final double cZ = tmpCz;
			final double aZ = tmpCz / tmpAz;
			for (int i = 0; i < mInputs.size(); ++i) {
				final MaterialLabel mmft = mInputs.get(i);
				final Material mat = mmft.getMaterial();
				final MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
				final int mfi = inputIndex(mft);
				// The material may or may not have the element...
				final int mmfi = inputIndex(mmft);
				final double mmf = point.getEntry(mmfi);
				if (mfi != -1) {
					rm.setEntry(iCz, mfi, mmf);
					final int awi = inputIndex(MaterialLabel.buildAtomicWeightTag(mat, elm));
					assert awi != -1;
					final double mf = point.getEntry(mfi);
					rm.setEntry(iCz, mmfi, mf);
					rm.setEntry(iCz, mfi, mmf);

					final double aw = point.getEntry(awi);
					final double kk = aZ * (1.0 - aZ / aw) / cZ;
					rm.setEntry(iAz, mfi, mmf * kk);
					rm.setEntry(iAz, awi, Math.pow(aZ / aw, 2.0) * (mf / cZ) * mmf);
					rm.setEntry(iAz, mmfi, mf * kk);
				}
			}
			rv.setEntry(iCz, cZ);
			rv.setEntry(iAz, aZ);

		}
		return Pair.create(rv, rm);
	}
}