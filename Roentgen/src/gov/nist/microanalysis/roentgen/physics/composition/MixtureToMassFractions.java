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
import gov.nist.microanalysis.roentgen.physics.Element;

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
			final Set<Material> mats //
	) {
		final List<Object> res = new ArrayList<>();
		for (final Material mat : mats) {
			res.addAll(MaterialLabel.massFractionTags(mat));
			res.addAll(MaterialLabel.atomWeightTags(mat));
			res.add(MaterialLabel.buildMaterialFractionTag(mat));
		}
		return res;
	}

	static private Set<Element> buildElements(final Set<Material> mats) {
		final Set<Element> elms = new HashSet<>();
		for (final Material mat : mats)
			elms.addAll(mat.getElementSet());
		return elms;
	}
	
	/**
	 * Converts a mixture of Composition objects into the mass fraction of the
	 * constituent elements.
	 * 
	 * @param htmlName String
	 * @param mats     Set&lt;Material&gt;
	 */
	public MixtureToMassFractions(final String htmlName, final Set<Material> mats) {
		super(buildInputs(mats), buildOutputs(new Material(htmlName, buildElements(mats)), mats));
		mNewMaterial = new Material(htmlName, buildElements(mats));
	}

	public Material getNewMaterial() {
		return mNewMaterial;
	}
	
	
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final List<MaterialMassFraction> mmfts = new ArrayList<>();
		for (final Object inTag : getInputLabels())
			if (inTag instanceof MaterialMassFraction)
				mmfts.add((MaterialMassFraction) inTag);
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MaterialLabel.MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final MaterialLabel.AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (final MaterialMassFraction mmft : mmfts) {
				final Material mat = mmft.getMaterial();
				final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
				// The material may or may not have the element...
				final double mmf = getValue(mmft, point);
				if (hasValue(mft)) {
					final MaterialLabel.AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mat, elm);
					assert hasValue(awt);
					final double mf = getValue(mft, point);
					final double aw = getValue(awt, point);
					tmpCz += mf * mmf;
					tmpAz += mf * mmf / aw;
				}
			}
			final double cZ = tmpCz;
			final double aZ = tmpCz / tmpAz;
			for (final MaterialMassFraction mmft : mmfts) {
				final Material mat = mmft.getMaterial();
				final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
				// The material may or may not have the element...
				final double mI = getValue(mmft, point);
				writeJacobian(iCz, mft, mI, rm);
				if (hasValue(mft)) {
					final MaterialLabel.AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mat, elm);
					assert hasValue(awt);
					final double cIZ = getValue(mft, point);
					writeJacobian(iCz, mmft, cIZ, rm);
					writeJacobian(iCz, mft, mI, rm);

					final double aIZ = getValue(awt, point);
					final double kk = aZ * (1.0 - aZ / aIZ) / cZ;
					writeJacobian(iAz, mft, mI * kk, rm);
					writeJacobian(iAz, awt, Math.pow(aZ / aIZ, 2.0) * (cIZ / cZ) * mI, rm);
					writeJacobian(iAz, mmft, cIZ * kk, rm);
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
		final List<MaterialMassFraction> mmfts = new ArrayList<>();
		for (final Object inTag : getInputLabels())
			if (inTag instanceof MaterialMassFraction)
				mmfts.add((MaterialMassFraction) inTag);
		for (final Element elm : mNewMaterial.getElementSet()) {
			double tmpCz = 0.0, tmpAz = 0.0;
			final MaterialLabel.MassFraction resCz = MaterialLabel.buildMassFractionTag(mNewMaterial, elm);
			final MaterialLabel.AtomicWeight resAz = MaterialLabel.buildAtomicWeightTag(mNewMaterial, elm);
			final int iCz = outputIndex(resCz);
			final int iAz = outputIndex(resAz);
			for (final MaterialMassFraction mmft : mmfts) {
				final Material mat = mmft.getMaterial();
				final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(mat, elm);
				// The material may or may not have the element...
				final double mmf = getValue(mmft, point);
				if (hasValue(mft)) {
					final MaterialLabel.AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mat, elm);
					assert hasValue(awt);
					final double mf = getValue(mft, point);
					final double aw = getValue(awt, point);
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
}