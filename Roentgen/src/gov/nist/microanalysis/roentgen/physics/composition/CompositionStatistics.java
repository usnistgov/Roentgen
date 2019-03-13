/**
 *
 */
package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * <p>
 * Compute compositional statistics
 * </p>
 * <ul>
 * <li>Mean Atomic Weight</li>
 * <li>Mean Atomic Number</li>
 * <li>Analytical Total</li>
 * </ul>
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class CompositionStatistics //
		extends LabeledMultivariateJacobianFunction //
		implements ILabeledMultivariateFunction {

	private static List<Object> buildOutputTags(final Material mat) {
		final List<Object> res = new ArrayList<>();
		res.add(MaterialLabel.buildAnalyticalTotalTag(mat));
		res.add(MaterialLabel.buildMeanAtomicNumberTag(mat));
		res.add(MaterialLabel.buildMeanAtomicWeighTag(mat));
		return res;
	}
	
	private static List<Object> buildInputs(final Material mat) {
		final List<Object> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildMassFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}


	private final Material mMaterial;

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public CompositionStatistics(final Material mat) {
		super(buildInputs(mat), buildOutputTags(mat));
		mMaterial = mat;
	}


	/*
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		assert point.getDimension() == getInputDimension();
		double meanA = 0.0, total = 0.0, meanZ = 0.0;
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), point.getDimension());
		final int totIdx = outputIndex(MaterialLabel.buildAnalyticalTotalTag(mMaterial));
		final int maIdx = outputIndex(MaterialLabel.buildMeanAtomicWeighTag(mMaterial));
		final int mzIdx = outputIndex(MaterialLabel.buildMeanAtomicNumberTag(mMaterial));
		for (final Element elm : mMaterial.getElementSet()) {
			final MassFraction mft = MaterialLabel.buildMassFractionTag(mMaterial, elm);
			final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
			final double mf = getValue(mft, point);
			final double aw = getValue(awt, point);
			total += mf;
			meanA += mf * aw;
			meanZ += mf * elm.getAtomicNumber();
			writeJacobian(totIdx, mft, 1.0, rm);
			writeJacobian(maIdx, awt, mf, rm);
			writeJacobian(maIdx, mft, aw, rm);
			writeJacobian(mzIdx, mft, elm.getAtomicNumber(), rm);
		}
		rv.setEntry(totIdx, total);
		rv.setEntry(maIdx, meanA);
		rv.setEntry(mzIdx, meanZ);
		return Pair.create(rv, rm);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction
	 * #optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(final RealVector point) {
		assert point.getDimension() == getInputDimension();
		double meanA = 0.0, total = 0.0, meanZ = 0.0;
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final int totIdx = outputIndex(MaterialLabel.buildAnalyticalTotalTag(mMaterial));
		final int maIdx = outputIndex(MaterialLabel.buildMeanAtomicWeighTag(mMaterial));
		final int mzIdx = outputIndex(MaterialLabel.buildMeanAtomicNumberTag(mMaterial));
		for (final Element elm : mMaterial.getElementSet()) {
			final MassFraction mft = MaterialLabel.buildMassFractionTag(mMaterial, elm);
			final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
			final double mf = getValue(mft, point);
			final double aw = getValue(awt, point);
			total += mf;
			meanA += mf * aw;
			meanZ += mf * elm.getAtomicNumber();
		}
		rv.setEntry(totIdx, total);
		rv.setEntry(maIdx, meanA);
		rv.setEntry(mzIdx, meanZ);
		return rv;
	}
	
	public String toString() {
		return "Composition Statistics";
	}


}
