/**
 *
 */
package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AnalyticalTotalTag;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.CompositionalStatisticTag;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MeanATag;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MeanZTag;

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
		extends ExplicitMeasurementModel<MaterialLabel, CompositionalStatisticTag> //
		implements ILabeledMultivariateFunction<MaterialLabel, CompositionalStatisticTag> {

	private static List<MaterialLabel> buildInputs(
			final Material mat
	) {
		final List<MaterialLabel> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildMassFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}

	private static List<CompositionalStatisticTag> buildOutputTags(
			final Material mat
	) {
		final List<CompositionalStatisticTag> res = new ArrayList<>();
		res.add(MaterialLabel.buildAnalyticalTotalTag(mat));
		res.add(MaterialLabel.buildMeanAtomicNumberTag(mat));
		res.add(MaterialLabel.buildMeanAtomicWeighTag(mat));
		return res;
	}

	private final Material mMaterial;

	/**
	 * @param inputLabels
	 * @param outputLabels
	 * @throws ArgumentException
	 */
	public CompositionStatistics(
			final Material mat
	) throws ArgumentException {
		super(buildInputs(mat), buildOutputTags(mat));
		mMaterial = mat;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction
	 * #optimized(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public RealVector optimized(
			final RealVector point
	) {
		final RealVector res = buildResult();
		double meanA = 0.0, total = 0.0, meanZ = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			final double mf = getArg(MaterialLabel.buildMassFractionTag(mMaterial, elm), point);
			final double aw = getArg(MaterialLabel.buildAtomicWeightTag(mMaterial, elm), point);
			total += mf;
			meanA += mf * aw;
			meanZ += mf * elm.getAtomicNumber();
		}
		setResult(MaterialLabel.buildAnalyticalTotalTag(mMaterial), res, total);
		setResult(MaterialLabel.buildMeanAtomicWeighTag(mMaterial), res, meanA);
		setResult(MaterialLabel.buildMeanAtomicNumberTag(mMaterial), res, meanZ);
		return res;
	}

	@Override
	public String toString() {
		return "Composition Statistics[" + mMaterial + "]";
	}

	/*
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		assert point.getDimension() == getInputDimension();
		double meanA = 0.0, total = 0.0, meanZ = 0.0;
		final RealVector res = buildResult();
		final RealMatrix jac = buildJacobian();
		final AnalyticalTotalTag totTag = MaterialLabel.buildAnalyticalTotalTag(mMaterial);
		final MeanATag maIdx = MaterialLabel.buildMeanAtomicWeighTag(mMaterial);
		final MeanZTag mzIdx = MaterialLabel.buildMeanAtomicNumberTag(mMaterial);
		for (final Element elm : mMaterial.getElementSet()) {
			final MassFraction mft = MaterialLabel.buildMassFractionTag(mMaterial, elm);
			final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
			final double mf = getArg(mft, point);
			final double aw = getArg(awt, point);
			total += mf;
			meanA += mf * aw;
			meanZ += mf * elm.getAtomicNumber();
			setJacobian(mft, totTag, jac, 1.0);
			setJacobian(awt, maIdx, jac, mf);
			setJacobian(mft, maIdx, jac, aw);
			setJacobian(mft, mzIdx, jac, elm.getAtomicNumber());
		}
		setResult(totTag, res, total);
		setResult(maIdx, res, meanA);
		setResult(mzIdx, res, meanZ);
		return Pair.create(res, jac);
	}

}
