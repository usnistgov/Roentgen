package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * Calculates the specified element using a set of valences to compute the mass
 * fraction of the element relative to the other elements in the material.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class ElementByStoichiometry //
		extends ExplicitMeasurementModel<MaterialLabel, MassFraction> //
		implements ILabeledMultivariateFunction<MaterialLabel, MassFraction> {

	static public ElementByStoichiometry buildDefaultOxygen(
			final Material mat
	) throws ArgumentException {
		return new ElementByStoichiometry(mat, Element.Oxygen, StoichiometeryRules.getOxygenDefaults());
	}

	private static List<MaterialLabel> buildInput(
			final Material mat, //
			final Element outputElm //
	) {
		final Set<Element> inputElms = new HashSet<>();
		inputElms.addAll(mat.getElementSet());
		inputElms.remove(outputElm);
		final List<MaterialLabel> res = new ArrayList<>();
		for (final Element elm : inputElms) {
			if (elm != outputElm) {
				res.add(MaterialLabel.buildMassFractionTag(mat, elm));
				res.add(MaterialLabel.buildAtomicWeightTag(mat, elm));
			}
		}
		res.add(MaterialLabel.buildAtomicWeightTag(mat, outputElm));
		return res;
	}

	private final Map<Element, Integer> mValences = new HashMap<>();

	private final Element mElement;

	private final Material mMaterial;

	/**
	 * @param inputLabels
	 * @param outputLabels
	 * @throws ArgumentException
	 */
	public ElementByStoichiometry(
			final Material mat, //
			final Element outputElm, //
			final Map<Element, Integer> valences //
	) throws ArgumentException {
		super(buildInput(mat, outputElm),
				Collections.singletonList(MaterialLabel.buildMassFractionTag(mat, outputElm)));
		mValences.putAll(valences);
		mElement = outputElm;
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
		final double ai = getArg(MaterialLabel.buildAtomicWeightTag(mMaterial, mElement), point);
		final double vi = mValences.get(mElement).doubleValue();
		double ci = 0.0;
		for (final Element elm : mMaterial.getElementSet()) {
			if (!elm.equals(mElement)) {
				final double cj = getArg(MaterialLabel.buildMassFractionTag(mMaterial, elm), point);
				assert cj >= 0.0 : elm.getAbbrev() + " = " + cj;
				final double aj = getArg(MaterialLabel.buildAtomicWeightTag(mMaterial, elm), point);
				final double vj = mValences.get(elm).doubleValue();
				ci -= (ai / aj) * (vj / vi) * cj;
			}
		}
		setResult(MaterialLabel.buildMassFractionTag(mMaterial, mElement), res, ci);
		return res;
	}

	@Override
	public String toString() {
		return mElement.getAbbrev() + "-by-stoichiometry";
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealVector res = buildResult();
		final RealMatrix jac = buildJacobian();
		final AtomicWeight awi = MaterialLabel.buildAtomicWeightTag(mMaterial, mElement);
		final double ai = getArg(awi, point);
		final double vi = mValences.get(mElement).doubleValue();
		double ci = 0.0;
		final MassFraction mfo = MaterialLabel.buildMassFractionTag(mMaterial, mElement);
		for (final Element elm : mMaterial.getElementSet()) {
			if (!elm.equals(mElement)) {
				final MassFraction mft = MaterialLabel.buildMassFractionTag(mMaterial, elm);
				final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mMaterial, elm);
				final double aj = getArg(awt, point), cj = getArg(mft, point);
				final double vj = mValences.get(elm).doubleValue();
				final double kkj = -(ai / aj) * (vj / vi);
				ci += kkj * cj;
				setJacobian(mft, mfo, jac, kkj);
				setJacobian(awt, mfo, jac, -(kkj * cj) / aj);
			}
		}
		setJacobian(awi, mfo, jac, ci / ai);
		setResult(mfo, res, ci);
		return Pair.create(res, jac);
	}
}
