package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class ElementByStoichiometry extends LabeledMultivariateJacobianFunction
		implements ILabeledMultivariateFunction {

	private final Map<Element, Integer> mValences = new HashMap<>();
	private final Element mElement;

	private static List<? extends Object> buildInput(//
			final Material mat, //
			final Element outputElm //
	) {
		Set<Element> inputElms = new HashSet<>();
		inputElms.addAll(mat.getElementSet());
		inputElms.remove(outputElm);
		final List<Object> res = new ArrayList<>();
		for (final Element elm : inputElms) {
			if (elm != outputElm) {
				res.add(MaterialLabel.buildMassFractionTag(mat, elm));
				res.add(MaterialLabel.buildAtomicWeightTag(mat, elm));
			}
		}
		res.add(MaterialLabel.buildAtomicWeightTag(mat, outputElm));
		return res;
	}

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public ElementByStoichiometry(final Material mat, final Element outputElm, final Map<Element, Integer> valences) {
		super(buildInput(mat, outputElm),
				Collections.singletonList(MaterialLabel.buildMassFractionTag(mat, outputElm)));
		mValences.putAll(valences);
		mElement = outputElm;
	}
	
	static public ElementByStoichiometry buildDefaultOxygen(Material mat) {
		return new ElementByStoichiometry(mat, Element.Oxygen, StoichiometeryRules.getOxygenDefaults());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see
	 * org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#
	 * value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(1);
		final RealMatrix rm = MatrixUtils.createRealMatrix(1, getInputDimension());
		double awi = Double.NaN;
		int awii = -1;
		double ci = 0.0, dcidawi = 0.0;
		for (int j = 0; j < getInputDimension(); ++j) {
			final Object lbl = getInputLabel(j);
			if (lbl instanceof MassFraction) {
				final MassFraction cjt = (MassFraction) lbl;
				final int awji = inputIndex(MaterialLabel.buildAtomicWeightTag(cjt.getMaterial(), cjt.getElement()));
				if (awii==-1) {
					awii = inputIndex(MaterialLabel.buildAtomicWeightTag(cjt.getMaterial(), mElement));
					awi = point.getEntry(awii);
				}
				final double cj = point.getEntry(j);
				final double awj = point.getEntry(awji);
				final double kkj = -(awi / awj)
						* (mValences.get(cjt.getElement()).doubleValue() / mValences.get(mElement).doubleValue());
				ci += kkj * cj;
				dcidawi += (kkj / awi) * cj;
				rm.setEntry(0, j, kkj);
				rm.setEntry(0, awji, -(kkj / awj) * cj);
			}
			rm.setEntry(0, awii, dcidawi);
		}
		rv.setEntry(0, ci);
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
		final RealVector rv = new ArrayRealVector(1);
		double ai = Double.NaN;
		double ci = 0.0;
		for (int j = 0; j < getInputDimension(); ++j) {
			final Object lbl = getInputLabel(j);
			if (lbl instanceof MassFraction) {
				final MassFraction mft = (MassFraction) lbl;
				final int awt = inputIndex(MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mft.getElement()));
				if (Double.isNaN(ai)) {
					final int awe = inputIndex(MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mElement));
					ai = point.getEntry(awe);
				}
				final double cj = point.getEntry(j);
				final double aj = point.getEntry(awt);
				ci -= (ai / aj)
						* (mValences.get(mft.getElement()).doubleValue() / mValences.get(mElement).doubleValue()) * cj;
			}
		}
		rv.setEntry(0, ci);
		return rv;
	}
	
	public String toString() {
		return mElement.getAbbrev()+"-by-stoichiometry";
	}

}
