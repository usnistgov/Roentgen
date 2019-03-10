package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
 * @author Nicholas W. M. Ritchie
 *
 */
public class ElementByStoichiometry extends LabeledMultivariateJacobianFunction
		implements ILabeledMultivariateFunction {

	private final Map<Element, Integer> mValences = new HashMap<>();
	private final Element mElement;

	private static List<? extends Object> buildInput(final Material mat, final List<Element> inputElms,
			final Element outputElm) {
		final List<Object> res = new ArrayList<>();
		for (final Element elm : inputElms) {
			res.add(MaterialLabel.buildMassFractionTag(mat, elm));
			res.add(MaterialLabel.buildAtomicWeightTag(mat, elm));
		}
		res.add(MaterialLabel.buildAtomicWeightTag(mat, outputElm));
		return res;
	}

	/**
	 * @param inputLabels
	 * @param outputLabels
	 */
	public ElementByStoichiometry(final Material mat, final List<Element> inputElms, final Element outputElm,
			final Map<Element, Integer> valences) {
		super(buildInput(mat, inputElms, outputElm),
				Collections.singletonList(MaterialLabel.buildMassFractionTag(mat, outputElm)));
		mValences.putAll(valences);
		mElement = outputElm;
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
		Object awit = null;
		double ci = 0.0, dcidawi = 0.0;
		for (int j = 0; j < getInputDimension(); ++j) {
			final Object lbl = getInputLabel(j);
			if (lbl instanceof MassFraction) {
				final MassFraction cjt = (MassFraction) lbl;
				final AtomicWeight awjt = MaterialLabel.buildAtomicWeightTag(cjt.getMaterial(), cjt.getElement());
				if (Double.isNaN(awi)) {
					awit = MaterialLabel.buildAtomicWeightTag(cjt.getMaterial(), mElement);
					awi = getValue(awit, point);
				}
				final double cj = getValue(cjt, point);
				final double awj = getValue(awjt, point);
				final double kkj = -(awi / awj) * (mValences.get(cjt.getElement()) / mValences.get(mElement));
				ci += kkj * cj;
				dcidawi += (kkj / awi) * cj;
				writeJacobian(0, cjt, kkj, rm);
				writeJacobian(0, awjt, -kkj / awj, rm);
			}
			writeJacobian(0, awit, dcidawi, rm);
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
				final AtomicWeight awt = MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mft.getElement());
				if (Double.isNaN(ai)) {
					final Object awe = MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mElement);
					ai = getValue(awe, point);
				}
				final double cj = getValue(mft, point);
				final double aj = getValue(awt, point);
				ci -= (ai / aj) * (mValences.get(mft.getElement()) / mValences.get(mElement)) * cj;
			}
		}
		rv.setEntry(0, ci);
		return rv;
	}

}
