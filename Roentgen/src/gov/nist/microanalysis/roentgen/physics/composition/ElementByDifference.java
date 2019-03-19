package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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
 * Computes a single element by assuming the mass total is unity and that all
 * the missing mass is found in the specified element.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class ElementByDifference //
		extends LabeledMultivariateJacobianFunction<MassFraction, MassFraction> //
		implements ILabeledMultivariateFunction<MassFraction, MassFraction> {

	private static List<MassFraction> buildInputs(final Material mat, final Element elmByDif) {
		final List<MassFraction> res = new ArrayList<>();
		for (final Element elm : mat.getElementSet())
			if (!elm.equals(elmByDif))
				res.add(MaterialLabel.buildMassFractionTag(mat, elm));
		return res;
	}

	public ElementByDifference(final Material mat, final Element elmByDiff) {
		super(buildInputs(mat, elmByDiff),
				Collections.singletonList(MaterialLabel.buildMassFractionTag(mat, elmByDiff)));
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		double sum = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			sum += point.getEntry(i);
		rv.setEntry(0, 1.0 - sum);
		return rv;
	}

	@Override
	public String toString() {
		return ((MassFraction) getOutputLabel(0)).getElement() + "-by-Stoichiometry";
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		double sum = 0.0;
		for (int i = 0; i < point.getDimension(); ++i)
			sum += point.getEntry(i);
		for (int iIdx = 0; iIdx < getInputDimension(); ++iIdx)
			rm.setEntry(0, iIdx, -1.0);
		rv.setEntry(0, 1.0 - sum);
		return Pair.create(rv, rm);

	}

}
