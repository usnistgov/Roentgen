package gov.nist.microanalysis.roentgen.physics.composition;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
/**
 * Used on composition data to convert mass fractions into normalized mass fractions.
 * 
 * @author nicholas
 *
 */
public class MassFractionToNormalized //
		extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	/**
	 * Constructs a AtomicFractionToMassFraction
	 *
	 * @param String Name of material
	 * @param Collection<Element> The elements present in the material
	 */
	public MassFractionToNormalized(final Material mat) {
		super(MaterialLabel.buildMassFractionTags(mat), MaterialLabel.buildNormMassFractionTags(mat));
	}

	private double denom(final RealVector point) {
		double res = 0.0;
		for (final Object tag : getInputLabels())
			if (tag instanceof MaterialLabel.MassFraction)
				res += getValue(tag, point);
		return res;
	}
	
	private double delta(final MaterialLabel a1, final MaterialLabel a2) {
		return a1.getElement().equals(a2.getElement()) ? 1.0 : 0.0;
	}


	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getInputDimension(), getOutputDimension());
		final double den = denom(point);
		for (int i1 = 0; i1 < getOutputDimension(); ++i1) {
			final MaterialLabel.NormalizedMassFraction nmft = (MaterialLabel.NormalizedMassFraction) getOutputLabel(i1);
			final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(nmft.getMaterial(), nmft.getElement());
			final double c1 = getValue(mft, point);
			final double n1 = c1 / den;
			vals.setEntry(i1, n1);
			for (int i2 = 0; i2 < getInputDimension(); ++i2) {
				final MaterialLabel.MassFraction mft2 = (MaterialLabel.MassFraction) getInputLabel(i2);
				jac.setEntry(i1, i2, (delta(nmft, mft2) - c1 / den)/den);
			}
		}
		return Pair.create(vals, jac);
	}
	
	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final double den = denom(point);
		for (int i1 = 0; i1 < getOutputDimension(); ++i1) {
			final MaterialLabel.NormalizedMassFraction nmft = (MaterialLabel.NormalizedMassFraction) getOutputLabel(i1);
			final MaterialLabel.MassFraction mft = MaterialLabel.buildMassFractionTag(nmft.getMaterial(), nmft.getElement());
			final double c1 = getValue(mft, point);
			vals.setEntry(i1, c1 / den);
		}
		return vals;
	}
}
