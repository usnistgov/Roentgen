package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomicWeight;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * Converts compositional data in atom fraction into mass fraction. Requires
 * also the atomic weights for each element in the material.
 *
 * @author nicholas
 *
 */
public class AtomFractionToMassFraction //
		extends LabeledMultivariateJacobianFunction<MaterialLabel, MassFraction> //
		implements ILabeledMultivariateFunction<MaterialLabel, MassFraction> {

	public final Material mMaterial;

	private static List<MaterialLabel> buildInputTags(final Material mat) {
		final List<MaterialLabel> res = new ArrayList<>();
		res.addAll(MaterialLabel.buildAtomFractionTags(mat));
		res.addAll(MaterialLabel.buildAtomicWeightTags(mat));
		return res;
	}

	/**
	 * Constructs a AtomicFractionToMassFraction instance.
	 *
	 * @param String html
	 * @param        Collection&lt;Element&gt; The elements present in the material.
	 */
	public AtomFractionToMassFraction(final Material mat) {
		super(buildInputTags(mat), MaterialLabel.buildMassFractionTags(mat));
		mMaterial = mat;
	}

	private double denom(final RealVector point) {
		final Map<Element, Number> frac = new HashMap<>();
		final Map<Element, Number> wgt = new HashMap<>();

		for (int i = 0; i < getInputDimension(); ++i) {
			final Object tag = getInputLabel(i);
			final double val = point.getEntry(i);
			if (tag instanceof MaterialLabel.AtomFraction)
				frac.put(((MaterialLabel.AtomFraction) tag).getElement(), val);
			else if (tag instanceof MaterialLabel.AtomicWeight)
				wgt.put(((MaterialLabel.AtomicWeight) tag).getElement(), val);
		}
		double res = 0.0;
		for (final Map.Entry<Element, Number> me : frac.entrySet())
			res += me.getValue().doubleValue() * wgt.get(me.getKey()).doubleValue();
		return res;
	}

	private double delta(final MaterialLabel a1, final MaterialLabel a2) {
		return a1.getElement().equals(a2.getElement()) ? 1.0 : 0.0;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final RealMatrix jac = NullableRealMatrix.build(getOutputDimension(), getInputDimension());
		final double den = denom(point);
		for (int o1 = 0; o1 < getOutputDimension(); ++o1) {
			final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) getOutputLabel(o1);
			final int awt1 = inputIndex(MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mft.getElement()));
			final int aft1 = inputIndex(MaterialLabel.buildAtomFractionTag(mft.getMaterial(), mft.getElement()));
			final double w1 = point.getEntry(awt1);
			final double a1 = point.getEntry(aft1);
			final double c1 = a1 * w1 / den;
			vals.setEntry(o1, c1);
			for (int i2 = 0; i2 < getInputDimension(); ++i2) {
				final Object inLabel = getInputLabel(i2);
				if (inLabel instanceof MaterialLabel.AtomFraction) {
					final MaterialLabel.AtomFraction aft = (MaterialLabel.AtomFraction) inLabel;
					final int awt2 = inputIndex(
							MaterialLabel.buildAtomicWeightTag(aft.getMaterial(), aft.getElement()));
					final double w2 = point.getEntry(awt2);
					jac.setEntry(o1, i2, (w1 / den) * (delta(mft, aft) - a1 * w2 / den));
				} else {
					assert inLabel instanceof MaterialLabel.AtomicWeight;
					final MaterialLabel.AtomicWeight awt = (MaterialLabel.AtomicWeight) inLabel;
					final AtomFraction aft2 = MaterialLabel.buildAtomFractionTag(awt.getMaterial(), awt.getElement());
					final int aft2i = inputIndex(aft2);
					final double a2 = point.getEntry(aft2i);
					jac.setEntry(o1, i2, (a1 / den) * (delta(mft, aft2) - w1 * a2 / den));
				}
			}
		}
		return Pair.create(vals, jac);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector vals = new ArrayRealVector(getOutputDimension());
		final double den = denom(point);
		for (int i1 = 0; i1 < getOutputDimension(); ++i1) {
			final MaterialLabel.MassFraction mft = (MaterialLabel.MassFraction) getOutputLabel(i1);
			final AtomicWeight awt1 = MaterialLabel.buildAtomicWeightTag(mft.getMaterial(), mft.getElement());
			final double w1 = point.getEntry(inputIndex(awt1));
			final AtomFraction aft1 = MaterialLabel.buildAtomFractionTag(mft.getMaterial(), mft.getElement());
			final double a1 = point.getEntry(inputIndex(aft1));
			final double c1 = a1 * w1 / den;
			vals.setEntry(i1, c1);
		}
		return vals;
	}

	@Override
	public String toString() {
		return "Atom fraction to Mass fraction[" + mMaterial + "]";
	}

}