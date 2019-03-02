package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

public abstract class UnmeasuredElementRule //
		extends LabeledMultivariateJacobianFunction {

	protected final SortedSet<Element> mElements;

	private static final List<MassFraction> buildOutputs(Material unk, Collection<Element> elms) {
		List<MassFraction> res = new ArrayList<>();
		for (Element elm : new HashSet<>(elms))
			res.add(MaterialLabel.buildMassFractionTag(unk, elm));
		return res;
	}

	protected UnmeasuredElementRule(Material unk, Collection<Element> elms, List<? extends Object> inputs) {
		super(inputs, buildOutputs(unk, elms));
		mElements = Collections.unmodifiableSortedSet(new TreeSet<>(elms));
	}

	public SortedSet<Element> getElements() {
		return mElements;
	}

	public MassFraction getMFOutput(int i) {
		return (MassFraction) getOutputLabel(i);
	}

	public static class ElementByDifference //
			extends UnmeasuredElementRule implements ILabeledMultivariateFunction {
		
		public ElementByDifference(Material unk, Element outElm) {
			super(unk, Collections.singletonList(outElm), buildInputs(unk, outElm));
		}
		
		private static List<MassFraction> buildInputs(Material unk, Element inElm) {
			List<MassFraction> res = new ArrayList<>(MaterialLabel.buildMassFractionTags(unk));
			assert unk.getElementSet().contains(inElm);
			assert res.contains(MaterialLabel.buildMassFractionTag(unk,inElm));
			res.remove(MaterialLabel.buildMassFractionTag(unk,inElm));
			return res;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
			RealVector rv = new ArrayRealVector(mElements.size());
			RealMatrix rm = MatrixUtils.createRealMatrix(mElements.size(), getInputDimension());
			double sum = 0.0;
			int oIdx = 0;
			for (Element elm : mElements) {
				for (int iIdx = 0; iIdx < getInputDimension(); ++iIdx) {
					MassFraction mf = (MassFraction) getInputLabel(iIdx);
					assert !mf.getElement().equals(elm);
					final double v = getValue(mf, point);
					if (v > 0.0) {
						sum += v;
						rm.setEntry(oIdx, iIdx, 1.0);
					}
				}
				if (sum < 1.0)
					rv.setEntry(oIdx, 1.0 - sum);
				else {
					rv.setEntry(oIdx, 0.0);
					for (int iIdx = 0; iIdx < getInputDimension(); ++iIdx)
						rm.setEntry(oIdx, iIdx, 0.0);
				}
				++oIdx;
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(RealVector point) {
			RealVector rv = new ArrayRealVector(mElements.size());
			double sum = 0.0;
			int i = 0;
			for (Element elm : mElements) {
				for (Object lbl : getInputLabels()) {
					assert lbl instanceof MassFraction;
					assert !((MassFraction) lbl).getElement().equals(elm);
					sum += Math.max(0.0, getValue(lbl, point));
				}
				rv.setEntry(i, Math.max(1.0 - sum, 0.0));
				++i;
			}
			return rv;
		}

	}

}