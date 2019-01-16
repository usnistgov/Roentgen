package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.FracTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;

final class CombineMassFractions extends LabeledMultivariateJacobianFunction {

	private final String mHTML;

	static private List<Object> buildOutputs(final String htmlName, final HashMap<Composition, Number> comps) {
		final Set<Element> elms = new HashSet<>();
		for (Composition comp : comps.keySet())
			elms.addAll(comp.getElementSet());
		final List<Object> outTags = new ArrayList<>();
		for (final Element elm : elms)
			outTags.add(Composition.buildMassFractionTag(htmlName, elm));
		return outTags;
	}

	static private List<Object> buildInputs(final HashMap<Composition, Number> comps) {
		final List<Object> inputTags = new ArrayList<>();
		for (Entry<Composition, Number> me : comps.entrySet()) {
			assert me.getKey().getNativeRepresentation() == Representation.MassFraction;
			inputTags.addAll(me.getKey().getLabels());
			inputTags.add(new FracTag(me.getKey()));
		}
		return inputTags;
	}

	@SafeVarargs
	static private HashMap<Composition, Number> convert(final Pair<Composition, Number>... comps) {
		HashMap<Composition, Number> res = new HashMap<>();
		for (Pair<Composition, Number> comp : comps)
			res.put(comp.getFirst(), comp.getSecond());
		return res;
	}

	@SafeVarargs
	public CombineMassFractions(final String htmlName, final Pair<Composition, Number>... comps) {
		super(buildInputs(convert(comps)), buildOutputs(htmlName, convert(comps)));
		mHTML = htmlName;
	}

	public CombineMassFractions(final String htmlName, final HashMap<Composition, Number> comps) {
		super(buildInputs(comps), buildOutputs(htmlName, comps));
		mHTML = htmlName;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final List<? extends Object> inTags = getInputLabels();
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		for (final Object inTag : inTags) {
			if (inTag instanceof FracTag) {
				final FracTag fracTag = (FracTag) inTag;
				final int fIdx = inputIndex(inTag);
				final double fracVal = getValue(inTag, point);
				for (final Object inTag2 : inTags) {
					if ((inTag2 instanceof MassFractionTag)
							&& (((MassFractionTag) inTag2).getHTML().equals(fracTag.getHTML()))) {
						final MassFractionTag mft = (MassFractionTag) inTag2;
						final int mfIdx = inputIndex(mft);
						final double mfVal = getValue(mft, point);
						final int outIdx = outputIndex(Composition.buildMassFractionTag(mHTML, mft.getElement()));
						rv.setEntry(outIdx, rv.getEntry(outIdx) + fracVal * mfVal);
						assert rm.getEntry(outIdx, mfIdx) == 0.0;
						rm.setEntry(outIdx, mfIdx, fracVal);
						rm.setEntry(outIdx, fIdx, rm.getEntry(outIdx, fIdx) + mfVal);
					}
				}
			}
		}
		return Pair.create(rv, rm);
	}
}