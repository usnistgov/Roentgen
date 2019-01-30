package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

class KRatioHModel2 //
		extends LabeledMultivariateJacobianFunction {

	private final Set<KRatioLabel> mKRatios;

	static List<? extends Object> buildOutputs(//
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		Composition unk = null;
		for (final KRatioLabel krl : kratios) {
			if (unk == null)
				unk = krl.getUnknown().getComposition();
			if (!unk.equals(krl.getUnknown().getComposition()))
				throw new ArgumentException("More than one unknown in KRatioHModel2");
		}
		return unk.massFractionTags();
	}

	static List<? extends Object> buildInputs(//
			final Set<KRatioLabel> kratios //
	) {
		final List<Object> res = new ArrayList<>();
		for (final KRatioLabel krl : kratios) {
			res.add(krl);
			res.add(Composition.buildMassFractionTag(krl.getUnknown().getComposition(), krl.getXRaySet().getElement()));
			res.add(Composition.buildMassFractionTag(krl.getStandard().getComposition(),
					krl.getXRaySet().getElement()));
			res.add(new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet()));
		}
		return res;
	}

	public KRatioHModel2(//
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		super(buildInputs(kratios), ImplicitMeasurementModel.buildHLabels(buildOutputs(kratios)));
		mKRatios = new HashSet<>(kratios);
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		for (final KRatioLabel kMeasTag : mKRatios) {
			final Element elm = kMeasTag.getXRaySet().getElement();
			final Object mfUnkTag = Composition.buildMassFractionTag(kMeasTag.getUnknown().getComposition(), elm);
			final Object mfStdTag = Composition.buildMassFractionTag(kMeasTag.getStandard().getComposition(), elm);
			final Object zafTag = new MatrixCorrectionLabel(kMeasTag.getUnknown(), kMeasTag.getStandard(),
					kMeasTag.getXRaySet());
			final Object hTag = new ImplicitMeasurementModel.HLabel(
					Composition.buildMassFractionTag(kMeasTag.getUnknown().getComposition(), elm));
			final int iKMeas = inputIndex(kMeasTag);
			final int iMFUnk = inputIndex(mfUnkTag);
			final int iMFStd = inputIndex(mfStdTag);
			final int iZAF = inputIndex(zafTag);
			final int oHTag = outputIndex(hTag);
			final double kMeas = point.getEntry(iKMeas);
			final double cUnk = point.getEntry(iMFUnk);
			final double cStd = point.getEntry(iMFStd);
			final double zaf = point.getEntry(iZAF);
			final double hi = kMeas - (cUnk / cStd) * zaf;
			rv.setEntry(oHTag, hi);
			rm.setEntry(oHTag, iKMeas, 1.0);
			rm.setEntry(oHTag, iMFUnk, (-1.0 / cStd) * zaf);
			rm.setEntry(oHTag, iMFStd, (cUnk / (cStd * cStd)) * zaf);
			rm.setEntry(oHTag, iZAF, -(cUnk / cStd));
		}
		return Pair.create(rv, rm);
	}

	@Override
	public String toString() {
		return "K-Ratio H-model";
	}
}