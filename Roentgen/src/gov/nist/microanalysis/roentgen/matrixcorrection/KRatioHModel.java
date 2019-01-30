package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

class KRatioHModel extends LabeledMultivariateJacobianFunction {

	private final UnknownMatrixCorrectionDatum mUnknown;
	private final Map<ElementXRaySet, StandardMatrixCorrectionDatum> mStandards;

	static List<? extends Object> buildOutputs(//
			final Composition unk //
	) {
		return unk.massFractionTags();
	}

	static List<? extends Object> buildInputs(//
			final UnknownMatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds //
	) {
		final List<Object> res = new ArrayList<>();
		for (final Map.Entry<ElementXRaySet, StandardMatrixCorrectionDatum> me : stds.entrySet()) {
			res.add(new KRatioLabel(unk, me.getValue(), me.getKey(), Method.Measured));
			res.add(Composition.buildMassFractionTag(unk.getComposition(), me.getKey().getElement()));
			res.add(Composition.buildMassFractionTag(me.getValue().getComposition(), me.getKey().getElement()));
			res.add(new MatrixCorrectionLabel(unk, me.getValue(), me.getKey()));
		}
		return res;
	}

	public KRatioHModel(//
			final UnknownMatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, StandardMatrixCorrectionDatum> stds //
	) {
		super(buildInputs(unk, stds), ImplicitMeasurementModel.buildHLabels(buildOutputs(unk.getComposition())));
		mUnknown = unk;
		mStandards = stds;
	}

	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector rv = new ArrayRealVector(getOutputDimension());
		final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
		for (final Map.Entry<ElementXRaySet, StandardMatrixCorrectionDatum> me : mStandards.entrySet()) {
			final Composition unk = mUnknown.getComposition();
			final KRatioLabel kMeasTag = new KRatioLabel(mUnknown, me.getValue(), me.getKey(), Method.Measured);
			final Element elm = me.getKey().getElement();
			final Object mfUnkTag = Composition.buildMassFractionTag(unk, elm);
			final Object mfStdTag = Composition.buildMassFractionTag(me.getValue().getComposition(), elm);
			final Object zafTag = new MatrixCorrectionLabel(mUnknown, me.getValue(), me.getKey());
			final Object hTag = new ImplicitMeasurementModel.HLabel(Composition.buildMassFractionTag(unk, elm));
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
}