package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel2;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioTag.Method;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * @author nicholas
 *
 */
public class KRatioCorrectionModel2 //
		extends ImplicitMeasurementModel2 {

	private static List<? extends Object> buildOutputs(//
			Composition unk //
	) {
		List<Object> res = new ArrayList<>();
		for (Element elm : unk.getElementSet())
			res.add(Composition.buildMassFractionTag(unk, elm));
		return res;
	}

	private static List<? extends Object> buildInputs(//
			final MatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, MatrixCorrectionDatum> stds //
	) {
		List<Object> res = new ArrayList<>();
		for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			res.add(new KRatioTag(unk, me.getValue(), me.getKey(), Method.Measured));
			res.add(Composition.buildMassFractionTag(unk.getComposition(), me.getKey().getElement()));
			res.add(Composition.buildMassFractionTag(me.getValue().getComposition(), me.getKey().getElement()));
			res.add(new MatrixCorrectionTag(unk, me.getValue(), me.getKey()));
		}
		return res;
	}

	private static class KRatioH extends NamedMultivariateJacobianFunction {

		private final MatrixCorrectionDatum mUnknown;
		private final Map<ElementXRaySet, MatrixCorrectionDatum> mStandards;

		public KRatioH(//
				final MatrixCorrectionDatum unk, //
				final Map<ElementXRaySet, MatrixCorrectionDatum> stds //
		) {
			super(buildInputs(unk, stds), buildHTags(buildOutputs(unk.getComposition())));
			mUnknown = unk;
			mStandards = stds;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
			RealVector rv = new ArrayRealVector(getOutputDimension());
			RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : mStandards.entrySet()) {
				final Composition unk = mUnknown.getComposition();
				final KRatioTag kMeasTag = new KRatioTag(mUnknown, me.getValue(), me.getKey(), Method.Measured);
				final Element elm = me.getKey().getElement();
				final Object mfUnkTag = Composition.buildMassFractionTag(unk, elm);
				final Object mfStdTag = Composition.buildMassFractionTag(me.getValue().getComposition(), elm);
				final Object zafTag = new MatrixCorrectionTag(mUnknown, me.getValue(), me.getKey());
				final Object hTag = new hTag(Composition.buildMassFractionTag(unk, elm));
				final int iKMeas = inputIndex(kMeasTag);
				final int iMFUnk = inputIndex(mfUnkTag);
				final int iMFStd = inputIndex(mfStdTag);
				final int iZAF = inputIndex(zafTag);
				final int oHTag = outputIndex(hTag);
				final double kMeas = point.getEntry(iKMeas);
				final double cUnk = point.getEntry(iMFUnk);
				final double cStd = point.getEntry(iMFStd);
				final double zaf = point.getEntry(iZAF);
				final double hi = kMeas - (cUnk/cStd)*zaf;
				rv.setEntry(oHTag, hi);
				rm.setEntry(oHTag, iKMeas, 1.0);
				rm.setEntry(oHTag, iMFUnk, (-1.0/cStd)*zaf);
				rm.setEntry(oHTag, iMFStd, (cUnk/(cStd*cStd))*zaf);
				rm.setEntry(oHTag, iZAF, - (cUnk/cStd));
			}
			return Pair.create(rv, rm);
		}
	}

	public KRatioCorrectionModel2(//
			final MatrixCorrectionDatum unk, //
			final Map<ElementXRaySet, MatrixCorrectionDatum> stds //
	) {
		super(new KRatioH(unk, stds), buildOutputs(unk.getComposition()));
	}
}
