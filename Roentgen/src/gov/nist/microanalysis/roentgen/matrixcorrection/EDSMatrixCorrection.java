/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * Takes the matrix corrections for the individual lines and sums them as appropriate
 * for k-ratios measured from multiple lines simultaneously as is the case in EDS.
 * 
 * Takes as inputs
 * <ol>
 * <li>Characteristic x-ray weights (XRayWeightTag)
 * <li>ZAF associated with single characteristic lines (XPPMatrixCorrection.zaTag)
 * </ol>
 * Returns as outputs
 * <ol>
 * <li>The effective ZAF for a set of characteristic x-ray lines from a single element
 * (MatrixCorrectionTag)
 * </ol>
 * @author nicholas
 *
 */
public class EDSMatrixCorrection extends NamedMultivariateJacobianFunction {

	public static class XRayWeightTag extends BaseTag<CharacteristicXRay, Object, Object> {

		public XRayWeightTag(CharacteristicXRay cxr) {
			super("W", cxr);
		}
	}

	private static List<? extends Object> buildInputTags(MatrixCorrectionDatum unk, //
			MatrixCorrectionDatum std, //
			List<ElementXRaySet> lexrs) {
		List<Object> res = new ArrayList<>();
		for (ElementXRaySet exrs : lexrs) {
			for (CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				res.add(new XRayWeightTag(cxr));
				res.add(XPPMatrixCorrection.zaTag(unk, std, cxr));
			}
		}
		return res;
	}

	private static List<? extends Object> buildOutputTags(MatrixCorrectionDatum unk, //
			MatrixCorrectionDatum std, //
			List<ElementXRaySet> lexrs) {
		List<Object> res = new ArrayList<>();
		for (ElementXRaySet exrs : lexrs)
			res.add(new MatrixCorrectionTag(unk, std, exrs));
		return res;
	}

	public EDSMatrixCorrection(MatrixCorrectionDatum unk, MatrixCorrectionDatum std, List<ElementXRaySet> lexrs) {
		super(buildInputTags(unk, std, lexrs), buildOutputTags(unk, std, lexrs));
	}

	@Override
	public Pair<RealVector, RealMatrix> value(RealVector point) {

		RealVector rv = new ArrayRealVector(getOutputDimension());
		RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

		for (int row = 0; row < getOutputDimension(); ++row) {
			Object tag = getOutputTags().get(row);
			assert tag instanceof MatrixCorrectionTag;
			MatrixCorrectionTag mct = (MatrixCorrectionTag) tag;
			ElementXRaySet exrs = mct.getElementXRaySet();
			MatrixCorrectionDatum std = mct.getStandard();
			MatrixCorrectionDatum unk = mct.getUnknown();
			double sumW = 0.0;
			Map<CharacteristicXRay, Pair<Integer, Integer>> idxMap = new HashMap<>();
			for (CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				int iw = inputIndex(new XRayWeightTag(cxr));
				int iz = inputIndex(XPPMatrixCorrection.zaTag(unk, std, cxr));
				idxMap.put(cxr, Pair.create(Integer.valueOf(iw), Integer.valueOf(iz)));
				sumW += point.getEntry(iw);
			}
			double zafEff = 0.0;
			for (CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				Pair<Integer, Integer> pwzi = idxMap.get(cxr);
				int iwi = pwzi.getFirst();
				int izi = pwzi.getSecond();
				final double wi = point.getEntry(iwi);
				final double zi = point.getEntry(izi);
				zafEff += wi * zi;
				rm.setEntry(row, iwi, zi / sumW - wi * zi / Math.pow(sumW, 2.0));
				rm.setEntry(row, izi, wi / sumW);
			}
			zafEff /= sumW;
			rv.setEntry(row, zafEff);
		}
		return Pair.create(rv, rm);
	}

}
