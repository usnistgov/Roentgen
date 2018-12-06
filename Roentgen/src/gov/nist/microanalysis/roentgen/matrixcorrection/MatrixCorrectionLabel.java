/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author nicholas
 *
 */
public class MatrixCorrectionLabel extends BaseLabel<UnknownMatrixCorrectionDatum, StandardMatrixCorrectionDatum, ElementXRaySet> {

	public MatrixCorrectionLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std, final ElementXRaySet exrs) {
		super("ZAF", unk, std, exrs);
	}
	
	public MatrixCorrectionLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std, final CharacteristicXRay cxr) {
		super("ZAF", unk, std, new ElementXRaySet(cxr));
	}
	
	
	public UnknownMatrixCorrectionDatum getUnknown() {
		return getObject1();
	}
	
	public StandardMatrixCorrectionDatum getStandard() {
		return getObject2();
	}

	public ElementXRaySet getElementXRaySet() {
		return getObject3();
	}
}
