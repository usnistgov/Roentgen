/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author nicho
 *
 */
public class MatrixCorrectionLabel extends BaseLabel<MatrixCorrectionDatum, MatrixCorrectionDatum, ElementXRaySet> {

	public MatrixCorrectionLabel(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final ElementXRaySet exrs) {
		super("ZAF", unk, std, exrs);
	}
	
	public MatrixCorrectionLabel(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final CharacteristicXRay cxr) {
		super("ZAF", unk, std, new ElementXRaySet(cxr));
	}
	
	
	public MatrixCorrectionDatum getUnknown() {
		return getObject1();
	}
	
	public MatrixCorrectionDatum getStandard() {
		return getObject2();
	}

	public ElementXRaySet getElementXRaySet() {
		return getObject3();
	}
}
