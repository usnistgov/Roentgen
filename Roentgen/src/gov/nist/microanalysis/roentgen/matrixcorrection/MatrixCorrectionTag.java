/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author nicho
 *
 */
public class MatrixCorrectionTag extends BaseTag<MatrixCorrectionDatum, MatrixCorrectionDatum, ElementXRaySet> {

	public MatrixCorrectionTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final ElementXRaySet exrs) {
		super("ZAF", unk, std, exrs);
	}
	
	public MatrixCorrectionTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final CharacteristicXRay cxr) {
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
