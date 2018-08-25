/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * @author nicho
 *
 */
public class MatrixCorrectionTag extends BaseTag<MatrixCorrectionDatum, MatrixCorrectionDatum, ElementXRaySet> {

	protected MatrixCorrectionTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final ElementXRaySet exrs) {
		super("ZAF", unk, std, exrs);
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
