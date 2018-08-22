/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;

/**
 * @author nicho
 *
 */
public class MatrixCorrectionTag extends BaseTag<MatrixCorrectionDatum, MatrixCorrectionDatum, CharacteristicXRay> {

	protected MatrixCorrectionTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final CharacteristicXRay cxr) {
		super("ZAF", unk, std, cxr);
	}
}
