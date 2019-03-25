package gov.nist.microanalysis.roentgen.physics.composition;

import gov.nist.juncertainty.models.Normalize;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.NormalizedMassFraction;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class MassFractionToNormalizedMassFraction extends Normalize<MassFraction, NormalizedMassFraction> {

	/**
	 * Constructs a StoichiometryToAtomFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 * @throws ArgumentException
	 */
	public MassFractionToNormalizedMassFraction(
			final Material mat
	) throws ArgumentException {
		super(MaterialLabel.buildMassFractionTags(mat), MaterialLabel.buildNormalizedMassFractionTags(mat));
	}

	@Override
	public String toString() {
		return "Mass Fraction-to-Normalized Mass Fraction";
	}
}
