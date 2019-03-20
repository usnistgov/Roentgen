package gov.nist.microanalysis.roentgen.physics.composition;

import gov.nist.microanalysis.roentgen.math.uncertainty.models.Normalize;
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
	 */
	public MassFractionToNormalizedMassFraction(final Material mat) {
		super(MaterialLabel.buildMassFractionTags(mat), MaterialLabel.buildNormalizedMassFractionTags(mat));
	}

	@Override
	public String toString() {
		return "Mass Fraction-to-Normalized Mass Fraction";
	}
}
