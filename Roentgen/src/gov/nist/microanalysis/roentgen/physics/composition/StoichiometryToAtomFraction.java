package gov.nist.microanalysis.roentgen.physics.composition;

import gov.nist.microanalysis.roentgen.math.uncertainty.models.Normalize;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.Stoichiometry;

/**
 * Converts from stoichiometry (Al2O3) to atom fraction f[Al]=0.4, f[O]=0.6
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public class StoichiometryToAtomFraction //
		extends Normalize<Stoichiometry, AtomFraction> {

	/**
	 * Constructs a StoichiometryToAtomFraction
	 *
	 * @param Composition   comp
	 * @param atomicWeights
	 */
	public StoichiometryToAtomFraction(final Material mat) {
		super(MaterialLabel.buildStoichiometryTags(mat), MaterialLabel.buildAtomFractionTags(mat));
	}

	@Override
	public String toString() {
		return "Stoichiometry-to-Atom Fraction";
	}
}