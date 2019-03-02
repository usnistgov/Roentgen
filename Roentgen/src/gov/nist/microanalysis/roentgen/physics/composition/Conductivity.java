package gov.nist.microanalysis.roentgen.physics.composition;

/**
 * A classification scheme for materials by electrical conductivity.
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
public enum Conductivity {

	/**
	 * Resistivity larger than approximately 10<sup>2</sup>&Omega;/m
	 */
	Insulator,

	/**
	 * Resistivity smaller than approximately 10<sup>2</sup>&Omega;/m but larger
	 * than a conductor
	 */
	Semiconductor,

	/**
	 * Resistivity smaller than approximately 10<sup>-3</sup>&Omega;/m Amorphous
	 * carbon is considered a conductor.
	 */
	Conductor,

	/**
	 * Resistivity of essentially zero &Omega;/m
	 */
	Superconductor,

}
