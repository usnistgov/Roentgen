package gov.nist.microanalysis.roentgen.physics.composition;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;

public final class MaterialMassFraction //
		extends BaseLabel<Material, Object, Object> {

	public MaterialMassFraction(final Material mat) {
		super("f", mat);
	}

	public Material getMaterial() {
		return getObject1();
	}
}