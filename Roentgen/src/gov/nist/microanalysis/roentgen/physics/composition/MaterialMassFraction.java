package gov.nist.microanalysis.roentgen.physics.composition;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;

public final class MaterialMassFraction //
		extends BaseLabel<String, Object, Object> {

	public MaterialMassFraction(final IMaterial comp) {
		this(comp.getHTMLName());
	}

	public MaterialMassFraction(String htmlName) {
		super("f", htmlName);
	}

	public String getHTML() {
		return getObject1();
	}
}