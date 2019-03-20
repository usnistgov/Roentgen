/**
 * 
 */
package gov.nist.microanalysis.roentgen.spectrum;

import gov.nist.microanalysis.roentgen.EPMALabel;

/**
 * @author nicho
 *
 */
public class SpectrumLabel extends EPMALabel.BaseLabel<Integer, Object, Object> {

	
	public static class Filtered extends SpectrumLabel {
		public Filtered(int index) {
			super("F<sub>S</sub>",Integer.valueOf(index));
		}
	}

	public static class Raw extends SpectrumLabel {
		public Raw(int index) {
			super("I", Integer.valueOf(index));
		}
	}
	
	public SpectrumLabel(String name, Integer index) {
		super(name, index);
	}

}
