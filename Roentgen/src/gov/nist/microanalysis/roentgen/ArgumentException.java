package gov.nist.microanalysis.roentgen;

/**
 * Thrown when there is a problem with the argument handed to a function.
 * 
 * @author Nicholas
 */
public class ArgumentException extends Exception {

	private static final long serialVersionUID = 3618920026934015298L;

	/**
	 * @param arg0
	 */
	public ArgumentException(final String arg0) {
		super(arg0);
	}

}
