package gov.nist.microanalysis.roentgen.utility;


/**
 * <p>
 * A centralized location for making choices about formating specific types of
 * numbers.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2017
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
public class FormatFactory {

	/**
	 * The default NumberFormat instance used to present mass fractions.
	 *
	 * @return NumberFormat
	 */
	public static BasicNumberFormat getMassFractionNumberFormat() {
		final BasicNumberFormat res = new BasicNumberFormat("0.0000");
		BasicNumberFormat.setSmallNumberFormat(1.0e-3, new HalfUpFormat("0.00E0"));
		return res;
	}

	/**
	 * The default NumberFormat instance used to present mass fractions.
	 *
	 * @return NumberFormat
	 */
	public static BasicNumberFormat getAtomFractionNumberFormat() {
		final BasicNumberFormat res = new BasicNumberFormat("0.0000");
		BasicNumberFormat.setSmallNumberFormat(1.0e-3, new HalfUpFormat("0.00E0"));
		return res;
	}

	/**
	 * The default NumberFormat instance used to present k-ratios fractions.
	 *
	 * @return NumberFormat
	 */
	public static BasicNumberFormat getKRatioNumberFormat() {
		final BasicNumberFormat res = new BasicNumberFormat("0.0000");
		BasicNumberFormat.setSmallNumberFormat(1.0e-3, new HalfUpFormat("0.00E0"));
		return res;
	}

}
