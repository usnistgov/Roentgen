package gov.nist.microanalysis.roentgen.physics;

import gov.nist.microanalysis.roentgen.physics.Shell.Principle;

/**
 * <p>
 * Implements a very simple mechanism for calculating the jump ratio.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
public class JumpRatio {

	/**
	 * Approximate jump ratios for the lines L1 and L2
	 */
	static private final double L_JUMPS[] = new double[] { 1.16, 1.41 };

	/**
	 * Approximate jump ratios for the lines M1, M2, M3, M4
	 */
	static private final double M_JUMPS[] = new double[] { 1.02, 1.04, 1.13, 1.38 };

	private final static JumpRatio INSTANCE = new JumpRatio();

	public static JumpRatio getInstance() {
		return INSTANCE;
	}

	/**
	 * Computes the total jump ratio for the specified element and line family.
	 *
	 * @param elm    An element
	 * @param family (One of AtomicShell.K_FAMILY, AtomicShell.L_FAMILY or
	 *               AtomicShell.M_FAMILY)
	 * @return The total jump ratio
	 */
	private double compute(final Element elm, final Principle family) {
		double tmp = 0.0, k = 1.0;
		final double z = elm.getAtomicNumber();
		switch (family) {
		case K:
			tmp = 0.924 - (0.00144 * z);
			break;
		case L: {
			k = 1.16 * 1.41;
			tmp = (0.548 - (0.00231 * z));
		}
		case M:
			return 3.56;
		default:
			return 1.0;
		}
		return k / (1.0 - (k * tmp));
	}

	public double compute(final AtomicShell sh) {
		switch (sh.getShell()) {
		case K:
			return compute(sh.getElement(), sh.getFamily());
		case L1:
			return L_JUMPS[0];
		case L2:
			return L_JUMPS[1];
		case L3:
			return compute(sh.getElement(), sh.getFamily()) / (L_JUMPS[0] * L_JUMPS[1]);
		case M1:
			return M_JUMPS[0];
		case M2:
			return M_JUMPS[1];
		case M3:
			return M_JUMPS[2];
		case M4:
			return M_JUMPS[3];
		case M5:
			return compute(sh.getElement(), sh.getFamily()) / (M_JUMPS[0] * M_JUMPS[1] * M_JUMPS[2] * M_JUMPS[3]);
		default:
			return 1.0;
		}
	}

}
