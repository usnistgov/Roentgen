package gov.nist.microanalysis.roentgen.physics;

/**
 * <p>
 * A Shell is an abstraction of the configuration in which an electron can
 * reside in an atom. The Shell is only the configuration, an AtomicShell
 * provides the additional Element information.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 280 $
 */
public enum Shell implements Comparable<Shell> {
	// K-shell
	K(Principle.K, Orbital.S, 1, "K"), //
	// L-shell
	L1(Principle.L, Orbital.S, 1, "L<sub>1</sub>"), //
	L2(Principle.L, Orbital.P, 1, "L<sub>2</sub>"), //
	L3(Principle.L, Orbital.P, 3, "L<sub>3</sub>"), //
	// M-shell
	M1(Principle.M, Orbital.S, 1, "M<sub>1</sub>"), //
	M2(Principle.M, Orbital.P, 1, "M<sub>2</sub>"), //
	M3(Principle.M, Orbital.P, 3, "M<sub>3</sub>"), //
	M4(Principle.M, Orbital.D, 3, "M<sub>4</sub>"), //
	M5(Principle.M, Orbital.D, 5, "M<sub>5</sub>"), //
	// N-shell
	N1(Principle.N, Orbital.S, 1, "N<sub>1</sub>"), //
	N2(Principle.N, Orbital.P, 1, "N<sub>2</sub>"), //
	N3(Principle.N, Orbital.P, 3, "N<sub>3</sub>"), //
	N4(Principle.N, Orbital.D, 3, "N<sub>4</sub>"), //
	N5(Principle.N, Orbital.D, 5, "N<sub>5</sub>"), //
	N6(Principle.N, Orbital.F, 5, "N<sub>6</sub>"), //
	N7(Principle.N, Orbital.F, 7, "N<sub>7</sub>"), //
	// O-shell
	O1(Principle.O, Orbital.S, 1, "O<sub>1</sub>"), //
	O2(Principle.O, Orbital.P, 1, "O<sub>2</sub>"), //
	O3(Principle.O, Orbital.P, 3, "O<sub>3</sub>"), //
	O4(Principle.O, Orbital.D, 3, "O<sub>4</sub>"), //
	O5(Principle.O, Orbital.D, 5, "O<sub>5</sub>"), //
	O6(Principle.O, Orbital.F, 5, "O<sub>6</sub>"), //
	O7(Principle.O, Orbital.F, 7, "O<sub>7</sub>"), //
	O8(Principle.O, Orbital.G, 7, "O<sub>8</sub>"), //
	O9(Principle.O, Orbital.G, 9, "O<sub>9</sub>"), //
	// P-shell
	P1(Principle.P, Orbital.S, 1, "P<sub>1</sub>"), //
	P2(Principle.P, Orbital.P, 1, "P<sub>2</sub>"), //
	P3(Principle.P, Orbital.P, 3, "P<sub>3</sub>"), //
	P4(Principle.P, Orbital.D, 3, "P<sub>4</sub>"), //
	P5(Principle.P, Orbital.D, 5, "P<sub>5</sub>"), //
	P6(Principle.P, Orbital.F, 5, "P<sub>6</sub>"), //
	P7(Principle.P, Orbital.F, 7, "P<sub>7</sub>"), //
	P8(Principle.P, Orbital.G, 7, "P<sub>8</sub>"), //
	P9(Principle.P, Orbital.G, 9, "P<sub>9</sub>"), //
	P10(Principle.P, Orbital.H, 9, "P<sub>10</sub>"), //
	P11(Principle.P, Orbital.H, 11, "P<sub>11</sub>"), //
	// Q-shell
	Q1(Principle.Q, Orbital.S, 1, "Q<sub>1</sub>"), //
	Q2(Principle.Q, Orbital.P, 1, "Q<sub>2</sub>"), //
	Q3(Principle.Q, Orbital.P, 3, "Q<sub>3</sub>"), //
	Q4(Principle.Q, Orbital.D, 3, "Q<sub>4</sub>"), //
	Q5(Principle.Q, Orbital.D, 5, "Q<sub>5</sub>"), //
	Q6(Principle.Q, Orbital.F, 5, "Q<sub>6</sub>"), //
	Q7(Principle.Q, Orbital.F, 7, "Q<sub>7</sub>"), //
	Q8(Principle.Q, Orbital.G, 7, "Q<sub>8</sub>"), //
	Q9(Principle.Q, Orbital.G, 9, "Q<sub>9</sub>"), //
	Q10(Principle.Q, Orbital.H, 9, "Q<sub>10</sub>"), //
	Q11(Principle.Q, Orbital.H, 11, "Q<sub>11</sub>"), //
	Q12(Principle.Q, Orbital.I, 11, "Q<sub>10</sub>"), //
	Q13(Principle.Q, Orbital.I, 13, "Q<sub>11</sub>");

	/**
	 * Principle quantum number
	 */
	private final Principle mPrinciple;
	/**
	 * Orbital angular momentum 0->S, 1->P, 2->D, 3->F
	 */
	private final Orbital mOrbital;
	/**
	 * Total angular momentum (times 2)
	 */
	private final int mTotal;

	private final String mHtml;

	public static final Shell[] K_FAMILY = new Shell[] { Shell.K };

	public static final Shell[] L_FAMILY = new Shell[] { Shell.L1, Shell.L2, Shell.L3 };

	public static final Shell[] M_FAMILY = new Shell[] { Shell.M1, Shell.M2, Shell.M3, Shell.M4, Shell.M5 };

	public static final Shell[] N_FAMILY = new Shell[] { Shell.N1, Shell.N2, Shell.N3, Shell.N4, Shell.N5, Shell.N6,
			Shell.N7 };

	public static final Shell[] O_FAMILY = new Shell[] { Shell.O1, Shell.O2, Shell.O3, Shell.O4, Shell.O5, Shell.O6,
			Shell.O7, Shell.O8, Shell.O9 };

	public static final Shell[] P_FAMILY = new Shell[] { Shell.P1, Shell.P2, Shell.P3, Shell.P4, Shell.P5, Shell.P6,
			Shell.P7, Shell.P8, Shell.P9, Shell.P10, Shell.P11 };

	public enum Orbital {
		S(0), P(1), D(2), F(3), G(4), H(5), I(6);
		private final int mL;

		Orbital(final int l) {
			mL = l;
		}

		public int getOrbitalAngularMomentum() {
			return mL;
		}

	};

	public enum Principle {
		K(1), L(2), M(3), N(4), O(5), P(6), Q(7), R(8);

		private final int mN;

		private Principle(final int n) {
			mN = n;
		}

		public int getPrincipleQuantumNumber() {
			return mN;
		}

		@Override
		public String toString() {
			return super.toString();
		}

	}

	private Shell(final Principle fam, final Orbital orb, final int j, final String html) {
		mPrinciple = fam;
		mOrbital = orb;
		mTotal = j;
		mHtml = html;
	}

	/**
	 * Returns the next shell in the same family as this (or null) L1.nextInFamily()
	 * -> L2 L3.nextInFamily() -> null
	 *
	 * @return Shell
	 */
	public Shell nextInFamily() {
		final int ord = ordinal();
		if (ord < Shell.values().length - 1) {
			final Shell next = Shell.values()[ord + 1];
			if (next.mPrinciple.equals(this.mPrinciple))
				return next;
		}
		return null;
	}

	/**
	 * Returns the next shell in the same family as this (or null) L1.nextInFamily()
	 * -> null L3.nextInFamily() -> L2
	 *
	 * @return Shell
	 */
	public Shell previousInFamily() {
		final int ord = ordinal();
		if (ord > 0) {
			final Shell next = Shell.values()[ord - 1];
			if (next.mPrinciple.equals(this.mPrinciple))
				return next;
		}
		return null;
	}

	public Shell.Orbital getOrbital() {
		return mOrbital;
	}

	public int getPrincipleQuantumNumber() {
		return mPrinciple.getPrincipleQuantumNumber();
	}

	public Shell.Principle getFamily() {
		return mPrinciple;
	}

	public int getOrbitalQuantumNumber() {
		return mOrbital.getOrbitalAngularMomentum();
	}

	public double getTotalAngularMomentum() {
		return mTotal * 0.5;
	}

	public int getCapacity() {
		return 1 + mTotal;
	}

	public static Shell parse(final String str) {
		for (final Shell sh : Shell.values())
			if (sh.name().equals(str))
				return sh;
		return null;
	}

	/**
	 * Returns a string representation in the standard IUPAC form K->"K", L1->
	 * "L<sub>1</sub>", L2->"L<sub>2</sub>"...
	 *
	 * @return String in HTML form.
	 */
	public String toHTML() {
		return mHtml;
	}

	/**
	 * Returns the shell name in standard atomic physics notation: K -> "1S",
	 * L1->"2S", L2->"2P<sub>1/2</sub>", L3->"2P<sub>3/2</sub>", ...
	 *
	 * @return String in HTML
	 */
	public String toAtomicNotation() {
		return Integer.toString(mPrinciple.getPrincipleQuantumNumber()) + mOrbital.toString()
				+ (mOrbital.getOrbitalAngularMomentum() == 0 ? "" : "<sub>" + Integer.toString(mTotal) + "/2</sub>");
	}

}
