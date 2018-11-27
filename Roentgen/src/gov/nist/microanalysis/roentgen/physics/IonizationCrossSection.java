package gov.nist.microanalysis.roentgen.physics;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * <p>
 * The documentation from the file <code>xion.f</code> as sent to NWMR by Cesc
 * Salvat in early September 2008.
 * </p>
 * <p>
 * This implementation is based on DBPW computations of the ionization cross
 * section. It is likely to be as good or better than anything available
 * elsewhere. It particular there is very little data for the L and M shells and
 * this computation is likely to be the best resource for this information.
 * </p>
 * <p>
 * This function delivers the total cross section for electron impact ionization
 * of K, L and M shells of neutral atoms of the elements from hydrogen (IZ=1) to
 * einsteinium (IZ=99). It uses a parameterization of numerical cross sections
 * computed with the distorted-wave and the plane-wave first-Born
 * approximations.
 * </p>
 * <h4>References:</h4>
 * <ul>
 * <li>D. Bote and F. Salvat, "Calculations of inner-shell ionization by
 * electron impact with the distorted-wave and plane-wave Born approximations",
 * Phys. Rev. A77, 042701 (2008).</li>
 * <li>D. Bote et al., "Cross sections for ionization of K, L and M shells of
 * atoms by impact of electrons and positrons with energies up to 1 GeV.
 * Analytical formulas", At. and Nucl. Data Tables (in preparation).</li>
 * </ul>
 * <h4>Input Arguments:</h4>
 * <ul>
 * <li>EEV ..... kinetic energy of the projectile electron (in eV).</li>
 * <li>IZ ...... atomic number of the target atom (IZ=1 to 99).</li>
 * <li>ISH ..... active target electron shell, 1=K, 2=L1, 3=L2, 4=L3, 5=M1, ...,
 * 9=M5.</li>
 * </ul>
 * <h4>Output value:</h4>
 * <ul>
 * <li>XIONE ... ionization cross section (in cm**2) of the ISH shell.</li>
 * </ul>
 * <p>
 * <i>D. Bote, F. Salvat, A. Jablonski and C.J. Powell (September 2008)</i>
 * </p>
 * <p>
 *
 * @author nicholas
 * @version 1.0
 */
public class IonizationCrossSection //
		implements UnivariateFunction {

	private static class Implementation {

		final static double FPIAB2 = 4.0 * Math.PI * Math.pow(PhysicalConstants.BohrRadiusCGS, 2.0);
		final static double REV = PhysicalConstants.ElectronRestMass_eV;
		private final double[][][] mA; // [100][9][5] - [Z][shell][5]
		private final double[][] mBe; // [100][9] - [Z][shell]
		private final double[][] mAnlj; // [100][9] - [Z][shell]
		private final double[][][] mG; // [100][9][4] - [Z][shell][4]

		private static double[] parse(final String line) {
			final String[] items = line.split(",");
			final double[] res = new double[items.length];
			for (int i = 0; i < items.length; ++i)
				res[i] = Double.parseDouble(items[i]);
			return res;
		}

		private boolean isSupported(final AtomicShell sh) {
			final int z = sh.getElement().getAtomicNumber();
			return (z < mA.length) && (sh.getShell().ordinal() < mA[z].length);

		}

		public Implementation() throws IOException, ParseException {
			{
				final InputStream is = getClass().getResourceAsStream("SalvatXionB.csv");
				final InputStreamReader isr = new InputStreamReader(is);
				final BufferedReader br = new BufferedReader(isr);
				mBe = new double[100][];
				mAnlj = new double[100][];
				mG = new double[100][][];
				final int B_BLOCK_SIZE = 6;
				for (int z = 1; z < 100; ++z) {
					final double[] line = parse(br.readLine());
					assert Math.round(line[0]) == z;
					final int len = line.length - 1;
					assert (len % B_BLOCK_SIZE) == 0;
					final int trCx = len / B_BLOCK_SIZE;
					mBe[z] = new double[trCx];
					mAnlj[z] = new double[trCx];
					mG[z] = new double[trCx][4];
					for (int i = 0; i < len; ++i) {
						final int nn = i / B_BLOCK_SIZE;
						switch (i % B_BLOCK_SIZE) {
						case 0:
							mBe[z][nn] = line[i + 1];
							break;
						case 1:
							mAnlj[z][nn] = line[i + 1];
							break;
						default:
							assert ((i % B_BLOCK_SIZE) - 2) >= 0;
							assert ((i % B_BLOCK_SIZE) - 2) < 4;
							mG[z][nn][(i % B_BLOCK_SIZE) - 2] = line[i + 1];
							break;
						}
					}
				}
			}
			{
				final InputStream is = getClass().getResourceAsStream("SalvatXionA.csv");
				final InputStreamReader isr = new InputStreamReader(is);
				final BufferedReader br = new BufferedReader(isr);
				final double a[][][] = new double[100][][];
				final int A_BLOCK_SIZE = 5;
				for (int z = 1; z < 100; ++z) {
					final double[] line = parse(br.readLine());
					assert Math.round(line[0]) == z;
					final int len = line.length - 1;
					assert (len % A_BLOCK_SIZE) == 0;
					final int trCx = len / A_BLOCK_SIZE;
					a[z] = new double[trCx][A_BLOCK_SIZE];
					for (int i = 0; i < len; ++i)
						a[z][i / A_BLOCK_SIZE][i % A_BLOCK_SIZE] = line[i + 1];
				}
				mA = a;
			}
		}

		/**
		 * Computes the ionization cross section for the specified AtomicShell and
		 * incident electron energy.
		 *
		 * @param atomicNumber Atomic number of the element (1 to 99)
		 * @param shellIdx     Shell to be ionized (K->0 to M5->8)
		 * @param edgeEnergy   Excitation edge energy in eV
		 * @param energy       Electron energy in eV
		 * @return The ionization cross section in cm<sup>2</sup>
		 */
		public double compute( //
				final int atomicNumber, //
				final int shellIdx, //
				final double edgeEnergy, //
				final double energy //
		) {
			double xione = 0.0;
			if ((edgeEnergy >= 1.0) && (energy > edgeEnergy)) {
				final double u = energy / edgeEnergy;
				if (u <= 16.0) {
					assert atomicNumber < mA.length : atomicNumber + ">=" + mA.length;
					assert shellIdx < mA[atomicNumber].length : shellIdx + ">=" + mA[atomicNumber].length;
					final double[] as = mA[atomicNumber][shellIdx];
					final double opo = 1.0 / (1.0 + u);
					final double opo2 = Math.pow(opo, 2.0);
					final double ffitlo = as[0] + (as[1] * u) + (opo * (as[2] + (opo2 * (as[3] + (opo2 * as[4])))));
					xione = (u - 1.0) * Math.pow(ffitlo / u, 2.0); // ok
				} else {
					assert atomicNumber < mBe.length : atomicNumber + " >= " + mBe.length;
					assert shellIdx < mBe[atomicNumber].length : shellIdx + " >= " + mBe[atomicNumber].length + " for "
							+ " for Z=" + atomicNumber + " shell=" + shellIdx;
					assert shellIdx < mG[atomicNumber].length : shellIdx + ">= " + mG[atomicNumber].length + " for "
							+ shellIdx;
					final double beta2 = (energy * (energy + (2.0 * REV))) / Math.pow(energy + REV, 2.0);
					final double x = Math.sqrt(energy * (energy + (2.0 * REV))) / REV;
					final double[] g = mG[atomicNumber][shellIdx];
					final double ffitup = (((2.0 * Math.log(x)) - beta2) * (1.0 + (g[0] / x))) + g[1]
							+ (g[2] * Math.pow((REV * REV) / ((energy + REV) * (energy + REV)), 0.25)) + (g[3] / x);
					final double factr = mAnlj[atomicNumber][shellIdx] / beta2;
					xione = ((factr * u) / (u + mBe[atomicNumber][shellIdx])) * ffitup;
				}
			}
			return FPIAB2 * xione;
		}

		public double compute(final AtomicShell shell, final double energy) {
			return compute(//
					shell.getElement().getAtomicNumber(), //
					shell.getShell().ordinal(), //
					shell.getEdgeEnergy(), //
					energy);
		}

		public TreeSet<AtomicShell> getSupportedAtomicShells(final Element elm) {
			final double[] be = mBe[elm.getAtomicNumber()];
			final TreeSet<AtomicShell> res = new TreeSet<>();
			for (int i = 0; i < be.length; ++i) {
				final AtomicShell sh = AtomicShell.find(elm, Shell.values()[i]);
				if (sh != null)
					res.add(sh);
			}
			return res;
		}

	}

	// Static variable initialization is thread safe on a per ClassLoader basis
	private static final Implementation IMPLEMENTATION = buildImplementation();

	private static Implementation buildImplementation() {
		try {
			return new Implementation();
		} catch (IOException | ParseException e) {
			e.printStackTrace();
		}
		return null;
	}

	private final AtomicShell mShell;

	public IonizationCrossSection(final AtomicShell shell) {
		final int iz = shell.getElement().getAtomicNumber();
		assert ((iz >= Element.Hydrogen.getAtomicNumber()) && (iz <= Element.Einsteinium
				.getAtomicNumber())) : "Unsupported element in BoteSalvatCrossSection.computeShell: "
						+ shell.getElement().getAbbrev();
		assert isSupported(shell) : "Unsupported shell in BoteSalvatCrossSection.computeShell: " + shell.toString();
		mShell = shell;
	}

	/**
	 * Computes the ionization cross section for the specified AtomicShell and
	 * incident electron energy.
	 *
	 * @param shell  AtomicShel
	 * @param energy in eV
	 * @return The ionization cross section in cm<sup>2</sup>
	 */
	public static double value(final AtomicShell shell, final double energy) {
		return IMPLEMENTATION.compute(shell, energy);
	}

	/**
	 * Computes the ionization cross section for the specified AtomicShell and
	 * incident electron energy.
	 *
	 * @param energy in eV
	 * @return The ionization cross section in cm<sup>2</sup>
	 */
	@Override
	public double value(final double energy) {
		return IMPLEMENTATION.compute(mShell, energy);
	}

	public static boolean isSupported(final AtomicShell sh) {
		return IMPLEMENTATION.isSupported(sh);
	}

	public AtomicShell getAtomicShell() {
		return mShell;
	}

	public static Set<AtomicShell> getSupportedAtomicShells(final Element elm) {
		return IMPLEMENTATION.getSupportedAtomicShells(elm);
	}

	@Override
	public String toString() {
		return "IonizationCrossSection[" + mShell.toString() + "]";
	}
}