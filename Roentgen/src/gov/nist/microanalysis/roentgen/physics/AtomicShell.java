package gov.nist.microanalysis.roentgen.physics;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.duckandcover.html.IToHTML;
import com.google.common.base.Objects;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.physics.Shell.Principle;

/**
 * <p>
 * An AtomicShell object describes the Shell location of an electron within an
 * atom of a specific Element. The AtomicShell objects only represent the filled
 * shells in a ground state atom.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 285 $
 */
public class AtomicShell implements Comparable<AtomicShell>, IToHTML {

	private final Element mElement;
	private final Shell mShell;
	private final double mEdgeEnergy;

	/**
	 * Chantler edge-energies are consistent with the MACs in the Chantler based
	 * implementation of {@link ElementalMAC}. The alternative is ENDF6 edge
	 * energies which is consistent with the ENDF6 database used for relaxation
	 * rates. It is probably better to use Chantler.
	 */

	static private Map<Element, TreeMap<Shell, AtomicShell>> mData = loadData();

	private static Map<Element, TreeMap<Shell, AtomicShell>> loadData() {
		try {
			final Map<Element, TreeMap<Shell, AtomicShell>> res = new HashMap<>();
			try (final InputStream is = AtomicShell.class.getResourceAsStream("Chantler2005EdgeData.txt")) {
				try (final BufferedReader br = new BufferedReader(new InputStreamReader(is))) {
					String line = br.readLine();
					while (line != null) {
						line = br.readLine();
						if (line != null) {
							final String[] items = line.split(" ");
							assert (items.length % 2) == 1;
							final Element elm = Element.parse(items[0]);
							final TreeMap<Shell, AtomicShell> elmData = new TreeMap<>();
							for (int i = 1; i < items.length; i += 2) {
								final Shell sh = Shell.valueOf(items[i]);
								final Double ev = Double.parseDouble(items[i + 1]);
								elmData.put(sh, new AtomicShell(elm, sh, ev));
							}
							res.put(elm, elmData);
						}
					}
				}
			}
			addExtra(res);
			return Collections.unmodifiableMap(res);
		} catch (final Exception e) {
			Globals.getLogger().error("Error loading AtomicShell data.", e);
			throw new Error(e);
		}
	}

	private static void addExtra(final Map<Element, TreeMap<Shell, AtomicShell>> res, final Element elm, final Shell sh,
			final double energy) {
		res.get(elm).put(sh, new AtomicShell(elm, sh, energy));

	}

	/**
	 * <p>
	 * There are a number of element's for whom there are no edge energies for close
	 * to valence shells in the Chantler database. For these, I've used the numbers
	 * in:
	 * </p>
	 * <p>
	 * http://www.chemistry.uoguelph.ca/educmat/atomdata/bindener/grp2num.htm
	 * </p>
	 * <p>
	 * The values are all small (near valence) so the exact value isn't critical and
	 * these would seem to be close enough for our purposes.
	 * </p>
	 *
	 * @param res
	 */
	private static void addExtra(final Map<Element, TreeMap<Shell, AtomicShell>> res) {
		addExtra(res, Element.Magnesium, Shell.M1, 7.64624); // ok
		addExtra(res, Element.Aluminum, Shell.M2, 10.62); // ok
		addExtra(res, Element.Aluminum, Shell.M3, 5.98577); // ok
		addExtra(res, Element.Silicon, Shell.M3, 8.15169); // ok
		addExtra(res, Element.Potassium, Shell.N1, 4.34066); // ok
		addExtra(res, Element.Calcium, Shell.N1, 6.11316); // ok
		addExtra(res, Element.Scandium, Shell.M4, 8.06667); // ok
		addExtra(res, Element.Scandium, Shell.M5, 8.06667); // ok
		addExtra(res, Element.Scandium, Shell.N1, 6.56144); // ok
		addExtra(res, Element.Titanium, Shell.M4, 8.1); // ok
		addExtra(res, Element.Titanium, Shell.M5, 8.1); // ok
		addExtra(res, Element.Titanium, Shell.N1, 6.8282); // ok
		addExtra(res, Element.Vanadium, Shell.M5, 8.13333); // ok
		addExtra(res, Element.Vanadium, Shell.N1, 6.7463); // ok
		addExtra(res, Element.Chromium, Shell.N1, 6.76664); // ok
		addExtra(res, Element.Manganese, Shell.N1, 7.43402); // ok
		addExtra(res, Element.Iron, Shell.N1, 7.9024); // ok
		addExtra(res, Element.Cobalt, Shell.N1, 7.881); // ok
		addExtra(res, Element.Nickel, Shell.N1, 7.6398); // ok
		addExtra(res, Element.Copper, Shell.N1, 7.7264); // ok
		addExtra(res, Element.Zinc, Shell.N1, 9.3941); // ok
		addExtra(res, Element.Gallium, Shell.N1, 11.); // ok
		addExtra(res, Element.Gallium, Shell.N2, 5.9993); // ok
		addExtra(res, Element.Gallium, Shell.N3, 5.9993); // ok
		addExtra(res, Element.Germanium, Shell.N1, 14.3); // ok
		addExtra(res, Element.Germanium, Shell.N2, 7.9); // ok
		addExtra(res, Element.Germanium, Shell.N3, 7.9); // ok
		addExtra(res, Element.Arsenic, Shell.N1, 17.); // ok
		addExtra(res, Element.Arsenic, Shell.N2, 9.8152); // ok
		addExtra(res, Element.Arsenic, Shell.N3, 9.8152); // ok
		addExtra(res, Element.Selenium, Shell.N1, 20.15); // ok
		addExtra(res, Element.Selenium, Shell.N2, 9.7524); // ok
		addExtra(res, Element.Selenium, Shell.N3, 9.7524); // ok
		addExtra(res, Element.Bromine, Shell.N1, 23.8); // ok
		addExtra(res, Element.Bromine, Shell.N2, 11.814); // ok
		addExtra(res, Element.Bromine, Shell.N3, 11.814); // ok
		addExtra(res, Element.Krypton, Shell.N1, 27.51); // ok
		addExtra(res, Element.Krypton, Shell.N2, 14); // ok
		addExtra(res, Element.Krypton, Shell.N3, 14); // ok
		addExtra(res, Element.Rubidium, Shell.N1, 32.4); // ok
		addExtra(res, Element.Rubidium, Shell.N2, 16.5167); // ok
		addExtra(res, Element.Rubidium, Shell.N3, 16.5167); // ok
		addExtra(res, Element.Yttrium, Shell.N4, 6.38); // ok
		addExtra(res, Element.Yttrium, Shell.N5, 6.38); // ok
		addExtra(res, Element.Zirconium, Shell.N5, 8.61); // ok
		addExtra(res, Element.Niobium, Shell.N5, 7.17); // ok
		addExtra(res, Element.Silver, Shell.O1, 7.5762); // ok
		addExtra(res, Element.Cadmium, Shell.O1, 8.9937); // ok
		addExtra(res, Element.Indium, Shell.O1, 10.0); // ok
		addExtra(res, Element.Indium, Shell.O2, 5.7864); // ok
		addExtra(res, Element.Indium, Shell.O3, 5.7864); // ok
		addExtra(res, Element.Tin, Shell.O1, 12.0); // ok
		addExtra(res, Element.Tin, Shell.O2, 7.3438); // ok
		addExtra(res, Element.Tin, Shell.O3, 7.3438); // ok
		addExtra(res, Element.Antimony, Shell.O1, 15.); // ok
		addExtra(res, Element.Antimony, Shell.O2, 8.64); // ok
		addExtra(res, Element.Antimony, Shell.O3, 8.64); // ok
		addExtra(res, Element.Tellurium, Shell.O1, 17.84); // ok
		addExtra(res, Element.Tellurium, Shell.O2, 9.0096); // ok
		addExtra(res, Element.Tellurium, Shell.O3, 9.0096);// ok
		addExtra(res, Element.Iodine, Shell.O1, 20.61); // ok
		addExtra(res, Element.Iodine, Shell.O2, 10.451); // ok
		addExtra(res, Element.Iodine, Shell.O3, 10.451); // ok
		addExtra(res, Element.Xenon, Shell.O1, 23.39); // ok
		addExtra(res, Element.Xenon, Shell.O2, 12.563); // ok
		addExtra(res, Element.Xenon, Shell.O3, 12.563); // ok
		addExtra(res, Element.Cerium, Shell.N7, 6.); // ok
		addExtra(res, Element.Praseodymium, Shell.N7, 6.); // ok
		addExtra(res, Element.Neodymium, Shell.N7, 6.); // ok
		addExtra(res, Element.Promethium, Shell.N7, 6.); // ok
		addExtra(res, Element.Samarium, Shell.N7, 6.); // ok
		addExtra(res, Element.Europium, Shell.N7, 6.); // ok
		addExtra(res, Element.Lutetium, Shell.O4, 6.6); // ok
		addExtra(res, Element.Lutetium, Shell.O5, 6.6); // ok
		addExtra(res, Element.Hafnium, Shell.O5, 7.); // ok
		addExtra(res, Element.Tantalum, Shell.O5, 8.3); // ok
		addExtra(res, Element.Tungsten, Shell.O5, 9.); // ok
		addExtra(res, Element.Platinum, Shell.P1, 9.); // ok
		addExtra(res, Element.Gold, Shell.P1, 9.2257); // ok
	}

	/**
	 * Returns true if AtomicShell data (edge energy) is available.
	 *
	 * @param elm
	 * @return true if data is available; false otherwise.
	 */
	public boolean dataAvailable(final Element elm) {
		return mData.containsKey(elm);
	}

	/**
	 * Returns Set of Element objects for which AtomicShell data is available.
	 *
	 * @return A Set&lt;Element&gt; for which AtomicShell data is available.
	 */
	public Set<Element> availableElements() {
		return mData.keySet();
	}

	private AtomicShell(final Element elm, final Shell shell, final double energy) {
		mElement = elm;
		mShell = shell;
		mEdgeEnergy = energy;
	}

	/**
	 * Returns the Shell associated with this AtomicShell
	 *
	 * @return Shell
	 */
	public Shell getShell() {
		return mShell;
	}

	/**
	 * Returns the Element associated with this AtomicShell
	 *
	 * @return Element
	 */
	public Element getElement() {
		return mElement;
	}

	/**
	 * Is this Shell non-empty (contains an electron) in the ground state for the
	 * specified Element?
	 *
	 * @param elm   Element
	 * @param shell Shell
	 * @return true if this shell contains at least one electron in the ground
	 *         state; false otherwise.
	 */
	public static boolean nonEmpty(final Element elm, final Shell shell) {
		return find(elm, shell) != null;
	}

	/**
	 * Returns the AtomicShell object representing the specific Element and Shell or
	 * null if none exists. Use this instead of the private AtomicShell constructor.
	 *
	 * @param elm   Element
	 * @param shell Shell
	 * @return AtomicShell or null
	 */
	public static AtomicShell find(final Element elm, final Shell shell) {
		assert mData.get(elm) != null : elm;
		return mData.get(elm).get(shell);
	}

	public static AtomicShell parse(final String str) {
		final String[] strs = str.split(" ");
		AtomicShell res = null;
		if (strs.length == 2) {
			final Element elm = Element.parse(strs[0]);
			final Shell sh = Shell.valueOf(strs[1]);
			res = find(elm, sh);
		}
		return res;
	}

	/**
	 * Returns an array of the AtomicShell objects associated with the specified
	 * Element.
	 *
	 * @param elm
	 * @return AtomicShell[]
	 */
	public static AtomicShell[] forElement(final Element elm) {
		final Collection<AtomicShell> values = mData.get(elm).values();
		return values != null ? values.toArray(new AtomicShell[values.size()]) : new AtomicShell[0];
	}

	/**
	 * Returns an array of the AtomicShell objects associated with the specified
	 * Element.
	 *
	 * @param elm
	 * @param eMax Maximum edge energy
	 * @return List&lt;AtomicShell&gt; sorted by edge energy (high energy first)
	 */
	public static List<AtomicShell> forElement(final Element elm, final double eMax) {
		final ArrayList<AtomicShell> res = new ArrayList<>();
		for (final AtomicShell as : mData.get(elm).values())
			if (as.mEdgeEnergy <= eMax)
				res.add(as);
		res.sort(new Comparator<AtomicShell>() {
			@Override
			public int compare(final AtomicShell o1, final AtomicShell o2) {
				return Double.compare(o2.mEdgeEnergy, o1.mEdgeEnergy);
			}
		});
		return Collections.unmodifiableList(res);
	}

	/**
	 * Returns the AtomicShell for the specified Element that is the closest in edge
	 * energy below the specified ionizationE.
	 *
	 * @param elm
	 * @param ionizationE
	 * @return AtomicShell
	 */
	public static AtomicShell closestIonizable(final Element elm, final double ionizationE) {
		AtomicShell res = null;
		double bestE = 0.0;
		for (final AtomicShell sh : mData.get(elm).values()) {
			final double ee = sh.mEdgeEnergy;
			if ((ee > bestE) && (ee < ionizationE)) {
				res = sh;
				bestE = ee;
			}
		}
		return res;
	}

	/**
	 * Returns a list of AtomicShell in fam below energy eMax (sorted high to low)
	 *
	 * @param elm  {@link Element}
	 * @param fam  {@link Principle}
	 * @param eMax in eV
	 * @return List<AtomicShell> List of AtomicShell in fam below energy eMax
	 *         (sorted high to low)
	 */
	public static List<AtomicShell> forElement(final Element elm, final Principle fam, final double eMax) {
		final ArrayList<AtomicShell> res = new ArrayList<>();
		for (final AtomicShell ash : mData.get(elm).values())
			if (ash.getFamily().equals(fam) && (ash.mEdgeEnergy < eMax))
				res.add(ash);
		return Collections.unmodifiableList(res);
	}

	/**
	 * Returns an array of the AtomicShell objects associated with the specified
	 * Element and the specified set of shells.
	 *
	 * @param elm
	 * @param shells
	 * @return AtomicShell[]
	 */
	public static AtomicShell[] forElement(final Element elm, final Shell[] shells) {
		final ArrayList<AtomicShell> res = new ArrayList<>();
		final TreeMap<Shell, AtomicShell> atSh = mData.get(elm);
		for (final Shell sh : shells)
			if (atSh.containsKey(sh))
				res.add(atSh.get(sh));
		return res.toArray(new AtomicShell[res.size()]);
	}

	/**
	 * Returns the edge energy for this AtomicShell
	 *
	 * @return double The edge energy in eV
	 */
	public double getEdgeEnergy() {
		return mEdgeEnergy;
	}

	/**
	 * Uses element and shell to create a hash code.
	 *
	 * @return int
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hashCode(mElement, mShell);
	}

	/**
	 * Compares element and shell.
	 *
	 * @param obj
	 * @return boolean
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final AtomicShell other = (AtomicShell) obj;
		return Objects.equal(mElement, other.mElement) && //
				Objects.equal(mShell, other.mShell);
	}

	/**
	 * Comparison is done by element then shell.
	 *
	 * @param o
	 * @return int
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(final AtomicShell o) {
		int res = mElement.compareTo(o.mElement);
		if (res == 0)
			res = mShell.compareTo(o.mShell);
		return res;
	}

	/**
	 * The standard IUPAC name of the shell in HTML.
	 *
	 * @return String in HTML
	 */
	@Override
	public String toHTML(final IToHTML.Mode mode) {
		switch (mode) {
		default:
		case TERSE:
		case NORMAL:
			return mElement.getAbbrev() + "&nbsp;" + mShell.toHTML();
		case VERBOSE: {
			final StringBuffer sb = new StringBuffer();
			sb.append("<table><tr><th>Property</th><th>Value</th>");
			sb.append("<tr><td>Name</td><td>" + toHTML(Mode.TERSE) + "</td>");
			sb.append("<tr><td>Family</td><td>" + mShell.getFamily().toString() + "</td>");
			sb.append("<tr><td>Notation</td><td>" + mShell.toAtomicNotation() + "</td>");
			sb.append("<tr><td>Capacity</td><td>" + mShell.getCapacity() + "</td>");
			// sb.append("<tr><td>Occupancy</td><td>" + mElement. + "</td>");
			sb.append("<tr><td>Orbital Angular<br/>Momentum</td><td>" + mShell.getOrbital().getOrbitalAngularMomentum()
					+ "</td>");
			sb.append("<tr><td>Total Angular<br/>Momentum</td><td>" + mShell.getTotalAngularMomentum() + "</td>");
			sb.append("</table>");
			return sb.toString();
		}
		}

	}

	public Principle getFamily() {
		return mShell.getFamily();
	}

	public static boolean sameFamily(final AtomicShell as1, final AtomicShell as2) {
		return as1.getFamily().equals(as2.getFamily());
	}

	public AtomicShell nextInFamily() {
		final Shell next = mShell.nextInFamily();
		return next != null ? find(mElement, next) : null;
	}

	public AtomicShell previousInFamily() {
		final Shell prev = mShell.previousInFamily();
		return prev != null ? find(mElement, prev) : null;
	}

	/**
	 * The name of the shell
	 *
	 * @return
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return mElement.getAbbrev() + " " + mShell.toString();
	}

	public double getJumpRatio() {
		return JumpRatio.getInstance().compute(this);
	}

	public double getIonizationFraction() {
		final double r = JumpRatio.getInstance().compute(this);
		final double res = r >= 1.0 ? (r - 1.0) / r : 0.0;
		assert (res >= 0.0) && (res <= 1.0) : res;
		return res;
	}

}
