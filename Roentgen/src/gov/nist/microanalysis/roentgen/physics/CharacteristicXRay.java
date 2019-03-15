package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.apache.commons.lang3.concurrent.ConcurrentException;
import org.apache.commons.lang3.concurrent.LazyInitializer;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.math.Utility;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Shell.Principle;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * Represents a characteristic x-ray, a type of generic x-ray. X-rays are
 * ordered by energy.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 307 $
 */
public class CharacteristicXRay extends XRay implements IToHTML, Comparable<XRay> {

	private final XRayTransition mTransition;
	private final Element mElement;
	private final LazyInitializer<Double> mWeight = new LazyInitializer<Double>() {

		@Override
		protected Double initialize() {
			final Principle family = mTransition.getInner().getFamily();
			final Set<AtomicShell> shells = XRayEmissionProbability.getAvailableShells(mElement);
			double e0 = 0.0;
			for (final AtomicShell shell : shells)
				if (shell.getFamily().equals(family))
					e0 = Math.max(e0, shell.getEdgeEnergy());
			final Map<CharacteristicXRay, Double> weights = XRayEmissionProbability.getWeights(mElement, 2.5 * e0);
			final Map<CharacteristicXRay, Double> all = new HashMap<>();
			for (final Map.Entry<CharacteristicXRay, Double> me : weights.entrySet()) {
				final CharacteristicXRay cxr = me.getKey();
				if (cxr.getInner().getFamily().equals(family))
					all.put(cxr, all.getOrDefault(cxr, 0.0) + me.getValue());
			}
			return (all.size() > 0) && (all.containsKey(CharacteristicXRay.this))
					? all.get(CharacteristicXRay.this) / Utility.max(all.values())
					: 0.0;
		}

	};

	private static double getEnergy(final Element elm, final XRayTransition tr) {
		return AtomicShell.find(elm, tr.getInner()).getEdgeEnergy()
				- AtomicShell.find(elm, tr.getOuter()).getEdgeEnergy();
	}

	public double getEdgeEnergy() {
		return (AtomicShell.find(mElement, mTransition.getInner())).getEdgeEnergy();
	}

	public AtomicShell getInner() {
		return AtomicShell.find(mElement, mTransition.getInner());
	}

	public AtomicShell getOuter() {
		return AtomicShell.find(mElement, mTransition.getOuter());
	}

	/**
	 * Returns the K, L, M etc family in which this transition is found.
	 *
	 * @return {@link com.duckandcover.roentgen.physics.Shell.Principle}
	 */
	public Principle getFamily() {
		return mTransition.getInner().getFamily();
	}

	public static CharacteristicXRay find(final CharacteristicXRay[] xrays, final XRayTransition tr) {
		for (final CharacteristicXRay xray : xrays)
			if (xray.mTransition.equals(tr))
				return xray;
		return null;
	}

	/**
	 * Returns the relative line weight (overvoltage = 2.5) for this line relative
	 * to the other lines in the same family.
	 *
	 * @return double
	 */
	public double getWeight() {
		try {
			return mWeight.get().doubleValue();
		} catch (final ConcurrentException e) {
			Globals.getLogger().catching(e);
			System.exit(1);
		}
		return 0.0;
	}

	/**
	 * Returns the relative line weight (overvoltage = 2.5) for this line relative
	 * to the other lines in the same family with estimated uncertainty.
	 *
	 * @return {@link UncertainValue}
	 */
	public UncertainValue getWeightUV() {
		final double w = getWeight();
		double dw;
		if (w > 0.5)
			dw = 0.05 * w;
		else if (w > 0.1)
			dw = 0.1 * w;
		else if (w > 0.01)
			dw = 0.2 * w;
		else
			dw = 0.5 * w;
		return new UncertainValue(w, dw);
	}

	/**
	 * Does this transition exist in the energy database?
	 *
	 * @param elm
	 * @param tr
	 * @return boolean
	 */
	public static boolean exists(final Element elm, final XRayTransition tr) {
		return (AtomicShell.find(elm, tr.getInner()) != null) && (AtomicShell.find(elm, tr.getOuter()) != null);
	}

	/**
	 * Returns all the CharacteristicXray objects which the specified element could
	 * emit based on the
	 *
	 * @param elm
	 * @return {@link CharacteristicXray} array
	 */
	public static CharacteristicXRay[] forElement(final Element elm) {
		final ArrayList<CharacteristicXRay> res = new ArrayList<>();
		for (final XRayTransition tr : XRayTransition.values())
			if (exists(elm, tr))
				res.add(CharacteristicXRay.create(elm, tr));
		// Special for Li & Be
		if ((elm == Element.Lithium) || (elm == Element.Beryllium))
			res.add(CharacteristicXRay.create(elm, XRayTransition.KL1));
		return res.toArray(new CharacteristicXRay[res.size()]);
	}

	/**
	 * Returns all the CharacteristicXray objects which the specified element could
	 * emit based on the
	 *
	 * @param elm
	 * @return {@link CharacteristicXray} array
	 */
	public static CharacteristicXRay[] forAtomicShell(final AtomicShell ashell) {
		final ArrayList<CharacteristicXRay> res = new ArrayList<>();
		final Element elm = ashell.getElement();
		for (final XRayTransition tr : XRayTransition.values())
			if ((tr.getInner().compareTo(ashell.getShell()) == 0) && exists(elm, tr))
				res.add(CharacteristicXRay.create(elm, tr));
		// Special for Li & Be
		if (((elm == Element.Lithium) || (elm == Element.Beryllium)) && (ashell.getShell() == Shell.K))
			res.add(CharacteristicXRay.create(elm, XRayTransition.KL1));
		return res.toArray(new CharacteristicXRay[res.size()]);
	}

	/**
	 * Returns all the characteristic x-rays for the specified element between the
	 * max and min energies specified.
	 *
	 * @param elm  Element
	 * @param minE double (eV)
	 * @param maxE double (eV)
	 * @return {@link CharacteristicXray} array
	 */
	public static CharacteristicXRay[] forElement(final Element elm, final double minE, final double maxE) {
		final ArrayList<CharacteristicXRay> res = new ArrayList<>();
		for (final XRayTransition tr : XRayTransition.values())
			if (exists(elm, tr)) {
				final CharacteristicXRay xr = CharacteristicXRay.create(elm, tr);
				if ((xr.getEnergy() >= minE) && (xr.getEnergy() <= maxE))
					res.add(xr);
			}
		// Special for Li & Be
		if ((res.size() == 0) && (exists(elm, XRayTransition.KL1))) {
			final CharacteristicXRay xr = CharacteristicXRay.create(elm, XRayTransition.KL1);
			if ((xr.getEnergy() >= minE) && (xr.getEnergy() <= maxE))
				res.add(xr);
		}
		return res.toArray(new CharacteristicXRay[res.size()]);
	}

	public static CharacteristicXRay create(final Element elm, final XRayTransition tr) {
		CharacteristicXRay res = null;
		try {
			final double energy = getEnergy(elm, tr);
			res = new CharacteristicXRay(elm, tr, energy);
		} catch (final Exception e) {
			// System.out.println("Energy unavailable for " + elm.getAbbrev() + " " + tr.toString());
		}
		return res;
	}

	private CharacteristicXRay(final Element elm, final XRayTransition tr, final double energy) {
		super(energy);
		assert exists(elm, tr) : "The transition " + elm + " " + tr + " does not exist.";
		mElement = elm;
		mTransition = tr;
	}

	/**
	 * Gets the current value assigned to element
	 *
	 * @return Returns the element.
	 */
	public Element getElement() {
		return mElement;
	}

	/**
	 * Returns the Transition which produced this x-ray.
	 *
	 * @return Transition
	 */
	public XRayTransition getTransition() {
		return mTransition;
	}

	/**
	 * A string describing this x-ray in a form like "Fe K-L3"
	 *
	 * @return String
	 * @see com.duckandcover.roentgen.physics.XRay#toString()
	 */
	@Override
	public String toString() {
		return mElement.getAbbrev() + " " + mTransition.toString();
	}

	static public CharacteristicXRay parse(final String str) {
		final int pos = str.indexOf(" ");
		try {
			final Element elm = Element.parse(str.substring(0, pos));
			final XRayTransition xrt = XRayTransition.parse(str.substring(pos + 1));
			return (elm != null) && (xrt != null) && CharacteristicXRay.exists(elm, xrt)
					? CharacteristicXRay.create(elm, xrt)
					: null;
		} catch (final Exception e) {
			return null;
		}
	}

	/**
	 * The name of this characteristic x-ray in HTML form.
	 *
	 * @return
	 * @see com.duckandcover.roentgen.physics.XRay#toHTML()
	 */
	@Override
	public String toHTML(final IToHTML.Mode mode) {
		switch (mode) {
		default:
		case TERSE:
			return mElement.getAbbrev() + "&nbsp;" + mTransition.toHTML();
		case NORMAL: {
			final Table t = new Table();
			final BasicNumberFormat bnf = new BasicNumberFormat("#,##0.0");
			t.addRow(Table.td(toHTML(Mode.TERSE) + "&nbsp;(" + getTransition().toSeigbahn() + ")"),
					Table.td(bnf.formatHTML(getEnergy()) + " eV"));
			return t.toHTML(Mode.VERBOSE);
		}
		case VERBOSE: {
			final Table t = new Table();
			final BasicNumberFormat bnf = new BasicNumberFormat("#,##0.0");
			t.addRow(Table.th("Name"), Table.th("Seigbahn"), Table.th("Inner"), Table.th("Outer"), Table.th("Weight"),
					Table.th("Energy"));
			t.addRow(Table.td(toHTML(Mode.TERSE)), //
					Table.td(mElement.getAbbrev() + " " + getTransition().toSeigbahn()), //
					Table.td(getInner().toHTML(Mode.TERSE)), //
					Table.td(getOuter().toHTML(Mode.TERSE)), //
					Table.td(bnf.formatHTML(getWeight() * 100.0) + " %"), //
					Table.td(bnf.formatHTML(getEnergy()) + " eV"));
			return t.toHTML(Mode.VERBOSE);
		}
		}
	}

	/**
	 * Hash code derived from element, transition, source and direction.
	 *
	 * @return int
	 * @see com.duckandcover.roentgen.physics.XRay#hashCode()
	 */
	@Override
	public int hashCode() {
		return Objects.hash(mElement, mTransition);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final CharacteristicXRay other = (CharacteristicXRay) obj;
		return Objects.equals(mElement, other.mElement) && //
				Objects.equals(mTransition, other.mTransition);
	}

	@Override
	public int compareTo(final XRay xr) {
		if (xr instanceof CharacteristicXRay) {
			final CharacteristicXRay cxr = (CharacteristicXRay) xr;
			int res = mElement.compareTo(cxr.mElement);
			if (res == 0)
				res = mTransition.compareTo(cxr.mTransition);
			assert (res != 0) || (mElement.equals(cxr.mElement) && mTransition.equals(cxr.mTransition));
			return res;
		} else
			return 1;
	}

}
