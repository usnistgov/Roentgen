package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.math3.util.ArithmeticUtils;

import com.duckandcover.html.IToHTML;
import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.physics.Element;

public class StoichiometeryRules implements IToHTML {

	final private Element mElement;
	final private Map<Element, Composition> mRules = new TreeMap<>();

	private StoichiometeryRules(final Element elm) {
		mElement = elm;
	}

	public void add(final Element elm, final Composition comp) {
		final Set<Element> elms = comp.getElementSet();
		Preconditions.checkArgument(elms.contains(elm));
		Preconditions.checkArgument(elms.contains(mElement));
		Preconditions.checkArgument(elms.size() == 2);
		mRules.put(elm, comp);
	}

	public Element getElement() {
		return mElement;
	}
	
	public Composition getRule(final Element elm) {
		return mRules.getOrDefault(elm, Composition.pureElement(elm));
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		default:
		case TERSE:
			return "Default Oxygen Stoichiometry Rules";
		case NORMAL:
		case VERBOSE: {
			final StringBuffer sb = new StringBuffer();
			sb.append("<h3>" + mElement.toString() + "-by-Stoichiometry</h3>");
			sb.append("<table>");
			sb.append("<tr><th>Element</th><th>Material</th><th>Composition</th></tr>");
			for (final Map.Entry<Element, Composition> me : mRules.entrySet()) {
				final Element elm = me.getKey();
				final Composition sc = me.getValue();
				sb.append("<tr>");
				sb.append("<td>" + elm.toString() + "</td>");
				sb.append("<td>" + sc.toHTML(Mode.TERSE) + "</td>");
				sb.append("<td>" + sc.toHTML(mode) + "</td>");
				sb.append("</tr>");
			}
			sb.append("</table>");
			return sb.toString();
		}
		}
	}

	public Map<Element, Composition> getRuleMap() {
		return Collections.unmodifiableMap(mRules);
	}
	
	public static Map<Element, Integer> getOxygenDefaults(){
		final int[] valences = { 0, 1, 0, 1, 2, 3, 4, 5, -2, 1, 0, 1, 2, 3, 4, 5, 6, 5, 0, 1, 2, 3, 4, 5, 2, 2, 2, 2, 2,
				2, 2, 3, 4, 3, 6, 5, 0, 1, 2, 3, 4, 5, 6, 2, 4, 4, 2, 1, 2, 3, 2, 3, 4, 5, 0, 1, 2, 3, 4, 3, 3, 3, 3, 3,
				3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 4, 4, 4, 4, 3, 2, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 4, 4, 4 };
		Map<Element, Integer> res=new HashMap<>();
		for(int i=1;i<valences.length;++i)
			res.put(Element.byAtomicNumber(i), valences[i]);
		return res;
	}
	
	

	public static StoichiometeryRules defaultOxygenByStoichiomery() //
			throws ArgumentException {
		final Element elmO = Element.Oxygen;
		final StoichiometeryRules res = new StoichiometeryRules(elmO);
		 Map<Element, Integer> valences = getOxygenDefaults();
		final int valenceO = valences.get(elmO);
		assert valenceO == -2;
		for (int z = 1; z <= 94; ++z)
			if (z != 8) {
				final Element other = Element.byAtomicNumber(z);
				final int valenceOther = valences.get(other);
				if (valenceOther > 0) {
					final int gcd = ArithmeticUtils.gcd(valenceOther, -valenceO);
					final Map<Element, Integer> stoic = new TreeMap<>();
					stoic.put(other, Integer.valueOf(-valenceO / gcd));
					stoic.put(elmO, Integer.valueOf(valenceOther / gcd));
					final StringBuffer html = new StringBuffer();
					html.append(other.getAbbrev());
					if (stoic.get(other).intValue() > 1)
						html.append("<sub>" + stoic.get(other).toString() + "</sub>");
					html.append(elmO.getAbbrev());
					if (stoic.get(elmO).intValue() > 1)
						html.append("<sub>" + stoic.get(elmO).toString() + "</sub>");
					res.add(other, Composition.stoichiometry(html.toString(), stoic));
				}
			}
		return res;
	}
}