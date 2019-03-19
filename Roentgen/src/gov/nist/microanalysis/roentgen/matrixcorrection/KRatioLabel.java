package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealVector;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.Shell.Principle;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;

/**
 * <p>
 * The KRatioTag class is intended for use within the {@link IUncertainValues}
 * framework to identify an input/output value as being a k-ratio. The KRatioTag
 * object contains the information about a single measurement intended to
 * measure the composition of a material. It contains information about the
 * unknown sample, the conditions under which the standard and unknown were
 * measured and the composition of the standard. Any information necessary to
 * convert an x-ray intensity measurement into a composition should be contained
 * within the KRatioSet.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 312 $
 */

public class KRatioLabel//
		extends EPMALabel.BaseLabel<UnknownMatrixCorrectionDatum, StandardMatrixCorrectionDatum, ElementXRaySet> //
		implements IToHTML, Comparable<KRatioLabel> {

	public enum Method {
		Measured, Calculated
	};

	private final Method mMethod;

	public KRatioLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final ElementXRaySet trans, final Method meth) {
		super("k", unk, std, trans);
		assert trans.size() >= 1;
		mMethod = meth;
	}

	public KRatioLabel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
			final CharacteristicXRay trans, final Method meth) {
		super("k", unk, std, new ElementXRaySet(trans));
		mMethod = meth;
	}

	public UnknownMatrixCorrectionDatum getUnknown() {
		return getObject1();
	}

	public StandardMatrixCorrectionDatum getStandard() {
		return getObject2();
	}

	public ElementXRaySet getXRaySet() {
		return getObject3();
	}
	
	public Element getElement() {
		return getObject3().getElement();
	}

	public Method getMethod() {
		return mMethod;
	}

	public KRatioLabel asMeasured() {
		if (isMeasured())
			return this;
		else
			return new KRatioLabel(getUnknown(), getStandard(), getXRaySet(), Method.Measured);
	}

	public KRatioLabel asCalculated() {
		if (isCalculated())
			return this;
		else
			return new KRatioLabel(getUnknown(), getStandard(), getXRaySet(), Method.Calculated);
	}

	public boolean isMeasured() {
		return mMethod == Method.Measured;
	}

	public boolean isCalculated() {
		return mMethod == Method.Calculated;
	}

	@Override
	public int hashCode() {
		return super.hashCode() + 31 * mMethod.hashCode();
	}

	static public boolean areAllSameUnknownElements(final Set<KRatioLabel> krs) {
		Set<Element> elms = null;
		for (final KRatioLabel krl : krs) {
			final UnknownMatrixCorrectionDatum tmp = krl.getUnknown();
			if (elms == null)
				elms = tmp.getElementSet();
			else if (!elms.equals(tmp.getElementSet()))
				return false;
		}
		return true;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		final KRatioLabel other = (KRatioLabel) obj;
		return mMethod == other.mMethod;
	}

	@Override
	public String toHTML(final Mode mode) {
		final String meth = mMethod == Method.Calculated ? "Calc" : "Meas";
		final String xrs = getXRaySet().toHTML(Mode.TERSE);
		if (mode == Mode.TERSE)
			return "k<sub>" + xrs + "," + meth + "</sub>";
		else if (mode == Mode.NORMAL)
			return "k<sub>" + getUnknown().getMaterial().getHTMLName() + "," //
					+ getStandard().getMaterial().getHTMLName() + "," + xrs + "," + meth + "</sub>";
		else {
			final Table table = new Table();

			table.addRow(Table.th("K-ratio", 2));
			table.addRow(Table.td("Unknown"), Table.td(getObject1().toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Standard"), Table.td(getObject2().toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Transitions"), Table.td(getObject3().toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Method"), Table.td(mMethod.toString()));
			return table.toHTML(Mode.VERBOSE);
		}
	}

	@Override
	public String toString() {
		final String meth = mMethod == Method.Calculated ? "Calc" : "Meas";
		final String xrts = HTML.stripTags(getXRaySet().toHTML(Mode.TERSE));
		return "k[" + xrts.substring(1, xrts.length() - 1) + "," + meth + "]";
	}

	@Override
	public int compareTo(final KRatioLabel o) {
		int c = getObject3().getElement().compareTo(o.getObject3().getElement());
		if (c == 0) {
			final Principle tp = getObject3().getBrightest().getFamily(),
					op = o.getObject3().getBrightest().getFamily();
			c = tp.compareTo(op);
		}
		if (c == 0) {
			final CharacteristicXRay tb = getObject3().getBrightest(), ob = o.getObject3().getBrightest();
			c = tb.compareTo(ob);
		}
		if (!equals(o)) {
			if (c == 0)
				c = getObject3().compareTo(o.getObject3());
			if (c == 0)
				c = getObject1().toHTML(Mode.VERBOSE).compareTo(o.getObject1().toHTML(Mode.VERBOSE).toString());
			if (c == 0)
				c = getObject2().toHTML(Mode.VERBOSE).compareTo(o.getObject2().toHTML(Mode.VERBOSE).toString());
			if (c == 0)
				c = hashCode() < o.hashCode() ? -1 : 1;
		}
		return c;
	}

	public static UncertainValues<KRatioLabel> extractKRatios(final RealVector res, final List<? extends Object> labels,
			final Method meth) {
		final Map<KRatioLabel, Number> vals = new HashMap<>();
		for (int i = 0; i < labels.size(); ++i) {
			final Object label = labels.get(i);
			if (label instanceof KRatioLabel) {
				final KRatioLabel calc = (KRatioLabel) label;
				assert calc.isCalculated();
				final KRatioLabel meas = new KRatioLabel(calc.getUnknown(), calc.getStandard(), calc.getXRaySet(),
						meth);
				vals.put(meas, res.getEntry(i));
			}
		}
		return new UncertainValues<KRatioLabel>(vals);
	}

}
