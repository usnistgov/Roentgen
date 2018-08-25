package gov.nist.microanalysis.roentgen.matrixcorrection;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
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
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 312 $
 */

public class KRatioTag extends BaseTag<MatrixCorrectionDatum, MatrixCorrectionDatum, ElementXRaySet> implements IToHTML, Comparable<KRatioTag> {

	public KRatioTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final ElementXRaySet trans) {
		super("K-ratio", unk, std, trans);
		assert trans.size() >= 1;
	}

	public KRatioTag(final MatrixCorrectionDatum unk, final MatrixCorrectionDatum std, final CharacteristicXRay trans) {
		super("K-ratio", unk, std, new ElementXRaySet(trans));
	}

	@Override
	public String toHTML(final Mode mode) {
		if (mode == Mode.TERSE)
			return "k<sub>"+getObject3().toHTML(Mode.TERSE)+"</sub>";
		else if (mode == Mode.NORMAL)
			return "k<sub>"+getObject3().toHTML(Mode.TERSE) + " using " + getObject1().toHTML(Mode.TERSE)+"</sub>";
		else {
			final Table table = new Table();
			table.addRow(Table.th("Item"), Table.th("Description"));
			table.addRow(Table.td("Unknown"), Table.td(getObject1().toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Standard"), Table.td(getObject2().toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Transitions"), Table.td(getObject3().toHTML(Mode.VERBOSE)));
			return table.toHTML(Mode.VERBOSE);
		}
	}

	@Override
	public int compareTo(final KRatioTag o) {
		int c = getObject3().getElement().compareTo(o.getObject3().getElement());
		if (c == 0) {
			final Principle tp = getObject3().getBrightest().getFamily(), op = o.getObject3().getBrightest().getFamily();
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

}
