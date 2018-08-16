package gov.nist.microanalysis.roentgen.microanalysis;

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

public class KRatioTag extends BaseTag implements IToHTML, Comparable<KRatioTag> {

	private final MeasurementDatum mUnknown;
	private final MeasurementDatum mStandard;
	private final ElementXRaySet mXRays;

	public KRatioTag(final MeasurementDatum std, final MeasurementDatum unk, final ElementXRaySet trans) {
		super("K-ratio", std, unk, trans);
		assert unk.isSuitableAsUnknown();
		assert std.isSuitableAsStandard(trans.getElement());
		assert trans.size() >= 1;
		mStandard = std;
		mUnknown = unk;
		mXRays = trans;
	}

	public KRatioTag(final MeasurementDatum std, final MeasurementDatum unk, final CharacteristicXRay trans) {
		super("K-ratio", std, unk, trans);
		assert unk.isSuitableAsUnknown();
		assert std.isSuitableAsStandard(trans.getElement());
		mStandard = std;
		mUnknown = unk;
		mXRays = new ElementXRaySet(trans);
	}

	@Override
	public String toHTML(final Mode mode) {
		if (mode == Mode.TERSE)
			return mXRays.toHTML(Mode.TERSE);
		else if (mode == Mode.NORMAL)
			return mXRays.toHTML(Mode.TERSE) + " using " + mStandard.toHTML(Mode.TERSE);
		else {
			final Table table = new Table();
			table.addRow(Table.th("Item"), Table.th("Description"));
			table.addRow(Table.td("Standard"), Table.td(mStandard.toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Unknown"), Table.td(mUnknown.toHTML(Mode.VERBOSE)));
			table.addRow(Table.td("Transitions"), Table.td(mXRays.toHTML(Mode.VERBOSE)));
			return table.toHTML(Mode.VERBOSE);
		}
	}

	@Override
	public int compareTo(final KRatioTag o) {
		int c = mXRays.getElement().compareTo(o.mXRays.getElement());
		if (c == 0) {
			final Principle tp = mXRays.getBrightest().getFamily(), op = o.mXRays.getBrightest().getFamily();
			c = tp.compareTo(op);
		}
		if (c == 0) {
			final CharacteristicXRay tb = mXRays.getBrightest(), ob = o.mXRays.getBrightest();
			c = tb.compareTo(ob);
		}
		if (!equals(o)) {
			if (c == 0)
				c = mXRays.compareTo(o.mXRays);
			if (c == 0)
				c = mStandard.toHTML(Mode.VERBOSE).compareTo(o.mStandard.toHTML(Mode.VERBOSE).toString());
			if (c == 0)
				c = mUnknown.toHTML(Mode.VERBOSE).compareTo(o.mUnknown.toHTML(Mode.VERBOSE).toString());
			if (c == 0)
				c = hashCode() < o.hashCode() ? -1 : 1;
		}
		return c;
	}

}
