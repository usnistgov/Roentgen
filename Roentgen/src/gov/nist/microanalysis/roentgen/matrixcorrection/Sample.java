package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Objects;
import java.util.Optional;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class Sample implements IToHTML {

	private static int mNextIndex = 0;

	public enum Conductivity {
		Insulator, Semiconductor, Conductor
	};

	private final int mIndex;
	private final String mName;
	private final Optional<Layer> mCoating;
	private final Optional<Composition> mComposition;
	private final Optional<Conductivity> mConductivity;

	@Override
	public int hashCode() {
		return Objects.hash(mCoating, mComposition, mConductivity, mName, mIndex);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final Sample other = (Sample) obj;
		return (mIndex == other.mIndex) && //
				Objects.equals(mCoating, other.mCoating) && //
				Objects.equals(mComposition, other.mComposition) && //
				Objects.equals(mName, other.mName) && //
				Objects.equals(mConductivity, other.mConductivity);
	}

	public Sample(final String name, final Layer coating, final Composition comp, final Conductivity conduct) {
		mIndex = (++mNextIndex);
		mName = name;
		mCoating = Optional.ofNullable(coating);
		mComposition = Optional.ofNullable(comp);
		mConductivity = Optional.ofNullable(conduct);
	}

	public boolean isSuitableAsStandard(final Element elm) {
		return mComposition.isPresent() && mComposition.get().contains(elm);
	}

	public Composition getComposition() {
		return mComposition.get();
	}

	public boolean isCoated() {
		return mCoating.isPresent();
	}

	public Layer getCoating() {
		return mCoating.get();
	}

	public Optional<Conductivity> getConductivity() {
		return mConductivity;
	}

	@Override
	public String toHTML(final Mode mode) {
		final Table t = new Table();
		t.addRow(Table.th("Index"), Table.td(Integer.toString(mIndex)));
		t.addRow(Table.th("Name"), Table.td(HTML.escape(mName)));
		if (mConductivity.isPresent())
			t.addRow(Table.th("Conductivity"), Table.td(mConductivity.toString()));
		else
			t.addRow(Table.th("Conductivity"), Table.td("Unknown"));
		if (isCoated())
			t.addRow(Table.th("Coated"), Table.td(HTML.toHTML(mCoating, Mode.TERSE)));
		else
			t.addRow(Table.th("Coated"), Table.td("-- None --"));
		if (mComposition.isPresent())
			t.addRow(Table.th("Composition"), Table.td(HTML.toHTML(mComposition, Mode.NORMAL)));
		else
			t.addRow(Table.th("Composition"), Table.td("Unknown"));
		return t.toHTML(mode);
	}

	public String getName() {
		return mName;
	}
}
