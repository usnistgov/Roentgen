/**
 * 
 */
package gov.nist.microanalysis.roentgen.microanalysis;

import java.util.Objects;
import java.util.Optional;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;

/**
 * @author Nicholas
 *
 */
public class Sample implements IToHTML {	
		
	private static int mNextIndex = 0;

	private final int mIndex;
	private final String mName;
	private final Optional<Layer> mCoating;
	private final Optional<Composition> mComposition;
	
	@Override
	public int hashCode() {
		return Objects.hash(mCoating, mComposition, mName, mIndex);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Sample other = (Sample) obj;
		return (mIndex == other.mIndex) && //
				Objects.equals(mCoating, other.mCoating) && //
				Objects.equals(mComposition, other.mComposition) && //
				Objects.equals(mName, other.mName);
	}

	public Sample(String name, Layer coating, Composition comp) {
		mIndex = (++mNextIndex);
		mName = name;
		mCoating = Optional.ofNullable(coating);
		mComposition = Optional.ofNullable(comp);
	}

	public boolean isSuitableAsStandard(Element elm) {
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

	@Override
	public String toHTML(Mode mode) {
		Table t = new Table();
		t.addRow(Table.th("Index"), Table.td(Integer.toString(mIndex)));
		t.addRow(Table.th("Name"), Table.td(HTML.escape(mName)));
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
