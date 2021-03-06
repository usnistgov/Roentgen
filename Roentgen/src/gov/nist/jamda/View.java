package gov.nist.jamda;

import java.io.File;
import java.io.IOException;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTMLExt;

import gov.nist.microanalysis.roentgen.ArgumentException;

/**
 * Implements a 2D window onto a HyperSheet
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class View implements IToHTML, IToHTMLExt {

	private final DataFrame mSheet;
	private final Extent mExtent;

	public View(final DataFrame sheet, final Extent extent) throws ArgumentException {
		if (extent.getDimensionality() > 2)
			throw new ArgumentException("The dimensionality of a view's extent must be 1 or 2.");
		mSheet = sheet;
		mExtent = extent;
	}

	final Extent getExtent() {
		return mExtent;
	}

	@Override
	public String toHTML(final Mode mode, final File base, final String dir) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String toHTML(final Mode mode) {

		if (mExtent.getDimensionality() == 1) {

		} else {

		}

		// TODO Auto-generated method stub
		return null;
	}

}
