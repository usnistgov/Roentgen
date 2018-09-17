package gov.nist.jamda;

import java.io.File;
import java.io.IOException;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTMLExt;

/**
 * Implements a 2D window onto a HyperSheet
 * 
 * @author nicho
 *
 */
public class View implements IToHTML, IToHTMLExt {
	
	private final DataFrame mSheet;
	private int mXAxis;
	private int mYAxis;
	
	
	public View(DataFrame sheet, int xAxis, int yAxis){
		mSheet = sheet;
		mXAxis = xAxis;
		mYAxis = yAxis;
	}
	
	final int getXDimension() {
		return mXAxis;
	}
	
	final int getYDimension() {
		return mYAxis;
	}

	@Override
	public String toHTML(Mode mode, File base, String dir) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String toHTML(Mode mode) {
		// TODO Auto-generated method stub
		return null;
	}
	
	
	
	
	

}
