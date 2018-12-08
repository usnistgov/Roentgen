package gov.nist.microanalysis.janedoe;

import java.util.ArrayList;
import java.util.List;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

/**
 * Implements a 2<sup>k-l</sup> design matrix.
 *
 * @author Nicholas
 */
public class FractionalFactorialDesign implements IToHTML {

	private final int mK;
	private final int[][] mLDesign;

	public FractionalFactorialDesign(final int k, final int[][] lDesign) {
		mK = k;
		mLDesign = lDesign.clone();
	}

	public int fractionalLength() {
		return mLDesign.length;
	}

	public int numberOfMeasurements() {
		return (int) Math.round(Math.pow(2.0, mK - mLDesign.length));
	}

	/**
	 * Calculates the design matrix in Yates order
	 *
	 * @return
	 */
	public int[][] designMatrix() {
		final int nMeas = numberOfMeasurements();
		final int[][] res = new int[nMeas][mK];
		int runLen = 1;
		for (int factor = 0; factor < mK - mLDesign.length; ++factor) {
			for (int j = 0; j < nMeas; ++j)
				res[j][factor] = -1 + 2 * ((j / (runLen * 2)) % 2);
			runLen *= 2;
		}
		for (int kp = 0; kp < mLDesign.length; ++kp) {
			final int factor = kp + mK - mLDesign.length;
			for (int j = 0; j < nMeas; ++j) {
				int val = 1;
				for (final int kk : mLDesign[kp])
					val *= res[j][kk];
				res[j][factor] = val;
			}
		}
		return res;
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE:
		case NORMAL:
			return "2<sup>" + mK + "-" + mLDesign.length + "</sup> Fractional Factorial Design";
		default: {
			final int[][] dm = designMatrix();
			final int nMeas = numberOfMeasurements();
			final Table t = new Table();
			{
				final List<Item> row = new ArrayList<>();
				row.add(Table.th("n"));
				for (int factor = 0; factor < mK; ++factor)
					row.add(Table.th("X<sub>" + factor + "</sub>"));
				t.addRow(row);
			}
			for (int j = 0; j < nMeas; ++j) {
				final List<Item> row = new ArrayList<>();
				row.add(Table.th(Integer.toString(j)));
				for (int factor = 0; factor < mK; ++factor)
					row.add(Table.td(Integer.toString(dm[j][factor])));
				t.addRow(row);
			}
			return t.toHTML(mode);
		}
		}
	}

}
