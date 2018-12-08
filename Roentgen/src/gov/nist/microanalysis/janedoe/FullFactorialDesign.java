package gov.nist.microanalysis.janedoe;

import java.util.ArrayList;
import java.util.List;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.math.Utility;

/**
 * Constructs the design matrix for a full factorial experimental design in
 * Yates ordering.
 *
 * @author Nicholas
 */
public class FullFactorialDesign implements IToHTML {

	private final int[] mLevels;

	public FullFactorialDesign(final int[] levels) {
		mLevels = levels.clone();
	}

	public int numberOfMeasurements() {
		return Utility.sum(mLevels);
	}

	/**
	 * Calculates the design matrix in Yates order
	 *
	 * @return
	 */
	public int[][] designMatrix() {
		final int nMeas = numberOfMeasurements();
		final int nFactors = mLevels.length;
		final int[][] res = new int[nMeas][nFactors];
		int runLen = 1;
		for (int factor = 0; factor < nFactors; ++factor) {
			for (int j = 0; j < nMeas; ++j)
				res[j][factor] = (j / (runLen * mLevels[factor])) % mLevels[factor];
			runLen *= mLevels[factor];
		}
		return res;
	}

	@Override
	public String toHTML(final Mode mode) {
		switch (mode) {
		case TERSE:
			return "Full factorial design with " + mLevels.length + " factors";
		case NORMAL: {
			final StringBuffer res = new StringBuffer();
			final int max = Utility.max(mLevels);
			final int[] c = new int[max + 1];
			for (int i = 0; i < mLevels.length; ++i)
				++c[mLevels[i]];
			for (int i = 0; i < max; ++i)
				if (c[i] > 0) {
					res.append(Integer.toString(i));
					res.append("<sup>");
					res.append(Integer.toString(c[i]));
					res.append("</sup>");
				}
			res.append(" Full Factorial Design");
			return res.toString();
		}

		default: {
			final int[][] dm = designMatrix();
			final int nMeas = numberOfMeasurements();
			final int nFactors = mLevels.length;
			final Table t = new Table();
			{
				final List<Item> row = new ArrayList<>();
				row.add(Table.th("n"));
				for (int factor = 0; factor < nFactors; ++factor)
					row.add(Table.th("X<sub>" + factor + "</sub>"));
				t.addRow(row);
			}
			for (int j = 0; j < nMeas; ++j) {
				final List<Item> row = new ArrayList<>();
				row.add(Table.th(Integer.toString(j)));
				for (int factor = 0; factor < nFactors; ++factor)
					row.add(Table.td(Integer.toString(dm[j][factor])));
				t.addRow(row);
			}
			return t.toHTML(mode);
		}
		}
	}
}
