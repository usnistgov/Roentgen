package gov.nist.juncertainty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.FastIndex;

/**
 * A very light wrapper around the Jacobian to make it safer and easier to
 * access entries using labels.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class Jacobian<G, H> //
		implements IToHTML {

	private final FastIndex<G> mInputs;
	private final FastIndex<H> mOutputs;
	private final RealMatrix mJacobian;

	/**
	 * Constructs a Jacobian object associated with a {@link RealMatrix} where the
	 * rows are labeled by <code>outputs</code> and the columns are labeled by
	 * <code>inputs</code>.
	 * 
	 * @param inputs A {@link List} of input variable labels
	 * @param outputs A {@link List} of output variable labels
	 * @param rm A matrix of dimensions outputs.size() rows x inputs.size() columns
	 */
	public Jacobian(
			final List<? extends G> inputs, //
			final List<? extends H> outputs, //
			final RealMatrix rm
	) {
		assert rm.getColumnDimension() == inputs.size();
		assert rm.getRowDimension() == outputs.size();
		mInputs = new FastIndex<>(inputs);
		mOutputs = new FastIndex<>(outputs);
		mJacobian = rm;
	}

	/**
	 * Returns the Jacobian entry associated with the specified input and output
	 * variable. Note the argument order which is different from the natural
	 * ordering of elements in the Jacobian matrix.
	 *
	 * @param input An input variable label
	 * @param output An output variable label
	 * @return double The partial derivative of output with respect to input
	 */
	public double getEntry(
			final G input, final H output
	) {
		return mJacobian.getEntry(mOutputs.indexOf(output), mInputs.indexOf(input));
	}

	/**
	 * @return A list of the input labels
	 */
	public List<G> getInputLabels() {
		return Collections.unmodifiableList(mInputs);
	}

	/**
	 * @return A list of output labels
	 */
	public List<H> getOutputLabels() {
		return Collections.unmodifiableList(mOutputs);
	}

	public String toHTML(
			final Mode mode, //
			final BasicNumberFormat bnf
	) {
		switch (mode) {
		case TERSE:
			return "Jacobian[N[Outputs]=" + mOutputs.size() + ",N[Inputs]=" + mInputs.size() + "]";
		case NORMAL:
		case VERBOSE:
		default: {
			final Table t = new Table();
			t.addHeaderRow("Inputs", getInputLabels());
			for (final H outLabel : getOutputLabels()) {
				final List<Item> items = new ArrayList<>();
				items.add(Table.th(HTML.toHTML(outLabel, Mode.TERSE)));
				for (final G inLabel : getInputLabels()) {
					final String num = mode == Mode.NORMAL ? bnf.format(getEntry(inLabel, outLabel))
							: bnf.formatHTML(getEntry(inLabel, outLabel));
					items.add(Table.td(num));
				}
				t.addRow(items);
			}
			return t.toHTML(Mode.NORMAL);
		}
		}
	}

	public String toHTML(
			final Mode mode
	) {
		return toHTML(mode, new BasicNumberFormat("0.00E0"));
	}
}
