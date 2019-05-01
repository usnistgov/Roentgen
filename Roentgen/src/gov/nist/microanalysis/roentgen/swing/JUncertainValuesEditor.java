package gov.nist.microanalysis.roentgen.swing;

import java.text.NumberFormat;

import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

import com.jgoodies.forms.factories.CC;
import com.jgoodies.forms.layout.FormLayout;

import gov.nist.juncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class JUncertainValuesEditor<H> //
		extends JPanel {

	private final UncertainValues<H> mInitial;
	// private UncertainValues mCurrent;

	private JTable mValuesTable;
	private ValuesTableModel mValuesModel;

	private JTable mCovariancesTable;
	private CovarTableModel mCovarModel;

	private NumberFormat mValFormat;

	/**
	 * 
	 */
	private static final long serialVersionUID = 5228516310659834676L;

	private class ValuesTableModel extends DefaultTableModel {

		private static final long serialVersionUID = -1103957341860868854L;

		public ValuesTableModel(int dim) {
			super(dim, 2);
		}
	};

	private class CovarTableModel extends DefaultTableModel {

		private static final long serialVersionUID = -1103957341860868854L;

		public CovarTableModel(int dim) {
			super(dim, dim);
		}
	};

	private void init() {
		final int size = mInitial != null ? mInitial.getDimension() : 3;
		final FormLayout fl = new FormLayout("pref, 4dlu, pref", "pref");
		this.setLayout(fl);
		mValuesModel = new ValuesTableModel(size);
		mCovarModel = new CovarTableModel(size);
		if (mInitial != null) {
			for (int r = 0; r < mInitial.getDimension(); ++r) {
				final H label = mInitial.getLabel(r);
				mValuesModel.setValueAt(label, r, 0);
				mValuesModel.setValueAt(mValFormat.format(mInitial.getUncertainValue(label)), r, 0);
				for (int c = 0; c < mInitial.getDimension(); ++c)
					mCovarModel.setValueAt(mValFormat.format(mInitial.getCovariance(r, c)), r, c);
			}
		}
		mValuesTable = new JTable(mValuesModel);
		mCovariancesTable = new JTable(mCovarModel);
		this.add(mValuesTable, CC.xy(0, 0));
		this.add(mCovariancesTable, CC.xy(2, 0));
		this.revalidate();
	}

	public void setValueFormat(NumberFormat bnf) {
		mValFormat = bnf;
		repaint();
	}

	public JUncertainValuesEditor(UncertainValues<H> initVals) {
		super(false);
		mInitial = initVals;
		mValFormat = new HalfUpFormat("0.00E0");
		init();
	}
	
	public JUncertainValuesEditor() {
		this(null);
	}


}
