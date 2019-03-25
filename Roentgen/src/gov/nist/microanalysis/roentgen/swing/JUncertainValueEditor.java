package gov.nist.microanalysis.roentgen.swing;

import java.awt.Color;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.text.NumberFormat;
import java.text.ParseException;

import javax.swing.JPanel;
import javax.swing.JTextField;

import com.jgoodies.forms.factories.CC;
import com.jgoodies.forms.layout.FormLayout;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class JUncertainValueEditor //
		extends JPanel {

	private static final long serialVersionUID = 6148891525568000402L;

	private UncertainValue mInitial;
	private UncertainValue mCurrent;

	private final JTextField jValueField;
	private final JTextField jUncertaintyField;

	private final NumberFormat mFormat;

	/**
	 * 
	 */
	public JUncertainValueEditor(UncertainValue init) {
		jValueField = new JTextField();
		jValueField.addFocusListener(new FocusListener() {

			@Override
			public void focusLost(FocusEvent e) {
				readout();
			}

			@Override
			public void focusGained(FocusEvent e) {
				jValueField.selectAll();
			}
		});
		jUncertaintyField = new JTextField();
		jUncertaintyField.addFocusListener(new FocusListener() {

			@Override
			public void focusLost(FocusEvent e) {
				readout();
			}

			@Override
			public void focusGained(FocusEvent e) {
				jUncertaintyField.selectAll();
			}
		});
		mFormat = new HalfUpFormat("0.00E0");
		init();
		setValue(init);
	}

	private void init() {
		FormLayout fl = new FormLayout("fill:def:grow, 5dlu, fill:def:grow", "pref");
		setLayout(fl);
		add(jValueField, CC.xy(0, 0));
		add(jUncertaintyField, CC.xy(2, 0));
		revalidate();
	}

	public void setValue(UncertainValue value) {
		mInitial = new UncertainValue(value);
		mCurrent = new UncertainValue(value);
		update();
	}
	
	public void reset() {
		mCurrent = new UncertainValue(mInitial);
		update();
	}

	private void update() {
		if (mCurrent != null) {
			jValueField.setText(mFormat.format(mCurrent.doubleValue()));
			jUncertaintyField.setText(mFormat.format(mCurrent.uncertainty()));
		} else {
			jValueField.setText("1.0");
			jUncertaintyField.setText("0.0");
		}
	}

	private void readout() {
		double val = Double.NaN, unc = Double.NaN;
		try {
			val = mFormat.parse(jValueField.getText()).doubleValue();
		} catch (ParseException e) {
			jValueField.setText(mFormat.format(mCurrent.doubleValue()));
			jValueField.setBackground(Color.pink);
		}
		try {
			unc = mFormat.parse(jValueField.getText()).doubleValue();
		} catch (ParseException e) {
			jUncertaintyField.setText(mFormat.format(mCurrent.doubleValue()));
			jUncertaintyField.setBackground(Color.pink);
		}
		if (!(Double.isNaN(val) || Double.isNaN(unc))) {
			mCurrent = new UncertainValue(val, unc);
			update();
		}
	}
}
