/**
 * 
 */
package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * @author nicho
 *
 */
public class StandardMatrixCorrectionDatum extends MatrixCorrectionDatum {

	private final Composition mComposition;

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public StandardMatrixCorrectionDatum(Composition comp, UncertainValue beamEnergy, UncertainValue takeOffAngle) {
		super(beamEnergy, takeOffAngle);
		mComposition = comp;
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public StandardMatrixCorrectionDatum(Composition comp, UncertainValue beamEnergy, UncertainValue takeOffAngle, double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mComposition=comp;
	}

	public Composition getComposition() {
		return mComposition;
	}
	
	
	@Override
	public String toHTML(Mode mode) {
		BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return mComposition.toHTML(Mode.TERSE) + " at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";
		else {
			Table t = new Table();
			t.addRow(Table.td("Composition"), Table.td(mComposition.toHTML(Mode.VERBOSE)));
			t.addRow(Table.td("Is standard?"), Table.td("True"));
			t.addRow(Table.td("Beam Energy"), Table.td(bnf.formatHTML(mBeamEnergy)));
			t.addRow(Table.td("Take-off angle"),
					Table.td(bnf.formatHTML(mTakeOffAngle.multiply(180.0 / Math.PI)) + "&deg;"));
			if (mRoughness.isPresent())
				t.addRow(Table.td("Roughness"), Table.td(bnf.formatHTML(mRoughness.get().doubleValue())));
			return t.toHTML(mode);
		}
	}

	public String toString() {
		DecimalFormat df = new DecimalFormat("0.0");
		return mComposition.toString() + " at " + df.format(mBeamEnergy.doubleValue()) + " keV";
	}

	
	
}
