package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Objects;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A package describing a set of arguments to the matrix correction algorithm.
 * 
 * @author nicho
 *
 */
public class MatrixCorrectionDatum implements IToHTML {
	
	private final Composition mComposition;
	private final boolean mIsStandard;
	private final UncertainValue mBeamEnergy;
	private final UncertainValue mTakeOffAngle;
	
	@Override
	public int hashCode() {
		return Objects.hash(mBeamEnergy, mComposition, mTakeOffAngle, Boolean.valueOf(mIsStandard));
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MatrixCorrectionDatum other = (MatrixCorrectionDatum) obj;
		return Objects.equals(mBeamEnergy, other.mBeamEnergy) && //
				Objects.equals(mComposition, other.mComposition) && //
				Objects.equals(mTakeOffAngle, other.mTakeOffAngle) &&
				(mIsStandard==other.mIsStandard);
	}

	
	public MatrixCorrectionDatum(Composition comp, boolean isStandard, UncertainValue beamEnergy, UncertainValue takeOffAngle){
		mComposition = comp;
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mIsStandard=isStandard;
	}
	
	public Composition getComposition() {
		return mComposition;
	}

	public UncertainValue getBeamEnergy() {
		return mBeamEnergy;
	}

	public UncertainValue getTakeOffAngle() {
		return mTakeOffAngle;
	}
	
	public boolean isStandard() {
		return mIsStandard;
	}

	@Override
	public String toHTML(Mode mode) {
		BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if(mode!=Mode.VERBOSE)
			return mComposition.toHTML(Mode.TERSE)+" at "+bnf.formatHTML(mBeamEnergy.doubleValue())+" keV";
		else {
			Table t = new Table();
			t.addRow(Table.td("Composition"), Table.td(mComposition.toHTML(Mode.VERBOSE)));
			t.addRow(Table.td("Is standard?"), Table.td(Boolean.toString(mIsStandard)));
			t.addRow(Table.td("Beam Energy"), Table.td(bnf.formatHTML(mBeamEnergy)));
			t.addRow(Table.td("Composition"), Table.td(bnf.formatHTML(mTakeOffAngle)));
			return t.toHTML(mode);
		}
	}
}
