package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.Objects;
import java.util.Optional;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A package describing a set of arguments to the matrix correction algorithm.
 * 
 * @author nicholas
 *
 */
public class MatrixCorrectionDatum implements IToHTML {

	private final Composition mComposition;
	private final boolean mIsStandard;
	private final UncertainValue mBeamEnergy;
	private final UncertainValue mTakeOffAngle;
	private final Optional<Double> mRoughness;

	@Override
	public int hashCode() {
		return Objects.hash(mBeamEnergy, mComposition, mTakeOffAngle, Boolean.valueOf(mIsStandard), mRoughness);
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
				Objects.equals(mTakeOffAngle, other.mTakeOffAngle) && //
				(mIsStandard == other.mIsStandard) && //
				Objects.equals(mRoughness, other.mRoughness);
	}

	/**
	 * @param comp
	 * @param isStandard
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 */
	public MatrixCorrectionDatum(Composition comp, boolean isStandard, UncertainValue beamEnergy,
			UncertainValue takeOffAngle, double roughness) {
		mComposition = comp;
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mIsStandard = isStandard;
		mRoughness = Optional.of(Double.valueOf(roughness));
	}

	/**
	 * Helper to simplify computing roughness.
	 * 
	 * @param dimension In nanometers
	 * @param density   g/cm<sup>3</sup>
	 * @return double cm (g/cm<sup>3</sup>)
	 */
	static public double roughness(double dimension, double density) {
		return 1.0e-7 * dimension * density;
	}

	public MatrixCorrectionDatum(Composition comp, boolean isStandard, UncertainValue beamEnergy,
			UncertainValue takeOffAngle) {
		mComposition = comp;
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mIsStandard = isStandard;
		mRoughness = Optional.empty();
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

	public double getRoughness() {
		return mRoughness.orElse(0.0);
	}
	
	public boolean hasRoughness() {
		return mRoughness.isPresent();
	}

	@Override
	public String toHTML(Mode mode) {
		BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return mComposition.toHTML(Mode.TERSE) + " at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";
		else {
			Table t = new Table();
			t.addRow(Table.td("Composition"), Table.td(mComposition.toHTML(Mode.VERBOSE)));
			t.addRow(Table.td("Is standard?"), Table.td(Boolean.toString(mIsStandard)));
			t.addRow(Table.td("Beam Energy"), Table.td(bnf.formatHTML(mBeamEnergy)));
			t.addRow(Table.td("Composition"), Table.td(bnf.formatHTML(mTakeOffAngle)));
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
