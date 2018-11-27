package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A package describing a set of arguments to the matrix correction algorithm.
 * 
 * {@link MatrixCorrectionDatum} come in three different flavors:
 * <ol>
 * <li>Standards - with a defined Composition and isStandard()==true</li>
 * <li>Unknown estimates - with a defined Composition and
 * isStandard()==false</li>
 * <li>True unknowns - with only an Element set</li>
 * </ol>
 * 
 * @author nicholas
 *
 */
public class MatrixCorrectionDatum //
		implements IToHTML {

	// Either mComposition or mElements is defined (never both!)
	private final Optional<Composition> mComposition;
	private final Optional<Set<Element>> mElements;
	private final boolean mIsStandard;
	private final UncertainValue mBeamEnergy;
	private final UncertainValue mTakeOffAngle;
	private final Optional<Double> mRoughness;

	@Override
	public int hashCode() {
		return Objects.hash(mBeamEnergy, mComposition, mElements, mIsStandard, mTakeOffAngle, mRoughness);
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
				Objects.equals(mElements, other.mElements) && //
				(mIsStandard == other.mIsStandard) && //
				Objects.equals(mTakeOffAngle, other.mTakeOffAngle) && //
				Objects.equals(mRoughness, other.mRoughness);
	}

	/**
	 * @param comp
	 * @param isStandard
	 * @param beamEnergy   keV
	 * @param takeOffAngle degrees
	 * @param roughness    in mass thickness cm * g/cm^3 or g/cm^2
	 */
	public MatrixCorrectionDatum(//
			Composition comp, //
			boolean isStandard, //
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle, //
			double roughness //
	) {
		mComposition = Optional.of(comp.asMassFraction());
		mElements = Optional.empty();
		mIsStandard = isStandard;
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mRoughness = Double.isNaN(roughness) ? Optional.empty() : Optional.of(Double.valueOf(roughness));
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

	public MatrixCorrectionDatum( //
			Composition comp, //
			boolean isStandard, //
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle //
	) {
		this(comp, isStandard, beamEnergy, takeOffAngle, Double.NaN);
	}

	public MatrixCorrectionDatum(//
			Set<Element> elms, //
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle, //
			double roughness //
	) {
		mComposition = Optional.empty();
		mElements = Optional.of(new HashSet<>(elms));
		mIsStandard = false;
		mBeamEnergy = beamEnergy;
		mTakeOffAngle = takeOffAngle;
		mRoughness = Double.isNaN(roughness) ? Optional.empty() : Optional.of(Double.valueOf(roughness));
	}

	public MatrixCorrectionDatum(//
			Set<Element> elms, //
			UncertainValue beamEnergy, //
			UncertainValue takeOffAngle //
	) {
		this(elms, beamEnergy, takeOffAngle, Double.NaN);
	}

	public Composition getComposition() {
		return mComposition.get();
	}

	public Set<Element> getElementSet() {
		return mComposition.isPresent() ? mComposition.get().getElementSet() : mElements.get();
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

	public boolean isUnknown() {
		return mElements.isPresent();
	}

	public double getRoughness() {
		return mRoughness.orElse(0.0);
	}

	public boolean hasRoughness() {
		return mRoughness.isPresent();
	}

	public String getNameAsHTML(Mode mode) {
		return mComposition.isPresent() ? mComposition.get().toHTML(mode) : "Unknown";
	}

	@Override
	public String toHTML(Mode mode) {
		BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return getNameAsHTML(Mode.TERSE) + " at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";
		else {
			Table t = new Table();
			t.addRow(Table.td("Composition"), Table.td(getNameAsHTML(Mode.VERBOSE)));
			t.addRow(Table.td("Is standard?"), Table.td(Boolean.toString(isStandard())));
			t.addRow(Table.td("Beam Energy"), Table.td(bnf.formatHTML(mBeamEnergy)));
			t.addRow(Table.td("Take-off angle"), Table.td(bnf.formatHTML(mTakeOffAngle.multiply(180.0/Math.PI))+"&deg;"));
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
