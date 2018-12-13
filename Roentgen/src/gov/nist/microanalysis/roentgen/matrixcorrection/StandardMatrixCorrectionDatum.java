package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.Objects;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A MatrixCorrectionDatum associated with a Standard (specifies the
 * Composition). Used by MatrixCorrection algorithms as input.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class StandardMatrixCorrectionDatum //
		extends MatrixCorrectionDatum {

	@Override
	public int hashCode() {
		return super.hashCode() ^ mComposition.hashCode();
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		final StandardMatrixCorrectionDatum other = (StandardMatrixCorrectionDatum) obj;
		return Objects.equals(mComposition, other.mComposition);
	}

	private final Composition mComposition;

	/**
	 * @param comp
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public StandardMatrixCorrectionDatum(final Composition comp, final UncertainValue beamEnergy,
			final UncertainValue takeOffAngle) {
		super(beamEnergy, takeOffAngle);
		mComposition = comp.asMassFraction();
	}

	/**
	 * /**
	 * 
	 * @param comp
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public StandardMatrixCorrectionDatum(final Composition comp, final UncertainValue beamEnergy,
			final UncertainValue takeOffAngle, final double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mComposition = comp.asMassFraction();
	}

	@Override
	public Composition getComposition() {
		assert mComposition.getNativeRepresentation() == Representation.MassFraction;
		return mComposition;
	}

	@Override
	public String toHTML(final Mode mode) {
		final BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return mComposition.toHTML(Mode.TERSE) + " at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";
		else {
			final Table t = new Table();
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

	@Override
	public String toString() {
		final DecimalFormat df = new DecimalFormat("0.0");
		return mComposition.toString() + " at " + df.format(mBeamEnergy.doubleValue()) + " keV";
	}
}
