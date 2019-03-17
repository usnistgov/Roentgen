package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Objects;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import gov.nist.microanalysis.roentgen.utility.HalfUpFormat;

/**
 * A MatrixCorrectionDatum associated with a Standard (specifies the
 * Composition). Used by MatrixCorrection algorithms as input.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class StandardMatrixCorrectionDatum //
		extends MatrixCorrectionDatum {

	private final Composition mComposition;

	@Override
	public int hashCode() {
		return 37 * super.hashCode() + mComposition.hashCode();
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

	/**
	 * @param comp
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public StandardMatrixCorrectionDatum( //
			final Composition comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle //
	) {
		super(beamEnergy, takeOffAngle);
		mComposition = comp;
	}

	/**
	 *
	 * @param comp
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public StandardMatrixCorrectionDatum( //
			final Composition comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness //
	) {
		super(beamEnergy, takeOffAngle, roughness);
		mComposition = comp;
	}

	/**
	 *
	 * @param comp
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public StandardMatrixCorrectionDatum( //
			final Composition comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness, final Layer coating //
	) {
		super(beamEnergy, takeOffAngle, roughness, coating);
		mComposition = comp;
	}

	@Override
	public String toHTML(final Mode mode) {
		final BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return "Standard[" + mComposition.toHTML(Mode.TERSE) + "]"; // + " at " +
																		// bnf.formatHTML(mBeamEnergy.doubleValue()) + "
																		// keV";
		else {
			final Table t = new Table();
			t.addRow(Table.td("Composition"), Table.td(mComposition.toHTML(Mode.VERBOSE)));
			t.addRow(Table.td("Is standard?"), Table.td("True"));
			t.addRow(Table.td("Beam Energy"), Table.td(bnf.formatHTML(mBeamEnergy)));
			t.addRow(Table.td("Take-off angle"),
					Table.td(bnf.formatHTML(mTakeOffAngle.multiply(180.0 / Math.PI)) + "&deg;"));
			if (mRoughness.isPresent())
				t.addRow(Table.td("Roughness"), Table.td(bnf.formatHTML(mRoughness.get().doubleValue())));
			if (mCoating.isPresent())
				t.addRow(Table.td("Coating"), Table.td(mCoating));
			else
				t.addRow(Table.td("Coating"), Table.td("Not coated"));
			return t.toHTML(mode);
		}
	}

	@Override
	public String toString() {
		final HalfUpFormat df = new HalfUpFormat("0.0");
		return mComposition.toString() + " at " + df.format(mBeamEnergy.doubleValue()) + " keV";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum#
	 * getMaterial()
	 */
	@Override
	public Material getMaterial() {
		return mComposition.getMaterial();
	}

	public Composition getComposition() {
		assert mComposition.hasRepresentation(Representation.MassFraction);
		return mComposition;
	}
}
