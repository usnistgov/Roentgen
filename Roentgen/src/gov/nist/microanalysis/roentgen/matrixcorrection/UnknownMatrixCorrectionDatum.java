package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A MatrixCorrectionDatum associated with an unknown. A list of candidate
 * elements must be provided. A estimated composition can optionally be
 * assigned. Used by MatrixCorrection algorithms as input.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class UnknownMatrixCorrectionDatum //
		extends MatrixCorrectionDatum {

	@Override
	public int hashCode() {
		return 37 * super.hashCode() + Objects.hash(mElements, mEstimate);
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		final UnknownMatrixCorrectionDatum other = (UnknownMatrixCorrectionDatum) obj;
		return Objects.equals(mElements, other.mElements) && Objects.equals(mEstimate, other.mEstimate);
	}

	private final Set<Element> mElements;
	private Optional<Composition> mEstimate;

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public UnknownMatrixCorrectionDatum(final Set<Element> elms, final UncertainValue beamEnergy,
			final UncertainValue takeOffAngle) {
		super(beamEnergy, takeOffAngle);
		mElements = elms;
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public UnknownMatrixCorrectionDatum(final Composition comp, final UncertainValue beamEnergy,
			final UncertainValue takeOffAngle) {
		this(comp, beamEnergy, takeOffAngle, Double.NaN);
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(final Set<Element> elms, final UncertainValue beamEnergy,
			final UncertainValue takeOffAngle, final double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mElements = elms;
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(//
			final Composition comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mElements = comp.getElementSet();
		mEstimate = Optional.of(comp);
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(//
			final Composition comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness, //
			Layer coating) {
		super(beamEnergy, takeOffAngle, roughness, coating);
		mElements = comp.getElementSet();
		mEstimate = Optional.of(comp);
	}

	public Set<Element> getElementSet() {
		return Collections.unmodifiableSet(mElements);
	}

	private boolean validate(final Composition comp) {
		for (final Element elm : comp.getElementSet())
			if (!mElements.contains(elm))
				return false;
		return true;
	}

	public void setEstimated(final Composition comp) {
		assert validate(comp);
		mEstimate = Optional.of(comp);
		mElements.clear();
		mElements.addAll(comp.getElementSet());
	}

	public void clearEstimate() {
		mEstimate = Optional.empty();
	}

	public Composition getEstimate() {
		assert mEstimate.get().hasRepresentation(Representation.MassFraction);
		return mEstimate.get();
	}

	@Override
	public Composition getComposition() {
		assert mEstimate.isPresent();
		assert mEstimate.get().hasRepresentation(Representation.MassFraction);
		return mEstimate.get();
	}

	public boolean hasEstimate() {
		return mEstimate.isPresent();
	}

	@Override
	public String toHTML(final Mode mode) {
		final BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE) {
			final String elms = getElementSet().stream().map((final Element elm) -> elm.getAbbrev())
					.collect(Collectors.joining(", "));
			return "Unknown[" + elms + "] at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";

		} else {
			final Table t = new Table();
			t.addRow(Table.td("Elements"), Table.td(mElements.toString()));
			if (mEstimate.isPresent())
				t.addRow(Table.td("Estimate"), Table.td(mEstimate.get().toHTML(Mode.NORMAL)));
			t.addRow(Table.td("Is standard?"), Table.td("False"));
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
		final DecimalFormat df = new DecimalFormat("0.0");
		final String elms = getElementSet().stream().map((final Element elm) -> elm.getAbbrev())
				.collect(Collectors.joining(", "));
		return "Unknown[" + elms + "] at " + df.format(mBeamEnergy.doubleValue()) + " keV";
	}

}
