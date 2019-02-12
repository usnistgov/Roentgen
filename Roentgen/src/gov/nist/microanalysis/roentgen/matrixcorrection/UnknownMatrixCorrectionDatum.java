package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
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
	
	private final Material mMaterial;

	@Override
	public int hashCode() {
		return 37 * super.hashCode() + Objects.hash(mMaterial);
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
		return Objects.equals(mMaterial, other.mMaterial);
	}

	/**
	 * @param comp Not null Material or Composition
	 * @param beamEnergy Not null keV
	 * @param takeOffAngle Not null radians
	 * @param roughness g/cm<sup>2</sup>
	 * @param coating Layer
	 */
	public UnknownMatrixCorrectionDatum(//
			final Material mat, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness, //
			Layer coating) {
		super(beamEnergy, takeOffAngle, roughness, coating);
		assert mat != null;
		mMaterial = mat;
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public UnknownMatrixCorrectionDatum(//
			Material estimate, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle//
	) {
		this(estimate, beamEnergy, takeOffAngle, Double.NaN, null);
	}

	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(//
			final Material comp, //
			final UncertainValue beamEnergy, //
			final UncertainValue takeOffAngle, //
			final double roughness) {
		this(comp, beamEnergy, takeOffAngle, roughness, null);
	}

	public Set<Element> getElementSet() {
		return Collections.unmodifiableSet(mMaterial.getElementSet());
	}

	public Material getMaterial() {
		return mMaterial;
	}
	
	
	@Override
	public String toHTML(final Mode mode) {
		final BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE) {
			final String elms = getElementSet().stream().map((final Element elm) -> elm.getAbbrev())
					.collect(Collectors.joining(", "));
			return "Unknown[" + elms + "]"; // " at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";

		} else {
			final Table t = new Table();
			t.addRow(Table.td("Estimate"), Table.td(mMaterial.toHTML(Mode.NORMAL)));
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
