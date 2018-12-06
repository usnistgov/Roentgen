package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Optional;
import java.util.Set;

import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * A {@link MatrixCorrectionDatum} representing an Unknown sample containing the specified elements.
 * 
 * @author nicholas
 *
 */
public class UnknownMatrixCorrectionDatum extends MatrixCorrectionDatum {

	
	private final Set<Element> mElements;
	private Optional<Composition> mEstimate;
	
	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public UnknownMatrixCorrectionDatum(Set<Element> elms, UncertainValue beamEnergy, UncertainValue takeOffAngle) {
		super(beamEnergy, takeOffAngle);
		mElements=elms;
	}
	
	
	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 */
	public UnknownMatrixCorrectionDatum(Composition comp, UncertainValue beamEnergy, UncertainValue takeOffAngle) {
		this(comp, beamEnergy, takeOffAngle, Double.NaN);
	}


	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(Set<Element> elms, UncertainValue beamEnergy, UncertainValue takeOffAngle, double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mElements=elms;
	}
	
	/**
	 * @param beamEnergy
	 * @param takeOffAngle
	 * @param roughness
	 */
	public UnknownMatrixCorrectionDatum(Composition comp, UncertainValue beamEnergy, UncertainValue takeOffAngle, double roughness) {
		super(beamEnergy, takeOffAngle, roughness);
		mElements=comp.getElementSet();
		mEstimate = Optional.of(comp);
	}

	
	
	Set<Element> getElementSet(){
		return Collections.unmodifiableSet(mElements);
	}
	
	public void setEstimated(Composition comp) {
		mEstimate=Optional.of(comp);
	}
	
	public void clearEstimate() {
		mEstimate=Optional.empty();
	}

	public Composition getEstimate() {
		return mEstimate.get();
	}
	
	public Composition getComposition() {
		assert mEstimate.isPresent();
		return mEstimate.get();
	}
	
	
	public boolean hasEstimate() {
		return mEstimate.isPresent();
	}
	
	@Override
	public String toHTML(Mode mode) {
		BasicNumberFormat bnf = new BasicNumberFormat("0.0");
		if (mode != Mode.VERBOSE)
			return "Unknown at " + bnf.formatHTML(mBeamEnergy.doubleValue()) + " keV";
		else {
			Table t = new Table();
			t.addRow(Table.td("Elements"), Table.td(mElements.toString()));
			if(mEstimate.isPresent())
				t.addRow(Table.td("Estimate"),Table.td(mEstimate.get().toHTML(Mode.NORMAL)));
			t.addRow(Table.td("Is standard?"), Table.td("False"));
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
		return "Unknown at " + df.format(mBeamEnergy.doubleValue()) + " keV";
	}

	
	
}
