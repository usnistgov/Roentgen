package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.d3x.morpheus.frame.DataFrame;
import com.d3x.morpheus.range.Range;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ParallelMeasurementModelBuilder;
import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.juncertainty.utility.FastIndex;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.EPMALabel.MaterialMAC;
import gov.nist.microanalysis.roentgen.math.Constraint;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.Layer;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction;
import gov.nist.microanalysis.roentgen.physics.XRay;
import gov.nist.microanalysis.roentgen.physics.XRayEmissionProbability;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.AtomFraction;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

/**
 * <p>
 * This class implements the XPP matrix correction algorithm as described in the
 * "Green Book" (Heinrich & Newbury 1991) by Pouchou and Pichoir. It also
 * implements all the Jacobians necessary to implement an uncertainty
 * calculation to propagate uncertainties in the input parameters into output
 * parameters. It breaks the calculation into a series of much simpler steps for
 * which the Jacobian can be readily calculated.
 * <p>
 *
 * <p>
 * Since it is derived from {@link ExplicitMeasurementModel}, it computes not
 * only the value (ie. the k-ratio) but also sensitivity matrix (Jacobian) that
 * maps uncertainty in the input parameters into uncertainty in the output
 * parameters.
 * <p>
 *
 * <p>
 * The following parameters can have associated uncertainties.
 * <ul>
 * <li>The composition of the standards (user specified)</li>
 * <li>The composition of the unknown (user specified)</li>
 * <li>The beam energy (user specified)</li>
 * <li>The take-off angle (user specified)</li>
 * <li>The mass absorption coefficient (computed using FFAST)</li>
 * <li>The mean ionization coefficient (computed using the equation in
 * PAP1991)</li>
 * <li>Surface roughness (user specified)</li>
 * </ul>
 *
 * <p>
 * This class calculates XPP assuming that all the data associate with the
 * unknown was collected at the same beam energy.
 * <p>
 *
 *
 * @author Nicholas
 */
public class XPPMatrixCorrection2 //
		extends MatrixCorrectionModel2 //
		implements IToHTML {

	private static final double MIN_EPS = 1.0e-6; // For steps stepa and stepAB

	private static class Stepa // Checked 16-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			return Collections.singletonList(MatrixCorrectionModel2.shellLabel("a", datum, shell));
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		public Stepa(final MatrixCorrectionDatum datum, //
				final AtomicShell shell) throws ArgumentException {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(iphi0, iF, iP, iRbar, ib);

			final double phi0 = point[iphi0];
			final double P = point[iP];
			final double Rbar = point[iRbar];
			final double F = point[iF];
			final double b = point[ib];

			final RealVector vals = buildResult();

			final int oa = outputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(oa);

			final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // C1
			final double eps = (a - b) / b;
			vals.setEntry(oa, Math.abs(eps) >= MIN_EPS ? a : (1.0 + MIN_EPS) * b);
			return vals;
		}

		@Override
		public String toString() {
			return "Stepa[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(iphi0, iF, iP, iRbar, ib);

			final double phi0 = point.getEntry(iphi0);
			final double P = point.getEntry(iP);
			final double Rbar = point.getEntry(iRbar);
			final double F = point.getEntry(iF);
			final double b = point.getEntry(ib);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oa = outputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(oa);

			final double b2 = Math.pow(b, 2.0);
			final double den = Math.pow(phi0 + b2 * F * Rbar - 2.0 * b * F, 2.0);
			final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // Ok
			final double eps = (a - b) / b;
			if (Math.abs(eps) < MIN_EPS) {
				vals.setEntry(oa, (1.0 + MIN_EPS) * b); // C1
				jac.setEntry(oa, ib, 1.0 + MIN_EPS); // C1
			} else {
				final double dadphi0 = (3.0 * b2 * F + P - 2.0 * b2 * b * F * Rbar) / den; // Ok
				final double dadP = 1.0 / (b * F * (2.0 - b * Rbar) - phi0); // Ok
				final double dadb = -2.0
						* (F * P + phi0 * phi0 - b * F * (phi0 + P * Rbar) + b2 * F * (F - phi0 * Rbar)) / den; // Ok
				final double dadF = b * (P * (-2.0 + b * Rbar) + b * phi0 * (2.0 * b * Rbar - 3.0)) / den; // Ok
				final double dadRbar = (b2 * F * (P + 2.0 * b * phi0 - b2 * F)) / den; // Ok

				vals.setEntry(oa, a); // C1
				jac.setEntry(oa, iP, dadP); // C1
				jac.setEntry(oa, ib, dadb); // C1
				jac.setEntry(oa, iphi0, dadphi0); // C1
				jac.setEntry(oa, iF, dadF); // C1
				jac.setEntry(oa, iRbar, dadRbar); // C1
			}
			return Pair.create(vals, jac);
		}
	}

	private static class StepAaBb //
			extends CompositeMeasurementModel<EPMALabel> //
	{

		private static List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps(
				final MatrixCorrectionDatum datum, final AtomicShell shell) throws ArgumentException {
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
			res.add(new StepRphi0(datum, shell));
			res.add(new StepFRbar(datum, shell));
			res.add(new StepPb(datum, shell));
			res.add(new Stepa(datum, shell));
			res.add(new StepAB(datum, shell));
			return res;
		}

		public StepAaBb(final MatrixCorrectionDatum datum, final AtomicShell shell) throws ArgumentException {
			super("StepAaBb", buildSteps(datum, shell));
		}

		@Override
		public String toString() {
			return "StepAaBb[]";
		}

	}

	private static class StepAB // Checked 17-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("a", datum, shell));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("A", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("B", datum, shell));
			return res;
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		public StepAB(final MatrixCorrectionDatum datum, final AtomicShell shell) throws ArgumentException {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(iphi0, iF, iP, ib, ia);

			final double phi0 = point[iphi0];
			final double F = point[iF];
			final double P = point[iP];
			final double b = point[ib];
			final double a = point[ia];

			final RealVector vals = buildResult();

			final int oA = outputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mShell));
			final int oB = outputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mShell));
			checkIndices(oA, oB);
			// Page 62 Plugging B[b,F,(a-b)/a,phi0,P] into the expression for
			// A[B,b,F,(a-b)/a,phi0,P] we get
			vals.setEntry(oA, a * (b * (b * F - 2.0 * phi0) - P) / Math.pow(a - b, 2.0)); // C2
			vals.setEntry(oB, b * (a * b * F - P - phi0 * (a + b)) / (a - b)); // C2
			return vals;
		}

		@Override
		public String toString() {
			return "StepAB[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(iphi0, iF, iP, ib, ia);

			final double phi0 = point.getEntry(iphi0);
			final double P = point.getEntry(iP);
			final double F = point.getEntry(iF);
			final double b = point.getEntry(ib);
			final double a = point.getEntry(ia);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oA = outputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mShell));
			final int oB = outputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mShell));
			checkIndices(oA, oB);
			// See Heinrich EPQ Page 62
			vals.setEntry(oA, a * (b * b * F - P - 2.0 * b * phi0) / Math.pow(a - b, 2.0)); // C4
			jac.setEntry(oA, iphi0, -2.0 * a * b / Math.pow(a - b, 2.0)); // C4
			jac.setEntry(oA, iP, -a / Math.pow(a - b, 2.0)); // C4
			jac.setEntry(oA, iF, a * b * b / Math.pow(a - b, 2.0)); // C4
			jac.setEntry(oA, ib, (2.0 * a * (a * b * F - P - (a + b) * phi0)) / Math.pow(a - b, 3.0)); // C4
			jac.setEntry(oA, ia, (a + b) * (P + 2 * b * phi0 - b * b * F) / Math.pow(a - b, 3.0)); // C3

			vals.setEntry(oB, b * (a * b * F - P - phi0 * (a + b)) / (a - b)); // C4
			jac.setEntry(oB, iphi0, -b * (a + b) / (a - b)); // C4
			jac.setEntry(oB, iP, -b / (a - b)); // C4
			jac.setEntry(oB, iF, (a * b * b) / (a - b)); // X
			jac.setEntry(oB, ib, (a * a * (2.0 * b * F - phi0) + b * b * phi0 - a * (b * (b * F + 2.0 * phi0) + P))
					/ Math.pow(a - b, 2.0)); // C4
			jac.setEntry(oB, ia, b * (P + 2.0 * b * phi0 - b * b * F) / Math.pow(a - b, 2.0)); // C2

			return Pair.create(vals, jac);
		}
	}

	private static class StepChi // Checked 16-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> { // C1

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(matMacLabel(datum.getMaterial(), cxr));
			res.add(MatrixCorrectionModel2.takeOffAngleLabel(datum));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.chiLabel(datum, cxr));
		}

		private final MatrixCorrectionDatum mDatum;

		private final CharacteristicXRay mXRay;

		public StepChi(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr) throws ArgumentException {
			super(buildInputs(datum, cxr), buildOutputs(datum, cxr));
			mDatum = datum;
			mXRay = cxr;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final RealVector vals = buildResult();

			final int ochi = outputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			checkIndices(ochi);

			final double toa = point[inputIndex(MatrixCorrectionModel2.takeOffAngleLabel(mDatum))];
			assert (toa > 0.0) && (toa < 0.5 * Math.PI) : "toa = " + toa;
			final double csc = 1.0 / Math.sin(toa);
			final double mac = point[inputIndex(MatrixCorrectionModel2.matMacLabel(mDatum.getMaterial(), mXRay))];
			final double chi = mac * csc; // C1
			vals.setEntry(ochi, chi);

			return vals;
		}

		@Override
		public String toString() {
			return "StepChi[" + mDatum + "," + mXRay + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final MatrixCorrectionDatumLabel toaT = MatrixCorrectionModel2.takeOffAngleLabel(mDatum);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int ochi = outputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			checkIndices(ochi);

			final int toai = inputIndex(toaT);
			final double toa = point.getEntry(toai);
			assert (toa > 0.0) && (toa < 0.5 * Math.PI) : "toa = " + toa;
			final double csc = 1.0 / Math.sin(toa);

			final MaterialMAC macT = MatrixCorrectionModel2.matMacLabel(mDatum.getMaterial(), mXRay);
			final int maci = inputIndex(macT);
			final double mac = point.getEntry(maci);
			final double chi = mac * csc;
			jac.setEntry(ochi, maci, csc);
			jac.setEntry(ochi, toai, -1.0 * chi / Math.tan(toa)); // C1
			vals.setEntry(ochi, chi);

			return Pair.create(vals, jac);
		}

	}

	private static class StepConductiveCoating //
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		private static final List<EPMALabel> buildInputLabels(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FofChiLabel(datum, cxr));
			if (datum.hasCoating()) {
				final Layer coating = datum.getCoating();
				res.add(MatrixCorrectionModel2.coatingMassThickness(coating));
				res.add(MatrixCorrectionModel2.matMacLabel(coating.getMaterial(), cxr));
				res.add(MatrixCorrectionModel2.takeOffAngleLabel(datum));
			}
			return res;

		}

		private static final List<EPMALabel> buildOutputLabels(
				//
				final MatrixCorrectionDatum datum, final CharacteristicXRay cxr //
		) {
			return Collections.singletonList(MatrixCorrectionModel2.FofChiReducedLabel(datum, cxr));
		}

		private final MatrixCorrectionDatum mDatum;

		private final CharacteristicXRay mXRay;

		public StepConductiveCoating(final MatrixCorrectionDatum mcd, //
				final CharacteristicXRay cxr) throws ArgumentException {
			super(buildInputLabels(mcd, cxr), buildOutputLabels(mcd, cxr));
			mDatum = mcd;
			mXRay = cxr;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			double coatTrans = 1.0;
			final double Fx = point[inputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay))];

			final RealVector vals = buildResult();

			if (mDatum.hasCoating()) {
				final Layer coating = mDatum.getCoating();
				final int iMassTh = inputIndex(MatrixCorrectionModel2.coatingMassThickness(coating));
				// A simple model for absorption of x-rays by a thin coating...
				final double massTh = point[iMassTh];
				final int itoa = inputIndex(MatrixCorrectionModel2.takeOffAngleLabel(mDatum));
				final int imu = inputIndex(MatrixCorrectionModel2.matMacLabel(coating.getMaterial(), mXRay));
				final double toa = point[itoa];
				final double csc = 1.0 / Math.sin(toa);
				final double mu = point[imu];
				coatTrans = Math.exp(-mu * massTh * csc);
				assert coatTrans > 0.1 : "coatTrans = " + coatTrans;
			}
			final int oFxRed = outputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, mXRay));
			vals.setEntry(oFxRed, coatTrans * Fx);
			return vals;
		}

		@Override
		public String toString() {
			final Layer coating = mDatum.getCoating();
			return "StepCoating[" + (coating != null ? "coating=" + coating.toString() : "Uncoated") + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int iFx = inputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			final double Fx = point.getEntry(iFx);
			final int oFxRed = outputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, mXRay));
			double coatTrans = 1.0;

			if (mDatum.hasCoating()) {
				final Layer coating = mDatum.getCoating();
				final int iMassTh = inputIndex(MatrixCorrectionModel2.coatingMassThickness(coating));
				// A simple model for absorption of x-rays by a thin coating...
				final double massTh = point.getEntry(iMassTh);
				final int itoa = inputIndex(MatrixCorrectionModel2.takeOffAngleLabel(mDatum));
				final int imu = inputIndex(MatrixCorrectionModel2.matMacLabel(coating.getMaterial(), mXRay));
				final double toa = point.getEntry(itoa);
				final double csc = 1.0 / Math.sin(toa);
				final double mu = point.getEntry(imu);
				coatTrans = Math.exp(-mu * massTh * csc);
				final double dcoatTransdmu = coatTrans * (-massTh * csc);
				final double dcoatTransdmassTh = coatTrans * (-mu * csc);
				final double dcoatTransdtoa = coatTrans * massTh * mu * csc * csc * Math.cos(toa);
				jac.setEntry(oFxRed, iMassTh, dcoatTransdmassTh * Fx);
				jac.setEntry(oFxRed, imu, dcoatTransdmu * Fx);
				jac.setEntry(oFxRed, itoa, dcoatTransdtoa * Fx);
				assert coatTrans > 0.1 : "coatTrans = " + coatTrans;
			}
			jac.setEntry(oFxRed, iFx, coatTrans);
			vals.setEntry(oFxRed, coatTrans * Fx);

			return Pair.create(vals, jac);
		}
	}

	private static class StepFRbar // Checked 15-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(
				//
				final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, datum.getMaterial()));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			return res;
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		public StepFRbar(final MatrixCorrectionDatum datum, //
				final AtomicShell shell) throws ArgumentException {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			final int iOneOverS = inputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int iQlaE0 = inputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iR = inputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			checkIndices(iZbarb, iOneOverS, iQlaE0, iPhi0, iR);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = point[inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum))] / Ea;
			final double Zbarb = point[iZbarb];
			final double oneOverS = point[iOneOverS];
			final double QlaE0 = point[iQlaE0];
			final double phi0 = point[iPhi0];
			final double R = point[iR];

			final double logZbarb = Math.log(Zbarb);
			final double X = 1.0 + 1.3 * logZbarb;
			final double Y = 0.2 + 0.005 * Zbarb;
			final double pu42 = Math.pow(u0, 0.42);
			final double FoRbar = 1.0 + (X * Math.log(1.0 + Y * (1.0 - 1.0 / pu42)) / Math.log(1.0 + Y));

			final RealVector vals = buildResult();

			final int oRbar = outputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int oF = outputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			checkIndices(oRbar, oF);

			// F and partials
			final double F = R * oneOverS / QlaE0; // C1
			vals.setEntry(oF, F);

			// Rbar and partials
			double Rbar = Double.NaN;
			if (FoRbar >= phi0)
				Rbar = F / FoRbar; // C1
			else
				Rbar = F / phi0;
			vals.setEntry(oRbar, Rbar);
			return vals;
		}

		@Override
		public String toString() {
			return "StepFRbar[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iE0 = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			final int iOneOverS = inputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int iQlaE0 = inputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iR = inputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			checkIndices(iZbarb, iOneOverS, iQlaE0, iPhi0, iR);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = point.getEntry(iE0) / Ea;
			final double Zbarb = point.getEntry(iZbarb);
			final double oneOverS = point.getEntry(iOneOverS);
			final double QlaE0 = point.getEntry(iQlaE0);
			final double phi0 = point.getEntry(iPhi0);
			final double R = point.getEntry(iR);

			final double logZbarb = Math.log(Zbarb);
			final double X = 1.0 + 1.3 * logZbarb; // Ok
			final double Y = 0.2 + 0.005 * Zbarb; // Ok
			final double dXdZbarb = 1.3 / Zbarb, dYdZbarb = 0.005; // Ok
			final double pu42 = Math.pow(u0, 0.42);
			final double log1pY = Math.log(1.0 + Y);
			final double arg = 1.0 + Y * (1.0 - 1.0 / pu42);
			final double logArg = Math.log(arg);
			final double FoRbar = 1.0 + X * logArg / log1pY; // Corrected 15-Jan-2019
			final double dFoRbardX = logArg / log1pY; // Ok
			final double dFoRbardY = X * ((pu42 - 1.0) * log1pY / (pu42 + Y * (pu42 - 1.0)) - logArg / (1.0 + Y)) //
					/ Math.pow(log1pY, 2.0); // Ok
			final double dFoRbardZbarb = dFoRbardX * dXdZbarb + dFoRbardY * dYdZbarb; // Ok
			final double dFoRbardu0 = (0.42 * X * Y) / (u0 * pu42 * log1pY * arg); // Ok

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oRbar = outputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int oF = outputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			checkIndices(oRbar, oF);

			// F and partials
			final double F = R * oneOverS / QlaE0; // C1
			vals.setEntry(oF, F);
			final double dFdR = oneOverS / QlaE0; // C1
			jac.setEntry(oF, iR, dFdR);
			final double dFdOneOverS = R / QlaE0; // C1
			jac.setEntry(oF, iOneOverS, dFdOneOverS);
			final double dFdQlaE0 = -1.0 * R * oneOverS / Math.pow(QlaE0, 2.0); // C1
			jac.setEntry(oF, iQlaE0, dFdQlaE0);

			// Rbar and partials
			double Rbar, dRbardF;
			if (FoRbar >= phi0) {
				Rbar = F / FoRbar; // C1
				dRbardF = 1.0 / FoRbar;
				final double dU0dE0 = 1.0 / Ea;
				final double dRbardFoRbar = -F / Math.pow(FoRbar, 2); // Ok
				jac.setEntry(oRbar, iE0, dRbardFoRbar * dFoRbardu0 * dU0dE0); // Ok
				jac.setEntry(oRbar, iZbarb, dRbardFoRbar * dFoRbardZbarb); // Ok
			} else {
				Rbar = F / phi0;
				dRbardF = 1.0 / phi0;
				jac.setEntry(oRbar, iPhi0, -F / Math.pow(phi0, 2));
				// dRbardE0 and dRbardZbarb = 0!
			}
			vals.setEntry(oRbar, Rbar);
			jac.setEntry(oRbar, iR, dRbardF * dFdR);
			jac.setEntry(oRbar, iOneOverS, dRbardF * dFdOneOverS);
			jac.setEntry(oRbar, iQlaE0, dRbardF * dFdQlaE0);
			return Pair.create(vals, jac);
		}

	}

	private static class StepFx // Checked 16-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr) {
			final List<EPMALabel> res = new ArrayList<>();
			final AtomicShell shell = cxr.getInner();
			res.add(MatrixCorrectionModel2.shellLabel("A", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("B", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("a", datum, shell));
			res.add(MatrixCorrectionModel2.chiLabel(datum, cxr));
			res.add(MatrixCorrectionModel2.roughnessLabel(datum));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.FofChiLabel(datum, cxr));
		}

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		public StepFx(final MatrixCorrectionDatum unk, //
				final CharacteristicXRay cxr) throws ArgumentException {
			super(buildInputs(unk, cxr), buildOutputs(unk, cxr));
			mDatum = unk;
			mXRay = cxr;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iA = inputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mXRay.getInner()));
			final int iB = inputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mXRay.getInner()));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mXRay.getInner()));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mXRay.getInner()));
			final int iChi = inputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mXRay.getInner()));
			checkIndices(iA, iB, ib, iPhi0, iChi, ia);

			final double A = point[iA];
			final double B = point[iB];
			final double b = point[ib];
			final double phi0 = point[iPhi0];
			final double chi = point[iChi];
			final double a = point[ia];
			final double dz = point[inputIndex(MatrixCorrectionModel2.roughnessLabel(mDatum))];

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector vals = buildResult();

			final double k0 = b + chi;
			final double eps = (a - b) / b;
			final double k1 = chi + b * (1 + eps);
			// Must include dz here since z = 0 +- dz != 0
			final double Fx = Math.exp(-chi * dz) * (B / k0 - (A * b * eps) / k1 + phi0) / k0;

			vals.setEntry(oFx, Fx); // C2
			return vals;
		}

		@Override
		public String toString() {
			return "StepFx[" + mDatum + "," + mXRay + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iA = inputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mXRay.getInner()));
			final int iB = inputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mXRay.getInner()));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mXRay.getInner()));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mXRay.getInner()));
			final int iChi = inputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mXRay.getInner()));
			final int itr = inputIndex(MatrixCorrectionModel2.roughnessLabel(mDatum));
			checkIndices(iA, iB, ib, iPhi0, iChi, ia, itr);

			final double A = point.getEntry(iA);
			final double B = point.getEntry(iB);
			final double b = point.getEntry(ib);
			final double phi0 = point.getEntry(iPhi0);
			final double chi = point.getEntry(iChi);
			final double a = point.getEntry(ia);
			final double dz = point.getEntry(itr);

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final double expDzChi = Math.exp(-chi * dz);

			final double eps = (a - b) / b;
			final double bpchi = b + chi;
			final double arg = chi + b * (1.0 + eps);
			final double powArg = Math.pow(arg, -2.0);
			final double powBpchi2 = Math.pow(bpchi, -2.0);

			final double Fx = expDzChi * (B / bpchi - (A * b * eps) / arg + phi0) / bpchi;
			final double dFxdA = -expDzChi * (b * eps) / (bpchi * arg);
			final double dFxdB = expDzChi * powBpchi2;
			final double dFxdphi0 = expDzChi / bpchi;
			final double dFxdchi = (expDzChi * ((A * b * eps) * powArg - B * powBpchi2) - (1 + bpchi * dz) * Fx)
					/ bpchi;
			final double dFxddz = -chi * Fx;

			final double dFxdb = expDzChi * ((A - phi0) * (b + chi) - 2 * B) / Math.pow(b + chi, 3.0);
			final double dFxda = -A * expDzChi / Math.pow(a + chi, 2.0);

			vals.setEntry(oFx, Fx); // C2
			jac.setEntry(oFx, iA, dFxdA); // C2
			jac.setEntry(oFx, iB, dFxdB); // C2
			jac.setEntry(oFx, ib, dFxdb); // C2
			jac.setEntry(oFx, iPhi0, dFxdphi0); // C2
			jac.setEntry(oFx, iChi, dFxdchi); // C2
			jac.setEntry(oFx, ia, dFxda); // C2
			jac.setEntry(oFx, itr, dFxddz);

			return Pair.create(vals, jac);
		}
	}

	private static class StepMJZBarb // Checked 14-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final Material comp, //
				final boolean isStandard //
		) {
			final List<EPMALabel> res = new ArrayList<>();
			for (final Element elm : comp.getElementSet()) {
				if (isStandard)
					res.add(MaterialLabel.buildMassFractionTag(comp, elm));
				res.add(MaterialLabel.buildAtomicWeightTag(comp, elm));
				res.add(MatrixCorrectionModel2.meanIonizationLabel(elm));
			}
			return res;
		}

		public static List<EPMALabel> buildOutputs(final Material mat) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel("M", mat));
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel("J", mat));
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mat));
			return res;
		}

		private final Material mMaterial;

		public StepMJZBarb(final Material mat, //
				final boolean isStandard //
		) throws ArgumentException {
			super(buildInputs(mat, isStandard), buildOutputs(mat));
			mMaterial = mat;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final List<Element> elms = new ArrayList<>(mMaterial.getElementSet());
			final double[] Ci = new double[elms.size()];
			final double[] ZoA = new double[elms.size()];
			final double[] Z = new double[elms.size()];
			final double[] Ji = new double[elms.size()];
			for (int i = 0; i < Ci.length; ++i) {
				final Element elm = elms.get(i);
				Ci[i] = getArg(MaterialLabel.buildMassFractionTag(mMaterial, elm), point);
				Z[i] = elm.getAtomicNumber();
				Ji[i] = point[inputIndex(MatrixCorrectionModel2.meanIonizationLabel(elm))];
				final double a = point[inputIndex(MaterialLabel.buildAtomicWeightTag(mMaterial, elm))];
				ZoA[i] = Z[i] / a;
			}

			final RealVector vals = buildResult();

			final int oM = outputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("M", mMaterial));
			final int oJ = outputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("J", mMaterial));
			final int oZbarb = outputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mMaterial));
			checkIndices(oM, oJ, oZbarb);

			// Calculate M
			double M = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				M += Ci[i] * ZoA[i];
			vals.setEntry(oM, M); // c1
			// Calculate J
			double lnJ = 0.0;
			for (int i = 0; i < Ji.length; ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * ZoA[i];
			lnJ /= M;
			vals.setEntry(oJ, Math.exp(lnJ)); // C2
			// Calculate Zbarb
			double Zbt = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // Ok
			assert Zbt > 0.0 : //
			Zbt;
			vals.setEntry(oZbarb, Zbt * Zbt);
			return vals;
		}

		@Override
		public String toString() {
			return "MJZBarb[" + mMaterial.toString() + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final ArrayList<Element> elms = new ArrayList<>(mMaterial.getElementSet());
			final MassFraction[] lCi = new MassFraction[elms.size()];
			final int[] iJi = new int[elms.size()];
			final int[] iAi = new int[elms.size()];
			final double[] Ci = new double[elms.size()];
			final double[] Ai = new double[elms.size()];
			final double[] Z = new double[elms.size()];
			final double[] Ji = new double[elms.size()];
			for (int i = 0; i < Ci.length; ++i) {
				final Element elm = elms.get(i);
				lCi[i] = MaterialLabel.buildMassFractionTag(mMaterial, elm);
				iJi[i] = inputIndex(MatrixCorrectionModel2.meanIonizationLabel(elm));
				iAi[i] = inputIndex(MaterialLabel.buildAtomicWeightTag(mMaterial, elm));
				Ci[i] = getArg(lCi[i], point);
				assert Ci[i] >= 0.0 : elm.getAbbrev() + "=" + Ci[i];
				Ji[i] = point.getEntry(iJi[i]);
				Ai[i] = point.getEntry(iAi[i]);
				Z[i] = elm.getAtomicNumber();
			}

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final Material mat = mMaterial;
			final MatrixCorrectionModel2.MaterialBasedLabel lM = new MatrixCorrectionModel2.MaterialBasedLabel("M",
					mat);
			final int oM = outputIndex(lM);
			final MatrixCorrectionModel2.MaterialBasedLabel lJ = new MatrixCorrectionModel2.MaterialBasedLabel("J",
					mat);
			final int oJ = outputIndex(lJ);
			final MatrixCorrectionModel2.MaterialBasedLabel lZbarb = new MatrixCorrectionModel2.MaterialBasedLabel(
					ZBARB, mat);
			final int oZbarb = outputIndex(lZbarb);
			checkIndices(oM, oJ, oZbarb);

			// Calculate M & partials
			double M = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				M += Ci[i] * Z[i] / Ai[i]; // Cx
			vals.setEntry(oM, M); // c1
			for (int i = 0; i < lCi.length; ++i) {
				setJacobian(lCi[i], lM, jac, (Z[i] / Ai[i])); // Cx
				jac.setEntry(oM, iAi[i], -Ci[i] * (Z[i] / (Ai[i] * Ai[i]))); // Cx
			}

			// Calculate J and partials
			double lnJ = 0.0;
			for (int i = 0; i < Ji.length; ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * Z[i] / Ai[i];
			lnJ /= M;
			final double J = Math.exp(lnJ); // keV
			vals.setEntry(oJ, J); // C2
			for (int i = 0; i < lCi.length; ++i) {
				final double lnJi = Math.log(Ji[i]);
				final double kk = (J / M) * (Z[i] / Ai[i]);
				jac.setEntry(oJ, iJi[i], kk * (Ci[i] / Ji[i])); // Cx
				setJacobian(lCi[i], lJ, jac, kk * (lnJi - lnJ)); // Cx
				jac.setEntry(oJ, iAi[i], kk * (Ci[i] / Ai[i]) * (lnJ - lnJi)); // Cx
			}
			// Calculate Zbarb and partials
			double Zbt = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // Cx
			vals.setEntry(oZbarb, Zbt * Zbt);
			for (int i = 0; i < elms.size(); ++i)
				setJacobian(lCi[i], lZbarb, jac, 2.0 * Math.sqrt(Z[i] * Zbt)); // Cx (fixed 26-Jul-2019)
			return Pair.create(vals, jac);
		}

	}

	private static class StepPb // Checked 16-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, datum.getMaterial()));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		StepPb(final MatrixCorrectionDatum comp, //
				final AtomicShell shell) throws ArgumentException {
			super(buildInputs(comp, shell), buildOutputs(comp, shell));
			mDatum = comp;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(iZbarb, iRbar, iF, iphi0);

			final double Zbarb = point[iZbarb];
			final double e0 = point[inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum))];
			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			final double Rbar = point[iRbar];
			final double F = point[iF];
			final double phi0 = point[iphi0];

			// P and partials
			final double expArg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
			final double g = 0.22 * Math.log(4.0 * Zbarb) * (1.0 - 2.0 * expArg);
			final double h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 2.0);
			final double b = Math.sqrt(2.0) * (1.0 + Math.sqrt(1.0 - Rbar * phi0 / F)) / Rbar;
			final double Rbar2 = Math.pow(Rbar, 2.0);

			// Two different ways to compute gh4
			final double gh4_1 = g * Math.pow(h, 4.0);
			final double gh4_2 = 0.9 * b * Rbar2 * (b - 2.0 * phi0 / F);

			final RealVector vals = buildResult();

			final int oP = outputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ob = outputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(oP, ob);

			if (gh4_1 < gh4_2) {
				// Depends on Zbarb, u0, F, Rbar
				final double P = gh4_1 * F / Rbar2;
				vals.setEntry(oP, P);
			} else {
				// Depends on F, Rbar, b, phi0
				final double P = gh4_2 * F / Rbar2;
				vals.setEntry(oP, P);
			}
			vals.setEntry(ob, b);

			return vals;
		}

		@Override
		public String toString() {
			return "StepP[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			final int iE0 = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(iZbarb, iRbar, iF, iphi0);

			final double Zbarb = point.getEntry(iZbarb);
			final double e0 = point.getEntry(iE0);
			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double Rbar = point.getEntry(iRbar);
			final double F = point.getEntry(iF);
			final double phi0 = point.getEntry(iphi0);

			final double u0 = e0 / Ea;
			final double du0de0 = 1.0 / Ea;
			final double k11_375 = 11.0 / 375.0;
			// P and partials
			final double kg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
			final double log4Zbarb = Math.log(4.0 * Zbarb);
			final double g = 0.22 * log4Zbarb * (1.0 - 2.0 * kg); // Ok
			final double dgdzbarb = g / (Zbarb * log4Zbarb) + k11_375 * kg * (u0 - 1.0) * log4Zbarb; // Ok
			final double dgdu0 = k11_375 * kg * Zbarb * log4Zbarb;

			final double h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 2.0); // Ok
			final double dhdzbarb = 20.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 3.0); // Ok
			final double dhdu0 = -Math.pow(((1 + u0 / 10) * Zbarb), -2.0);

			final double b = Math.sqrt(2.0) * (1.0 + Math.sqrt(1.0 - Rbar * phi0 / F)) / Rbar;
			final double dbdrbar = -b / Rbar - phi0 / (F * Rbar * Math.sqrt(2.0 * (1.0 - phi0 * Rbar / F)));
			final double dbdphi0 = -1.0 / (F * Math.sqrt(2.0 * (1.0 - phi0 * Rbar / F)));
			final double dbdF = phi0 / (F * F * Math.sqrt(2 * (1.0 - phi0 * Rbar / F)));

			// Two different ways to compute gh4
			final double gh4_1 = g * Math.pow(h, 4.0);
			final double gh4_2 = 0.9 * b * Math.pow(Rbar, 2.0) * (b - 2.0 * phi0 / F);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oP = outputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ob = outputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(oP, ob);

			if (gh4_1 < gh4_2) {
				// Use gh4_1. Depends on Zbarb, u0, F, Rbar
				final double P = gh4_1 * F / Math.pow(Rbar, 2.0);
				final double dPdzbarb = P * ((1.0 / g) * dgdzbarb + (4.0 / h) * dhdzbarb);
				final double dPdzu0 = P * ((1 / g) * dgdu0 + (4.0 / h) * dhdu0);

				vals.setEntry(oP, P);
				jac.setEntry(oP, iZbarb, dPdzbarb);
				jac.setEntry(oP, iE0, dPdzu0 * du0de0);
				jac.setEntry(oP, iF, P / F);
				jac.setEntry(oP, iRbar, -2.0 * P / Rbar);
			} else {
				// Depends on F, Rbar, b, phi0
				final double P = 0.9 * b * (b * F - 2.0 * phi0);
				final double dPdb = 1.8 * (b * F - phi0);
				final double dPdF = dPdb * dbdF + 0.9 * b * b;
				final double dPdphi0 = -1.8 * b;
				vals.setEntry(oP, P);
				jac.setEntry(oP, iphi0, dPdphi0 + dPdb * dbdphi0);
				jac.setEntry(oP, iF, dPdF);
				jac.setEntry(oP, iRbar, dPdb * dbdrbar);
				jac.setEntry(oP, iRbar, 0.0);
			}

			final double k3 = Math.sqrt(1. - (phi0 * Rbar) / F);

			vals.setEntry(ob, b);
			jac.setEntry(ob, iRbar,
					(-1. * (-0.5 * phi0 * Rbar + F * (1. + k3)) * Math.sqrt(2.0)) / (F * Math.pow(Rbar, 2.0) * k3));
			jac.setEntry(ob, iphi0, -Math.sqrt(2.0) / (2. * F * k3));
			jac.setEntry(ob, iF, (phi0 * Math.sqrt(2.0)) / (2. * Math.pow(F, 2) * k3));

			return Pair.create(vals, jac);
		}

	}

	private static class StepQlaE0OneOverS // Checked 14-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(
				//
				final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel("M", datum.getMaterial()));
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel("J", datum.getMaterial()));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(new MatrixCorrectionModel2.IonizationExponentLabel(shell));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			return res;
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		public StepQlaE0OneOverS(final MatrixCorrectionDatum datum, //
				final AtomicShell shell) throws ArgumentException {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int tagE0 = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final int tagm = inputIndex(new MatrixCorrectionModel2.IonizationExponentLabel(mShell));
			final int iJ = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("J", mDatum.getMaterial()));
			final int iM = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("M", mDatum.getMaterial()));
			checkIndices(iJ, iM);

			final double e0 = point[tagE0];
			final double m = point[tagm];
			final double J = point[iJ];
			final double M = point[iM];

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			assert u0 > 1.0 : e0 + " over " + Ea + " for " + mDatum + " - " + mShell;
			final double logU0 = Math.log(u0);
			// QlaE0 and partials
			final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
			vals.setEntry(oQlaE0, QlaE0); // C1
			jac.setEntry(oQlaE0, tagE0, Math.pow(u0, -1.0 - m) * (1.0 - m * logU0) / Math.pow(Ea, 3.0)); // C1
			jac.setEntry(oQlaE0, tagm, -QlaE0 * logU0);

			final double[] D = { 6.6e-6, 1.12e-5 * (1.35 - 0.45 * J * J), 2.2e-6 / J };
			final double[] P = { 0.78, 0.1, 0.25 * J - 0.5 };
			final double[] T = { 1.0 + P[0] - m, 1.0 + P[1] - m, 1.0 + P[2] - m };

			final double v0ou0 = Ea / J; // c1

			double OoS = 0.0;
			final double kk = J / (Ea * M);
			for (int k = 0; k < 3; ++k) {
				final double u0tk = Math.pow(u0, T[k]);
				final double v0ou0_pk = Math.pow(v0ou0, P[k]);
				OoS += kk * D[k] * v0ou0_pk * (u0tk * (T[k] * logU0 - 1.0) + 1.0) / (T[k] * T[k]); // d
			}
			vals.setEntry(oOneOverS, OoS); // C2
			return vals;
		}

		@Override
		public String toString() {
			return "StepQlaE0OneOverS[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int tagE0 = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final int tagm = inputIndex(new MatrixCorrectionModel2.IonizationExponentLabel(mShell));
			final int iJ = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("J", mDatum.getMaterial()));
			final int iM = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel("M", mDatum.getMaterial()));
			checkIndices(iJ, iM);

			final double e0 = point.getEntry(tagE0);
			final double m = point.getEntry(tagm);
			final double J = point.getEntry(iJ);
			final double M = point.getEntry(iM);

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final double Ea = eVtokeV(mShell.getEdgeEnergy());

			final double u0 = e0 / Ea;
			final double du0de0 = 1 / Ea;
			assert u0 > 1.0 : e0 + " over " + Ea;
			final double logU0 = Math.log(u0), logU0_2 = logU0 * logU0;

			final double v0 = e0 / J;
			final double dv0dJ = -e0 / (J * J);
			final double v0ou0 = Ea / J; // V0/U0 = (E0/J)/(E0/Ea) = Ea/J
			// QlaE0 and partials
			final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
			vals.setEntry(oQlaE0, QlaE0); // C1
			jac.setEntry(oQlaE0, tagE0, (1.0 - m * logU0) / (Math.pow(u0, 1.0 + m) * Math.pow(Ea, 3.0))); // C1
			jac.setEntry(oQlaE0, tagm, -QlaE0 * logU0); // Ok

			final double[] D = { 6.6e-6, 1.12e-5 * (1.35 - 0.45 * J * J), 2.2e-6 / J }; // Ok!
			final double[] P = { 0.78, 0.1, 0.25 * J - 0.5 }; // Ok!
			final double[] T = { 1.0 + P[0] - m, 1.0 + P[1] - m, 1.0 + P[2] - m }; // Ok!

			final double[] dDdJ = { 0.0, 2.0 * 1.12e-5 * (-0.45 * J), -2.2e-6 / (J * J) }; // Ok!
			final double[] dPdJ = { 0.0, 0.0, 0.25 }; // Ok!
			final double[] dTdJ = dPdJ; // Ok!
			final double dTdmk = -1.0; // For all 3

			// Compute the sum of the k-terms
			double OoS = 0.0, dOoSdJ = 0.0, dOoSdU0 = 0.0, dOoSdm = 0.0;
			final double kk = J / (Ea * M); // J dependence in sum
			for (int k = 0; k < 3; ++k) {
				final double u0tk = Math.pow(u0, T[k]);
				final double v0ou0_pk = Math.pow(v0ou0, P[k]);
				final double OoSk = kk * D[k] * v0ou0_pk * (u0tk * (T[k] * logU0 - 1.0) + 1.0) / (T[k] * T[k]); // d
				OoS += OoSk;
				dOoSdm += (kk * u0tk * v0ou0_pk * D[k] * logU0_2 - 2.0 * OoSk) * (dTdmk / T[k]); // d
				dOoSdJ += (1.0 / J + dDdJ[k] / D[k] + Math.log(v0ou0) * dPdJ[k] + (P[k] / v0) * dv0dJ
						- 2.0 * dTdJ[k] / T[k]) * OoSk + //
						kk * u0tk * v0ou0_pk * D[k] * logU0_2 * dTdJ[k] / T[k]; // d
				dOoSdU0 += (kk * u0tk * v0ou0_pk * D[k] * logU0) / u0; // d
			}
			// E0/J0 = J/Ea
			vals.setEntry(oOneOverS, OoS); // C2
			jac.setEntry(oOneOverS, iM, (-1.0 / M) * OoS); // C2
			jac.setEntry(oOneOverS, iJ, dOoSdJ); // C2
			jac.setEntry(oOneOverS, tagE0, dOoSdU0 * du0de0); // C2
			jac.setEntry(oOneOverS, tagm, dOoSdm);

			return Pair.create(vals, jac);
		}
	}

	private static class StepRphi0 // Checked 15-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final MatrixCorrectionDatum datum) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, datum.getMaterial()));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		private final MatrixCorrectionDatum mDatum;

		private final AtomicShell mShell;

		public StepRphi0(final MatrixCorrectionDatum datum, //
				final AtomicShell shell) throws ArgumentException {
			super(buildInputs(datum), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			checkIndices(iZbarb);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double e0 = point[inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum))];
			final double u0 = e0 / Ea;
			final double Zbarb = point[iZbarb];

			final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3))); // Ok!
			final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55); // Ok!
			final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar); // Ok!
			final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0); // Ok!
			final double Gu0 = (u0 * (1.0 + q) - (2.0 + q) + Math.pow(u0, -1.0 - q)) / ((1.0 + q) * (2.0 + q) * Ju0);
			final double R = 1.0 - etabar * Wbar * (1.0 - Gu0);

			final double u0_mr = Math.pow(u0, 2.3 * etabar - 2.0);
			final double phi0 = 1.0 + 3.3 * Math.pow(etabar, 1.2) * (1.0 - u0_mr); // Ok!

			final RealVector vals = buildResult();

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			vals.setEntry(oR, R);

			vals.setEntry(oPhi0, phi0);
			return vals;
		}

		@Override
		public String toString() {
			return "StepRPhi0[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.MaterialBasedLabel(ZBARB, mDatum.getMaterial()));
			final int tagE0 = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			checkIndices(iZbarb);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double Zbarb = point.getEntry(iZbarb);
			final double e0 = point.getEntry(tagE0);

			final double u0 = e0 / Ea;
			final double du0de0 = 1.0 / Ea;

			final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3))); // Ok!
			final double detabardZbarb = 1.75e-3 //
					+ 7.215e-3 * Math.pow(Zbarb, 0.3) * Math.exp(-0.015 * Math.pow(Zbarb, 1.3)); // Ok!

			final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55); // Ok!
			final double dWbardeta = 1.0 / 3.7 + 4.55 * Math.pow(etabar, 3.55); // Ok!
			final double dWbardZbarb = dWbardeta * detabardZbarb; // Ok!

			final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar); // Ok!
			final double dqdWbar = Math.pow(Wbar - 1.0, -2.0); // Ok!
			final double dqdZbarb = dqdWbar * dWbardZbarb; // Ok

			final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0); // Ok!
			final double dJu0du0 = Math.log(u0); // Ok!

			final double opq = 1.0 + q;
			final double Gu0 = (u0 * opq - (2.0 + q) + Math.pow(u0, -1.0 - q)) / //
					(opq * (2.0 + q) * Ju0); // Ok!
			// assert Math.abs(Gu0 - ((u0 - 1.0 - (1.0 - Math.pow(u0, -1.0 - q)) / opq) /
			// ((2.0 + q) * Ju0))) < Math.abs(1.0e-6 * Gu0);
			final double dGu0du0 = (1.0 - Math.pow(u0, -2.0 - q)) / ((2.0 + q) * Ju0) //
					- ((dJu0du0 / Ju0) * Gu0); // Ok!
			final double dGu0de0 = dGu0du0 * du0de0; // Ok!
			final double dGu0dq = (((Math.pow(u0, opq) - 1.0) - opq * Math.log(u0)) //
					/ (Ju0 * Math.pow(u0, opq) * Math.pow(opq, 2.0)) - Gu0) //
					/ (2.0 + q); // Ok!
			final double dGu0dZbarb = dGu0dq * dqdZbarb;

			final double R = 1.0 - etabar * Wbar * (1.0 - Gu0); // Ok!
			final double dRde0 = etabar * Wbar * dGu0de0; // Ok!
			final double dRdZbarb = etabar * (Wbar * dGu0dZbarb - (1.0 - Gu0) * dWbardZbarb) - //
					(1.0 - Gu0) * Wbar * detabardZbarb; // Ok!

			final double u0_mr = Math.pow(u0, 2.3 * etabar - 2.0);
			final double phi0 = 1.0 + 3.3 * Math.pow(etabar, 1.2) * (1.0 - u0_mr); // Ok!
			final double dphi0detabar = Math.pow(etabar, 0.2) //
					* (3.96 * (1.0 - u0_mr) - 7.59 * etabar * u0_mr * Math.log(u0)); // Ok!
			final double dphi0dZbarb = dphi0detabar * detabardZbarb; // Ok!
			final double dphi0du0 = -3.3 * Math.pow(etabar, 1.2) * (2.3 * etabar - 2.0) * u0_mr / u0; // Ok!
			final double dphi0de0 = dphi0du0 * du0de0; // Ok!

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			vals.setEntry(oR, R);
			jac.setEntry(oR, iZbarb, dRdZbarb);
			jac.setEntry(oR, tagE0, dRde0);

			vals.setEntry(oPhi0, phi0);
			jac.setEntry(oPhi0, iZbarb, dphi0dZbarb);
			jac.setEntry(oPhi0, tagE0, dphi0de0);

			return Pair.create(vals, jac);
		}

	}

	/**
	 * This class performs the full XPP calculation for a single composition and the
	 * associated set of characteristic x-rays. It breaks the calculation into a
	 * series of steps which involve 1) only the composition; 2) the composition and
	 * the atomic shells; and 3) the composition and the characteristic x-rays. It
	 * is designed to break the calculation down into small independent blocks each
	 * of which is relatively fast to calculate. This minimizes and postpones the
	 * need to build and copy large Jacobian matrices and is thus slightly more
	 * computationally efficient.
	 *
	 * @author Nicholas W. M. Ritchie
	 */
	private static final class StepXPP //
			extends CompositeMeasurementModel<EPMALabel> //
	{

		private static List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps(
				final MatrixCorrectionDatum datum, final CharacteristicXRaySet exrs) throws ArgumentException {
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
			res.add(new StepMJZBarb(datum.getMaterial(), datum instanceof StandardMatrixCorrectionDatum));
			final Set<AtomicShell> shells = exrs.getSetOfInnerAtomicShells();
			{
				final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
				for (final AtomicShell shell : shells)
					step.add(new StepQlaE0OneOverS(datum, shell));
				res.add(ParallelMeasurementModelBuilder.join("QlaOoS", step));
			}
			{
				final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
				for (final AtomicShell shell : shells)
					step.add(new StepAaBb(datum, shell));
				res.add(ParallelMeasurementModelBuilder.join("AaBb", step));
			}
			{
				final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepChi(datum, cxr));
				res.add(ParallelMeasurementModelBuilder.join("Chi", step));
			}
			{
				final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepFx(datum, cxr));
				res.add(ParallelMeasurementModelBuilder.join("Fx", step));
			}
			{
				final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepConductiveCoating(datum, cxr));
				res.add(ParallelMeasurementModelBuilder.join("Coating", step));
			}
			return res;
		}

		public StepXPP(final MatrixCorrectionDatum datum, final CharacteristicXRaySet cxrs,
				final List<EPMALabel> outputs) throws ArgumentException {
			super("XPP[" + datum + "]", buildSteps(datum, cxrs), outputs);
		}

	};

	private static class StepZA // Checked 16-Jan-2019
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> {

		public static List<EPMALabel> buildInputs(final UnknownMatrixCorrectionDatum unk, //
				final MatrixCorrectionDatum std, //
				final CharacteristicXRay cxr //
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FofChiReducedLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", unk, cxr.getInner()));
			res.add(MatrixCorrectionModel2.FofChiReducedLabel(std, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", std, cxr.getInner()));
			return res;
		}

		public static List<EPMALabel> buildOutputs(final UnknownMatrixCorrectionDatum unk,
				final StandardMatrixCorrectionDatum std, final CharacteristicXRay cxr) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FxFLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.FxFLabel(std, cxr));
			res.add(MatrixCorrectionModel2.zLabel(unk, std, cxr));
			res.add(MatrixCorrectionModel2.aLabel(unk, std, cxr));
			res.add(MatrixCorrectionModel2.zaLabel(unk, std, cxr));
			return res;
		}

		private final UnknownMatrixCorrectionDatum mUnknown;

		private final StandardMatrixCorrectionDatum mStandard;

		private final CharacteristicXRay mXRay;

		public StepZA(final UnknownMatrixCorrectionDatum unk, //
				final StandardMatrixCorrectionDatum std, //
				final CharacteristicXRay cxr) throws ArgumentException {
			super(buildInputs(unk, std, cxr), buildOutputs(unk, std, cxr));
			mUnknown = unk;
			mStandard = std;
			mXRay = cxr;
		}

		@Override
		public RealVector computeValue(final double[] point) {
			final int iFxu = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mUnknown, mXRay));
			final int iFxs = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mStandard, mXRay));
			final int iFu = inputIndex(MatrixCorrectionModel2.shellLabel("F", mUnknown, mXRay.getInner()));
			final int iFs = inputIndex(MatrixCorrectionModel2.shellLabel("F", mStandard, mXRay.getInner()));
			checkIndices(iFxu, iFxs, iFu, iFs);

			final double Fxu = point[iFxu];
			final double Fxs = point[iFxs];
			final double Fu = point[iFu];
			final double Fs = point[iFs];

			final RealVector vals = buildResult();

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			final int oZA = outputIndex(MatrixCorrectionModel2.zaLabel(mUnknown, mStandard, mXRay));
			checkIndices(oFxFu, oFxFs, oZ, oA);

			vals.setEntry(oA, (Fs * Fxu) / (Fu * Fxs)); // C2
			vals.setEntry(oZ, (Fu / Fs)); // C2
			vals.setEntry(oZA, Fxu / Fxs);
			vals.setEntry(oFxFu, Fxu / Fu); // C2
			vals.setEntry(oFxFs, Fxs / Fs); // C2
			return vals;
		}

		@Override
		public String toString() {
			return "StepZA[" + mStandard + "," + mUnknown + "," + mXRay + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iFxu = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mUnknown, mXRay));
			final int iFxs = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mStandard, mXRay));
			final int iFu = inputIndex(MatrixCorrectionModel2.shellLabel("F", mUnknown, mXRay.getInner()));
			final int iFs = inputIndex(MatrixCorrectionModel2.shellLabel("F", mStandard, mXRay.getInner()));

			checkIndices(iFxu, iFxs, iFu, iFs);

			final double Fxu = point.getEntry(iFxu);
			final double Fxs = point.getEntry(iFxs);
			final double Fu = point.getEntry(iFu);
			final double Fs = point.getEntry(iFs);

			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			final int oZA = outputIndex(MatrixCorrectionModel2.zaLabel(mUnknown, mStandard, mXRay));
			checkIndices(oFxFu, oFxFs, oZ, oA);

			final double a = (Fs * Fxu) / (Fu * Fxs);
			vals.setEntry(oA, a); // C2
			jac.setEntry(oA, iFxu, a / Fxu); // C2
			jac.setEntry(oA, iFxs, -a / Fxs); // C2
			jac.setEntry(oA, iFs, a / Fs); // C2
			jac.setEntry(oA, iFu, -a / Fu); // C2

			final double z = Fu / Fs;
			vals.setEntry(oZ, z); // C2
			jac.setEntry(oZ, iFs, -z / Fs); // C2
			jac.setEntry(oZ, iFu, 1.0 / Fs); // C2

			final double za = Fxu / Fxs;
			vals.setEntry(oZA, za);
			jac.setEntry(oZA, iFxu, 1.0 / Fxs);
			jac.setEntry(oZA, iFxs, -za / Fxs);

			final double FxFu = Fxu / Fu;
			vals.setEntry(oFxFu, FxFu); // C2
			jac.setEntry(oFxFu, iFu, -FxFu / Fu); // C2
			jac.setEntry(oFxFu, iFxu, 1.0 / Fu); // C2

			final double FxFs = Fxs / Fs;
			vals.setEntry(oFxFs, FxFs); // C2
			jac.setEntry(oFxFs, iFs, -FxFs / Fs); // C2
			jac.setEntry(oFxFs, iFxs, 1.0 / Fs); // C2

			return Pair.create(vals, jac);
		}
	}

	private static final String ONE_OVER_S = "<sup>1</sup>/<sub>S</sub>";

	private static final String QLA = "<html>Q<sub>l</sub><sup>a</sup>";

	private static final String ZBARB = "Z<sub>barb</sub>";

	private static final String RBAR = "R<sub>bar</sub>";

	/**
	 * Builds a {@link UncertainValuesCalculator} around an instance of the
	 * {@link XPPMatrixCorrection2} algorithm that returns a basic set of default
	 * outputs.
	 *
	 * @param kratios
	 * @param estUnknown
	 * @return {@link UncertainValuesCalculator}
	 * @throws ArgumentException
	 */
	public static UncertainValuesCalculator<EPMALabel> buildAnalytical(final Set<KRatioLabel> kratios, //
			final Map<MassFraction, Double> estUnknown, //
			final boolean allOutputs) throws ArgumentException {
		Material unkMat = null;
		for (final MassFraction mf : estUnknown.keySet()) {
			if (unkMat == null)
				unkMat = mf.getMaterial();
			if (!unkMat.equals(mf.getMaterial()))
				throw new ArgumentException("Two different materials [" + mf.getMaterial() + " and " + unkMat
						+ " are present in the estimated unknown.");
		}
		final List<EPMALabel> outputs = allOutputs ? Collections.emptyList() : buildDefaultOutputs(kratios);
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(kratios, outputs);
		final UncertainValues<EPMALabel> inputs = xpp.buildInput(unkMat);
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(estUnknown);
		return new UncertainValuesCalculator<EPMALabel>(xpp, inputs);
	}

	/**
	 * Returns a very basic set of output that includes the ZAF, Z, A terms and the
	 * C<sub>Z</sub> in the standard for each k-ratio in the set.
	 *
	 *
	 * @param kratios Set&lt;EPMALabel&gt;
	 * @return List<EPMALabel>
	 */
	static public List<EPMALabel> buildDefaultOutputs(final Set<KRatioLabel> kratios) {
		final FastIndex<EPMALabel> res = new FastIndex<>();
		for (final KRatioLabel krl : kratios) {
			final UnknownMatrixCorrectionDatum unk = krl.getUnknown();
			final StandardMatrixCorrectionDatum std = krl.getStandard();
			res.add(MatrixCorrectionModel2.zafLabel(krl));
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				res.add(MatrixCorrectionModel2.zLabel(unk, std, cxr));
				res.add(MatrixCorrectionModel2.aLabel(unk, std, cxr));
				res.add(MatrixCorrectionModel2.zaLabel(unk, std, cxr));
			}
			res.add(krl.as(Method.Calculated));
		}
		return res;
	};

	/**
	 * Builds a {@link UncertainValuesCalculator} around an instance of the
	 * {@link XPPMatrixCorrection2} algorithm that returns a basic set of default
	 * outputs.
	 *
	 * @param kratios
	 * @param estUnknown
	 * @return {@link UncertainValuesCalculator}
	 * @throws ArgumentException
	 */
	public static UncertainValuesCalculator<EPMALabel> buildFiniteDifference(final Set<KRatioLabel> kratios, //
			final Map<MassFraction, Double> estUnknown, //
			final double frac, //
			final boolean allOutputs) throws ArgumentException {
		Material unkMat = null;
		for (final MassFraction mf : estUnknown.keySet()) {
			if (unkMat == null)
				unkMat = mf.getMaterial();
			if (!unkMat.equals(mf.getMaterial()))
				throw new ArgumentException("Two different materials [" + mf.getMaterial() + " and " + unkMat
						+ " are present in the estimated unknown.");
		}
		final List<EPMALabel> outputs = allOutputs ? Collections.emptyList() : buildDefaultOutputs(kratios);
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(kratios, outputs);
		final UncertainValues<EPMALabel> inputs = xpp.buildInput(unkMat);
		xpp.addConstraints(xpp.buildConstraints(inputs));
		xpp.addAdditionalInputs(estUnknown);
		final UncertainValuesCalculator<EPMALabel> res = new UncertainValuesCalculator<EPMALabel>(xpp, inputs);
		res.setCalculator(res.new FiniteDifference(xpp.delta(res.getInputLabels(), inputs, frac)));
		return res;
	}

	/**
	 * Builds a {@link UncertainValuesCalculator} around an instance of the
	 * {@link XPPMatrixCorrection2} algorithm that returns a basic set of default
	 * outputs.
	 *
	 * @param kratios
	 * @param estUnknown
	 * @return {@link UncertainValuesCalculator}
	 * @throws ArgumentException
	 */
	public static UncertainValuesCalculator<EPMALabel> buildMonteCarlo(final Set<KRatioLabel> kratios, //
			final Map<MassFraction, Double> estUnknown, //
			final int nEvals, //
			final boolean allOutputs) throws ArgumentException {
		final UncertainValuesCalculator<EPMALabel> res = buildAnalytical(kratios, estUnknown, allOutputs);
		res.setCalculator(res.new MonteCarlo(nEvals));
		return res;
	}

	/**
	 * Computes the ionization exponent ('m') in the ionization cross section model
	 * Qla(U) = log(U)/(U<sup>m</sup> El<sup>2</sup>).
	 *
	 * @param sh   AtomicShell
	 * @param frac A scale for the relative uncertainty
	 * @return UncertainValue
	 */
	public static UncertainValue computeIonizationExponent(final AtomicShell sh, final double frac) {
		double m, dm;
		switch (sh.getFamily()) {
		case K:
			m = 0.86 + 0.12 * Math.exp(-Math.pow(sh.getElement().getAtomicNumber() / 5, 2.0));
			dm = m * 0.01;
			break;
		case L:
			m = 0.82;
			dm = 0.02 * m;
			break;
		case M:
		case N:
		default:
			m = 0.78;
			dm = 0.05 * m;
			break;
		}
		return new UncertainValue(m, frac * dm);
	}

	/**
	 * This method joins together the individual Composition computations into a
	 * single final step that outputs the final matrix corrections.
	 *
	 * @param unk
	 * @param stds
	 * @return List&lt;NamedMultivariateJacobianFunction&gt;
	 * @throws ArgumentException
	 */
	private static List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps(
			final Set<KRatioLabel> kratios, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		final Map<UnknownMatrixCorrectionDatum, CharacteristicXRaySet> unks = new HashMap<>();
		final Map<StandardMatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		final Set<Composition> allComps = new HashSet<>();
		for (final KRatioLabel krl : kratios) {
			assert krl.getMethod().equals(Method.Measured) : krl + " not measured.";
			final UnknownMatrixCorrectionDatum unk = krl.getUnknown();
			final StandardMatrixCorrectionDatum std = krl.getStandard();
			final ElementXRaySet exrs = krl.getXRaySet();
			if (!unks.containsKey(unk))
				unks.put(unk, new CharacteristicXRaySet());
			unks.get(unk).addAll(exrs);
			if (!stds.containsKey(std))
				stds.put(std, new CharacteristicXRaySet());
			stds.get(std).addAll(exrs);
			allComps.add(std.getComposition());
			if (std.hasCoating())
				allComps.add(std.getCoating().getComposition());
			if (unk.hasCoating())
				allComps.add(unk.getCoating().getComposition());
		}

		// Build the LMJF list in reverse to determine the input required by later
		// functions first
		final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
		final FastIndex<EPMALabel> reqInputs = new FastIndex<>(outputLabels);
		// Need to do it this way to ensure that certain items aren't double
		// calculated...
		final MultiE0MultiLineModel multiE0 = new MultiE0MultiLineModel(kratios);
		res.add(multiE0);
		reqInputs.addMissing(multiE0.getInputLabels());
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay())
					step.add(new StepZA(krl.getUnknown(), krl.getStandard(), cxr));
			final ExplicitMeasurementModel<EPMALabel, EPMALabel> za = ParallelMeasurementModelBuilder.join("ZA", step);
			reqInputs.addMissing(za.getInputLabels());
			res.add(za);
		}
		{
			if (outputLabels.size() == 0)
				reqInputs.clear();
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> step = new ArrayList<>();
			final CharacteristicXRaySet cxrs = new CharacteristicXRaySet();
			for (final Map.Entry<StandardMatrixCorrectionDatum, CharacteristicXRaySet> me : stds.entrySet()) {
				step.add(new StepXPP(me.getKey(), me.getValue(), reqInputs));
				cxrs.addAll(me.getValue());
			}
			for (final Map.Entry<UnknownMatrixCorrectionDatum, CharacteristicXRaySet> me : unks.entrySet())
				step.add(new StepXPP(me.getKey(), me.getValue(), reqInputs));
			res.add(ParallelMeasurementModelBuilder.join("XPP", step));
		}
		for (final UnknownMatrixCorrectionDatum unk : unks.keySet()) {
			res.add(MaterialMACFunction.buildMaterialMACFunctions(unk.getMaterial(), kratios));
			if (unk.hasCoating())
				res.add(MaterialMACFunction.buildMaterialMACFunctions(unk.getCoating().getMaterial(), kratios));
		}
		final ParallelMeasurementModelBuilder<MaterialLabel, MaterialLabel> builder = //
				new ParallelMeasurementModelBuilder<>("Compositions");
		for (final Composition comp : allComps)
			builder.add(comp.getFunction());
		res.add(builder.build());
		Collections.reverse(res);
		return res;
	}

	/**
	 * <p>
	 * Do some very basic checks on the input and output indices.
	 * <ul>
	 * <li>Not duplicated..</li>
	 * <li>Not missing</li>
	 * </ul>
	 * </p>
	 *
	 * @param idxs
	 */
	static private void checkIndices(final int... idxs) {
		final int len = idxs.length;
		for (int i = 0; i < len; ++i) {
			assert idxs[i] >= 0 : "checkIndices[" + i + "] = " + idxs[i];

			for (int j = i + 1; j < len; ++j)
				assert idxs[i] != idxs[j] : "idx[" + i + "] == idx[" + j + "]";
		}
	}

	private static final double eVtokeV(final double eV) {
		return 0.001 * eV;
	}

	/**
	 * @param kratios      Set of KRatioLabel objects.
	 * @param outputLabels
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2(final Set<KRatioLabel> kratios, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		super("XPP Matrix Correction", kratios, buildSteps(kratios, outputLabels), outputLabels);
	}

	public Map<EPMALabel, Constraint> buildConstraints(final UncertainValuesBase<EPMALabel> inputs) {
		final Map<EPMALabel, Constraint> res = new HashMap<EPMALabel, Constraint>();
		// Add constraints for the standards
		for (final MassFraction mf : inputs.getLabels(MassFraction.class))
			res.put(mf, new Constraint.Range(0.0, 1.0e6));
		for (final AtomFraction mf : inputs.getLabels(AtomFraction.class))
			res.put(mf, new Constraint.Range(0.0, 1.0e6));
		// Add constraints for the unknown(s)
		final Set<Material> unkMats = new HashSet<>();
		for (final KRatioLabel krl : inputs.getLabels(KRatioLabel.class))
			unkMats.add(krl.getUnknown().getMaterial());
		for (final Material unkMat : unkMats)
			for (final Element elm : unkMat.getElementSet())
				res.put(MaterialLabel.buildMassFractionTag(unkMat, elm), new Constraint.Range(0.0, 1.0e6));
		// Add elemental MAC constraints for the MC model
		for (final ElementalMAC.ElementMAC em : inputs.getLabels(ElementalMAC.ElementMAC.class))
			res.put(em, new Constraint.Range(0.1 * inputs.getEntry(em), 10.0 * inputs.getEntry(em)));
		// Add elemental MAC constraints for the MC model
		for (final CoatingThickness th : inputs.getLabels(CoatingThickness.class)) {
			final double thickness = Math.abs(inputs.getEntry(th));
			res.put(th, new Constraint.Range(0.1 * thickness, 10.0 * thickness));
		}
		return res;
	}

	/**
	 * This is a quick way to build many of the input parameters. Many of the input
	 * parameters are computed or tabulated. Some are input as experimental
	 * conditions like beam energy which are provided in the {@link KRatioLabel}
	 * objects. Using this information, it is possible to calculate all the
	 * necessary inputs to this {@link ExplicitMeasurementModel}.
	 *
	 * @param estUnknown
	 * @return {@link UncertainValues}&lt;EPMALabel&gt;
	 * @throws ArgumentException
	 */
	@Override
	public UncertainValues<EPMALabel> buildInput(final Material unkMat) throws ArgumentException {
		final List<UncertainValuesBase<? extends EPMALabel>> results = new ArrayList<>();
		// Make sure that there are no replicated Compositions
		results.add(unkMat.getAtomicWeights());
		{
			final Set<Composition> allComps = new HashSet<>();
			for (final KRatioLabel krl : mKRatios) {
				allComps.add(krl.getStandard().getComposition());
				for (final MatrixCorrectionDatum mcd : Arrays.asList(krl.getStandard(), krl.getUnknown())) {
					if (mcd.hasCoating())
						allComps.add(mcd.getCoating().getComposition());
				}
			}
			// Add the mass fraction and atomic weight tags
			for (final Composition comp : allComps)
				results.add(comp.getInputs());
		}
		results.add(buildMeanIonizationPotentials());
		results.add(buildElementalMACs());
		{
			final Map<EPMALabel, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				{
					final UnknownMatrixCorrectionDatum mcd = krl.getUnknown();
					vals.put(MatrixCorrectionModel2.beamEnergyLabel(mcd), mcd.getBeamEnergy());
				}
				{
					final StandardMatrixCorrectionDatum mcd = krl.getStandard();
					vals.put(MatrixCorrectionModel2.beamEnergyLabel(mcd), mcd.getBeamEnergy());
				}
			}
			results.add(new UncertainValues<>(vals));
		}
		final Set<Layer> coatings = new HashSet<>();
		for (final KRatioLabel krl : mKRatios) {
			if (krl.getStandard().hasCoating())
				coatings.add(krl.getStandard().getCoating());
			if (krl.getUnknown().hasCoating())
				coatings.add(krl.getUnknown().getCoating());
		}
		if (coatings.size() > 0) {
			final Map<EPMALabel, Number> vals = new HashMap<>();
			for (final Layer coating : coatings)
				vals.put(MatrixCorrectionModel2.coatingMassThickness(coating), coating.getMassThickness());
			results.add(new UncertainValues<>(vals));
		}
		{
			final Map<EPMALabel, Number> vals = new HashMap<>();
			final Set<AtomicShell> allShells = new HashSet<>();
			for (final KRatioLabel krl : mKRatios)
				for (final XRay cxr : krl.getXRaySet())
					allShells.add(((CharacteristicXRay) cxr).getInner());
			for (final AtomicShell sh : allShells)
				vals.put(new IonizationExponentLabel(sh), computeIonizationExponent(sh, 1.0));
			assert !vals.isEmpty();
			results.add(new UncertainValues<>(vals));
		}
		{
			final Map<EPMALabel, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getUnknown()),
						krl.getUnknown().getTakeOffAngle());
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getStandard()),
						krl.getStandard().getTakeOffAngle());
			}
			results.add(new UncertainValues<>(vals));
		}
		{
			final Map<EPMALabel, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getUnknown()),
						new UncertainValue(0.0, krl.getUnknown().getRoughness()));
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getStandard()),
						new UncertainValue(0.0, krl.getStandard().getRoughness()));
			}
			if (!vals.isEmpty())
				results.add(new UncertainValues<>(vals));
		}
		final CharacteristicXRaySet allCxr = new CharacteristicXRaySet();
		for (final KRatioLabel krl : mKRatios)
			allCxr.addAll(krl.getXRaySet());
		final int sz = allCxr.size();
		if (sz > 0) {
			final Map<EPMALabel, Number> weightT = new HashMap<>();
			for (final CharacteristicXRay cxr : allCxr.getSetOfCharacteristicXRay()) {
				final XRayEmissionProbability xrep = new XRayEmissionProbability(cxr.getInner());
				weightT.put(new MatrixCorrectionModel2.XRayWeightLabel(cxr), xrep.getWeightU(cxr));
			}
			results.add(new UncertainValues<>(weightT));
		}
		final Map<EPMALabel, Number> mon = new HashMap<>();
		for (final KRatioLabel krl : mKRatios) {
			for (final MatrixCorrectionDatum mcd : Arrays.asList(krl.getStandard(), krl.getUnknown()))
				for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
					final EPMALabel sfLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mcd,
							cxr.getInner());
					if (!mon.containsKey(sfLbl))
						mon.put(sfLbl, new UncertainValue(1.0, 0.01));
				}
		}
		results.add(new UncertainValues<>(mon));
		return UncertainValues.<EPMALabel>asUncertainValues(UncertainValuesBase.combine(results, true));
	}

	/**
	 * Mean ionization potential
	 *
	 * @param elm
	 * @return {@link UncertainValue}
	 */
	public UncertainValue computeJi(final Element elm) {
		final double z = elm.getAtomicNumber();
		final double j = 1.0e-3 * z * (10.04 + 8.25 * Math.exp(-z / 11.22));
		return new UncertainValue(j, 0.03 * j);
	}

	public DataFrame<Double, String> computePhiRhoZCurve(final Map<EPMALabel, Double> outputs, //
			final double rhoZmax, //
			final double dRhoZ, //
			final double minWeight //
	) {
		Range<Double> rz = Range.of(0.0, rhoZmax, Math.max(rhoZmax / 1000.0, dRhoZ));
		DataFrame<Double, String> res = DataFrame.of(rz, String.class, columns -> {
		});
		for (final KRatioLabel krl : mKRatios) {
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				if (cxr.getWeight() > minWeight) {
					final StandardMatrixCorrectionDatum std = krl.getStandard();
					final double aStd = outputs.get(MatrixCorrectionModel2.shellLabel("a", std, cxr.getInner()))
							.doubleValue();
					final double bStd = outputs.get(MatrixCorrectionModel2.shellLabel("b", std, cxr.getInner()))
							.doubleValue();
					final double AStd = outputs.get(MatrixCorrectionModel2.shellLabel("A", std, cxr.getInner()))
							.doubleValue();
					final double BStd = outputs.get(MatrixCorrectionModel2.shellLabel("B", std, cxr.getInner()))
							.doubleValue();
					final double chiStd = outputs.get(new MatrixCorrectionModel2.ChiLabel(std, cxr)).doubleValue();
					final double phi0Std = outputs.get(new MatrixCorrectionModel2.Phi0Label(std, cxr.getInner()))
							.doubleValue();
					final UnknownMatrixCorrectionDatum unk = krl.getUnknown();
					final double aUnk = outputs.get(MatrixCorrectionModel2.shellLabel("a", unk, cxr.getInner()))
							.doubleValue();
					final double bUnk = outputs.get(MatrixCorrectionModel2.shellLabel("b", unk, cxr.getInner()))
							.doubleValue();
					final double AUnk = outputs.get(MatrixCorrectionModel2.shellLabel("A", unk, cxr.getInner()))
							.doubleValue();
					final double BUnk = outputs.get(MatrixCorrectionModel2.shellLabel("B", unk, cxr.getInner()))
							.doubleValue();
					final double chiUnk = outputs.get(new MatrixCorrectionModel2.ChiLabel(unk, cxr)).doubleValue();
					final double phi0Unk = outputs.get(new MatrixCorrectionModel2.Phi0Label(unk, cxr.getInner()))
							.doubleValue();
					res.cols().add(std.toString() + "[" + cxr.toString() + ", gen]",
							rz.map(rhoZ -> phiRhoZ(rhoZ, aStd, bStd, AStd, BStd, phi0Std)));
					res.cols().add(std.toString() + "[" + cxr.toString() + ", emit]",
							rz.map(rhoZ -> phiRhoZ(rhoZ, aStd, bStd, AStd, BStd, phi0Std) * Math.exp(-chiStd * rhoZ)));
					res.cols().add(unk.toString() + "[" + cxr.toString() + ", gen]",
							rz.map(rhoZ -> phiRhoZ(rhoZ, aUnk, bUnk, AUnk, BUnk, phi0Unk)));
					res.cols().add(unk.toString() + "[" + cxr.toString() + ", emit]",
							rz.map(rhoZ -> phiRhoZ(rhoZ, aUnk, bUnk, AUnk, BUnk, phi0Unk) * Math.exp(-chiUnk * rhoZ)));
				}
			}
		}
		return res;
	}

	/**
	 * Useful for calculating the
	 *
	 * @param inputs
	 * @param frac
	 * @return
	 */
	public RealVector delta(final List<? extends EPMALabel> list, //
			final UncertainValuesBase<EPMALabel> inputs, //
			final double frac) {
		final RealVector res = inputs.getValues(list).copy();
		for (final EPMALabel lbl : inputs.getLabels(MatrixCorrectionModel2.RoughnessLabel.class))
			res.setEntry(list.indexOf(lbl), 1.0e-7);
		for (final EPMALabel lbl : inputs.getLabels(MatrixCorrectionModel2.CoatingThickness.class))
			res.setEntry(list.indexOf(lbl), 1.0e-7);
		return res.mapMultiply(frac);
	}

	public Set<Element> getElements() {
		final Set<Element> elms = new TreeSet<>();
		for (final KRatioLabel krl : mKRatios)
			elms.add(krl.getElement());
		return elms;
	}

	@Override
	public String toHTML(final Mode mode) {
		final Report r = new Report("XPP");
		r.addHeader("XPP Matrix Correction");
		final Table tbl = new Table();
		tbl.addRow(Table.th("Unknown"), Table.th("Standard"), Table.th("Transitions"));
		for (final KRatioLabel krl : mKRatios)
			tbl.addRow(Table.td(krl.getUnknown()), Table.td(krl.getStandard()), Table.td(krl.getXRaySet()));
		r.add(tbl);
		r.addSubHeader("Inputs / Outputs");
		r.addHTML(super.toHTML(mode));
		return r.toHTML(mode);
	}

	@Override
	public String toString() {
		return "XPPMatrixCorrection2";
	}

	private UncertainValues<EPMALabel> buildElementalMACs() throws ArgumentException {
		final Set<Element> elms = new HashSet<>();
		final Set<CharacteristicXRay> scxr = new HashSet<>();
		{
			for (final KRatioLabel krl : mKRatios) {
				for (final MatrixCorrectionDatum mcd : Arrays.asList(krl.getStandard(), krl.getUnknown())) {
					elms.addAll(mcd.getMaterial().getElementSet());
					if (mcd.hasCoating())
						elms.addAll(mcd.getCoating().getMaterial().getElementSet());
				}
				scxr.addAll(krl.getXRaySet().getSetOfCharacteristicXRay());
			}
		}
		final Map<EPMALabel, Number> macInps = new HashMap<>();
		final ElementalMAC em = new ElementalMAC();
		for (final Element elm : elms)
			for (final CharacteristicXRay cxr : scxr) {
				final ElementalMAC.ElementMAC label = new ElementalMAC.ElementMAC(elm, cxr);
				final UncertainValue mac = em.compute(elm, cxr);
				macInps.put(label, mac);
			}
		return new UncertainValues<EPMALabel>(macInps);
	}

	private UncertainValues<EPMALabel> buildMeanIonizationPotentials() {
		final Set<Element> elms = new HashSet<>();
		for (final KRatioLabel krl : mKRatios) {
			elms.addAll(krl.getStandard().getMaterial().getElementSet());
			elms.addAll(krl.getUnknown().getMaterial().getElementSet());
		}
		final RealVector vals = new ArrayRealVector(elms.size());
		final RealVector var = new ArrayRealVector(vals.getDimension());
		final List<EPMALabel> tags = new ArrayList<>();
		int i = 0;
		for (final Element elm : elms) {
			tags.add(MatrixCorrectionModel2.meanIonizationLabel(elm));
			final UncertainValue j = computeJi(elm);
			vals.setEntry(i, j.doubleValue());
			var.setEntry(i, j.variance());
			++i;
		}
		assert tags.size() == vals.getDimension();
		final UncertainValues<EPMALabel> mip = new UncertainValues<EPMALabel>(tags, vals, var);
		return mip;
	}

	private double phiRhoZ(final double rhoZ, final double a, final double b, final double A, final double B,
			final double phi0) {
		return A * Math.exp(-a * rhoZ) + (B * rhoZ + phi0 - A) * Math.exp(-b * rhoZ);
	}
}
