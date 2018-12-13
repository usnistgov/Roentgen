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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunctionBuilder;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction.MaterialMAC;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import joinery.DataFrame;

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
 * Since it is derived from {@link LabeledMultivariateJacobianFunction}, it
 * computes not only the value (ie. the k-ratio) but also sensitivity matrix
 * (Jacobian) that maps uncertainty in the input parameters into uncertainty in
 * the output parameters.
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
		implements ILabeledMultivariateFunction, IToHTML {

	// The types of variables to compute Jacobian elements.
	private final Set<MatrixCorrectionModel2.Variates> mVariates = new HashSet<>();

	// The material MACs computed in getMaterialMACs or updated by setMaterialMACs()
	private UncertainValues mMaterialMACS = null;

	public static final String EPS = "&epsilon;";
	public static final String ONE_OVER_S = "<sup>1</sup>/<sub>S</sub>";
	public static final String QLA = "<html>Q<sub>l</sub><sup>a</sup>";
	public static final String ZBARB = "Z<sub>barb</sub>";
	public static final String RBAR = "R<sub>bar</sub>";
	public static final String E_0 = "E<sub>0</sub>";

	private static final double eVtokeV(final double eV) {
		return 0.001 * eV;
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
	static void checkIndices(final int... idxs) {
		final int len = idxs.length;
		for (int i = 0; i < len; ++i) {
			assert idxs[i] >= 0;
			for (int j = i + 1; j < len; ++j)
				assert idxs[i] != idxs[j];
		}
	}

	private static class StepMJZBarb // C1
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final Composition mComposition;

		public StepMJZBarb(final Composition comp, final Set<MatrixCorrectionModel2.Variates> variates,
				final boolean isStandard) {
			super(buildInputs(comp, variates, isStandard), buildOutputs(comp));
			mComposition = comp;
		}

		public static List<? extends Object> buildOutputs(final Composition comp) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel("M", comp));
			res.add(new MatrixCorrectionModel2.CompositionLabel("J", comp));
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, comp));
			return res;
		}

		public static List<? extends Object> buildInputs(final Composition comp,
				final Set<MatrixCorrectionModel2.Variates> variates, final boolean isStandard) {
			final List<Object> res = new ArrayList<>();
			final MatrixCorrectionModel2.Variates varType = isStandard
					? MatrixCorrectionModel2.Variates.StandardComposition
					: MatrixCorrectionModel2.Variates.UnknownComposition;
			for (final Element elm : comp.getElementSet()) {
				if (variates.contains(varType))
					res.add(Composition.buildMassFractionTag(comp, elm));
				if (variates.contains(MatrixCorrectionModel2.Variates.MeanIonizationPotential))
					res.add(MatrixCorrectionModel2.meanIonizationLabel(elm));
			}
			return res;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Composition comp = mComposition;
			final List<Element> elms = new ArrayList<>(comp.getElementSet());
			final Object[] tagCi = new Object[elms.size()];
			final Object[] tagJi = new Object[elms.size()];
			final double[] Ci = new double[elms.size()];
			final double[] ZoA = new double[elms.size()];
			final double[] Z = new double[elms.size()];
			final double[] Ji = new double[elms.size()];
			for (int i = 0; i < Ci.length; ++i) {
				final Element elm = elms.get(i);
				tagCi[i] = Composition.buildMassFractionTag(comp, elm);
				tagJi[i] = MatrixCorrectionModel2.meanIonizationLabel(elm);
				Ci[i] = getValue(tagCi[i], point);
				Z[i] = elm.getAtomicNumber();
				Ji[i] = getValue(tagJi[i], point);
				ZoA[i] = Z[i] / elm.getAtomicWeight();
			}

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oM = outputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mComposition));
			final int oJ = outputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mComposition));
			final int oZbarb = outputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mComposition));
			checkIndices(oM, oJ, oZbarb);

			// Calculate M & partials
			double M = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				M += Ci[i] * ZoA[i];
			rv.setEntry(oM, M); // c1
			for (int i = 0; i < tagCi.length; ++i)
				writeJacobian(oM, tagCi[i], ZoA[i], rm); // c1

			// Calculate J and partials
			double lnJ = 0.0;
			for (int i = 0; i < elms.size(); ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * ZoA[i];
			lnJ /= M;
			final double J = Math.exp(lnJ); // keV
			rv.setEntry(oJ, J); // C2
			for (int i = 0; i < tagCi.length; ++i) {
				writeJacobian(oJ, tagJi[i], Ci[i] * J * ZoA[i] / (M * Ji[i]), rm); // C3
				double k1 = 0.0, k2 = 0.0;
				for (int j = 0; j < tagCi.length; ++j) {
					if (j != i) {
						k1 += Ci[j] * ZoA[j];
						k2 -= Ci[j] * ZoA[j] * Math.log(Ji[j]);
					}
				}
				writeJacobian(oJ, tagCi[i], J * ZoA[i] * (k1 * Math.log(Ji[i]) + k2) / (M * M), rm); // C3
			}
			// Calculate Zbarb and partials
			double Zbt = 0.0;
			for (int i = 0; i < elms.size(); ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // C2
			final double Zbarb = Math.pow(Zbt, 2.0);
			rv.setEntry(oZbarb, Zbarb);
			for (int i = 0; i < elms.size(); ++i)
				writeJacobian(oZbarb, tagCi[i], 2.0 * Math.sqrt(Z[i]) * Zbt, rm); // Ci
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final Composition comp = mComposition;
			final List<Element> elms = new ArrayList<>(comp.getElementSet());
			final double[] Ci = new double[elms.size()];
			final double[] ZoA = new double[elms.size()];
			final double[] Z = new double[elms.size()];
			final double[] Ji = new double[elms.size()];
			for (int i = 0; i < Ci.length; ++i) {
				final Element elm = elms.get(i);
				Ci[i] = getValue(Composition.buildMassFractionTag(comp, elm), point);
				Z[i] = elm.getAtomicNumber();
				Ji[i] = getValue(MatrixCorrectionModel2.meanIonizationLabel(elm), point);
				ZoA[i] = Z[i] / elm.getAtomicWeight();
			}

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oM = outputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mComposition));
			final int oJ = outputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mComposition));
			final int oZbarb = outputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mComposition));
			checkIndices(oM, oJ, oZbarb);

			// Calculate M & partials
			double M = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				M += Ci[i] * ZoA[i];
			rv.setEntry(oM, M); // c1

			// Calculate J and partials
			double lnJ = 0.0;
			for (int i = 0; i < elms.size(); ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * ZoA[i];
			lnJ /= M;
			final double J = Math.exp(lnJ); // keV
			rv.setEntry(oJ, J); // C2
			// Calculate Zbarb and partials
			double Zbt = 0.0;
			for (int i = 0; i < elms.size(); ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // C2
			final double Zbarb = Math.pow(Zbt, 2.0);
			rv.setEntry(oZbarb, Zbarb);
			return rv;
		}

		@Override
		public String toString() {
			return "MJZBarb[" + mComposition.toString() + "]";
		}

	}

	private static class StepQlaE0OneOverS // C1
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel("M", datum.getComposition()));
			res.add(new MatrixCorrectionModel2.CompositionLabel("J", datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variates.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			return res;
		}

		public StepQlaE0OneOverS(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(datum, variates), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final int iJ = inputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mDatum.getComposition()));
			final int iM = inputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mDatum.getComposition()));
			checkIndices(iJ, iM);

			final double e0 = getValue(tagE0, point);
			final double J = point.getEntry(iJ);
			final double M = point.getEntry(iM);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			assert u0 > 1.0 : e0 + " over " + Ea;

			double m = 1.0;
			final double Za = mShell.getElement().getAtomicNumber();
			switch (mShell.getFamily().getPrincipleQuantumNumber()) {
			case 1: // K
				m = 0.86 + 0.12 * Math.exp(-1.0 * Math.pow(Za / 5, 2.0));
				break;
			case 2: // L
				m = 0.82;
				break;
			case 3: // M
				m = 0.78;
				break;
			}

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final double logU0 = Math.log(u0);
			// QlaE0 and partials
			final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
			rv.setEntry(oQlaE0, QlaE0); // C1
			writeJacobian(oQlaE0, tagE0, Math.pow(u0, -1.0 - m) * (1.0 - m * logU0) / Math.pow(Ea, 3.0), rm); // C1

			// oneOverS and partials
			final double D1 = 6.6e-6, D2 = 1.12e-5 * (1.35 - 0.45 * J), D3 = 2.2e-6 / J;
			final double P1 = 0.78, P2 = 0.1, P3 = 0.25 * J - 0.5;
			final double T1 = 1.0 + P1 - m, T2 = 1.0 + P2 - m, T3 = 1.0 + P3 - m;
			final double dD2dJ = 1.12e-5 * (-0.45), dD3dJ = -1.0 * 2.2e-6 / (J * J), dP3dJ = 0.25, dT3dJ = 0.25;
			final double v0ou0 = Ea / J; // c1

			final double T1_2 = Math.pow(T1, 2), T2_2 = Math.pow(T2, 2), T3_2 = Math.pow(T3, 2), T3_3 = Math.pow(T3, 3);
			final double v0ou0_P1 = Math.pow(v0ou0, P1), v0ou0_P2 = Math.pow(v0ou0, P2), v0ou0_P3 = Math.pow(v0ou0, P3);
			final double u0_T1 = Math.pow(u0, T1), u0_T2 = Math.pow(u0, T2), u0_T3 = Math.pow(u0, T3);
			final double logV0oU0 = Math.log(v0ou0);
			final double dOoSdJ = (T3
					* (T3_2 * (-(D1 * v0ou0_P1 * (-1.0 + P1) * T2_2 * (1 - u0_T1 + u0_T1 * T1 * logU0))
							- v0ou0_P2 * (-1.0 + P2) * T1_2 * D2 * (1 - u0_T2 + u0_T2 * T2 * logU0)
							+ v0ou0_P2 * J * T1_2 * (1 - u0_T2 + u0_T2 * T2 * logU0) * dD2dJ) - (-1.0 + u0_T3)
									* v0ou0_P3 * J * T1_2 * T2_2 * dD3dJ
							+ u0_T3 * v0ou0_P3 * J * T1_2 * T2_2 * logU0 * T3 * dD3dJ)
					- v0ou0_P3 * T1_2 * T2_2 * D3
							* (-2.0 * (-1.0 + u0_T3) * J * dT3dJ
									- u0_T3 * logU0 * T3_2 * (1 - P3 + J * logV0oU0 * dP3dJ + J * logU0 * dT3dJ)
									+ T3 * (-1.0 + u0_T3 - (-1.0 + u0_T3) * P3 + (-1.0 + u0_T3) * J * logV0oU0 * dP3dJ
											+ 2.0 * u0_T3 * J * logU0 * dT3dJ)))
					/ (Ea * M * T1_2 * T2_2 * T3_3);

			final double dOoSdE0 = (J * (D1 * u0_T1 * v0ou0_P1 + u0_T2 * v0ou0_P2 * D2 + u0_T3 * v0ou0_P3 * D3) * logU0)
					/ (e0 * Ea * M);

			final double dOoSdM = -((J * ((D1 * v0ou0_P1 * (1 - u0_T1 + u0_T1 * T1 * logU0)) / T1_2
					+ (v0ou0_P2 * D2 * (1.0 - u0_T2 + u0_T2 * T2 * logU0)) / T2_2
					+ (v0ou0_P3 * D3 * (1.0 - u0_T3 + u0_T3 * logU0 * T3)) / T3_2)) / (Ea * Math.pow(M, 2)));

			final double oneOverS = ((J / (Ea * M)) //
					* ((D1 * v0ou0_P1 * (1.0 - u0_T1 + T1 * u0_T1 * logU0)) / T1_2
							+ (D2 * v0ou0_P2 * (1.0 - u0_T2 + T2 * u0_T2 * logU0)) / T2_2 //
							+ (D3 * v0ou0_P3 * (1.0 - u0_T3 + T3 * u0_T3 * logU0)) / T3_2));
			rv.setEntry(oOneOverS, oneOverS); // C2
			rm.setEntry(oOneOverS, iM, dOoSdM); // C2
			rm.setEntry(oOneOverS, iJ, dOoSdJ); // C2
			writeJacobian(oOneOverS, tagE0, dOoSdE0, rm); // C2

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iJ = inputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mDatum.getComposition()));
			final int iM = inputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mDatum.getComposition()));
			checkIndices(iJ, iM);

			final double e0 = getValue(MatrixCorrectionModel2.beamEnergyLabel(mDatum), point);
			final double J = point.getEntry(iJ);
			final double M = point.getEntry(iM);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			assert u0 > 1.0;

			double m = 1.0;
			final double Za = mShell.getElement().getAtomicNumber();
			switch (mShell.getFamily().getPrincipleQuantumNumber()) {
			case 1: // K
				m = 0.86 + 0.12 * Math.exp(-1.0 * Math.pow(Za / 5, 2.0));
				break;
			case 2: // L
				m = 0.82;
				break;
			case 3: // M
				m = 0.78;
				break;
			}

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final double logU0 = Math.log(u0);
			// QlaE0 and partials
			final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
			rv.setEntry(oQlaE0, QlaE0); // C1

			// oneOverS and partials
			final double D1 = 6.6e-6, D2 = 1.12e-5 * (1.35 - 0.45 * J), D3 = 2.2e-6 / J;
			final double P1 = 0.78, P2 = 0.1, P3 = 0.25 * J - 0.5;
			final double T1 = 1.0 + P1 - m, T2 = 1.0 + P2 - m, T3 = 1.0 + P3 - m;
			final double v0ou0 = Ea / J; // c1

			final double T1_2 = Math.pow(T1, 2), T2_2 = Math.pow(T2, 2), T3_2 = Math.pow(T3, 2);
			final double v0ou0_P1 = Math.pow(v0ou0, P1), v0ou0_P2 = Math.pow(v0ou0, P2), v0ou0_P3 = Math.pow(v0ou0, P3);
			final double u0_T1 = Math.pow(u0, T1), u0_T2 = Math.pow(u0, T2), u0_T3 = Math.pow(u0, T3);

			final double oneOverS = ((J / (Ea * M)) //
					* ((D1 * v0ou0_P1 * (1.0 - u0_T1 + T1 * u0_T1 * logU0)) / T1_2
							+ (D2 * v0ou0_P2 * (1.0 - u0_T2 + T2 * u0_T2 * logU0)) / T2_2 //
							+ (D3 * v0ou0_P3 * (1.0 - u0_T3 + T3 * u0_T3 * logU0)) / T3_2));
			rv.setEntry(oOneOverS, oneOverS); // C2

			return rv;
		}

		@Override
		public String toString() {
			return "StepQlaE0OneOverS[" + mDatum + "," + mShell + "]";
		}
	}

	private static class StepRphi0 extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variates.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			return res;
		}

		public StepRphi0(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(datum, variates), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			checkIndices(iZbarb);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double e0 = getValue(tagE0, point);
			final double u0 = e0 / Ea;
			final double Zbarb = point.getEntry(iZbarb);

			final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3)));
			final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55);
			final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar);
			final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0);
			final double Gu0 = (u0 - 1.0 - (1.0 - 1.0 / Math.pow(u0, 1.0 + q)) / (1.0 + q)) / ((2.0 + q) * Ju0);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			final double R = 1.0 - etabar * Wbar * (1.0 - Gu0);
			rv.setEntry(oR, R);

			final double detabardZbarb = 0.00175
					+ (0.007215 * Math.pow(Zbarb, 0.3)) / Math.exp(0.015 * Math.pow(Zbarb, 1.3));

			final double dWbardZbarb = (0.27027027027027023 + 4.55 * Math.pow(etabar, 3.55)) * detabardZbarb;

			final double dqdZbarb = dWbardZbarb / Math.pow(-1.0 + Wbar, 2.0);

			final double dGu0dZbarb = -((Math.pow(u0, -1.0 - q)
					* (3.0 - 4.0 * Math.pow(u0, 1.0 + q) + Math.pow(u0, 2.0 + q) + 2.0 * Math.log(u0)
							+ (2.0 - 4.0 * Math.pow(u0, 1.0 + q) + 2.0 * Math.pow(u0, 2.0 + q) + 3.0 * Math.log(u0)) * q
							+ ((-1.0 + u0) * Math.pow(u0, 1.0 + q) + Math.log(u0)) * Math.pow(q, 2.0))
					* dqdZbarb) / (Ju0 * Math.pow(1.0 + q, 2.0) * Math.pow(2.0 + q, 2.0)));

			final double dRdZbarb = -((1 - Gu0) * Wbar * detabardZbarb) + etabar * Wbar * dGu0dZbarb
					- etabar * (1.0 - Gu0) * dWbardZbarb;

			final double dJu0de0 = Math.log(u0) / Ea;

			final double dGu0de0 = ((1.0 / Ea - Ea / (Math.pow(e0, 2.0) * Math.pow(u0, q))) * Ju0
					+ (1.0 - e0 / Ea + 1.0 / (1.0 + q) - Ea / (Math.pow(u0, q) * (e0 + e0 * q))) * dJu0de0)
					/ ((2.0 + q) * Math.pow(Ju0, 2.0));

			final double dRde0 = etabar * Wbar * dGu0de0;

			rm.setEntry(oR, iZbarb, dRdZbarb);
			writeJacobian(oR, tagE0, dRde0, rm);

			final double r = 2.0 - 2.3 * etabar;
			final double phi0 = 1.0 + 3.3 * (1.0 - 1.0 / Math.pow(u0, r)) * Math.pow(etabar, 1.2);
			rv.setEntry(oPhi0, phi0);
			rm.setEntry(oPhi0, iZbarb,
					detabardZbarb * (Math.pow(etabar, 0.2) * (3.96 - 3.96 * Math.pow(u0, -2 + 2.3 * etabar))
							- 7.59 * Math.pow(etabar, 1.2) * Math.pow(u0, -2 + 2.3 * etabar) * Math.log(u0)));
			writeJacobian(oPhi0, tagE0, (3.3 * Math.pow(etabar, 1.2) * r * Math.pow(u0, -1.0 - r)) / Ea, rm);
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			checkIndices(iZbarb);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double e0 = getValue(MatrixCorrectionModel2.beamEnergyLabel(mDatum), point);
			final double u0 = e0 / Ea;
			final double Zbarb = point.getEntry(iZbarb);

			final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3)));
			final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55);
			final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar);
			final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0);
			final double Gu0 = (u0 - 1.0 - (1.0 - 1.0 / Math.pow(u0, 1.0 + q)) / (1.0 + q)) / ((2.0 + q) * Ju0);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			final double R = 1.0 - etabar * Wbar * (1.0 - Gu0);
			rv.setEntry(oR, R);

			final double r = 2.0 - 2.3 * etabar;
			final double phi0 = 1.0 + 3.3 * (1.0 - 1.0 / Math.pow(u0, r)) * Math.pow(etabar, 1.2);
			rv.setEntry(oPhi0, phi0);
			return rv;
		}

		@Override
		public String toString() {
			return "StepRPhi0[" + mDatum + "," + mShell + "]";
		}

	}

	private static class StepFRbar extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variates.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			return res;
		}

		public StepFRbar(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(datum, shell, variates), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			final int iOneOverS = inputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int iQlaE0 = inputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iR = inputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			checkIndices(iZbarb, iOneOverS, iQlaE0, iPhi0, iR);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = getValue(tagE0, point) / Ea;
			final double Zbarb = point.getEntry(iZbarb);
			final double oneOverS = point.getEntry(iOneOverS);
			final double QlaE0 = point.getEntry(iQlaE0);
			final double phi0 = point.getEntry(iPhi0);
			final double R = point.getEntry(iR);

			final double logZbarb = Math.log(Zbarb);
			final double X = 1.0 + 1.3 * logZbarb;
			final double Y = 0.2 + 0.005 * Zbarb;
			final double pu42 = Math.pow(u0, 0.42);
			final double FoRbar = 1.0 + (X * Math.log(1.0 + Y * (1.0 - 1.0 / pu42)) / Math.log(1.0 + Y));

			final double Zbarb240 = 240. + Zbarb;
			final double kk = Math.log(Zbarb240 / 200.0);
			final double dFoRbardZbarb = (((-200.0 + 200.0 * pu42 + (-260.0 + 260.0 * pu42) * logZbarb) * kk)
					/ (-40.0 - Zbarb + pu42 * Zbarb240)
					+ Math.log(1.2 + (-0.2 - Zbarb / 200.0) / pu42 + Zbarb / 200.0)
							* ((-200.0 - 260.0 * logZbarb) / Zbarb240 + (260.0 * kk) / Zbarb))
					/ (200.0 * kk * kk);
			final double dFoRbardu0 = (0.546 * (40.0 + Zbarb) * (0.7692307692307692 + logZbarb))
					/ (u0 * (-40.0 - Zbarb + pu42 * Zbarb240) * kk);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oRbar = outputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int oF = outputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			checkIndices(oRbar, oF);

			// F and partials
			final double F = R * oneOverS / QlaE0; // C1
			rv.setEntry(oF, F);
			final double dFdR = oneOverS / QlaE0; // C1
			rm.setEntry(oF, iR, dFdR);
			final double dFdOneOverS = R / QlaE0; // C1
			rm.setEntry(oF, iOneOverS, dFdOneOverS);
			final double dFdQlaE0 = -1.0 * R * oneOverS / Math.pow(QlaE0, 2.0); // C1
			rm.setEntry(oF, iQlaE0, dFdQlaE0);

			// Rbar and partials
			double Rbar = Double.NaN;
			double dRbardF = Double.NaN;
			if (FoRbar >= phi0) {
				Rbar = F / FoRbar; // C1
				dRbardF = 1.0 / FoRbar;
				final double dRbardFoRbar = -F / Math.pow(FoRbar, 2);
				writeJacobian(oRbar, tagE0, dRbardFoRbar * dFoRbardu0 / Ea, rm);
				rm.setEntry(oRbar, iZbarb, dRbardFoRbar * dFoRbardZbarb);
			} else {
				Rbar = F / phi0;
				dRbardF = 1.0 / phi0;
				rm.setEntry(oRbar, iPhi0, -F / Math.pow(phi0, 2));
			}
			rv.setEntry(oRbar, Rbar);
			rm.setEntry(oRbar, iR, dRbardF * dFdR);
			rm.setEntry(oRbar, iOneOverS, dRbardF * dFdOneOverS);
			rm.setEntry(oRbar, iQlaE0, dRbardF * dFdQlaE0);

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			final int iOneOverS = inputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int iQlaE0 = inputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iR = inputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			checkIndices(iZbarb, iOneOverS, iQlaE0, iPhi0, iR);

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = getValue(MatrixCorrectionModel2.beamEnergyLabel(mDatum), point) / Ea;
			final double Zbarb = point.getEntry(iZbarb);
			final double oneOverS = point.getEntry(iOneOverS);
			final double QlaE0 = point.getEntry(iQlaE0);
			final double phi0 = point.getEntry(iPhi0);
			final double R = point.getEntry(iR);

			final double logZbarb = Math.log(Zbarb);
			final double X = 1.0 + 1.3 * logZbarb;
			final double Y = 0.2 + 0.005 * Zbarb;
			final double pu42 = Math.pow(u0, 0.42);
			final double FoRbar = 1.0 + (X * Math.log(1.0 + Y * (1.0 - 1.0 / pu42)) / Math.log(1.0 + Y));

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oRbar = outputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int oF = outputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			checkIndices(oRbar, oF);

			// F and partials
			final double F = R * oneOverS / QlaE0; // C1
			rv.setEntry(oF, F);

			// Rbar and partials
			double Rbar = Double.NaN;
			if (FoRbar >= phi0)
				Rbar = F / FoRbar; // C1
			else
				Rbar = F / phi0;
			rv.setEntry(oRbar, Rbar);
			return rv;
		}

		@Override
		public String toString() {
			return "StepFRbar[" + mDatum + "," + mShell + "]";
		}

	}

	private static class StepP extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variates.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		StepP(final MatrixCorrectionDatum comp, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(comp, shell, variates), buildOutputs(comp, shell));
			mDatum = comp;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(iZbarb, iRbar, iF, iphi0);

			final double Zbarb = point.getEntry(iZbarb);
			final double e0 = getValue(tagE0, point);
			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			final double Rbar = point.getEntry(iRbar);
			final double F = point.getEntry(iF);
			final double phi0 = point.getEntry(iphi0);

			// P and partials
			final double kg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
			final double g = 0.22 * Math.log(4.0 * Zbarb) * (1.0 - 2.0 * kg);
			final double h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 2.0);
			final double b = Math.sqrt(2.0) * (1.0 + Math.sqrt(1.0 - Rbar * phi0 / F)) / Rbar;
			final double sqrt2 = Math.sqrt(2.0), Rbar2 = Math.pow(Rbar, 2.0);

			// Two different ways to compute gh4
			final double gh4_1 = g * Math.pow(h, 4.0);
			final double gh4_2 = 0.9 * b * Rbar2 * (b - 2.0 * phi0 / F);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oP = outputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ob = outputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(oP, ob);

			if (gh4_1 < gh4_2) {
				// Depends on Zbarb, u0, F, Rbar
				final double P = gh4_1 * F / Rbar2;
				rv.setEntry(oP, P);
				final double k1 = Math.exp(((-1.0 + u0) * Zbarb) / 15.0);
				final double Zbarb2 = Math.pow(Zbarb, 2), Zbarb3 = Math.pow(Zbarb, 3.0);
				final double k2 = 10.0 * Zbarb2 + u0 * (-10.0 + Zbarb2);
				final double u02 = Math.pow(u0, 2.0);
				rm.setEntry(oP, iZbarb,
						(11.0 * F * Math.pow(k2, 3.0)
								* (15.0 * (-2.0 + k1) * k2 + 2.0
										* (-10.0 * Zbarb3 + u02 * Zbarb * (-10.0 + Zbarb2)
												+ u0 * (-1200.0 + 600.0 * k1 + 10.0 * Zbarb + 9.0 * Zbarb3))
										* Math.log(4.0 * Zbarb)))
								/ (750.0 * k1 * Rbar2 * Math.pow(10.0 + u0, 4.0) * Math.pow(Zbarb, 9.0)));
				writeJacobian(oP, tagE0, //
						((11.0 * F * Math.pow(k2, 3.0)
								* (-3000.0 * k1 + u02 * Zbarb * (-10.0 + Zbarb2) + 20.0 * u0 * Zbarb * (-5.0 + Zbarb2)
										+ 100.0 * (60.0 + Zbarb3))
								* Math.log(4.0 * Zbarb))
								/ (375.0 * k1 * Rbar2 * Math.pow(10.0 + u0, 5.0) * Math.pow(Zbarb, 8.0))) / Ea,
						rm);
				rm.setEntry(oP, iF, P / F);
				rm.setEntry(oP, iRbar, -2.0 * P / Rbar);
			} else {
				// Depends on F, Rbar, b, phi0
				final double P = gh4_2 * F / Rbar2;
				rv.setEntry(oP, P);
				final double k2 = Math.sqrt(1.0 - (phi0 * Rbar) / F);
				rm.setEntry(oP, iphi0, (-0.9 * sqrt2 * (-3.0 * phi0 * Math.pow(Rbar, 3.0) + F * Rbar2 * (2.0 + 2.0 * k2)
						+ Math.pow(F, 2.0) * (1. + k2) * sqrt2)) / (F * Math.pow(Rbar, 3.0) * k2));
				rm.setEntry(oP, iF, (-0.9 * sqrt2
						* (Math.pow(phi0, 2.0) * Math.pow(Rbar, 4.0) + Math.pow(F, 3.0) * (-4.0 - 4.0 * k2) * sqrt2
								+ Math.pow(F, 2.0) * phi0 * Rbar * (3.0 + k2) * sqrt2))
						/ (Math.pow(F, 2.0) * Math.pow(Rbar, 4.0) * k2));
				rm.setEntry(oP, iRbar,
						(-7.2 * sqrt2
								* (0.125 * Math.pow(phi0, 2.0) * Math.pow(Rbar, 4.0)
										+ F * phi0 * Math.pow(Rbar, 3.0) * (-0.25 - 0.25 * k2)
										+ Math.pow(F, 2.0) * phi0 * Rbar * (-0.875 - 0.375 * k2) * sqrt2
										+ Math.pow(F, 3.0) * (1.0 + 1.0 * k2) * sqrt2))
								/ (F * Math.pow(Rbar, 5.0) * k2));
			}

			final double k3 = Math.sqrt(1. - (phi0 * Rbar) / F);

			rv.setEntry(ob, b);
			rm.setEntry(ob, iRbar, (-1. * (-0.5 * phi0 * Rbar + F * (1. + k3)) * sqrt2) / (F * Rbar2 * k3));
			rm.setEntry(ob, iphi0, -sqrt2 / (2. * F * k3));
			rm.setEntry(ob, iF, (phi0 * sqrt2) / (2. * Math.pow(F, 2) * k3));

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iZbarb = inputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mDatum.getComposition()));
			final int iRbar = inputIndex(MatrixCorrectionModel2.shellLabel(RBAR, mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(iZbarb, iRbar, iF, iphi0);

			final double Zbarb = point.getEntry(iZbarb);
			final double e0 = getValue(MatrixCorrectionModel2.beamEnergyLabel(mDatum), point);
			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			final double Rbar = point.getEntry(iRbar);
			final double F = point.getEntry(iF);
			final double phi0 = point.getEntry(iphi0);

			// P and partials
			final double kg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
			final double g = 0.22 * Math.log(4.0 * Zbarb) * (1.0 - 2.0 * kg);
			final double h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 2.0);
			final double b = Math.sqrt(2.0) * (1.0 + Math.sqrt(1.0 - Rbar * phi0 / F)) / Rbar;
			final double Rbar2 = Math.pow(Rbar, 2.0);

			// Two different ways to compute gh4
			final double gh4_1 = g * Math.pow(h, 4.0);
			final double gh4_2 = 0.9 * b * Rbar2 * (b - 2.0 * phi0 / F);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oP = outputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ob = outputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(oP, ob);

			if (gh4_1 < gh4_2) {
				// Depends on Zbarb, u0, F, Rbar
				final double P = gh4_1 * F / Rbar2;
				rv.setEntry(oP, P);
			} else {
				// Depends on F, Rbar, b, phi0
				final double P = gh4_2 * F / Rbar2;
				rv.setEntry(oP, P);
			}
			rv.setEntry(ob, b);

			return rv;
		}

		@Override
		public String toString() {
			return "StepP[" + mDatum + "," + mShell + "]";
		}

	}

	private static class Stepa // C1
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			return Collections.singletonList(MatrixCorrectionModel2.shellLabel("a", datum, shell));
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		public Stepa(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
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

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oa = outputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(oa);

			final double b2 = Math.pow(b, 2.0);
			final double den = Math.pow(phi0 + b2 * F * Rbar - 2.0 * b * F, 2.0);
			final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // C1

			rv.setEntry(oa, a); // C1
			rm.setEntry(oa, iP, 1.0 / (b * F * (2.0 - b * Rbar) - phi0)); // C1
			rm.setEntry(oa, ib,
					(-2.0 * (F * P + phi0 * phi0 - b * F * (phi0 + P * Rbar) + b2 * F * (F - phi0 * Rbar))) / den); // C1
			rm.setEntry(oa, iphi0, (3.0 * b2 * F + P - 2.0 * b2 * b * F * Rbar) / den); // C1
			rm.setEntry(oa, iF, (b * (P * (-2.0 + b * Rbar) + b * phi0 * (2.0 * b * Rbar - 3.0))) / den); // C1
			rm.setEntry(oa, iRbar, (b2 * F * (P + 2.0 * b * phi0 - b2 * F)) / den); // C1

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
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

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oa = outputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			checkIndices(oa);

			final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // C1

			rv.setEntry(oa, a); // C1
			return rv;
		}

		@Override
		public String toString() {
			return "Stepa[" + mDatum + "," + mShell + "]";
		}
	};

	static private class StepEps // C1
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private static final double MIN_EPS = 1.0e-6;

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			return Collections.singletonList(MatrixCorrectionModel2.shellLabel(EPS, datum, shell));
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("a", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		public StepEps(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(ia, ib);

			final double a = point.getEntry(ia);
			final double b = point.getEntry(ib);

			final int oeps = outputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mShell));
			checkIndices(oeps);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final double eps = (a - b) / b;
			if (Math.abs(eps) >= MIN_EPS) {
				rv.setEntry(oeps, eps); // C1
				rm.setEntry(oeps, ia, 1.0 / b); // C1
				rm.setEntry(oeps, ib, -a / (b * b)); // C1
			} else {
				rv.setEntry(oeps, MIN_EPS); // C1
				rm.setEntry(oeps, ia, 0.0); // C1
				rm.setEntry(oeps, ib, 0.0); // C1
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int ia = inputIndex(MatrixCorrectionModel2.shellLabel("a", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(ia, ib);

			final double a = point.getEntry(ia);
			final double b = point.getEntry(ib);

			final int oeps = outputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mShell));
			checkIndices(oeps);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final double eps = (a - b) / b;
			if (Math.abs(eps) >= MIN_EPS) {
				rv.setEntry(oeps, eps); // C1
			} else {
				rv.setEntry(oeps, MIN_EPS); // C1
			}
			return rv;
		}

		@Override
		public String toString() {
			return "StepEPS[" + mDatum + "," + mShell + "]";
		}
	}

	private static class StepAB // C1
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("A", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("B", datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(EPS, datum, shell));
			return res;
		}

		public StepAB(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			super(buildInputs(datum, shell), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			final int ieps = inputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mShell));
			checkIndices(iphi0, iF, iP, ib, ieps);

			final double phi0 = point.getEntry(iphi0);
			final double P = point.getEntry(iP);
			final double F = point.getEntry(iF);
			final double b = point.getEntry(ib);
			final double eps = point.getEntry(ieps);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			final int oA = outputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mShell));
			final int oB = outputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mShell));
			checkIndices(oA, oB);
			// Page 62
			final double B = (b * b * F * (1.0 + eps) - P - phi0 * b * (2.0 + eps)) / eps; // C2-Ok
			rv.setEntry(oB, B);
			rm.setEntry(oB, iphi0, (-b * (2.0 + eps)) / eps); // C2-Ok
			rm.setEntry(oB, iP, -1.0 / eps); // C2-Ok
			rm.setEntry(oB, iF, (b * b * (1.0 + eps)) / eps); // C2-Ok
			rm.setEntry(oB, ib, (2.0 * b * F * (1.0 + eps) - phi0 * (2.0 + eps)) / eps); // C2-Ok
			rm.setEntry(oB, ieps, (P + 2.0 * b * phi0 - b * b * F) / (eps * eps)); // C2-Ok
			final double k1 = (1.0 + eps) / (eps * eps); // C1-Ok
			// Plugging B[b,F,eps,phi0,P] into the expression for A[B,b,F,eps,phi0,P] we get
			// A[b,F,eps,phi0,P].
			rv.setEntry(oA, k1 * (b * (b * F - 2.0 * phi0) - P) / b); // C1-Ok
			rm.setEntry(oA, iphi0, -2.0 * k1); // C1-Ok
			rm.setEntry(oA, iP, -k1 / b); // C1-Ok
			rm.setEntry(oA, iF, k1 * b); // C1-Ok
			rm.setEntry(oA, ib, k1 * (F + P / (b * b))); // C1-Ok
			rm.setEntry(oA, ieps, (((2.0 + eps) * (P + 2.0 * b * phi0 - b * b * F)) / (b * Math.pow(eps, 3.0)))); // C1-Ok
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iphi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			final int iF = inputIndex(MatrixCorrectionModel2.shellLabel("F", mDatum, mShell));
			final int iP = inputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			final int ieps = inputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mShell));
			checkIndices(iphi0, iF, iP, ib, ieps);

			final double phi0 = point.getEntry(iphi0);
			final double P = point.getEntry(iP);
			final double F = point.getEntry(iF);
			final double b = point.getEntry(ib);
			final double eps = point.getEntry(ieps);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final int oA = outputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mShell));
			final int oB = outputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mShell));
			checkIndices(oA, oB);
			// Page 62
			final double B = (b * b * F * (1.0 + eps) - P - phi0 * b * (2.0 + eps)) / eps; // C2-Ok
			rv.setEntry(oB, B);
			final double k1 = (1.0 + eps) / (eps * eps); // C1-Ok
			// Plugging B[b,F,eps,phi0,P] into the expression for A[B,b,F,eps,phi0,P] we get
			// A[b,F,eps,phi0,P].
			rv.setEntry(oA, k1 * (b * (b * F - 2.0 * phi0) - P) / b); // C1-Ok
			return rv;
		}

		@Override
		public String toString() {
			return "StepAB[" + mDatum + "," + mShell + "]";
		}
	}

	private static class StepAaBb extends SerialLabeledMultivariateJacobianFunction
			implements ILabeledMultivariateFunction {

		private static List<LabeledMultivariateJacobianFunction> buildSteps(final MatrixCorrectionDatum datum,
				final AtomicShell shell, final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
			res.add(new StepRphi0(datum, shell, variates));
			res.add(new StepFRbar(datum, shell, variates));
			res.add(new StepP(datum, shell, variates));
			res.add(new Stepa(datum, shell));
			res.add(new StepEps(datum, shell));
			res.add(new StepAB(datum, shell));
			return res;
		}

		public StepAaBb(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variates> variates) throws ArgumentException {
			super("StepAaBb", buildSteps(datum, shell, variates));
		}

		@Override
		public String toString() {
			return "StepAaBb[]";
		}

	}

	private static class StepChi extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction { // C1

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		public StepChi(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(datum, cxr, variates), buildOutputs(datum, cxr));
			mDatum = datum;
			mXRay = cxr;
		}

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.chiLabel(datum, cxr));
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			if (variates.contains(MatrixCorrectionModel2.Variates.MassAbsorptionCofficient))
				res.add(matMacLabel(datum.getComposition(), cxr));
			if (variates.contains(MatrixCorrectionModel2.Variates.TakeOffAngle))
				res.add(MatrixCorrectionModel2.takeOffAngleLabel(datum));
			return res;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Object toaT = MatrixCorrectionModel2.takeOffAngleLabel(mDatum);
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = NullableRealMatrix.build(getOutputDimension(), getInputDimension());

			final int ochi = outputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			checkIndices(ochi);

			final double toa = getValue(toaT, point);
			assert (toa > 0.0) && (toa < 0.5 * Math.PI);
			final double csc = 1.0 / Math.sin(toa);

			final Object macT = MatrixCorrectionModel2.matMacLabel(mDatum.getComposition(), mXRay);
			final double mac = getValue(macT, point);
			final double chi = mac * csc;
			writeJacobian(ochi, macT, csc, rm);
			writeJacobian(ochi, toaT, -1.0 * chi / Math.tan(toa), rm); // C1
			rv.setEntry(ochi, chi);

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int ochi = outputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			checkIndices(ochi);

			final double toa = getValue(MatrixCorrectionModel2.takeOffAngleLabel(mDatum), point);
			assert (toa > 0.0) && (toa < 0.5 * Math.PI);
			final double csc = 1.0 / Math.sin(toa);
			final double mac = getValue(MatrixCorrectionModel2.matMacLabel(mDatum.getComposition(), mXRay), point);
			final double chi = mac * csc; // C1
			rv.setEntry(ochi, chi);

			return rv;
		}

		@Override
		public String toString() {
			return "StepChi[" + mDatum + "," + mXRay + "]";
		}

	}

	private static class StepFx extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			final AtomicShell shell = cxr.getInner();
			res.add(MatrixCorrectionModel2.shellLabel("A", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("B", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(EPS, datum, shell));
			res.add(MatrixCorrectionModel2.chiLabel(datum, cxr));
			if (variates.contains(MatrixCorrectionModel2.Variates.SurfaceRoughness))
				res.add(MatrixCorrectionModel2.roughnessLabel(datum));
			return res;
		}

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.FofChiLabel(datum, cxr));
		}

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		public StepFx(final MatrixCorrectionDatum unk, final CharacteristicXRay cxr,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(unk, cxr, variates), buildOutputs(unk, cxr));
			mDatum = unk;
			mXRay = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iA = inputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mXRay.getInner()));
			final int iB = inputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mXRay.getInner()));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mXRay.getInner()));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mXRay.getInner()));
			final int iChi = inputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			final int ieps = inputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mXRay.getInner()));
			checkIndices(iA, iB, ib, iPhi0, iChi, ieps);

			final double A = point.getEntry(iA);
			final double B = point.getEntry(iB);
			final double b = point.getEntry(ib);
			final double phi0 = point.getEntry(iPhi0);
			final double chi = point.getEntry(iChi);
			final double eps = point.getEntry(ieps);
			final MatrixCorrectionModel2.RoughnessLabel tagRoughness = MatrixCorrectionModel2.roughnessLabel(mDatum);
			final double sr = getValue(tagRoughness, point);

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final double k0 = b + chi;
			final double k1 = chi + b * (1 + eps);

			final double Fx = Math.exp(-chi * sr) * (B / k0 - (A * b * eps) / k1 + phi0) / k0;
			rv.setEntry(oFx, Fx); // C2
			rm.setEntry(oFx, iA, -(1.0 / k0) + 1.0 / k1); // C2
			rm.setEntry(oFx, iB, Math.pow(k0, -2.0)); // C2
			rm.setEntry(oFx, ib,
					(-(B / k0) + (A * b * eps) / k1
							+ k0 * (-(B / Math.pow(k0, 2.0)) - (A * chi * eps) / Math.pow(k1, 2.0)) - phi0)
							/ Math.pow(k0, 2.0)); // C2
			rm.setEntry(oFx, iPhi0, 1.0 / k0); // C2
			rm.setEntry(oFx, iChi, (-2.0 * B + k0 * ((A * b * eps * (2.0 * k0 + b * eps)) / Math.pow(k1, 2.0) - phi0))
					/ Math.pow(k0, 3.0)); // C2
			rm.setEntry(oFx, ieps, -((A * b) / Math.pow(k1, 2.0))); // C2
			writeJacobian(oFx, tagRoughness, -chi * Math.exp(-chi * sr) * Fx, rm);
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iA = inputIndex(MatrixCorrectionModel2.shellLabel("A", mDatum, mXRay.getInner()));
			final int iB = inputIndex(MatrixCorrectionModel2.shellLabel("B", mDatum, mXRay.getInner()));
			final int ib = inputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mXRay.getInner()));
			final int iPhi0 = inputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mXRay.getInner()));
			final int iChi = inputIndex(MatrixCorrectionModel2.chiLabel(mDatum, mXRay));
			final int ieps = inputIndex(MatrixCorrectionModel2.shellLabel(EPS, mDatum, mXRay.getInner()));
			checkIndices(iA, iB, ib, iPhi0, iChi, ieps);

			final double A = point.getEntry(iA);
			final double B = point.getEntry(iB);
			final double b = point.getEntry(ib);
			final double phi0 = point.getEntry(iPhi0);
			final double chi = point.getEntry(iChi);
			final double eps = point.getEntry(ieps);
			final double sr = getValue(MatrixCorrectionModel2.roughnessLabel(mDatum), point);

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final double k0 = b + chi;
			final double k1 = chi + b * (1 + eps);
			rv.setEntry(oFx, Math.exp(-chi * sr) * (B / k0 - (A * b * eps) / k1 + phi0) / k0); // C2
			return rv;
		}

		@Override
		public String toString() {
			return "StepFx[" + mDatum + "," + mXRay + "]";
		}

	}

	private static class StepZA extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final UnknownMatrixCorrectionDatum mUnknown;
		private final StandardMatrixCorrectionDatum mStandard;
		private final CharacteristicXRay mXRay;

		public static List<? extends Object> buildInputs( //
				final UnknownMatrixCorrectionDatum unk, //
				final MatrixCorrectionDatum std, //
				final CharacteristicXRay cxr, //
				final Set<MatrixCorrectionModel2.Variates> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FofChiLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", unk, cxr.getInner()));
			res.add(MatrixCorrectionModel2.FofChiLabel(std, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", std, cxr.getInner()));
			// if (variates.contains(MatrixCorrectionModel2.Variates.UnknownComposition))
			// res.add(Composition.buildMassFractionTag(unk.getEstimate(),
			// cxr.getElement()));
			// if (variates.contains(MatrixCorrectionModel2.Variates.StandardComposition))
			// res.add(Composition.buildMassFractionTag(std.getComposition(),
			// cxr.getElement()));
			return res;
		}

		public static List<? extends Object> buildOutputs(final UnknownMatrixCorrectionDatum unk,
				final StandardMatrixCorrectionDatum std, final CharacteristicXRay cxr) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FxFLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.FxFLabel(std, cxr));
			res.add(MatrixCorrectionModel2.zLabel(unk, std, cxr));
			res.add(MatrixCorrectionModel2.aLabel(unk, std, cxr));
			// res.add(MatrixCorrectionModel2.zafLabel(unk, std, cxr));
			// res.add(new KRatioLabel(unk, std, cxr, Method.Calculated));
			return res;
		}

		public StepZA(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputs(unk, std, cxr, variates), buildOutputs(unk, std, cxr));
			mUnknown = unk;
			mStandard = std;
			mXRay = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iFxu = inputIndex(MatrixCorrectionModel2.FofChiLabel(mUnknown, mXRay));
			final int iFxs = inputIndex(MatrixCorrectionModel2.FofChiLabel(mStandard, mXRay));
			final int iFu = inputIndex(MatrixCorrectionModel2.shellLabel("F", mUnknown, mXRay.getInner()));
			final int iFs = inputIndex(MatrixCorrectionModel2.shellLabel("F", mStandard, mXRay.getInner()));
			// final MassFractionTag cut =
			// Composition.buildMassFractionTag(mUnknown.getEstimate(), mXRay.getElement());
			// final MassFractionTag cst =
			// Composition.buildMassFractionTag(mStandard.getComposition(),
			// mXRay.getElement());

			checkIndices(iFxu, iFxs, iFu, iFs);

			final double Fxu = point.getEntry(iFxu);
			final double Fxs = point.getEntry(iFxs);
			final double Fu = point.getEntry(iFu);
			final double Fs = point.getEntry(iFs);
			// final double Cs = getValue(cst, point);
			// final double Cu = getValue(cut, point);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			// final int oZA = outputIndex(MatrixCorrectionModel2.zafLabel(mUnknown,
			// mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			// final int oK = outputIndex(new KRatioLabel(mUnknown, mStandard, mXRay,
			// Method.Calculated));
			checkIndices(oFxFu, oFxFs, oZ, oA);

			// rv.setEntry(oZA, Fxu / Fxs); // C2
			// rm.setEntry(oZA, iFxu, 1.0 / Fxs); // C2
			// rm.setEntry(oZA, iFxs, -1.0 * Fxu / (Fxs * Fxs)); // C2

			// final double k = (Fxu * Cu) / (Fxs * Cs);
			// rv.setEntry(oK, k);
			// rm.setEntry(oK, iFxu, k / Fxu);
			// rm.setEntry(oK, iFxs, -k / Fxs);
			// writeJacobian(oK, cut, k / Cu, rm);
			// writeJacobian(oK, cst, -k / Cs, rm);

			final double a = (Fs * Fxu) / (Fu * Fxs);
			rv.setEntry(oA, a); // C2
			rm.setEntry(oA, iFxu, a / Fxu); // C2
			rm.setEntry(oA, iFxs, -a / Fxs); // C2
			rm.setEntry(oA, iFs, a / Fs); // C2
			rm.setEntry(oA, iFu, -a / Fu); // C2

			final double z = Fu / Fs;
			rv.setEntry(oZ, z); // C2
			rm.setEntry(oZ, iFs, -z / Fs); // C2
			rm.setEntry(oZ, iFu, 1.0 / Fs); // C2

			final double FxFu = Fxu / Fu;
			rv.setEntry(oFxFu, FxFu); // C2
			rm.setEntry(oFxFu, iFu, -FxFu / Fu); // C2
			rm.setEntry(oFxFu, iFxu, 1.0 / Fu); // C2

			final double FxFs = Fxs / Fs;
			rv.setEntry(oFxFs, FxFs); // C2
			rm.setEntry(oFxFs, iFs, -FxFs / Fs); // C2
			rm.setEntry(oFxFs, iFxs, 1.0 / Fs); // C2

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iFxu = inputIndex(MatrixCorrectionModel2.FofChiLabel(mUnknown, mXRay));
			final int iFxs = inputIndex(MatrixCorrectionModel2.FofChiLabel(mStandard, mXRay));
			final int iFu = inputIndex(MatrixCorrectionModel2.shellLabel("F", mUnknown, mXRay.getInner()));
			final int iFs = inputIndex(MatrixCorrectionModel2.shellLabel("F", mStandard, mXRay.getInner()));
			// final int iCu =
			// inputIndex(Composition.buildMassFractionTag(mUnknown.getComposition(),
			// mXRay.getElement()));
			// final int iCs =
			// inputIndex(Composition.buildMassFractionTag(mStandard.getComposition(),
			// mXRay.getElement()));
			checkIndices(iFxu, iFxs, iFu, iFs); // , iCu, iCs);

			final double Fxu = point.getEntry(iFxu);
			final double Fxs = point.getEntry(iFxs);
			final double Fu = point.getEntry(iFu);
			final double Fs = point.getEntry(iFs);
			// final double Cs = point.getEntry(iCs);
			// final double Cu = point.getEntry(iCu);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			// final int oZA = outputIndex(MatrixCorrectionModel2.zafLabel(mUnknown,
			// mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			// final int oK = outputIndex(new KRatioLabel(mUnknown, mStandard, mXRay,
			// Method.Calculated));
			checkIndices(oFxFu, oFxFs, oZ, oA); // , oZA, oK);

			// rv.setEntry(oK, (Fxu * Cu) / (Fxs * Cs));
			// rv.setEntry(oZA, Fxu / Fxs); // C2
			rv.setEntry(oA, (Fs * Fxu) / (Fu * Fxs)); // C2
			rv.setEntry(oZ, (Fu / Fs)); // C2
			rv.setEntry(oFxFu, Fxu / Fu); // C2
			rv.setEntry(oFxFs, Fxs / Fs); // C2
			return rv;
		}

		@Override
		public String toString() {
			return "StepZA[" + mStandard + "," + mUnknown + "," + mXRay + "]";
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
	private static final class StepXPP extends SerialLabeledMultivariateJacobianFunction
			implements ILabeledMultivariateFunction {

		private static List<LabeledMultivariateJacobianFunction> buildSteps(final MatrixCorrectionDatum datum,
				final CharacteristicXRaySet exrs, final Set<MatrixCorrectionModel2.Variates> variates)
				throws ArgumentException {
			final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
			res.add(new StepMJZBarb(datum.getComposition(), variates, datum instanceof StandardMatrixCorrectionDatum));
			final Set<AtomicShell> shells = exrs.getSetOfInnerAtomicShells();
			{
				final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
				for (final AtomicShell shell : shells)
					step.add(new StepQlaE0OneOverS(datum, shell, variates));
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("QlaOoS", step));
			}
			{
				final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
				for (final AtomicShell shell : shells)
					step.add(new StepAaBb(datum, shell, variates));
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("AaBb", step));
			}
			{
				final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepChi(datum, cxr, variates));
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("Chi", step));
			}
			{
				final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepFx(datum, cxr, variates));
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("Fx", step));
			}
			return res;
		}

		public StepXPP(final MatrixCorrectionDatum datum, final CharacteristicXRaySet cxrs,
				final Set<MatrixCorrectionModel2.Variates> variates) throws ArgumentException {
			super("XPP[" + datum + "]", buildSteps(datum, cxrs, variates));
		}

		@Override
		public String toString() {
			return "StepXPP[]";
		}

	};

	/**
	 * This method joins together the individual Composition computations into a
	 * single final step that outputs the final matrix corrections.
	 *
	 * @param unk
	 * @param stds
	 * @return List&lt;NamedMultivariateJacobianFunction&gt;
	 * @throws ArgumentException
	 */
	private static List<LabeledMultivariateJacobianFunction> buildSteps( //
			final Set<KRatioLabel> kratios, //
			final Set<MatrixCorrectionModel2.Variates> variates) throws ArgumentException {
		final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();

		final Map<UnknownMatrixCorrectionDatum, CharacteristicXRaySet> unks = new HashMap<>();
		final Map<StandardMatrixCorrectionDatum, CharacteristicXRaySet> stds = new HashMap<>();
		for (final KRatioLabel krl : kratios) {
			assert krl.getMethod().equals(Method.Measured);
			final UnknownMatrixCorrectionDatum unk = krl.getUnknown();
			final StandardMatrixCorrectionDatum std = krl.getStandard();
			final ElementXRaySet exrs = krl.getXRaySet();
			if (!unks.containsKey(unk))
				unks.put(unk, new CharacteristicXRaySet());
			unks.get(unk).addAll(exrs);
			if (!stds.containsKey(std))
				stds.put(std, new CharacteristicXRaySet());
			stds.get(std).addAll(exrs);
		}
		{
			final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
			final CharacteristicXRaySet cxrs = new CharacteristicXRaySet();
			for (final Map.Entry<StandardMatrixCorrectionDatum, CharacteristicXRaySet> me : stds.entrySet()) {
				step.add(new StepXPP(me.getKey(), me.getValue(), variates));
				cxrs.addAll(me.getValue());
			}
			for (final Map.Entry<UnknownMatrixCorrectionDatum, CharacteristicXRaySet> me : unks.entrySet())
				step.add(new StepXPP(me.getKey(), me.getValue(), variates));
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("XPP", step));
		}
		if (true) {
			final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay())
					step.add(new StepZA(krl.getUnknown(), krl.getStandard(), cxr, variates));
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("ZA", step));
		}
		{
			final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				final UnknownMatrixCorrectionDatum unk = krl.getUnknown();
				final StandardMatrixCorrectionDatum std = krl.getStandard();
				final ElementXRaySet exrs = krl.getXRaySet();
				step.add(new MultiE0MultiLineModel(unk, std, exrs, variates));
			}
			if (step.size() > 0)
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("Combine", step));
		}
		return res;
	}

	public static Set<MatrixCorrectionModel2.Variates> defaultVariates() {
		final Set<MatrixCorrectionModel2.Variates> res = allVariates();
		res.remove(MatrixCorrectionModel2.Variates.SurfaceRoughness);
		return res;
	}

	public static Set<MatrixCorrectionModel2.Variates> allVariates() {
		final Set<MatrixCorrectionModel2.Variates> res = new HashSet<>(
				Arrays.asList(MatrixCorrectionModel2.Variates.values()));
		return res;
	}

	public static Set<MatrixCorrectionModel2.Variates> minimalVariates() {
		final Set<MatrixCorrectionModel2.Variates> res = new HashSet<>();
		res.add(MatrixCorrectionModel2.Variates.StandardComposition);
		res.add(MatrixCorrectionModel2.Variates.UnknownComposition);
		return res;
	}

	/**
	 * @param kratios  Set of KRatioLabel objects.
	 * @param variates Which inputs to consider when calculating uncertainties.
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2(//
			final Set<KRatioLabel> kratios, //
			final Set<MatrixCorrectionModel2.Variates> variates //
	) throws ArgumentException {
		super("XPP Matrix Correction", kratios, buildSteps(kratios, variates));
		mVariates.addAll(variates);
		initializeConstants();
	}

	/**
	 * @param kratios Set of KRatioLabel objects.
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2( //
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		this(kratios, defaultVariates());
	}

	/**
	 * @param krl  KRatioLabel A single {@link KRatioLabel}
	 * @param std  Composition
	 * @param scxr A {@link CharacteristicXRay} object. It is important that the
	 *             edge energies must be below the beam energy.
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2( //
			final KRatioLabel krl, //
			final Set<MatrixCorrectionModel2.Variates> variates //
	) throws ArgumentException {
		this(Collections.singleton(krl), variates);
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
		return new UncertainValue(j, "J[" + elm.getAbbrev() + "]", 0.03 * j);
	}

	/**
	 * Computes the material MACs as required for the standard and unknowns relative
	 * to all the {@link CharacteristicXRay}s that take part in the measurement.
	 *
	 * @param elmMacs
	 * @return UncertainValues containing {@link MaterialMAC} tags
	 * @throws ArgumentException
	 */
	public final UncertainValues computeMaterialMACs(final UncertainValues elmMacs) //
			throws ArgumentException {
		final Map<CharacteristicXRay, List<Composition>> mclc = new HashMap<>();
		final Set<Composition> comps = new HashSet<>();
		for (final KRatioLabel krl : mKRatios) {
			final Set<CharacteristicXRay> scxrs = krl.getXRaySet().getSetOfCharacteristicXRay();
			for (final CharacteristicXRay cxr : scxrs) {
				List<Composition> lc = mclc.get(cxr);
				if (lc == null) {
					lc = new ArrayList<>();
					mclc.put(cxr, lc);
				}
				final Composition stdComp = krl.getStandard().getComposition();
				comps.add(stdComp);
				if (!lc.contains(stdComp))
					lc.add(stdComp);
				final Composition unkComp = krl.getUnknown().getComposition();
				if (!lc.contains(unkComp))
					lc.add(unkComp);
			}
		}
		final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
		for (final Map.Entry<CharacteristicXRay, List<Composition>> me : mclc.entrySet())
			step.add(new MaterialMACFunction(me.getValue(), me.getKey()));
		final LabeledMultivariateJacobianFunction macs = LabeledMultivariateJacobianFunctionBuilder.join("MACs", step);
		final UncertainValues[] items = new UncertainValues[1 + comps.size()];
		comps.toArray(items);
		items[items.length - 1] = elmMacs;
		final UncertainValues all = UncertainValues.combine(items);
		return UncertainValues.propagate(macs, all);
	}

	/**
	 * This is a quick way to build many of the input parameters. Many of the input
	 * parameters are computed or tabulated. Some are input as experimental
	 * conditions like beam energy which are provided in the {@link KRatioLabel}
	 * objects. Using this information, it is possible to calculate all the
	 * necessary inputs to this {@link LabeledMultivariateJacobianFunction}.
	 *
	 * @param unk  MatrixCorrectionDatum associated with the unknown
	 * @param stds Set&lt;MatrixCorrectionDatum&gt; A set of
	 *             {@link MatrixCorrectionDatum} associated with the standards
	 * @return
	 * @throws ArgumentException
	 */
	public UncertainValues buildInput() //
			throws ArgumentException {
		return UncertainValues.extract(getInputLabels(), buildParameters(true));
	}

	/**
	 * Initialize the parameters that are being held constant.
	 *
	 * @param unk  MatrixCorrectionDatum associated with the unknown
	 * @param stds Set&lt;MatrixCorrectionDatum&gt; A set of
	 *             {@link MatrixCorrectionDatum} associated with the standards
	 * @throws ArgumentException
	 */
	private void initializeConstants() //
			throws ArgumentException {
		final UncertainValues constants = buildParameters(false);
		if (constants.getDimension() > 0) {
			final Map<Object, Double> mod = new HashMap<>();
			for (final Object tag : constants.getLabels())
				mod.put(tag, constants.getEntry(tag));
			initializeConstants(mod);
		}
	}

	/**
	 * Many of the input parameters are computed or tabulated. Some are input as
	 * experimental conditions like beam energy.
	 *
	 * @param withUnc true to return parameters with uncertainties or false for
	 *                those without
	 * @return UncertainValues object containing either the variable quantities or
	 *         the constant quantities
	 * @throws ArgumentException
	 */
	private UncertainValues buildParameters(final boolean withUnc) //
			throws ArgumentException {
		final List<UncertainValues> results = new ArrayList<>();
		if (mVariates.contains(MatrixCorrectionModel2.Variates.MeanIonizationPotential) == withUnc)
			results.add(buildMeanIonizationPotentials());
		if (mVariates.contains(MatrixCorrectionModel2.Variates.MassAbsorptionCofficient) == withUnc)
			results.add(getMaterialMACs());
		if (mVariates.contains(MatrixCorrectionModel2.Variates.BeamEnergy) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.beamEnergyLabel(krl.getUnknown()), krl.getUnknown().getBeamEnergy());
				vals.put(MatrixCorrectionModel2.beamEnergyLabel(krl.getStandard()), krl.getStandard().getBeamEnergy());
			}
			results.add(new UncertainValues(vals));
		}
		if (mVariates.contains(MatrixCorrectionModel2.Variates.TakeOffAngle) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getUnknown()),
						krl.getUnknown().getTakeOffAngle());
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getStandard()),
						krl.getStandard().getTakeOffAngle());
			}
			results.add(new UncertainValues(vals));
		}
		if (mVariates.contains(MatrixCorrectionModel2.Variates.SurfaceRoughness) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getUnknown()),
						new UncertainValue(0.0, krl.getUnknown().getRoughness()));
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getStandard()),
						new UncertainValue(0.0, krl.getStandard().getRoughness()));
			}
			results.add(new UncertainValues(vals));
		}
		if (mVariates.contains(MatrixCorrectionModel2.Variates.WeightsOfLines) == withUnc) {
			final CharacteristicXRaySet allCxr = new CharacteristicXRaySet();
			for (final KRatioLabel krl : mKRatios)
				allCxr.addAll(krl.getXRaySet());
			final int sz = allCxr.size();
			if (sz > 0) {
				final Map<Object, Number> weightT = new HashMap<>();
				for (final CharacteristicXRay cxr : allCxr.getSetOfCharacteristicXRay())
					weightT.put(new MatrixCorrectionModel2.XRayWeightLabel(cxr), cxr.getWeightUV());
				final UncertainValues uvs = new UncertainValues(weightT);
				results.add(uvs);
			}
		}
		if (mVariates.contains(MatrixCorrectionModel2.Variates.IonizationExponent) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			final Set<AtomicShell> shells = new HashSet<>();
			for (final KRatioLabel krl : mKRatios)
				shells.addAll(krl.getXRaySet().getSetOfInnerAtomicShells());
			for (final AtomicShell sh : shells)
				vals.put(new MatrixCorrectionModel2.IonizationExponentLabel(sh), MatrixCorrectionModel2.computeM(sh));
			results.add(new UncertainValues(vals));
		}

		{
			final Map<Object, Number> mon = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				for (final MatrixCorrectionDatum mcd : Arrays.asList(krl.getStandard(), krl.getUnknown()))
					for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
						final Object sfLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mcd,
								cxr.getInner());
						if (!mon.containsKey(sfLbl))
							mon.put(sfLbl, new UncertainValue(1.0, 0.01));
					}
			}
			results.add(new UncertainValues(mon));
		}
		// Make sure that there are no replicated Compositions
		final Set<Composition> allComps = new HashSet<>();
		for (final KRatioLabel krl : mKRatios) {
			if (mVariates.contains(MatrixCorrectionModel2.Variates.StandardComposition) == withUnc)
				allComps.add(krl.getStandard().getComposition());
			if (mVariates.contains(MatrixCorrectionModel2.Variates.UnknownComposition) == withUnc)
				allComps.add(krl.getUnknown().getComposition());
		}
		results.addAll(allComps);
		return UncertainValues.combine(results.toArray(new UncertainValues[results.size()]));
	}

	public Set<Element> getElements() {
		final Set<Element> elms = new TreeSet<>();
		for (final KRatioLabel krl : mKRatios)
			elms.add(krl.getXRaySet().getElement());
		return elms;
	}

	public void setMaterialMACs(final UncertainValues uvs) {
		mMaterialMACS = uvs;
	}

	public UncertainValues getMaterialMACs() //
			throws ArgumentException {
		if (mMaterialMACS == null)
			mMaterialMACS = buildMaterialMACs();
		return mMaterialMACS;
	}

	private UncertainValues buildMaterialMACs() //
			throws ArgumentException {
		final List<ElementalMAC.ElementMAC> elmMacs = new ArrayList<>();
		final Map<Composition, CharacteristicXRaySet> comps = new HashMap<>();
		final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
		{
			final Set<CharacteristicXRay> allCxr = new HashSet<>();
			for (final KRatioLabel krl : mKRatios) {
				{
					final Composition comp = krl.getUnknown().getComposition();
					if (!comps.containsKey(comp))
						comps.put(comp, new CharacteristicXRaySet());
					comps.get(comp).addAll(krl.getXRaySet());
				}
				{
					final Composition comp = krl.getStandard().getComposition();
					if (!comps.containsKey(comp))
						comps.put(comp, new CharacteristicXRaySet());
					comps.get(comp).addAll(krl.getXRaySet());
				}
				allCxr.addAll(krl.getXRaySet().getSetOfCharacteristicXRay());
			}
			for (final CharacteristicXRay cxr : allCxr) {
				final Set<Composition> mats = new HashSet<>();
				final Set<Element> elms = new HashSet<>();
				for (final Map.Entry<Composition, CharacteristicXRaySet> me : comps.entrySet())
					if (me.getValue().contains(cxr)) {
						mats.add(me.getKey());
						elms.addAll(me.getKey().getElementSet());
					}
				assert mats.size() > 1;
				funcs.add(new MaterialMACFunction(new ArrayList<>(mats), cxr));
				for (final Element elm : elms)
					elmMacs.add(new ElementalMAC.ElementMAC(elm, cxr));
			}
		}
		final LabeledMultivariateJacobianFunction all = //
				LabeledMultivariateJacobianFunctionBuilder.join("MaterialMACs", funcs);
		final UncertainValues macInps = new UncertainValues(elmMacs);
		final ElementalMAC em = new ElementalMAC();
		for (final ElementalMAC.ElementMAC emac : elmMacs)
			macInps.set(emac, em.compute(emac.getElement(), emac.getXRay()));
		final List<UncertainValues> inputs = new ArrayList<>();
		inputs.add(macInps);
		for (final Composition comp : comps.keySet())
			inputs.add(comp);
		return UncertainValues.propagate(all, UncertainValues.combine(inputs));
	}

	private UncertainValues buildMeanIonizationPotentials() {
		final Set<Element> elms = new HashSet<>();
		for (final KRatioLabel krl : mKRatios) {
			elms.addAll(krl.getStandard().getComposition().getElementSet());
			elms.addAll(krl.getUnknown().getComposition().getElementSet());
		}
		final RealVector vals = new ArrayRealVector(elms.size());
		final RealVector var = new ArrayRealVector(vals.getDimension());
		final List<Object> tags = new ArrayList<>();
		int i = 0;
		for (final Element elm : elms) {
			tags.add(MatrixCorrectionModel2.meanIonizationLabel(elm));
			final UncertainValue j = computeJi(elm);
			vals.setEntry(i, j.doubleValue());
			var.setEntry(i, j.variance());
			++i;
		}
		assert tags.size() == vals.getDimension();
		final UncertainValues mip = new UncertainValues(tags, vals, var);
		return mip;
	}

	public static LabeledMultivariateJacobianFunction extractUnknown(//
			final LabeledMultivariateJacobian jac, //
			final Composition unknown //
	) throws ArgumentException {
		final List<Object> inpOut = new ArrayList<>();
		final Composition mf = unknown;
		for (final Object mft : mf.getLabels()) {
			assert mft instanceof MassFractionTag;
			inpOut.add(mft);
		}
		return jac.extract(inpOut, inpOut);
	}

	public static LabeledMultivariateJacobianFunction extractRemainder( //
			final LabeledMultivariateJacobian jac, //
			final Composition unknown //
	) throws ArgumentException {
		final Composition mf = unknown;
		final List<Object> outTags = new ArrayList<>(jac.getOutputLabels());
		outTags.removeAll(mf.getLabels());
		final List<Object> inTags = new ArrayList<>(jac.getInputLabels());
		inTags.removeAll(mf.getLabels());
		return jac.extract(inTags, outTags);
	}

	private double phiRhoZ(final double rhoZ, final double a, final double b, final double A, final double B,
			final double phi0) {
		return A * Math.exp(-a * rhoZ) + (B * rhoZ + phi0 - A) * Math.exp(-b * rhoZ);
	}

	public DataFrame<Double> computePhiRhoZCurve( //
			final Map<Object, Double> outputs, //
			final double rhoZmax, //
			final double dRhoZ, //
			final double minWeight //
	) {
		final DataFrame<Double> res = new DataFrame<>();
		{
			final List<Double> vals = new ArrayList<>();
			for (double rhoZ = 0.0; rhoZ <= rhoZmax; rhoZ += dRhoZ)
				vals.add(rhoZ);
			res.add("rhoZ", vals);
		}
		for (final KRatioLabel krl : mKRatios) {
			for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay()) {
				if (cxr.getWeight() > minWeight) {
					{
						final MatrixCorrectionDatum std = krl.getStandard();
						final double a = outputs.get(MatrixCorrectionModel2.shellLabel("a", std, cxr.getInner()))
								.doubleValue();
						final double b = outputs.get(MatrixCorrectionModel2.shellLabel("b", std, cxr.getInner()))
								.doubleValue();
						final double A = outputs.get(MatrixCorrectionModel2.shellLabel("A", std, cxr.getInner()))
								.doubleValue();
						final double B = outputs.get(MatrixCorrectionModel2.shellLabel("B", std, cxr.getInner()))
								.doubleValue();
						final double chi = outputs.get(new MatrixCorrectionModel2.ChiLabel(std, cxr)).doubleValue();
						final double phi0 = outputs.get(new MatrixCorrectionModel2.Phi0Label(std, cxr.getInner()))
								.doubleValue();
						final List<Double> vals1 = new ArrayList<>();
						final List<Double> vals2 = new ArrayList<>();
						for (double rhoZ = 0.0; rhoZ <= rhoZmax; rhoZ += dRhoZ) {
							final double prz = phiRhoZ(rhoZ, a, b, A, B, phi0);
							vals1.add(prz);
							vals2.add(prz * Math.exp(-chi * rhoZ));
						}
						res.add(std.getComposition().toString() + "[" + cxr.toString() + ", gen]", vals1);
						res.add(std.getComposition().toString() + "[" + cxr.toString() + ", emit]", vals2);
					}
					{
						final MatrixCorrectionDatum unk = krl.getUnknown();
						final double a = outputs.get(MatrixCorrectionModel2.shellLabel("a", unk, cxr.getInner()))
								.doubleValue();
						final double b = outputs.get(MatrixCorrectionModel2.shellLabel("b", unk, cxr.getInner()))
								.doubleValue();
						final double A = outputs.get(MatrixCorrectionModel2.shellLabel("A", unk, cxr.getInner()))
								.doubleValue();
						final double B = outputs.get(MatrixCorrectionModel2.shellLabel("B", unk, cxr.getInner()))
								.doubleValue();
						final double chi = outputs.get(new MatrixCorrectionModel2.ChiLabel(unk, cxr)).doubleValue();
						final double phi0 = outputs.get(new MatrixCorrectionModel2.Phi0Label(unk, cxr.getInner()))
								.doubleValue();
						final List<Double> vals1 = new ArrayList<>();
						final List<Double> vals2 = new ArrayList<>();
						for (double rhoZ = 0.0; rhoZ <= rhoZmax; rhoZ += dRhoZ) {
							final double prz = phiRhoZ(rhoZ, a, b, A, B, phi0);
							vals1.add(prz);
							vals2.add(prz * Math.exp(-chi * rhoZ));
						}
						res.add(unk.getComposition().toString() + "[" + cxr.toString() + ", gen]", vals1);
						res.add(unk.getComposition().toString() + "[" + cxr.toString() + ", emit]", vals2);
					}
				}
			}
		}
		return res;
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
}