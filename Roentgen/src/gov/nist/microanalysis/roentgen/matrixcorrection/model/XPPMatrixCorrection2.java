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
import gov.nist.microanalysis.roentgen.physics.XRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import gov.nist.microanalysis.roentgen.physics.composition.Layer;
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

	private static final String EPS = "&epsilon;";
	private static final String ONE_OVER_S = "<sup>1</sup>/<sub>S</sub>";
	private static final String QLA = "<html>Q<sub>l</sub><sup>a</sup>";
	private static final String ZBARB = "Z<sub>barb</sub>";
	private static final String RBAR = "R<sub>bar</sub>";

	private static final boolean VALIDATE = false;

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
	static private void checkIndices(final int... idxs) {
		final int len = idxs.length;
		for (int i = 0; i < len; ++i) {
			assert idxs[i] >= 0;
			for (int j = i + 1; j < len; ++j)
				assert idxs[i] != idxs[j];
		}
	}

	/**
	 * Ensure that the optimized and compute value are identical.
	 * 
	 * @param lmjf     {@link ILabeledMultivariateFunction}
	 * @param point    Evaluation point
	 * @param computed Comparison value
	 */
	static private void checkOptimized(ILabeledMultivariateFunction lmjf, RealVector inp, RealVector computed) {
		final RealVector optimized = lmjf.optimized(inp);
		for (int i = 0; i < optimized.getDimension(); ++i) {
			final double opt = optimized.getEntry(i);
			final double val = computed.getEntry(i);
			assert Math.abs(opt - val) < 1.0e-6 * Math.abs(val + 1.0e-6);
		}
	}

	/**
	 * 
	 * 
	 * @param lmjf
	 * @param inp
	 * @param dinp
	 * @param computed
	 */
	static private void checkPartials(LabeledMultivariateJacobianFunction lmjf, RealVector inp, RealVector dinp,
			RealMatrix computed) {
		final RealMatrix delta = lmjf.computeDelta(inp, dinp);
		for (int r = 0; r < delta.getRowDimension(); r++)
			for (int c = 0; c < delta.getColumnDimension(); ++c)
				assert Math.abs(delta.getEntry(r, c) - computed.getEntry(r, c)) < 1.0e-3
						* (Math.abs(computed.getEntry(r, c)) + 1.0e-6) : //
				lmjf.getOutputLabel(r) + " " + lmjf.getInputLabel(c) + " <=> " + delta.getEntry(r, c) + "!="
						+ computed.getEntry(r, c);
	}

	private static class StepMJZBarb // Checked 14-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final Composition mComposition;

		public static List<? extends Object> buildInputs(final Composition comp,
				final Set<MatrixCorrectionModel2.Variate> variates, final boolean isStandard) {
			final List<Object> res = new ArrayList<>();
			final MatrixCorrectionModel2.Variate varType = isStandard
					? MatrixCorrectionModel2.Variate.StandardComposition
					: MatrixCorrectionModel2.Variate.UnknownComposition;
			for (final Element elm : comp.getElementSet()) {
				if (variates.contains(varType)) {
					res.add(Composition.buildMassFractionTag(comp, elm));
					res.add(Composition.buildAtomicWeightTag(comp, elm));
				}
				if (variates.contains(MatrixCorrectionModel2.Variate.MeanIonizationPotential))
					res.add(MatrixCorrectionModel2.meanIonizationLabel(elm));
			}
			return res;
		}

		public static List<? extends Object> buildOutputs(final Composition comp) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel("M", comp));
			res.add(new MatrixCorrectionModel2.CompositionLabel("J", comp));
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, comp));
			return res;
		}

		public StepMJZBarb(final Composition comp, final Set<MatrixCorrectionModel2.Variate> variates,
				final boolean isStandard) {
			super(buildInputs(comp, variates, isStandard), buildOutputs(comp));
			mComposition = comp;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Composition comp = mComposition;
			final List<Element> elms = new ArrayList<>(comp.getElementSet());
			final Object[] tagCi = new Object[elms.size()];
			final Object[] tagJi = new Object[elms.size()];
			final Object[] tagAi = new Object[elms.size()];
			final double[] Ci = new double[elms.size()];
			final double[] Ai = new double[elms.size()];
			final double[] Z = new double[elms.size()];
			final double[] Ji = new double[elms.size()];
			for (int i = 0; i < Ci.length; ++i) {
				final Element elm = elms.get(i);
				tagCi[i] = Composition.buildMassFractionTag(comp, elm);
				tagJi[i] = MatrixCorrectionModel2.meanIonizationLabel(elm);
				Ci[i] = getValue(tagCi[i], point);
				Ji[i] = getValue(tagJi[i], point);
				tagAi[i] = Composition.buildAtomicWeightTag(comp, elm);
				Ai[i] = getValue(tagAi[i], point);
				Z[i] = elm.getAtomicNumber();
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
				M += Ci[i] * Z[i] / Ai[i];
			rv.setEntry(oM, M); // c1
			for (int i = 0; i < tagCi.length; ++i) {
				writeJacobian(oM, tagCi[i], (Z[i] / Ai[i]), rm); // c1
				writeJacobian(oM, tagAi[i], -Ci[i] * (Z[i] / (Ai[i] * Ai[i])), rm); //
			}

			// Calculate J and partials
			double lnJ = 0.0;
			for (int i = 0; i < Ji.length; ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * Z[i] / Ai[i];
			lnJ /= M;
			final double J = Math.exp(lnJ); // keV
			rv.setEntry(oJ, J); // C2
			for (int i = 0; i < tagCi.length; ++i) {
				final double tmp = (J / M) * (Z[i] / Ai[i]);
				final double logJi = Math.log(Ji[i]);
				writeJacobian(oJ, tagJi[i], tmp * (Ci[i] / Ji[i]), rm); // Ok!
				writeJacobian(oJ, tagCi[i], tmp * (logJi - lnJ), rm); // Ok!
				writeJacobian(oJ, tagAi[i], tmp * (Ci[i] / Ai[i]) * (lnJ - logJi), rm); // Ok!
			}
			// Calculate Zbarb and partials
			double Zbt = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // Ok
			rv.setEntry(oZbarb, Zbt * Zbt);
			for (int i = 0; i < elms.size(); ++i)
				writeJacobian(oZbarb, tagCi[i], 2.0 * Math.sqrt(Z[i]) * Zbt, rm); // Ci
			if (VALIDATE) {
				checkOptimized(this, point, rv);
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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
				final double a = getValue(Composition.buildAtomicWeightTag(comp, elm), point);
				ZoA[i] = Z[i] / a;
			}

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oM = outputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mComposition));
			final int oJ = outputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mComposition));
			final int oZbarb = outputIndex(new MatrixCorrectionModel2.CompositionLabel(ZBARB, mComposition));
			checkIndices(oM, oJ, oZbarb);

			// Calculate M
			double M = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				M += Ci[i] * ZoA[i];
			rv.setEntry(oM, M); // c1
			// Calculate J
			double lnJ = 0.0;
			for (int i = 0; i < Ji.length; ++i)
				lnJ += Math.log(Ji[i]) * Ci[i] * ZoA[i];
			lnJ /= M;
			rv.setEntry(oJ, Math.exp(lnJ)); // C2
			// Calculate Zbarb
			double Zbt = 0.0;
			for (int i = 0; i < Ci.length; ++i)
				Zbt += Ci[i] * Math.sqrt(Z[i]); // Ok
			rv.setEntry(oZbarb, Zbt * Zbt);
			return rv;
		}

		@Override
		public String toString() {
			return "MJZBarb[" + mComposition.toString() + "]";
		}

	}

	private static class StepQlaE0OneOverS // Checked 14-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel("M", datum.getComposition()));
			res.add(new MatrixCorrectionModel2.CompositionLabel("J", datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variate.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			if (variates.contains(MatrixCorrectionModel2.Variate.IonizationExponent))
				res.add(new MatrixCorrectionModel2.IonizationExponentLabel(shell));
			return res;
		}

		public StepQlaE0OneOverS(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			super(buildInputs(datum, shell, variates), buildOutputs(datum, shell));
			mDatum = datum;
			mShell = shell;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final Object tagm = new MatrixCorrectionModel2.IonizationExponentLabel(mShell);
			final int iJ = inputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mDatum.getComposition()));
			final int iM = inputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mDatum.getComposition()));
			checkIndices(iJ, iM);

			final double e0 = getValue(tagE0, point);
			final double m = getValue(tagm, point);
			final double J = point.getEntry(iJ);
			final double M = point.getEntry(iM);

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

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
			rv.setEntry(oQlaE0, QlaE0); // C1
			writeJacobian(oQlaE0, tagE0, (1.0 - m * logU0) / (Math.pow(u0, 1.0 + m) * Math.pow(Ea, 3.0)), rm); // C1
			writeJacobian(oQlaE0, tagm, -QlaE0 * logU0, rm); // Ok

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
				double OoSk = kk * D[k] * v0ou0_pk * (u0tk * (T[k] * logU0 - 1.0) + 1.0) / (T[k] * T[k]); // d
				OoS += OoSk;
				dOoSdm += (kk * u0tk * v0ou0_pk * D[k] * logU0_2 - 2.0 * OoSk) * (dTdmk / T[k]); // d
				dOoSdJ += (1.0 / J + dDdJ[k] / D[k] + Math.log(v0ou0) * dPdJ[k] + (P[k] / v0) * dv0dJ
						- 2.0 * dTdJ[k] / T[k]) * OoSk + //
						kk * u0tk * v0ou0_pk * D[k] * logU0_2 * dTdJ[k] / T[k]; // d
				dOoSdU0 += (kk * u0tk * v0ou0_pk * D[k] * logU0) / u0; // d
			}
			// E0/J0 = J/Ea
			rv.setEntry(oOneOverS, OoS); // C2
			rm.setEntry(oOneOverS, iM, (-1.0 / M) * OoS); // C2
			rm.setEntry(oOneOverS, iJ, dOoSdJ); // C2
			writeJacobian(oOneOverS, tagE0, dOoSdU0 * du0de0, rm); // C2
			writeJacobian(oOneOverS, tagm, dOoSdm, rm);

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final Object tagE0 = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final Object tagm = new MatrixCorrectionModel2.IonizationExponentLabel(mShell);
			final int iJ = inputIndex(new MatrixCorrectionModel2.CompositionLabel("J", mDatum.getComposition()));
			final int iM = inputIndex(new MatrixCorrectionModel2.CompositionLabel("M", mDatum.getComposition()));
			checkIndices(iJ, iM);

			final double e0 = getValue(tagE0, point);
			final double m = getValue(tagm, point);
			final double J = point.getEntry(iJ);
			final double M = point.getEntry(iM);

			final int oOneOverS = outputIndex(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, mDatum, mShell));
			final int oQlaE0 = outputIndex(MatrixCorrectionModel2.shellLabel(QLA, mDatum, mShell));
			checkIndices(oOneOverS, oQlaE0);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final double Ea = eVtokeV(mShell.getEdgeEnergy());
			final double u0 = e0 / Ea;
			assert u0 > 1.0 : e0 + " over " + Ea;
			final double logU0 = Math.log(u0);
			// QlaE0 and partials
			final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
			rv.setEntry(oQlaE0, QlaE0); // C1
			writeJacobian(oQlaE0, tagE0, Math.pow(u0, -1.0 - m) * (1.0 - m * logU0) / Math.pow(Ea, 3.0), rm); // C1
			writeJacobian(oQlaE0, tagm, -QlaE0 * logU0, rm);

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
			rv.setEntry(oOneOverS, OoS); // C2
			return rv;
		}

		@Override
		public String toString() {
			return "StepQlaE0OneOverS[" + mDatum + "," + mShell + "]";
		}
	}

	private static class StepRphi0 // Checked 15-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variate.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			return res;
		}

		public StepRphi0(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
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
			final double Zbarb = point.getEntry(iZbarb);
			final double e0 = getValue(tagE0, point);

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
			assert Math.abs(Gu0 - ((u0 - 1.0 - (1.0 - Math.pow(u0, -1.0 - q)) / opq) / ((2.0 + q) * Ju0))) //
			< Math.abs(1.0e-6 * Gu0);
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

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			rv.setEntry(oR, R);
			rm.setEntry(oR, iZbarb, dRdZbarb);
			writeJacobian(oR, tagE0, dRde0, rm);

			rv.setEntry(oPhi0, phi0);
			rm.setEntry(oPhi0, iZbarb, dphi0dZbarb);
			writeJacobian(oPhi0, tagE0, dphi0de0, rm);

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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

			final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3))); // Ok!
			final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55); // Ok!
			final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar); // Ok!
			final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0); // Ok!
			final double Gu0 = (u0 * (1.0 + q) - (2.0 + q) + Math.pow(u0, -1.0 - q)) / ((1.0 + q) * (2.0 + q) * Ju0);
			final double R = 1.0 - etabar * Wbar * (1.0 - Gu0);

			final double u0_mr = Math.pow(u0, 2.3 * etabar - 2.0);
			final double phi0 = 1.0 + 3.3 * Math.pow(etabar, 1.2) * (1.0 - u0_mr); // Ok!

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oR = outputIndex(MatrixCorrectionModel2.shellLabel("R", mDatum, mShell));
			final int oPhi0 = outputIndex(MatrixCorrectionModel2.phi0Label(mDatum, mShell));
			checkIndices(oR, oPhi0);

			rv.setEntry(oR, R);

			rv.setEntry(oPhi0, phi0);
			return rv;
		}

		@Override
		public String toString() {
			return "StepRPhi0[" + mDatum + "," + mShell + "]";
		}

	}

	private static class StepFRbar // Checked 15-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variate.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.shellLabel(ONE_OVER_S, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(QLA, datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("R", datum, shell));
			return res;
		}

		public StepFRbar(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
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
				final double dU0dE0 = 1.0 / Ea;
				final double dRbardFoRbar = -F / Math.pow(FoRbar, 2); // Ok
				writeJacobian(oRbar, tagE0, dRbardFoRbar * dFoRbardu0 * dU0dE0, rm); // Ok
				rm.setEntry(oRbar, iZbarb, dRbardFoRbar * dFoRbardZbarb); // Ok
			} else {
				Rbar = F / phi0;
				dRbardF = 1.0 / phi0;
				rm.setEntry(oRbar, iPhi0, -F / Math.pow(phi0, 2));
				// dRbardE0 and dRbardZbarb = 0!
			}
			rv.setEntry(oRbar, Rbar);
			rm.setEntry(oRbar, iR, dRbardF * dFdR);
			rm.setEntry(oRbar, iOneOverS, dRbardF * dFdOneOverS);
			rm.setEntry(oRbar, iQlaE0, dRbardF * dFdQlaE0);
			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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

	private static class StepPb // Checked 16-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final AtomicShell mShell;

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum, final AtomicShell shell) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel("P", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			return res;
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.shellLabel(RBAR, datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("F", datum, shell));
			res.add(new MatrixCorrectionModel2.CompositionLabel(ZBARB, datum.getComposition()));
			if (variates.contains(MatrixCorrectionModel2.Variate.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(datum));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			return res;
		}

		StepPb(final MatrixCorrectionDatum comp, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) {
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

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oP = outputIndex(MatrixCorrectionModel2.shellLabel("P", mDatum, mShell));
			final int ob = outputIndex(MatrixCorrectionModel2.shellLabel("b", mDatum, mShell));
			checkIndices(oP, ob);

			if (gh4_1 < gh4_2) {
				// Use gh4_1. Depends on Zbarb, u0, F, Rbar
				final double P = gh4_1 * F / Math.pow(Rbar, 2.0);
				final double dPdzbarb = P * ((1.0 / g) * dgdzbarb + (4.0 / h) * dhdzbarb);
				final double dPdzu0 = P * ((1 / g) * dgdu0 + (4.0 / h) * dhdu0);

				rv.setEntry(oP, P);
				rm.setEntry(oP, iZbarb, dPdzbarb);
				writeJacobian(oP, tagE0, dPdzu0 * du0de0, rm);
				rm.setEntry(oP, iF, P / F);
				rm.setEntry(oP, iRbar, -2.0 * P / Rbar);
			} else {
				// Depends on F, Rbar, b, phi0
				final double P = 0.9 * b * (b * F - 2.0 * phi0);
				final double dPdb = 1.8 * (b * F - phi0);
				final double dPdF = dPdb * dbdF + 0.9 * b * b;
				final double dPdphi0 = -1.8 * b;
				rv.setEntry(oP, P);
				rm.setEntry(oP, iphi0, dPdphi0 + dPdb * dbdphi0);
				rm.setEntry(oP, iF, dPdF);
				rm.setEntry(oP, iRbar, dPdb * dbdrbar);
				rm.setEntry(oP, iRbar, 0.0);
			}

			final double k3 = Math.sqrt(1. - (phi0 * Rbar) / F);

			rv.setEntry(ob, b);
			rm.setEntry(ob, iRbar,
					(-1. * (-0.5 * phi0 * Rbar + F * (1. + k3)) * Math.sqrt(2.0)) / (F * Math.pow(Rbar, 2.0) * k3));
			rm.setEntry(ob, iphi0, -Math.sqrt(2.0) / (2. * F * k3));
			rm.setEntry(ob, iF, (phi0 * Math.sqrt(2.0)) / (2. * Math.pow(F, 2) * k3));

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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
			final double expArg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
			final double g = 0.22 * Math.log(4.0 * Zbarb) * (1.0 - 2.0 * expArg);
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

	private static class Stepa // Checked 16-Jan-2019
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
			final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // Ok
			final double dadphi0 = (3.0 * b2 * F + P - 2.0 * b2 * b * F * Rbar) / den; // Ok
			final double dadP = 1.0 / (b * F * (2.0 - b * Rbar) - phi0); // Ok
			final double dadb = -2.0 * (F * P + phi0 * phi0 - b * F * (phi0 + P * Rbar) + b2 * F * (F - phi0 * Rbar))
					/ den; // Ok
			final double dadF = b * (P * (-2.0 + b * Rbar) + b * phi0 * (2.0 * b * Rbar - 3.0)) / den; // Ok
			final double dadRbar = (b2 * F * (P + 2.0 * b * phi0 - b2 * F)) / den; // Ok

			rv.setEntry(oa, a); // C1
			rm.setEntry(oa, iP, dadP); // C1
			rm.setEntry(oa, ib, dadb); // C1
			rm.setEntry(oa, iphi0, dadphi0); // C1
			rm.setEntry(oa, iF, dadF); // C1
			rm.setEntry(oa, iRbar, dadRbar); // C1

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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

	static private class StepEps // Checked 16-Jan-2019
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
				// rm.setEntry(oeps, ia, 0.0); // C1
				// rm.setEntry(oeps, ib, 0.0); // C1
			}
			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
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

	private static class StepAB // Checked 17-Jan-2019
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
			final double B = (b * b * F * (1.0 + eps) - P - phi0 * b * (2.0 + eps)) / eps; // Ok
			final double dBdphi0 = -b * (2.0 + eps) / eps; // Ok
			final double dBdP = -1.0 / eps; // Ok
			final double dBdF = (b * b * (1.0 + eps)) / eps; // Ok
			final double dBdb = (2.0 * b * F * (1.0 + eps) - phi0 * (2.0 + eps)) / eps; // Ok
			final double dBdeps = (P + 2.0 * b * phi0 - b * b * F) / (eps * eps); // Ok

			final double kk = (1.0 + eps) / eps;
			final double A = (B / b + phi0 - b * F) * kk; // Ok
			final double dAdb = (dBdb / b - B / (b * b) - F) * kk; // Ok
			final double dAdphi0 = ((b + dBdphi0) / b) * kk; // Ok
			final double dAdF = ((dBdF - b * b) / b) * kk; // Ok
			final double dAdeps = -A / (eps * (1.0 + eps)) + (dBdeps / b) * kk; // Ok
			final double dAdP = kk * dBdP / b; // Ok

			rv.setEntry(oB, B);
			rm.setEntry(oB, iphi0, dBdphi0); // C2-Ok
			rm.setEntry(oB, iP, dBdP); // C2-Ok
			rm.setEntry(oB, iF, dBdF); // C2-Ok
			rm.setEntry(oB, ib, dBdb); // C2-Ok
			rm.setEntry(oB, ieps, dBdeps); // C2-Ok

			rv.setEntry(oA, A); // C1-Ok
			rm.setEntry(oA, iphi0, dAdphi0); // C1-Ok
			rm.setEntry(oA, iP, dAdP); // C1-Ok
			rm.setEntry(oA, iF, dAdF); // C1-Ok
			rm.setEntry(oA, ib, dAdb); // C1-Ok
			rm.setEntry(oA, ieps, dAdeps); // C1-Ok

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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

	private static class StepAaBb //
			extends SerialLabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private static List<LabeledMultivariateJacobianFunction> buildSteps(final MatrixCorrectionDatum datum,
				final AtomicShell shell, final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
			res.add(new StepRphi0(datum, shell, variates));
			res.add(new StepFRbar(datum, shell, variates));
			res.add(new StepPb(datum, shell, variates));
			res.add(new Stepa(datum, shell));
			res.add(new StepEps(datum, shell));
			res.add(new StepAB(datum, shell));
			return res;
		}

		public StepAaBb(final MatrixCorrectionDatum datum, final AtomicShell shell,
				final Set<MatrixCorrectionModel2.Variate> variates) throws ArgumentException {
			super("StepAaBb", buildSteps(datum, shell, variates));
		}

		@Override
		public String toString() {
			return "StepAaBb[]";
		}

	}

	private static class StepChi // Checked 16-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction { // C1

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		public StepChi(final MatrixCorrectionDatum datum, final CharacteristicXRay cxr,
				final Set<MatrixCorrectionModel2.Variate> variates) {
			super(buildInputs(datum, cxr, variates), buildOutputs(datum, cxr));
			mDatum = datum;
			mXRay = cxr;
		}

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.chiLabel(datum, cxr));
		}

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			if (variates.contains(MatrixCorrectionModel2.Variate.MassAbsorptionCofficient))
				res.add(matMacLabel(datum.getComposition(), cxr));
			if (variates.contains(MatrixCorrectionModel2.Variate.TakeOffAngle))
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

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
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

	private static class StepFx // Checked 16-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		public static List<? extends Object> buildInputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			final AtomicShell shell = cxr.getInner();
			res.add(MatrixCorrectionModel2.shellLabel("A", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("B", datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel("b", datum, shell));
			res.add(MatrixCorrectionModel2.phi0Label(datum, shell));
			res.add(MatrixCorrectionModel2.shellLabel(EPS, datum, shell));
			res.add(MatrixCorrectionModel2.chiLabel(datum, cxr));
			if (variates.contains(MatrixCorrectionModel2.Variate.SurfaceRoughness))
				res.add(MatrixCorrectionModel2.roughnessLabel(datum));
			if (variates.contains(MatrixCorrectionModel2.Variate.Coating) && (datum.hasCoating())) {
				res.add(MatrixCorrectionModel2.coatingMassThickness(datum));
				res.add(MatrixCorrectionModel2.takeOffAngleLabel(datum));
				res.add(MatrixCorrectionModel2.matMacLabel(datum.getCoating().getComposition(), cxr));
			}
			return res;
		}

		public static List<? extends Object> buildOutputs(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr) {
			return Collections.singletonList(MatrixCorrectionModel2.FofChiLabel(datum, cxr));
		}

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		public StepFx(final MatrixCorrectionDatum unk, final CharacteristicXRay cxr,
				final Set<MatrixCorrectionModel2.Variate> variates) {
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
			final MatrixCorrectionModel2.RoughnessLabel tagRoughness = MatrixCorrectionModel2.roughnessLabel(mDatum);
			checkIndices(iA, iB, ib, iPhi0, iChi, ieps);

			final double A = point.getEntry(iA);
			final double B = point.getEntry(iB);
			final double b = point.getEntry(ib);
			final double phi0 = point.getEntry(iPhi0);
			final double chi = point.getEntry(iChi);
			final double eps = point.getEntry(ieps);
			final double dz = getValue(tagRoughness, point);

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final double expDzChi = Math.exp(-chi * dz);

			final double bpchi = b + chi;
			final double arg = chi + b * (1.0 + eps);
			final double powArg = Math.pow(arg, -2.0);
			final double powBpchi2 = Math.pow(bpchi, -2.0);

			final double Fx = expDzChi * (B / bpchi - (A * b * eps) / arg + phi0) / bpchi;
			final double dFxdA = -expDzChi * (b * eps) / (bpchi * arg);
			final double dFxdB = expDzChi * powBpchi2;
			final double dFxdb = (-expDzChi * (B * powBpchi2 + A * chi * eps * powArg) - Fx) / bpchi;
			final double dFxdphi0 = expDzChi / bpchi;
			final double dFxdchi = (expDzChi * ((A * b * eps) * powArg - B * powBpchi2) - (1 + bpchi * dz) * Fx)
					/ bpchi;
			final double dFxddz = -chi * Fx;
			final double dFxdeps = -expDzChi * A * b * powArg;

			rv.setEntry(oFx, Fx); // C2
			rm.setEntry(oFx, iA, dFxdA); // C2
			rm.setEntry(oFx, iB, dFxdB); // C2
			rm.setEntry(oFx, ib, dFxdb); // C2
			rm.setEntry(oFx, iPhi0, dFxdphi0); // C2
			rm.setEntry(oFx, iChi, dFxdchi); // C2
			rm.setEntry(oFx, ieps, dFxdeps); // C2
			writeJacobian(oFx, tagRoughness, dFxddz, rm);
			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				RealVector dpt = point.mapMultiply(1.0e-6);
				final int itr = inputIndex(tagRoughness);
				if (itr >= 0)
					dpt.setEntry(itr, 1.0e-12);
				checkPartials(this, point, dpt, rm);
			}
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
			final double dz = getValue(MatrixCorrectionModel2.roughnessLabel(mDatum), point);

			final int oFx = outputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			checkIndices(oFx);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final double k0 = b + chi;
			final double k1 = chi + b * (1 + eps);
			// Must include dz here since z = 0 +- dz != 0
			final double Fx = Math.exp(-chi * dz) * (B / k0 - (A * b * eps) / k1 + phi0) / k0;

			rv.setEntry(oFx, Fx); // C2
			return rv;
		}

		@Override
		public String toString() {
			return "StepFx[" + mDatum + "," + mXRay + "]";
		}
	}

	private static class StepConductiveCoating //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final CharacteristicXRay mXRay;

		private static final List<? extends Object> buildInputLabels(final MatrixCorrectionDatum datum,
				final CharacteristicXRay cxr, final Set<Variate> variates) {
			final List<Object> res = new ArrayList<>();
			if (variates.contains(MatrixCorrectionModel2.Variate.Coating) && (datum.hasCoating())) {
				res.add(MatrixCorrectionModel2.coatingMassThickness(datum));
				res.add(MatrixCorrectionModel2.takeOffAngleLabel(datum));
				res.add(MatrixCorrectionModel2.matMacLabel(datum.getCoating().getComposition(), cxr));
			}
			res.add(MatrixCorrectionModel2.FofChiLabel(datum, cxr));
			return res;

		}

		private static final List<? extends Object> buildOutputLabels(//
				final MatrixCorrectionDatum datum, final CharacteristicXRay cxr //
		) {
			return Collections.singletonList(MatrixCorrectionModel2.FofChiReducedLabel(datum, cxr));
		}

		public StepConductiveCoating(final MatrixCorrectionDatum mcd, final CharacteristicXRay cxr,
				final Set<Variate> variates) {
			super(buildInputLabels(mcd, cxr, variates), buildOutputLabels(mcd, cxr));
			mDatum = mcd;
			mXRay = cxr;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final int iMassTh = inputIndex(MatrixCorrectionModel2.coatingMassThickness(mDatum));
			final int iFx = inputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			final int oFxRed = outputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, mXRay));

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			double coatTrans = 1.0;
			final double Fx = point.getEntry(iFx);
			// A simple model for absorption of x-rays by a thin coating...
			if (iMassTh != -1) {
				final double massTh = point.getEntry(iMassTh);
				final int itoa = inputIndex(MatrixCorrectionModel2.takeOffAngleLabel(mDatum));
				final int imu = inputIndex(
						MatrixCorrectionModel2.matMacLabel(mDatum.getCoating().getComposition(), mXRay));
				final double toa = point.getEntry(itoa);
				final double csc = 1.0 / Math.sin(toa);
				final double mu = point.getEntry(imu);
				coatTrans = Math.exp(-mu * massTh * csc);
				final double dcoatTransdmu = coatTrans * (-massTh * csc);
				final double dcoatTransdmassTh = coatTrans * (-mu * csc);
				final double dcoatTransdtoa = coatTrans * massTh * mu * csc * csc * Math.cos(toa);
				rm.setEntry(oFxRed, iMassTh, dcoatTransdmassTh * Fx);
				rm.setEntry(oFxRed, imu, dcoatTransdmu * Fx);
				rm.setEntry(oFxRed, itoa, dcoatTransdtoa * Fx);
				assert coatTrans > 0.9 : coatTrans;
			}
			rm.setEntry(oFxRed, iFx, coatTrans);
			rv.setEntry(oFxRed, coatTrans * Fx);

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iMassTh = inputIndex(MatrixCorrectionModel2.coatingMassThickness(mDatum));
			final int iFx = inputIndex(MatrixCorrectionModel2.FofChiLabel(mDatum, mXRay));
			final int oFxRed = outputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, mXRay));

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			double coatTrans = 1.0;
			final double Fx = point.getEntry(iFx);
			// A simple model for absorption of x-rays by a thin coating...
			if (iMassTh != -1) {
				final double massTh = point.getEntry(iMassTh);
				final int itoa = inputIndex(MatrixCorrectionModel2.takeOffAngleLabel(mDatum));
				final int imu = inputIndex(
						MatrixCorrectionModel2.matMacLabel(mDatum.getCoating().getComposition(), mXRay));
				final double toa = point.getEntry(itoa);
				final double csc = 1.0 / Math.sin(toa);
				final double mu = point.getEntry(imu);
				coatTrans = Math.exp(-mu * massTh * csc);
				assert coatTrans > 0.9 : coatTrans;
			}
			rv.setEntry(oFxRed, coatTrans * Fx);
			return rv;
		}

		@Override
		public String toString() {
			final Layer coating = mDatum.getCoating();
			return "StepCoating[" + (coating != null ? "coating=" + coating.toString() : "Uncoated") + "]";
		}
	}

	private static class StepZA // Checked 16-Jan-2019
			extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

		private final UnknownMatrixCorrectionDatum mUnknown;
		private final StandardMatrixCorrectionDatum mStandard;
		private final CharacteristicXRay mXRay;

		public static List<? extends Object> buildInputs( //
				final UnknownMatrixCorrectionDatum unk, //
				final MatrixCorrectionDatum std, //
				final CharacteristicXRay cxr, //
				final Set<MatrixCorrectionModel2.Variate> variates) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FofChiReducedLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", unk, cxr.getInner()));
			res.add(MatrixCorrectionModel2.FofChiReducedLabel(std, cxr));
			res.add(MatrixCorrectionModel2.shellLabel("F", std, cxr.getInner()));
			return res;
		}

		public static List<? extends Object> buildOutputs(final UnknownMatrixCorrectionDatum unk,
				final StandardMatrixCorrectionDatum std, final CharacteristicXRay cxr) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.FxFLabel(unk, cxr));
			res.add(MatrixCorrectionModel2.FxFLabel(std, cxr));
			res.add(MatrixCorrectionModel2.zLabel(unk, std, cxr));
			res.add(MatrixCorrectionModel2.aLabel(unk, std, cxr));
			return res;
		}

		public StepZA(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final CharacteristicXRay cxr, final Set<MatrixCorrectionModel2.Variate> variates) {
			super(buildInputs(unk, std, cxr, variates), buildOutputs(unk, std, cxr));
			mUnknown = unk;
			mStandard = std;
			mXRay = cxr;
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

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			checkIndices(oFxFu, oFxFs, oZ, oA);

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

			if (VALIDATE) {
				checkOptimized(this, point, rv);
				;
				checkPartials(this, point, point.mapMultiply(1.0e-6), rm);
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final int iFxu = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mUnknown, mXRay));
			final int iFxs = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mStandard, mXRay));
			final int iFu = inputIndex(MatrixCorrectionModel2.shellLabel("F", mUnknown, mXRay.getInner()));
			final int iFs = inputIndex(MatrixCorrectionModel2.shellLabel("F", mStandard, mXRay.getInner()));
			checkIndices(iFxu, iFxs, iFu, iFs);

			final double Fxu = point.getEntry(iFxu);
			final double Fxs = point.getEntry(iFxs);
			final double Fu = point.getEntry(iFu);
			final double Fs = point.getEntry(iFs);

			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final int oFxFu = outputIndex(MatrixCorrectionModel2.FxFLabel(mUnknown, mXRay));
			final int oFxFs = outputIndex(MatrixCorrectionModel2.FxFLabel(mStandard, mXRay));
			final int oZ = outputIndex(MatrixCorrectionModel2.zLabel(mUnknown, mStandard, mXRay));
			final int oA = outputIndex(MatrixCorrectionModel2.aLabel(mUnknown, mStandard, mXRay));
			checkIndices(oFxFu, oFxFs, oZ, oA);

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
				final CharacteristicXRaySet exrs, final Set<MatrixCorrectionModel2.Variate> variates)
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
			{
				final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
				for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
					step.add(new StepConductiveCoating(datum, cxr, variates));
				res.add(LabeledMultivariateJacobianFunctionBuilder.join("Coating", step));
			}
			return res;
		}

		public StepXPP(final MatrixCorrectionDatum datum, final CharacteristicXRaySet cxrs,
				final Set<MatrixCorrectionModel2.Variate> variates) throws ArgumentException {
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
			final Set<MatrixCorrectionModel2.Variate> variates) throws ArgumentException {
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
		{
			final List<LabeledMultivariateJacobianFunction> step = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				for (final CharacteristicXRay cxr : krl.getXRaySet().getSetOfCharacteristicXRay())
					step.add(new StepZA(krl.getUnknown(), krl.getStandard(), cxr, variates));
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("ZA", step));
		}
		// Need to do it this way to ensure that certain items aren't double
		// calculated...
		res.add(new MultiE0MultiLineModel(kratios, variates));
		return res;
	}

	/**
	 * @param kratios  Set of KRatioLabel objects.
	 * @param variates Which inputs to consider when calculating uncertainties.
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2(//
			final Set<KRatioLabel> kratios, //
			final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		super("XPP Matrix Correction", kratios, buildSteps(kratios, variates), variates);
		initializeConstants();
	}

	/**
	 * @param kratios Set of KRatioLabel objects.
	 * @throws ArgumentException
	 */
	public XPPMatrixCorrection2( //
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		this(kratios, MatrixCorrectionModel2.defaultVariates());
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
			final Set<MatrixCorrectionModel2.Variate> variates //
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
	 * Computes the material MACs as required for the standard and unknowns relative
	 * to all the {@link CharacteristicXRay}s that take part in the measurement.
	 *
	 * @param elmMacs
	 * @return UncertainValues containing {@link MaterialMAC} tags
	 * @throws ArgumentException
	 */
	public UncertainValues computeMaterialMACs(final UncertainValues elmMacs) //
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
		if (isSet(Variate.MeanIonizationPotential) == withUnc)
			results.add(buildMeanIonizationPotentials());
		if (isSet(Variate.MassAbsorptionCofficient) == withUnc)
			results.add(buildMaterialMACs(withUnc));
		if (isSet(Variate.BeamEnergy) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
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
			results.add(new UncertainValues(vals));
		}
		if (isSet(Variate.Coating) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				if (krl.getStandard().hasCoating()) {
					final MatrixCorrectionDatum mcd = krl.getStandard();
					vals.put(MatrixCorrectionModel2.coatingMassThickness(mcd), mcd.getCoating().getMassThickness());
				}
				if (krl.getUnknown().hasCoating()) {
					final MatrixCorrectionDatum mcd = krl.getUnknown();
					vals.put(MatrixCorrectionModel2.coatingMassThickness(mcd), mcd.getCoating().getMassThickness());
				}
			}
			if (!vals.isEmpty())
				results.add(new UncertainValues(vals));
		}
		if (isSet(Variate.IonizationExponent) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			final Set<AtomicShell> allShells = new HashSet<>();
			for (final KRatioLabel krl : mKRatios)
				for (final XRay cxr : krl.getXRaySet())
					allShells.add(((CharacteristicXRay) cxr).getInner());
			for (final AtomicShell sh : allShells)
				vals.put(new IonizationExponentLabel(sh), computeIonizationExponent(sh, 1.0));
			assert !vals.isEmpty();
			results.add(new UncertainValues(vals));
		}
		if (isSet(Variate.TakeOffAngle) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getUnknown()),
						krl.getUnknown().getTakeOffAngle());
				vals.put(MatrixCorrectionModel2.takeOffAngleLabel(krl.getStandard()),
						krl.getStandard().getTakeOffAngle());
			}
			results.add(new UncertainValues(vals));
		}
		if (isSet(Variate.SurfaceRoughness) == withUnc) {
			final Map<Object, Number> vals = new HashMap<>();
			for (final KRatioLabel krl : mKRatios) {
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getUnknown()),
						new UncertainValue(0.0, krl.getUnknown().getRoughness()));
				vals.put(MatrixCorrectionModel2.roughnessLabel(krl.getStandard()),
						new UncertainValue(0.0, krl.getStandard().getRoughness()));
			}
			if (!vals.isEmpty())
				results.add(new UncertainValues(vals));
		}
		if (isSet(Variate.WeightsOfLines) == withUnc) {
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
		if (isSet(Variate.SecondaryFluorescence) == withUnc) {
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
			if (isSet(Variate.StandardComposition) == withUnc)
				allComps.add(krl.getStandard().getComposition());
			if (isSet(Variate.UnknownComposition) == withUnc)
				allComps.add(krl.getUnknown().getComposition());
		}
		for (Composition comp : allComps) {
			results.add(comp.toMassFraction());
			results.add(comp.getAtomicWeights());
		}
		return UncertainValues.combine(results);
	}

	public Set<Element> getElements() {
		final Set<Element> elms = new TreeSet<>();
		for (final KRatioLabel krl : mKRatios)
			elms.add(krl.getXRaySet().getElement());
		return elms;
	}

	private UncertainValues buildMaterialMACs(final boolean withUnc) //
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
				if (krl.getUnknown().hasCoating() && (isSet(Variate.Coating) == withUnc)) {
					final Composition comp = krl.getUnknown().getCoating().getComposition();
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
				if (krl.getStandard().hasCoating() && (isSet(Variate.Coating) == withUnc)) {
					final Composition comp = krl.getStandard().getCoating().getComposition();
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
				// assert mats.size() > 1;
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
			inputs.add(comp.toMassFraction());
		return inputs.size() > 0 ? UncertainValues.propagate(all, UncertainValues.combine(inputs))
				: UncertainValues.NULL;
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
