package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRayEmissionProbability;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * Takes the matrix corrections for the individual lines and sums them as
 * appropriate for k-ratios measured from multiple lines simultaneously as is
 * the case in EDS. It also accounts for a difference in beam energy between the
 * standard and the unknown.
 *
 * Takes as inputs
 * <ol>
 * <li>Characteristic x-ray weights (XRayWeightTag)
 * <li>ZAF associated with single characteristic lines
 * (MatrixCorrectionModel2.zaLabel)</li>
 * <li>The ionization exponent 'm'
 * (MatrixCorrectionModel2.IonizationExponentLabel)</li>
 * </ol>
 * Returns as outputs
 * <ol>
 * <li>The effective ZAF for a set of characteristic x-ray lines from a single
 * element (MatrixCorrectionTag)</li>
 * <li>The calculated k-ratio (KRatioLabel)</li>
 * </ol>
 *
 * @author Nicholas W. M. Ritchie
 *
 */
class MultiE0MultiLineModel //
		extends SerialLabeledMultivariateJacobianFunction {

	private final ElementXRaySet mXRaySet;
	private final Set<MatrixCorrectionModel2.Variates> mVariates; //

	private static List<LabeledMultivariateJacobianFunction> buildSteps( //
			final UnknownMatrixCorrectionDatum unk, //
			final StandardMatrixCorrectionDatum std, //
			final ElementXRaySet exrs, //
			final Set<MatrixCorrectionModel2.Variates> variates //
	) {
		final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
		for (final AtomicShell shell : exrs.getSetOfInnerAtomicShells()) {
			res.add(new IonizationCrossSection(unk, shell));
			res.add(new IonizationCrossSection(std, shell));
		}
		res.add(new IntensityModel(unk, exrs, variates));
		res.add(new IntensityModel(std, exrs, variates));
		res.add(new kRatioModel(unk, std, exrs));
		return res;
	}

	public MultiE0MultiLineModel( //
			final UnknownMatrixCorrectionDatum unk, //
			final StandardMatrixCorrectionDatum std, //
			final ElementXRaySet exrs, //
			final Set<MatrixCorrectionModel2.Variates> variates //
	) throws ArgumentException {
		super("Multi-E0", buildSteps(unk, std, exrs, variates));
		mXRaySet = exrs;
		mVariates = variates;
	}

	private static Object intensityLabel(final MatrixCorrectionDatum mcd, final ElementXRaySet exrs) {
		return new BaseLabel<MatrixCorrectionDatum, ElementXRaySet, Object>("Intensity", mcd, exrs);
	}

	public static Object buildICXLabel(final MatrixCorrectionDatum mcd, final AtomicShell sh) {
		return new BaseLabel<MatrixCorrectionDatum, AtomicShell, Object>("ICX", mcd, sh);
	}

	private static class IntensityModel //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final MatrixCorrectionDatum mDatum;
		private final ElementXRaySet mXRaySet;

		private static List<? extends Object> buildInputTags(//
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs, //
				final Set<MatrixCorrectionModel2.Variates> variates //
		) {
			final List<Object> res = new ArrayList<>();
			for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				res.add(MatrixCorrectionModel2.FofChiLabel(mcd, cxr));
				if (variates.contains(MatrixCorrectionModel2.Variates.WeightsOfLines))
					res.add(new MatrixCorrectionModel2.XRayWeightLabel(cxr));
			}
			if (variates.contains(MatrixCorrectionModel2.Variates.IonizationExponent))
				for (final AtomicShell sh : exrs.getSetOfInnerAtomicShells()) {
					res.add(new MatrixCorrectionModel2.IonizationExponentLabel(sh));
					res.add(new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mcd, sh));
				}
			return res;
		}

		private static List<? extends Object> buildOutputTags(//
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs) {
			return Collections.singletonList(intensityLabel(mcd, exrs));
		}

		public IntensityModel(final MatrixCorrectionDatum mcd, final ElementXRaySet exrs,
				final Set<MatrixCorrectionModel2.Variates> variates) {
			super(buildInputTags(mcd, exrs, variates), buildOutputTags(mcd, exrs));
			mDatum = mcd;
			mXRaySet = exrs;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {

			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final int intIdx = outputIndex(intensityLabel(mDatum, mXRaySet));
			assert intIdx == 0;
			final Map<AtomicShell, Set<CharacteristicXRay>> shells = new HashMap<>();
			// Split the characteristic x-rays by shells
			for (final CharacteristicXRay cxr : mXRaySet.getSetOfCharacteristicXRay()) {
				final AtomicShell shell = cxr.getInner();
				if (!shells.containsKey(shell))
					shells.put(shell, new HashSet<>());
				shells.get(shell).add(cxr);
			}
			double res = 0.0;
			for (final Map.Entry<AtomicShell, Set<CharacteristicXRay>> me : shells.entrySet()) {
				final AtomicShell shell = me.getKey();
				final Object ionLbl = new MatrixCorrectionModel2.IonizationExponentLabel(shell);
				final Object secLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell);
				final double ionVal = getValue(ionLbl, point);
				final double secVal = getValue(secLbl, point);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final Object wLbl = new MatrixCorrectionModel2.XRayWeightLabel(cxr);
					final Object fxLbl = MatrixCorrectionModel2.FofChiLabel(mDatum, cxr);
					final double wVal = getValue(wLbl, point);
					final double fxVal = getValue(fxLbl, point);
					inner += wVal * fxVal;
					writeJacobian(intIdx, fxLbl, wVal * outer, rm);
					writeJacobian(intIdx, wLbl, fxVal * outer, rm);
				}
				writeJacobian(intIdx, ionLbl, secVal * inner, rm);
				writeJacobian(intIdx, secLbl, ionVal * inner, rm);
				res += outer * inner;
			}
			rv.setEntry(intIdx, res);
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final int intIdx = outputIndex(intensityLabel(mDatum, mXRaySet));
			final Map<AtomicShell, Set<CharacteristicXRay>> shells = new HashMap<>();
			// Split the characteristic x-rays by shells
			for (final CharacteristicXRay cxr : mXRaySet.getSetOfCharacteristicXRay()) {
				final AtomicShell shell = cxr.getInner();
				if (!shells.containsKey(shell))
					shells.put(shell, new HashSet<>());
				shells.get(shell).add(cxr);
			}
			double res = 0.0;
			for (final Map.Entry<AtomicShell, Set<CharacteristicXRay>> me : shells.entrySet()) {
				final AtomicShell shell = me.getKey();
				final Object ionLbl = new MatrixCorrectionModel2.IonizationExponentLabel(shell);
				final Object secLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell);
				final double ionVal = getValue(ionLbl, point);
				final double secVal = getValue(secLbl, point);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final Object wLbl = new MatrixCorrectionModel2.XRayWeightLabel(cxr);
					final Object fxLbl = MatrixCorrectionModel2.FofChiLabel(mDatum, cxr);
					final double wVal = getValue(wLbl, point);
					final double fxVal = getValue(fxLbl, point);
					inner += wVal * fxVal;
				}
				res += outer * inner;
			}
			rv.setEntry(intIdx, res);
			return rv;
		}

		@Override
		public String toString() {
			return "Intensity[" + mDatum.toString() + "," + mXRaySet.toString() + "]";
		}
	}

	private static class kRatioModel //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final UnknownMatrixCorrectionDatum mUnknown;
		private final StandardMatrixCorrectionDatum mStandard;
		private final ElementXRaySet mXRaySet;

		private static List<? extends Object> buildInputLabels( //
				final UnknownMatrixCorrectionDatum unk, //
				final StandardMatrixCorrectionDatum std, //
				final ElementXRaySet exrs) {
			final List<Object> res = new ArrayList<>();
			res.add(intensityLabel(unk, exrs));
			res.add(intensityLabel(std, exrs));
			res.add(Composition.buildMassFractionTag(unk.getComposition(), exrs.getElement()));
			res.add(Composition.buildMassFractionTag(std.getComposition(), exrs.getElement()));
			return res;
		}

		private static List<? extends Object> buildOutputLabels(//
				final UnknownMatrixCorrectionDatum unk, //
				final StandardMatrixCorrectionDatum std, //
				final ElementXRaySet exrs //
		) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.zafLabel(unk, std, exrs));
			res.add(new KRatioLabel(unk, std, exrs, Method.Calculated));
			return res;
		}

		public kRatioModel(final UnknownMatrixCorrectionDatum unk, final StandardMatrixCorrectionDatum std,
				final ElementXRaySet exrs) {
			super(buildInputLabels(unk, std, exrs), buildOutputLabels(unk, std, exrs));
			mUnknown = unk;
			mStandard = std;
			mXRaySet = exrs;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

			final Object zafLabel = MatrixCorrectionModel2.zafLabel(mUnknown, mStandard, mXRaySet);
			final Object kRatioLabel = new KRatioLabel(mUnknown, mStandard, mXRaySet, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			XPPMatrixCorrection2.checkIndices(kRow);
			final int zafRow = outputIndex(zafLabel);
			XPPMatrixCorrection2.checkIndices(zafRow);

			final Object intUnkLbl = intensityLabel(mUnknown, mXRaySet);
			final Object intStdLbl = intensityLabel(mStandard, mXRaySet);
			final Object cUnkLbl = Composition.buildMassFractionTag(mUnknown.getComposition(), mXRaySet.getElement());
			final Object cStdLbl = Composition.buildMassFractionTag(mStandard.getComposition(), mXRaySet.getElement());

			final double intUnk = getValue(intUnkLbl, point);
			final double intStd = getValue(intStdLbl, point);
			final double cUnk = getValue(cUnkLbl, point);
			final double cStd = getValue(cStdLbl, point);

			final double k = (cUnk * intUnk) / (cStd * intStd);
			rv.setEntry(kRow, k);
			writeJacobian(kRow, intUnkLbl, k / intUnk, rm);
			writeJacobian(kRow, intStdLbl, -k / intStd, rm);
			writeJacobian(kRow, cUnkLbl, k / cUnk, rm);
			writeJacobian(kRow, cStdLbl, -k / cStd, rm);

			rv.setEntry(zafRow, intUnk / intStd);
			writeJacobian(zafRow, intUnkLbl, 1.0 / intStd, rm);
			writeJacobian(zafRow, intStdLbl, -1.0 * intUnk / (intStd * intStd), rm);

			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());

			final Object zafLabel = MatrixCorrectionModel2.zafLabel(mUnknown, mStandard, mXRaySet);
			final Object kRatioLabel = new KRatioLabel(mUnknown, mStandard, mXRaySet, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			XPPMatrixCorrection2.checkIndices(kRow);
			final int zafRow = outputIndex(zafLabel);
			XPPMatrixCorrection2.checkIndices(zafRow);

			final Object intUnkLbl = intensityLabel(mUnknown, mXRaySet);
			final Object intStdLbl = intensityLabel(mStandard, mXRaySet);
			final Object cUnkLbl = Composition.buildMassFractionTag(mUnknown.getComposition(), mXRaySet.getElement());
			final Object cStdLbl = Composition.buildMassFractionTag(mStandard.getComposition(), mXRaySet.getElement());

			final double intUnk = getValue(intUnkLbl, point);
			final double intStd = getValue(intStdLbl, point);
			final double cUnk = getValue(cUnkLbl, point);
			final double cStd = getValue(cStdLbl, point);

			rv.setEntry(kRow, (cUnk * intUnk) / (cStd * intStd));
			rv.setEntry(zafRow, intUnk / intStd);

			return rv;
		}
	}

	private static class IonizationCrossSection //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final AtomicShell mShell;
		private final MatrixCorrectionDatum mDatum;

		static private List<Object> buildInputLabels(final MatrixCorrectionDatum mcd, final AtomicShell sh) {
			final List<Object> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.IonizationExponentLabel(sh));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(mcd));
			return res;
		}

		static private List<Object> buildOutputLabels(final MatrixCorrectionDatum mcd, final AtomicShell sh) {
			return Collections.singletonList(buildICXLabel(mcd, sh));
		}

		public IonizationCrossSection(final MatrixCorrectionDatum mcd, final AtomicShell sh) {
			super(buildInputLabels(mcd, sh), buildOutputLabels(mcd, sh));
			mShell = sh;
			mDatum = mcd;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			assert point.getDimension() == getInputDimension();
			assert getOutputDimension() == 1;
			final double eL = 1.0e-3 * mShell.getEdgeEnergy();
			final Object e0Lbl = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final double e0 = getValue(e0Lbl, point);
			final Object mLbl = new MatrixCorrectionModel2.IonizationExponentLabel(mShell);
			final double m = getValue(mLbl, point);

			final double u = e0 / eL;
			assert u > 1.0 : "The overvoltage for " + mDatum.toString() + " is less than unity.";
			final RealVector rv = new ArrayRealVector(1);
			final double logU = Math.log(u);
			final double icx = logU / (Math.pow(u, m) * eL * eL);
			rv.setEntry(0, icx);
			final RealMatrix rm = NullableRealMatrix.build(1, getInputDimension());
			writeJacobian(0, e0Lbl, (1.0 - m * logU) / (Math.pow(u, m + 1.0) * Math.pow(eL, 3.0)), rm);
			writeJacobian(0, mLbl, -icx * logU, rm);
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			assert point.getDimension() == getInputDimension();
			final double eL = 1.0e-3 * mShell.getEdgeEnergy();
			final Object e0Lbl = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			final double e0 = getValue(e0Lbl, point);
			final Object mLbl = new MatrixCorrectionModel2.IonizationExponentLabel(mShell);
			final double m = getValue(mLbl, point);
			final double u = e0 / eL;
			final RealVector rv = new ArrayRealVector(1);
			final double icx = Math.log(u) / (Math.pow(u, m) * eL * eL);
			rv.setEntry(0, icx);
			return rv;
		}
	}

	/**
	 * <p>
	 * Builds those inputs which are unlikely to have been calculated or provided
	 * previously. This includes the weights-of-lines and the exponent 'm' in the
	 * expression for the ionization cross section dependence on E0.
	 * </p>
	 * <p>
	 * Does not build massFractionLabel items
	 * </p>
	 *
	 *
	 *
	 * @param variables
	 * @return {@link UncertainValues}
	 */
	private UncertainValues buildInputs(final boolean variables) {
		final Map<Object, Number> vals = new HashMap<>();
		XRayEmissionProbability xrep = null;
		for (final CharacteristicXRay cxr : mXRaySet.getSetOfCharacteristicXRay()) {
			// MatrixCorrectionModel2.FofChiLabel(mcd, cxr) => Calculated elsewhere...
			if (mVariates.contains(MatrixCorrectionModel2.Variates.WeightsOfLines) == variables) {
				if ((xrep == null) || (xrep.getIonized() != cxr.getInner()))
					xrep = new XRayEmissionProbability(cxr.getInner());
				vals.put(new MatrixCorrectionModel2.XRayWeightLabel(cxr), xrep.getWeightU(cxr));
			}
		}
		if (mVariates.contains(MatrixCorrectionModel2.Variates.IonizationExponent) == variables)
			for (final AtomicShell sh : mXRaySet.getSetOfInnerAtomicShells())
				vals.put(new MatrixCorrectionModel2.IonizationExponentLabel(sh), MatrixCorrectionModel2.computeM(sh));
		return new UncertainValues(vals);
	}

	public UncertainValues buildInputs() {
		return buildInputs(true);
	}

	public UncertainValues buildConstants() {
		return buildInputs(false);
	}

}