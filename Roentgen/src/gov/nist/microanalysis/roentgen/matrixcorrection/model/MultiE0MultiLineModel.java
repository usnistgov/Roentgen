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
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunctionBuilder;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.Variate;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRayEmissionProbability;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;

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

	private final Set<KRatioLabel> mKRatios;
	private final Set<MatrixCorrectionModel2.Variate> mVariates; //

	private static List<LabeledMultivariateJacobianFunction> buildSteps( //
			final Set<KRatioLabel> kratios, //
			final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		Map<MatrixCorrectionDatum, Set<AtomicShell>> allMcd = new HashMap<>();

		for (KRatioLabel krl : kratios) {
			final Set<AtomicShell> shells = krl.getXRaySet().getSetOfInnerAtomicShells();
			{
				final MatrixCorrectionDatum mcd = krl.getStandard();
				if (!allMcd.containsKey(mcd))
					allMcd.put(mcd, new HashSet<>());
				allMcd.get(mcd).addAll(shells);
			}
			{
				final MatrixCorrectionDatum mcd = krl.getUnknown();
				if (!allMcd.containsKey(mcd))
					allMcd.put(mcd, new HashSet<>());
				allMcd.get(mcd).addAll(shells);
			}
		}

		final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
		{
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			for (Map.Entry<MatrixCorrectionDatum, Set<AtomicShell>> me : allMcd.entrySet())
				for (AtomicShell sh : me.getValue())
					funcs.add(new IonizationCrossSection(me.getKey(), sh, variates));
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("ICXs", funcs));
		}
		{
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				final ElementXRaySet cxrs = krl.getXRaySet();
				funcs.add(new IntensityModel(krl.getStandard(), cxrs, variates));
				funcs.add(new IntensityModel(krl.getUnknown(), cxrs, variates));
			}
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("Intensities", funcs));
		}
		{
			final List<LabeledMultivariateJacobianFunction> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				funcs.add(new kRatioModel(krl));
			res.add(LabeledMultivariateJacobianFunctionBuilder.join("K-ratios", funcs));
		}
		return res;
	}

	public MultiE0MultiLineModel( //
			final Set<KRatioLabel> kratios, final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		super("Multi-E0", buildSteps(kratios, variates));
		mKRatios = kratios;
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
				final Set<MatrixCorrectionModel2.Variate> variates //
		) {
			final List<Object> res = new ArrayList<>();
			for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				res.add(MatrixCorrectionModel2.FofChiReducedLabel(mcd, cxr));
				if (variates.contains(Variate.WeightsOfLines))
					res.add(new MatrixCorrectionModel2.XRayWeightLabel(cxr));
			}
			for (final AtomicShell sh : exrs.getSetOfInnerAtomicShells()) {
				res.add(buildICXLabel(mcd, sh));
				if (variates.contains(Variate.SecondaryFluorescence))
					res.add(new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mcd, sh));
			}
			return res;
		}

		private static List<? extends Object> buildOutputTags(//
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs) {
			return Collections.singletonList(intensityLabel(mcd, exrs));
		}

		public IntensityModel( //
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs, //
				final Set<Variate> variates //
		) {
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
				final Object ionLbl = buildICXLabel(mDatum, shell);
				final Object secLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell);
				final double ionVal = getValue(ionLbl, point);
				final double secVal = getValue(secLbl, point);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final Object wLbl = new MatrixCorrectionModel2.XRayWeightLabel(cxr);
					final Object fxLbl = MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr);
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
				final Object ionLbl = buildICXLabel(mDatum, shell);
				final Object secLbl = new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell);
				final double ionVal = getValue(ionLbl, point);
				final double secVal = getValue(secLbl, point);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final Object wLbl = new MatrixCorrectionModel2.XRayWeightLabel(cxr);
					final Object fxLbl = MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr);
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

		private final KRatioLabel mKRatio;

		private static List<? extends Object> buildInputLabels(final KRatioLabel krl) {
			final List<Object> res = new ArrayList<>();
			res.add(intensityLabel(krl.getUnknown(), krl.getXRaySet()));
			res.add(intensityLabel(krl.getStandard(), krl.getXRaySet()));
			res.add(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(), krl.getXRaySet().getElement()));
			res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(),
					krl.getXRaySet().getElement()));
			return res;
		}

		private static List<? extends Object> buildOutputLabels(final KRatioLabel krl) {
			final List<Object> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.zafLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet()));
			res.add(new KRatioLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet(), Method.Calculated));
			return res;
		}

		public kRatioModel(final KRatioLabel krl) {
			super(buildInputLabels(krl), buildOutputLabels(krl));
			mKRatio = krl;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			final UnknownMatrixCorrectionDatum unk = mKRatio.getUnknown();
			final StandardMatrixCorrectionDatum std = mKRatio.getStandard();
			final ElementXRaySet exrs = mKRatio.getXRaySet();
			final Object zafLabel = MatrixCorrectionModel2.zafLabel(mKRatio);
			final Object kRatioLabel = new KRatioLabel(unk, std, exrs, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			final int zafRow = outputIndex(zafLabel);

			final Object intUnkLbl = intensityLabel(unk, exrs);
			final Object intStdLbl = intensityLabel(std, exrs);
			final Object cUnkLbl = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final Object cStdLbl = MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement());

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
			final UnknownMatrixCorrectionDatum unk = mKRatio.getUnknown();
			final StandardMatrixCorrectionDatum std = mKRatio.getStandard();
			final ElementXRaySet exrs = mKRatio.getXRaySet();
			final Object zafLabel = MatrixCorrectionModel2.zafLabel(unk, std, exrs);
			final Object kRatioLabel = new KRatioLabel(unk, std, exrs, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			final int zafRow = outputIndex(zafLabel);

			final Object intUnkLbl = intensityLabel(unk, exrs);
			final Object intStdLbl = intensityLabel(std, exrs);
			final Object cUnkLbl = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final Object cStdLbl = MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement());

			final double intUnk = getValue(intUnkLbl, point);
			final double intStd = getValue(intStdLbl, point);
			final double cUnk = getValue(cUnkLbl, point);
			final double cStd = getValue(cStdLbl, point);

			rv.setEntry(kRow, (cUnk * intUnk) / (cStd * intStd));
			rv.setEntry(zafRow, intUnk / intStd);

			return rv;
		}

		public String toString() {
			return "kRatio[" + mKRatio + "]";
		}

	}

	private static class IonizationCrossSection //
			extends LabeledMultivariateJacobianFunction //
			implements ILabeledMultivariateFunction {

		private final AtomicShell mShell;
		private final MatrixCorrectionDatum mDatum;

		static private List<Object> buildInputLabels(final MatrixCorrectionDatum mcd, final AtomicShell sh,
				Set<Variate> variates) {
			final List<Object> res = new ArrayList<>();
			if (variates.contains(Variate.IonizationExponent))
				res.add(new MatrixCorrectionModel2.IonizationExponentLabel(sh));
			if (variates.contains(Variate.BeamEnergy))
				res.add(MatrixCorrectionModel2.beamEnergyLabel(mcd));
			return res;
		}

		static private List<Object> buildOutputLabels(final MatrixCorrectionDatum mcd, final AtomicShell sh) {
			return Collections.singletonList(buildICXLabel(mcd, sh));
		}

		public IonizationCrossSection(final MatrixCorrectionDatum mcd, final AtomicShell sh, Set<Variate> variates) {
			super(buildInputLabels(mcd, sh, variates), buildOutputLabels(mcd, sh));
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
			assert outputIndex(buildICXLabel(mDatum, mShell)) == 0;

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
			assert outputIndex(buildICXLabel(mDatum, mShell)) == 0;
			rv.setEntry(0, icx);
			return rv;
		}

		public String toString() {
			return "ICX2[" + mDatum + "," + mShell + "]";
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
	 * @param withUnc
	 * @return {@link UncertainValues}
	 */
	private UncertainValues buildInputs(final boolean withUnc) {
		final Map<Object, Number> vals = new HashMap<>();
		XRayEmissionProbability xrep = null;
		for (final KRatioLabel krl : mKRatios) {
			final ElementXRaySet exrs = krl.getXRaySet();
			for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				// MatrixCorrectionModel2.FofChiReducedLabel(mcd, cxr) => Calculated
				// elsewhere...
				if (mVariates.contains(MatrixCorrectionModel2.Variate.WeightsOfLines) == withUnc) {
					if ((xrep == null) || (xrep.getIonized() != cxr.getInner()))
						xrep = new XRayEmissionProbability(cxr.getInner());
					vals.put(new MatrixCorrectionModel2.XRayWeightLabel(cxr), xrep.getWeightU(cxr));
				}
			}
			if (mVariates.contains(MatrixCorrectionModel2.Variate.IonizationExponent) == withUnc)
				for (final AtomicShell sh : exrs.getSetOfInnerAtomicShells())
					vals.put(new MatrixCorrectionModel2.IonizationExponentLabel(sh),
							XPPMatrixCorrection2.computeIonizationExponent(sh, 1.0));
		}
		return new UncertainValues(vals);
	}

	public UncertainValues buildInputs() {
		return buildInputs(true);
	}

	public UncertainValues buildConstants() {
		return buildInputs(false);
	}

}