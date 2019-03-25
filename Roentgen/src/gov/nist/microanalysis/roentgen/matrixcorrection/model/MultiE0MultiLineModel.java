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

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ILabeledMultivariateFunction;
import gov.nist.juncertainty.ParallelMeasurementModelBuilder;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.ZAFMultiLineLabel;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

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
		extends CompositeMeasurementModel<EPMALabel> {

	private static class IntensityModel //
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> {

		private static List<EPMALabel> buildInputTags(
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs //
		) {
			final List<EPMALabel> res = new ArrayList<>();
			for (final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay()) {
				res.add(MatrixCorrectionModel2.FofChiReducedLabel(mcd, cxr));
				res.add(new MatrixCorrectionModel2.XRayWeightLabel(cxr));
			}
			for (final AtomicShell sh : exrs.getSetOfInnerAtomicShells()) {
				res.add(buildICXLabel(mcd, sh));
				res.add(new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mcd, sh));
			}
			return res;
		}

		private static List<EPMALabel> buildOutputTags(
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs
		) {
			return Collections.singletonList(intensityLabel(mcd, exrs));
		}

		private final MatrixCorrectionDatum mDatum;

		private final ElementXRaySet mXRaySet;

		public IntensityModel(
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs //
		) throws ArgumentException {
			super(buildInputTags(mcd, exrs), buildOutputTags(mcd, exrs));
			mDatum = mcd;
			mXRaySet = exrs;
		}

		@Override
		public RealVector optimized(
				final RealVector point
		) {
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
				final int ionLbl = inputIndex(buildICXLabel(mDatum, shell));
				final int secLbl = inputIndex( //
						new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell));
				final double ionVal = point.getEntry(ionLbl);
				final double secVal = point.getEntry(secLbl);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final int wLbl = inputIndex(new MatrixCorrectionModel2.XRayWeightLabel(cxr));
					final int fxLbl = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr));
					final double wVal = point.getEntry(wLbl);
					final double fxVal = point.getEntry(fxLbl);
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

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {

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
				final int ionIdx = inputIndex(buildICXLabel(mDatum, shell));
				final int secIdx = inputIndex( //
						new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell));
				final double ionVal = point.getEntry(ionIdx);
				final double secVal = point.getEntry(secIdx);
				final double outer = ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final int wLbl = inputIndex(new MatrixCorrectionModel2.XRayWeightLabel(cxr));
					final int fxLbl = inputIndex(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr));
					final double wVal = point.getEntry(wLbl);
					final double fxVal = point.getEntry(fxLbl);
					inner += wVal * fxVal;
					rm.setEntry(intIdx, fxLbl, wVal * outer);
					rm.setEntry(intIdx, wLbl, fxVal * outer);
				}
				rm.setEntry(intIdx, ionIdx, secVal * inner);
				rm.setEntry(intIdx, secIdx, ionVal * inner);
				res += outer * inner;
			}
			rv.setEntry(intIdx, res);
			return Pair.create(rv, rm);
		}
	}

	private static class IonizationCrossSection //
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> {

		static private List<EPMALabel> buildInputLabels(
				final MatrixCorrectionDatum mcd, final AtomicShell sh //
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(new MatrixCorrectionModel2.IonizationExponentLabel(sh));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(mcd));
			return res;
		}

		static private List<EPMALabel> buildOutputLabels(
				final MatrixCorrectionDatum mcd, final AtomicShell sh
		) {
			return Collections.singletonList(buildICXLabel(mcd, sh));
		}

		private final AtomicShell mShell;

		private final MatrixCorrectionDatum mDatum;

		public IonizationCrossSection(
				final MatrixCorrectionDatum mcd, //
				final AtomicShell sh
		) throws ArgumentException {
			super(buildInputLabels(mcd, sh), buildOutputLabels(mcd, sh));
			mShell = sh;
			mDatum = mcd;
		}

		@Override
		public RealVector optimized(
				final RealVector point
		) {
			assert point.getDimension() == getInputDimension();
			final double eL = 1.0e-3 * mShell.getEdgeEnergy();
			final int e0Idx = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final double e0 = point.getEntry(e0Idx);
			final int mIdx = inputIndex(new MatrixCorrectionModel2.IonizationExponentLabel(mShell));
			final double m = point.getEntry(mIdx);
			final double u = e0 / eL;
			final RealVector rv = new ArrayRealVector(1);
			final double icx = Math.log(u) / (Math.pow(u, m) * eL * eL);
			assert outputIndex(buildICXLabel(mDatum, mShell)) == 0;
			rv.setEntry(0, icx);
			return rv;
		}

		@Override
		public String toString() {
			return "ICX[" + mDatum + "," + mShell + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			assert point.getDimension() == getInputDimension();
			assert getOutputDimension() == 1;
			final double eL = 1.0e-3 * mShell.getEdgeEnergy();
			final int e0Idx = inputIndex(MatrixCorrectionModel2.beamEnergyLabel(mDatum));
			final double e0 = point.getEntry(e0Idx);
			final int mIdx = inputIndex(new MatrixCorrectionModel2.IonizationExponentLabel(mShell));
			final double m = point.getEntry(mIdx);
			assert outputIndex(buildICXLabel(mDatum, mShell)) == 0;

			final double u = e0 / eL;
			assert u > 1.0 : "The overvoltage for " + mDatum.toString() + " is less than unity.";
			final RealVector rv = new ArrayRealVector(1);
			final double logU = Math.log(u);
			final double icx = logU / (Math.pow(u, m) * eL * eL);
			rv.setEntry(0, icx);
			final RealMatrix rm = NullableRealMatrix.build(1, getInputDimension());
			rm.setEntry(0, e0Idx, (1.0 - m * logU) / (Math.pow(u, m + 1.0) * Math.pow(eL, 3.0)));
			rm.setEntry(0, mIdx, -icx * logU);
			return Pair.create(rv, rm);
		}

	}

	private static class KRatioZAFModel //
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> {

		private static List<EPMALabel> buildInputLabels(
				final KRatioLabel krl
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(intensityLabel(krl.getUnknown(), krl.getXRaySet()));
			res.add(intensityLabel(krl.getStandard(), krl.getXRaySet()));
			res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
			return res;
		}

		private static List<EPMALabel> buildOutputLabels(
				final KRatioLabel krl
		) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.zafLabel(krl));
			res.add(krl.asCalculated());
			return res;
		}

		private final KRatioLabel mKRatio;

		public KRatioZAFModel(
				final KRatioLabel krl
		) throws ArgumentException {
			super(buildInputLabels(krl), buildOutputLabels(krl));
			mKRatio = krl;
		}

		@Override
		public RealVector optimized(
				final RealVector point
		) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final UnknownMatrixCorrectionDatum unk = mKRatio.getUnknown();
			final StandardMatrixCorrectionDatum std = mKRatio.getStandard();
			final ElementXRaySet exrs = mKRatio.getXRaySet();
			final Object zafLabel = MatrixCorrectionModel2.zafLabel(mKRatio);
			final Object kRatioLabel = new KRatioLabel(unk, std, exrs, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			final int zafRow = outputIndex(zafLabel);

			final int intUnkIdx = inputIndex(intensityLabel(unk, exrs));
			final int intStdIdx = inputIndex(intensityLabel(std, exrs));
			final MassFraction lUnk = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final int cStdIdx = inputIndex( //
					MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement()));

			final double intUnk = point.getEntry(intUnkIdx);
			final double intStd = point.getEntry(intStdIdx);
			final double cUnk = getArg(lUnk, point);
			final double cStd = point.getEntry(cStdIdx);

			rv.setEntry(kRow, (cUnk * intUnk) / (cStd * intStd));
			rv.setEntry(zafRow, intUnk / intStd);

			return rv;
		}

		@Override
		public String toString() {
			return "kRatio/ZAF[" + mKRatio + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final RealVector vals = buildResult();
			final RealMatrix jac = buildJacobian();
			final UnknownMatrixCorrectionDatum unk = mKRatio.getUnknown();
			final StandardMatrixCorrectionDatum std = mKRatio.getStandard();
			final ElementXRaySet exrs = mKRatio.getXRaySet();
			final ZAFMultiLineLabel zafLabel = MatrixCorrectionModel2.zafLabel(mKRatio);
			final KRatioLabel kRatioLabel = new KRatioLabel(unk, std, exrs, Method.Calculated);
			final int kRow = outputIndex(kRatioLabel);
			final int zafRow = outputIndex(zafLabel);

			final int intUnkIdx = inputIndex(intensityLabel(unk, exrs));
			final int intStdIdx = inputIndex(intensityLabel(std, exrs));
			final MassFraction lUnk = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final int cStdIdx = inputIndex( //
					MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement()));

			final double intUnk = point.getEntry(intUnkIdx);
			final double intStd = point.getEntry(intStdIdx);
			final double cUnk = getArg(lUnk, point);
			final double cStd = point.getEntry(cStdIdx);

			final double k = (cUnk * intUnk) / (cStd * intStd);
			vals.setEntry(kRow, k);
			jac.setEntry(kRow, intUnkIdx, k / intUnk);
			jac.setEntry(kRow, intStdIdx, -k / intStd);
			setJacobian(lUnk, kRatioLabel, jac, k / cUnk);
			jac.setEntry(kRow, cStdIdx, -k / cStd);

			vals.setEntry(zafRow, intUnk / intStd);
			jac.setEntry(zafRow, intUnkIdx, 1.0 / intStd);
			jac.setEntry(zafRow, intStdIdx, -1.0 * intUnk / (intStd * intStd));

			return Pair.create(vals, jac);
		}

	}

	public static EPMALabel buildICXLabel(
			final MatrixCorrectionDatum mcd, final AtomicShell sh
	) {
		return new EPMALabel.BaseLabel<MatrixCorrectionDatum, AtomicShell, Object>("ICX", mcd, sh);
	}

	private static List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps(
			//
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		final Map<MatrixCorrectionDatum, Set<AtomicShell>> allMcd = new HashMap<>();

		for (final KRatioLabel krl : kratios) {
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

		final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> funcs = new ArrayList<>();
			for (final Map.Entry<MatrixCorrectionDatum, Set<AtomicShell>> me : allMcd.entrySet())
				for (final AtomicShell sh : me.getValue())
					funcs.add(new IonizationCrossSection(me.getKey(), sh));
			res.add(ParallelMeasurementModelBuilder.join("ICXs", funcs));
		}
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				final ElementXRaySet cxrs = krl.getXRaySet();
				funcs.add(new IntensityModel(krl.getStandard(), cxrs));
				funcs.add(new IntensityModel(krl.getUnknown(), cxrs));
			}
			res.add(ParallelMeasurementModelBuilder.join("Intensities", funcs));
		}
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				funcs.add(new KRatioZAFModel(krl));
			res.add(ParallelMeasurementModelBuilder.join("K-ratios", funcs));
		}
		return res;
	}

	private static EPMALabel intensityLabel(
			final MatrixCorrectionDatum mcd, final ElementXRaySet exrs
	) {
		return new EPMALabel.BaseLabel<MatrixCorrectionDatum, ElementXRaySet, Object>("Intensity", mcd, exrs);
	}

	public MultiE0MultiLineModel(
			//
			final Set<KRatioLabel> kratios //
	) throws ArgumentException {
		super("Multi-E0", buildSteps(kratios));
	}
}