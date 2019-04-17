package gov.nist.microanalysis.roentgen.matrixcorrection.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ILabeledMultivariateFunction;
import gov.nist.juncertainty.ParallelMeasurementModelBuilder;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.MatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.StandardMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.UnknownMatrixCorrectionDatum;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.FofChiReducedLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.EmittedIntensityLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.MatrixCorrectionDatumLabel;
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

	/**
	 * Computes the emitted intensity for the specified material and set of
	 * characteristic X-rays.
	 * 
	 * @author Nicholas W. M. Ritchie
	 *
	 */
	private static class EmittedIntensityModel //
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
			if (mcd instanceof StandardMatrixCorrectionDatum)
				res.add(MaterialLabel.buildMassFractionTag(mcd.getMaterial(), exrs.getElement()));
			return res;
		}

		private static List<EPMALabel> buildOutputTags(
				final MatrixCorrectionDatum mcd, //
				final ElementXRaySet exrs
		) {
			return Collections.singletonList(MatrixCorrectionModel2.buildEmittedIntensityLabel(mcd, exrs));
		}

		private final MatrixCorrectionDatum mDatum;

		private final ElementXRaySet mXRaySet;

		public EmittedIntensityModel(
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
			final RealVector rv = buildResult();
			final EmittedIntensityLabel intIdx = MatrixCorrectionModel2.buildEmittedIntensityLabel(mDatum, mXRaySet);
			final MassFraction mfLbl = MaterialLabel.buildMassFractionTag(mDatum.getMaterial(), mXRaySet.getElement());
			final Map<AtomicShell, Set<CharacteristicXRay>> shells = new HashMap<>();
			// Split the characteristic x-rays by shells
			for (final CharacteristicXRay cxr : mXRaySet.getSetOfCharacteristicXRay()) {
				final AtomicShell shell = cxr.getInner();
				if (!shells.containsKey(shell))
					shells.put(shell, new HashSet<>());
				shells.get(shell).add(cxr);
			}
			final double c = getArg(mfLbl, point);
			double res = 0.0;
			for (final Map.Entry<AtomicShell, Set<CharacteristicXRay>> me : shells.entrySet()) {
				final AtomicShell shell = me.getKey();
				final double ionVal = getArg(buildICXLabel(mDatum, shell), point);
				final double secVal = getArg(new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell),
						point);
				final double outer = ionVal * secVal * c;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final double wVal = getArg(new MatrixCorrectionModel2.XRayWeightLabel(cxr), point);
					final double fxVal = getArg(MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr), point);
					inner += wVal * fxVal;
				}
				res += outer * inner;
			}
			setResult(intIdx, rv, res);
			return rv;
		}

		@Override
		public String toString() {
			return "I<sub>emitted</sub>[" + mDatum.toString() + "," + mXRaySet.toString() + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final RealVector rv = buildResult();
			final RealMatrix rm = buildJacobian();

			final EmittedIntensityLabel intIdx = MatrixCorrectionModel2.buildEmittedIntensityLabel(mDatum, mXRaySet);
			final Map<AtomicShell, Set<CharacteristicXRay>> shells = new HashMap<>();
			final MassFraction mfLbl = MaterialLabel.buildMassFractionTag(mDatum.getMaterial(), mXRaySet.getElement());
			// Split the characteristic x-rays by shells
			for (final CharacteristicXRay cxr : mXRaySet.getSetOfCharacteristicXRay()) {
				final AtomicShell shell = cxr.getInner();
				if (!shells.containsKey(shell))
					shells.put(shell, new HashSet<>());
				shells.get(shell).add(cxr);
			}
			final double c = getArg(mfLbl, point);
			double res = 0.0;
			for (final Map.Entry<AtomicShell, Set<CharacteristicXRay>> me : shells.entrySet()) {
				final AtomicShell shell = me.getKey();
				final EPMALabel ionIdx = buildICXLabel(mDatum, shell);
				final SecondaryFluorescenceModel.SecondaryFluorescenceLabel secIdx = //
						new SecondaryFluorescenceModel.SecondaryFluorescenceLabel(mDatum, shell);
				final double ionVal = getArg(ionIdx, point);
				final double secVal = getArg(secIdx, point);
				final double outer = c * ionVal * secVal;
				double inner = 0.0;
				for (final CharacteristicXRay cxr : me.getValue()) {
					final MatrixCorrectionModel2.XRayWeightLabel wLbl = new MatrixCorrectionModel2.XRayWeightLabel(cxr);
					final FofChiReducedLabel fxLbl = MatrixCorrectionModel2.FofChiReducedLabel(mDatum, cxr);
					final double wVal = getArg(wLbl, point);
					final double fxVal = getArg(fxLbl, point);
					inner += wVal * fxVal;
					setJacobian(fxLbl, intIdx, rm, wVal * outer);
					setJacobian(wLbl, intIdx, rm, fxVal * outer);
				}
				setJacobian(ionIdx, intIdx, rm, secVal * inner);
				setJacobian(secIdx, intIdx, rm, ionVal * inner);
				res += outer * inner;
			}
			setResult(intIdx, rv, res);
			return Pair.create(rv, rm);
		}
	}

	private static class IonizationCrossSection //
			extends ExplicitMeasurementModel<EPMALabel, EPMALabel> //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> {

		static private List<EPMALabel> buildInputLabels(
				final MatrixCorrectionDatum mcd, //
				final Set<AtomicShell> shs //
		) {
			final List<EPMALabel> res = new ArrayList<>();
			for (AtomicShell sh : shs)
				res.add(new MatrixCorrectionModel2.IonizationExponentLabel(sh));
			res.add(MatrixCorrectionModel2.beamEnergyLabel(mcd));
			return res;
		}

		static private List<EPMALabel> buildOutputLabels(
				final MatrixCorrectionDatum mcd, //
				final Set<AtomicShell> shs
		) {
			List<EPMALabel> res = new ArrayList<>();
			for (AtomicShell sh : shs)
				res.add(buildICXLabel(mcd, sh));
			return res;
		}

		private final Set<AtomicShell> mShells;
		private final MatrixCorrectionDatum mDatum;

		public IonizationCrossSection(
				final MatrixCorrectionDatum mcd, //
				final Set<AtomicShell> ssh
		) throws ArgumentException {
			super(buildInputLabels(mcd, ssh), buildOutputLabels(mcd, ssh));
			mShells = ssh;
			mDatum = mcd;
		}

		@Override
		public RealVector optimized(
				final RealVector point
		) {
			final RealVector rv = buildResult();
			final MatrixCorrectionDatumLabel e0Idx = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			for (AtomicShell sh : mShells) {
				assert point.getDimension() == getInputDimension();
				final double eL = 1.0e-3 * sh.getEdgeEnergy();
				final double e0 = getArg(e0Idx, point);
				final MatrixCorrectionModel2.IonizationExponentLabel mIdx = new MatrixCorrectionModel2.IonizationExponentLabel(
						sh);
				final double m = getArg(mIdx, point);
				final double u = e0 / eL;
				final double icx = Math.log(u) / (Math.pow(u, m) * eL * eL);
				setResult(buildICXLabel(mDatum, sh), rv, icx);
			}
			return rv;
		}

		@Override
		public String toString() {
			return "ICX[" + mDatum + "," + mShells + "]";
		}

		@Override
		public Pair<RealVector, RealMatrix> value(
				final RealVector point
		) {
			final RealVector rv = buildResult();
			final RealMatrix rm = buildJacobian();
			final MatrixCorrectionDatumLabel e0Idx = MatrixCorrectionModel2.beamEnergyLabel(mDatum);
			for (AtomicShell sh : mShells) {
				final double eL = 1.0e-3 * sh.getEdgeEnergy();
				final double e0 = getArg(e0Idx, point);
				final MatrixCorrectionModel2.IonizationExponentLabel mIdx = new MatrixCorrectionModel2.IonizationExponentLabel(
						sh);
				final double m = getArg(mIdx, point);
				final EPMALabel icxLbl = buildICXLabel(mDatum, sh);
				final double u = e0 / eL;
				assert u > 1.0 : "The overvoltage for " + mDatum.toString() + " is less than unity.";
				final double logU = Math.log(u);
				final double icx = logU / (Math.pow(u, m) * eL * eL);
				setResult(icxLbl, rv, icx);
				setJacobian(e0Idx, icxLbl, rm, (1.0 - m * logU) / (Math.pow(u, m + 1.0) * Math.pow(eL, 3.0)));
				setJacobian(mIdx, icxLbl, rm, -icx * logU);
			}
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
			res.add(MatrixCorrectionModel2.buildEmittedIntensityLabel(krl.getUnknown(), krl.getXRaySet()));
			res.add(MatrixCorrectionModel2.buildEmittedIntensityLabel(krl.getStandard(), krl.getXRaySet()));
			res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
			// res.add(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(),
			// krl.getElement()));
			return res;
		}

		private static List<EPMALabel> buildOutputLabels(
				final KRatioLabel krl
				) {
			final List<EPMALabel> res = new ArrayList<>();
			res.add(MatrixCorrectionModel2.zafLabel(krl));
			res.add(krl.as(Method.Calculated));
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
			final RealVector rv = buildResult();
			final UnknownMatrixCorrectionDatum unk = mKRatio.getUnknown();
			final StandardMatrixCorrectionDatum std = mKRatio.getStandard();
			final ElementXRaySet exrs = mKRatio.getXRaySet();

			final ZAFMultiLineLabel zafLabel = MatrixCorrectionModel2.zafLabel(mKRatio);
			final KRatioLabel kRatioLabel = mKRatio.as(Method.Calculated);

			final EmittedIntensityLabel intUnkLbl = MatrixCorrectionModel2.buildEmittedIntensityLabel(unk, exrs);
			final EmittedIntensityLabel intStdLbl = MatrixCorrectionModel2.buildEmittedIntensityLabel(std, exrs);
			final MassFraction cUnkLbl = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final MassFraction cStdLbl = MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement());

			final double intUnk = getArg(intUnkLbl, point);
			final double intStd = getArg(intStdLbl, point);
			final double cUnk = getArg(cUnkLbl, point);
			final double cStd = getArg(cStdLbl, point);

			setResult(zafLabel, rv, (cStd * intUnk) / (cUnk * intStd));
			setResult(kRatioLabel, rv, intUnk / intStd);

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
			final KRatioLabel kRatioLabel = mKRatio.as(Method.Calculated);

			final EmittedIntensityLabel intUnkLbl = MatrixCorrectionModel2.buildEmittedIntensityLabel(unk, exrs);
			final EmittedIntensityLabel intStdLbl = MatrixCorrectionModel2.buildEmittedIntensityLabel(std, exrs);
			final MassFraction cUnkLbl = MaterialLabel.buildMassFractionTag(unk.getMaterial(), exrs.getElement());
			final MassFraction cStdLbl = MaterialLabel.buildMassFractionTag(std.getMaterial(), exrs.getElement());

			final double intUnk = getArg(intUnkLbl, point);
			final double intStd = getArg(intStdLbl, point);
			final double cUnk = getArg(cUnkLbl, point);
			final double cStd = getArg(cStdLbl, point);

			final double k = intUnk / intStd;
			setResult(kRatioLabel, vals, k);
			setJacobian(intUnkLbl, kRatioLabel, jac, 1.0 / intStd);
			setJacobian(intStdLbl, kRatioLabel, jac, -intUnk / (intStd * intStd));

			final double zaf = (cStd * intUnk) / (cUnk * intStd);
			setResult(zafLabel, vals, zaf);
			setJacobian(intUnkLbl, zafLabel, jac, zaf / intUnk);
			setJacobian(intStdLbl, zafLabel, jac, -zaf / intStd);
			setJacobian(cUnkLbl, zafLabel, jac, -zaf / cUnk);
			setJacobian(cStdLbl, zafLabel, jac, zaf / cStd);

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
				funcs.add(new IonizationCrossSection(me.getKey(), me.getValue()));
			res.add(ParallelMeasurementModelBuilder.join("ICXs", funcs));
		}
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				final ElementXRaySet cxrs = krl.getXRaySet();
				funcs.add(new EmittedIntensityModel(krl.getStandard(), cxrs));
				funcs.add(new EmittedIntensityModel(krl.getUnknown(), cxrs));
			}
			res.add(ParallelMeasurementModelBuilder.join("I<sub>emitted</sub>", funcs));
		}
		{
			final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> funcs = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				funcs.add(new KRatioZAFModel(krl));
			res.add(ParallelMeasurementModelBuilder.join("K-ratios", funcs));
		}
		return res;
	}

	public MultiE0MultiLineModel(
			//
			final Set<KRatioLabel> kratios //
			) throws ArgumentException {
		super("Multi-E0", buildSteps(kratios));
	}
}