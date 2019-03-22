package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.CompositeLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2.ZAFMultiLineLabel;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;
import gov.nist.microanalysis.roentgen.utility.FastIndex;

/**
 * <p>
 * Implements the implicit measurement model that defines the fundamental model
 * behind k-ratio-based compositional measurement - <i>k<sub>meas,z</sub> -
 * k<sub>calc,z</sub>(<b>C</b>)</i> = 0 for all <i>z</i>.
 * </p>
 *
 * <p>
 * Since the <i>k<sub>calc,z</sub>(<b>C</b>)</i> are not solvable for
 * <b><i>C</i></b>, it is necessary to implement a scheme for determining the
 * <i>C<sub>z</sub></i> that minimizes <i>k<sub>meas,z</sub> -
 * k<sub>calc,z</sub>(<b>C</b>)</i>. This non-linear optimization can be done
 * various different ways - classically the process is called "iteration" in
 * microanalysis and there are a handful of "iteration algorithms." Including
 * "simple iteration", "alpha-factor", "PAP iteration", "Wegstein iteration"
 * although it is also possible to use other non-linear optimization algorithms.
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class KRatioCorrectionModel2 //
		extends CompositeLabeledMultivariateJacobianFunction<EPMALabel> {

	static private final int MAX_ITERATIONS = 100;
	static private final double THRESH = 1.0e-4;

	private final MatrixCorrectionModel2 mModel;

	private final List<KRatioLabel> mKRatioSet;

	private final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> mExtraElementRules;

	private static class KR2HModel //
			extends ImplicitMeasurementModel.HModel<EPMALabel> //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> //
	{

		private final List<KRatioLabel> mKRatios;

		static List<EPMALabel> buildOutputs(//
				final List<KRatioLabel> kratios //
		) throws ArgumentException {
			Material unk = null;
			final Set<Element> elms = new TreeSet<Element>();
			for (final KRatioLabel krl : kratios) {
				if (unk == null)
					unk = krl.getUnknown().getMaterial();
				if (!unk.equals(krl.getUnknown().getMaterial()))
					throw new ArgumentException("More than one unknown in KRatioCorrectionModel2.HModel");
				elms.add(krl.getElement());
			}
			final List<EPMALabel> res = new ArrayList<>();
			for (final Element elm : elms)
				res.add(MaterialLabel.buildMassFractionTag(unk, elm));
			return res;
		}

		static List<EPMALabel> buildInputs(//
				final List<KRatioLabel> kratios //
		) {
			final Set<EPMALabel> res = new HashSet<>();
			for (final KRatioLabel krl : kratios) {
				res.add(krl);
				res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
				res.add(MatrixCorrectionModel2.zafLabel(krl));
			}
			return new ArrayList<>(res);
		}

		static List<EPMALabel> buildHLabels(//
				final List<KRatioLabel> kratios //
		) {
			final Set<EPMALabel> res = new HashSet<>();
			for (int i = 0; i < kratios.size(); ++i)
				res.add(new EPMALabel("H[" + i + "]"));
			return new ArrayList<>(res);
		}

		/**
		 * <p>
		 * Builds a k-ratio H model in which the inputs are the standard's mass
		 * fraction, matrix correction factors and the k-ratio for each measured element
		 * and the outputs are the unknown's mass fraction for each element.
		 * </p>
		 * <p>
		 * Because this is an H-model for an {@link ImplicitMeasurementModel}, the
		 * unknown's mass fractions are input as constants.
		 * </p>
		 * <p>
		 * The outputs of the {@link LabeledMultivariateJacobianFunction} can only be
		 * MassFraction values. (Might be able to relax this but why???)
		 * </p>
		 *
		 *
		 * @param kratios
		 * @param otherElms A model for calculating unmeasured elements
		 * @throws ArgumentException
		 */
		private KR2HModel(final List<KRatioLabel> kratios) //
				throws ArgumentException {
			super(buildInputs(kratios), buildOutputs(kratios), buildHLabels(kratios));
			mKRatios = new ArrayList<>(kratios);
		}

		/**
		 * <p>
		 * Builds a k-ratio H model in which the inputs are the standard's mass
		 * fraction, matrix correction factors and the k-ratio for each measured element
		 * and the outputs are the unknown's mass fraction for each element.
		 * </p>
		 * <p>
		 * Because this is an H-model for an {@link ImplicitMeasurementModel}, the
		 * unknown's mass fractions are input as constants.
		 * </p>
		 *
		 * @param kratios
		 * @throws ArgumentException
		 */
		/*
		 * private KR2HModel(final Set<KRatioLabel> kratios) // throws ArgumentException
		 * { super(buildInputs(kratios, null), buildOutputs(kratios, null)); mKRatios =
		 * new HashSet<>(kratios); mOtherElements = null; }
		 */

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			for (int oHTag = 0; oHTag < mKRatios.size(); ++oHTag) {
				final KRatioLabel kMeasTag = mKRatios.get(oHTag);
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final Material unkMat = unknown.getMaterial();
				final Material stdMat = standard.getMaterial();
				final MaterialLabel mfUnkTag = MaterialLabel.buildMassFractionTag(unkMat, elm);
				final MaterialLabel mfStdTag = MaterialLabel.buildMassFractionTag(stdMat, elm);
				final ZAFMultiLineLabel zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				final int iKMeas = inputIndex(kMeasTag);
				final int iMFUnk = inputIndex(mfUnkTag);
				final int iMFStd = inputIndex(mfStdTag);
				final int iZAF = inputIndex(zafTag);
				final double kMeas = point.getEntry(iKMeas);
				final double cUnk = point.getEntry(iMFUnk);
				final double cStd = point.getEntry(iMFStd);
				final double zaf = point.getEntry(iZAF);
				final double hi = kMeas - (cUnk / cStd) * zaf;
				// assert Math.abs(hi * kMeas) < 1.0e-4 : "h[" + kMeasTag + "] = " + hi;
				rv.setEntry(oHTag, hi);
				rm.setEntry(oHTag, iKMeas, 1.0);
				rm.setEntry(oHTag, iMFUnk, (-1.0 / cStd) * zaf);
				rm.setEntry(oHTag, iMFStd, (cUnk / (cStd * cStd)) * zaf);
				rm.setEntry(oHTag, iZAF, -(cUnk / cStd));
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			for (int oHTag = 0; oHTag < mKRatios.size(); ++oHTag) {
				final KRatioLabel kMeasTag = mKRatios.get(oHTag);
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final Material unkMat = unknown.getMaterial();
				final Material stdMat = standard.getMaterial();
				final MassFraction mfUnkTag = MaterialLabel.buildMassFractionTag(unkMat, elm);
				final MassFraction mfStdTag = MaterialLabel.buildMassFractionTag(stdMat, elm);
				final ZAFMultiLineLabel zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				final int iKMeas = inputIndex(kMeasTag);
				final int iMFUnk = inputIndex(mfUnkTag);
				final int iMFStd = inputIndex(mfStdTag);
				final int iZAF = inputIndex(zafTag);
				final double kMeas = point.getEntry(iKMeas);
				final double cUnk = point.getEntry(iMFUnk);
				final double cStd = point.getEntry(iMFStd);
				final double zaf = point.getEntry(iZAF);
				final double hi = kMeas - (cUnk / cStd) * zaf;
				// assert Math.abs(hi * kMeas) < 1.0e-4 : "h[" + kMeasTag + "] = " + hi;
				rv.setEntry(oHTag, hi);
			}
			return rv;
		}

		@Override
		public String toString() {
			return "h[k[meas]-k[calc]]";
		}
	}

	private static class AltComputeCUnknown //
			implements ILabeledMultivariateFunction<EPMALabel, EPMALabel> {

		private final Set<KRatioLabel> mKRatios;
		private final List<EPMALabel> mInputLabels;
		private final List<EPMALabel> mOutputLabels;

		/**
		 * <p>
		 * An explicit alternative to the implicit {@link KR2HModel}.
		 * </p>
		 *
		 *
		 * @param kratios
		 * @param otherElms A model for calculating unmeasured elements
		 * @throws ArgumentException
		 */
		private AltComputeCUnknown(final List<KRatioLabel> kratios) //
				throws ArgumentException {
			mKRatios = new HashSet<>(kratios);
			mInputLabels = KR2HModel.buildInputs(kratios);
			mOutputLabels = KR2HModel.buildOutputs(kratios);
		}

		@Override
		public List<EPMALabel> getInputLabels() {
			return mInputLabels;
		}

		@Override
		public List<EPMALabel> getOutputLabels() {
			return mOutputLabels;
		}

		@Override
		public RealVector optimized(final RealVector point) {
			final RealVector rv = new ArrayRealVector(mOutputLabels.size());
			for (final KRatioLabel kMeasTag : mKRatios) {
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final Material unkMat = unknown.getMaterial();
				final Material stdMat = standard.getMaterial();
				final Object mfUnkTag = MaterialLabel.buildMassFractionTag(unkMat, elm);
				final Object mfStdTag = MaterialLabel.buildMassFractionTag(stdMat, elm);
				final Object zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				final int iKMeas = mInputLabels.indexOf(kMeasTag);
				final int iMFStd = mInputLabels.indexOf(mfStdTag);
				final int iZAF = mInputLabels.indexOf(zafTag);
				final int iMFUnk = mOutputLabels.indexOf(mfUnkTag);
				final double kMeas = point.getEntry(iKMeas);
				final double cStd = point.getEntry(iMFStd);
				final double zaf = point.getEntry(iZAF);
				final double cUnk = cStd * (kMeas / zaf);
				rv.setEntry(iMFUnk, cUnk);
			}
			return rv;
		}

		@Override
		public String toString() {
			return "Alt[C[unknown]]";
		}

	}

	private static List<LabeledMultivariateJacobianFunction<? extends EPMALabel, ? extends EPMALabel>> buildSteps(//
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 mcm, //
			final LabeledMultivariateJacobianFunction<MaterialLabel, MaterialLabel> extraElms //
	) throws ArgumentException {
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		final Material unkMat = unk.getMaterial();
		final List<MassFraction> outputs = new ArrayList<>();
		for (final KRatioLabel krl : krs)
			outputs.add(MaterialLabel.buildMassFractionTag(unkMat, krl.getElement()));
		final List<LabeledMultivariateJacobianFunction<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
		// Perform the matrix correction model to calculate Ci from ki
		res.add(mcm);
		// Perform the implicit model to propagate uncertainty in the Ci
		final List<KRatioLabel> krl = new ArrayList<>(krs);
		final KR2HModel hModel = new KR2HModel(krl);
		res.add(new ImplicitMeasurementModel<EPMALabel>(hModel, outputs, new AltComputeCUnknown(krl)));
		res.add(new Composition.MassFractionToComposition(unkMat, extraElms));
		return res;
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 model, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		super("K-Ratio Model[" + model.toString() + "]", buildSteps(krs, model, null), outputLabels);
		mModel = model;
		mKRatioSet = new ArrayList<>(krs);
		mExtraElementRules = extraElms;
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		this(krs, new XPPMatrixCorrection2(krs, outputLabels), extraElms, outputLabels);
	}

	public List<KRatioLabel> getKRatios() {
		return Collections.unmodifiableList(mKRatioSet);
	}

	/**
	 * Returns an instance of the {@link MatrixCorrectionModel2} used to build this
	 * {@link KRatioCorrectionModel2} instance.
	 *
	 * @return {@link MatrixCorrectionModel2}
	 */
	public MatrixCorrectionModel2 getModel() {
		return mModel;
	}

	/**
	 * Builds the inputs to the {@link MatrixCorrectionModel2} instance.
	 *
	 * @param estUnknown  An estimate of the composition of the unknown
	 * @param measKratios The measured k-ratios
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public UncertainValuesBase<EPMALabel> buildInput(//
			final UncertainValues<MassFraction> estUnknown, //
			final UncertainValuesBase<KRatioLabel> measKratios //
	) throws ArgumentException {
		List<UncertainValuesBase<? extends EPMALabel>> list = //
				Arrays.asList(mModel.buildInput(estUnknown), measKratios);
		return UncertainValuesBase.<EPMALabel>combine(list, false);
	}

	/**
	 * Builds a {@link LabeledMultivariateJacobianFunction} which sequentially
	 * combines the explicit {@link XPPMatrixCorrection2} model with the implicit
	 * model {@link KRatioCorrectionModel2}. It returns a KRatioCorrectionModel2.
	 *
	 * @param keySet
	 * @param extraElms
	 * @param outputLabels
	 * @return {@link KRatioCorrectionModel2}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel( //
			final Set<KRatioLabel> keySet, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		List<EPMALabel> outputs = new ArrayList<>(outputLabels);
		if (outputs.size() > 0)
			for (EPMALabel lbl : buildDefaultOutputs(keySet))
				if (!outputs.contains(lbl))
					outputs.add(lbl);
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(keySet, outputs);
		return new KRatioCorrectionModel2(keySet, xpp, extraElms, outputs);
	}

	/**
	 * Builds a {@link LabeledMultivariateJacobianFunction} which sequentially
	 * combines the explicit {@link XPPMatrixCorrection2} model with the implicit
	 * model {@link KRatioCorrectionModel2}. It returns a KRatioCorrectionModel2.
	 *
	 * The model returns the values specified by buildDefaultOutputs(...).
	 *
	 * @param keySet
	 * @param extraElms
	 * @return {@link KRatioCorrectionModel2}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel( //
			final Set<KRatioLabel> keySet, //
			final LabeledMultivariateJacobianFunction<? extends MaterialLabel, ? extends MaterialLabel> extraElms //
	) throws ArgumentException {
		return buildXPPModel(keySet, extraElms, buildDefaultOutputs(keySet));
	}

	public static List<EPMALabel> buildDefaultOutputs(final Set<KRatioLabel> keySet) {
		final FastIndex<EPMALabel> outputs = new FastIndex<>(XPPMatrixCorrection2.buildDefaultOutputs(keySet));
		Material unk = null;
		for (final KRatioLabel krl : keySet) {
			if(unk==null)
				unk=krl.getUnknown().getMaterial();
			assert unk.equals(krl.getUnknown().getMaterial());
			outputs.addIfMissing(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(), krl.getElement()));
			outputs.addIfMissing(krl.asCalculated());
		}
		outputs.addMissing(MaterialLabel.buildAtomicWeightTags(unk));
		return outputs;
	}

	public Material getUnknownMaterial() {
		return mModel.getUnknownMaterial();
	}

	/**
	 * Determine the composition that minimizes the difference between the computed
	 * and measured k-ratios.
	 *
	 *
	 * @param kratios
	 * @return Composition
	 * @throws ArgumentException
	 */
	public Composition levenbergMarquardtOptimize(final UncertainValues<KRatioLabel> kratios) //
			throws ArgumentException {
		final Material unkMat = getUnknownMaterial();
		final Composition estComp = computeFirstEstimate(unkMat, kratios);
		// Figure out which elements are in the unknown
		final List<MassFraction> compInp = MaterialLabel.buildMassFractionTags(unkMat);
		final RealVector goal = new ArrayRealVector(kratios.getDimension());
		final List<KRatioLabel> kCalc = new ArrayList<KRatioLabel>();
		int j = 0;
		for (final KRatioLabel krl : kratios.getLabels()) {
			kCalc.add(krl.asCalculated());
			goal.setEntry(j, kratios.getEntry(krl));
			++j;
		}
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction<EPMALabel, EPMALabel> trimmed = //
				new TrimmedNamedMultivariateJacobianFunction<EPMALabel, EPMALabel>(this, compInp, kCalc);
		trimmed.setBaseInputs(buildInput(estComp.toMassFraction(), kratios));
		final RealVector start = new ArrayRealVector(trimmed.getInputDimension());
		final List<Element> elms = new ArrayList<>();
		for (int i = 0; i < trimmed.getInputDimension(); ++i) {
			final MassFraction mf = (MassFraction) trimmed.getInputLabel(i);
			elms.add(mf.getElement());
			start.setEntry(i, estComp.getEntry(mf));
		}
		final ConvergenceChecker<Evaluation> checker = new EvaluationRmsChecker(1.0e-5);
		final LeastSquaresProblem lsm = LeastSquaresFactory.create( //
				trimmed, // The trimmed model
				goal, // Goal of zero
				start, // Starting estimate
				checker, 100, 100);
		final LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		final Optimum res = lmo.optimize(lsm);
		return Composition.massFraction(unkMat, new ArrayList<>(elms), res.getPoint(), res.getCovariances(1.0e-10));
	}

	/**
	 * Calculate the first estimate of the composition and build a
	 * {@link MatrixCorrectionDatum} from the data.
	 *
	 * @param measKratios
	 * @return {@link Composition}
	 * @throws ArgumentException
	 */
	private Composition computeFirstEstimate( //
			final Material unkMat, //
			final UncertainValuesBase<KRatioLabel> measKratios //
	) throws ArgumentException {
		final Map<KRatioLabel, Double> krm = measKratios.getValueMap();
		final Map<Element, Number> est = new HashMap<>();
		for (final Map.Entry<KRatioLabel, Double> me : krm.entrySet()) {
			final KRatioLabel tag = me.getKey();
			assert tag.isMeasured();
			final ElementXRaySet xrs = tag.getXRaySet();
			final Composition std = tag.getStandard().getComposition();
			final MaterialLabel.MassFraction smfTag = MaterialLabel.buildMassFractionTag(std.getMaterial(),
					xrs.getElement());
			final double cStd = std.getEntry(smfTag);
			final double kR = me.getValue().doubleValue();
			est.put(xrs.getElement(), kR * cStd);
		}
		return Composition.massFraction(unkMat, est, mExtraElementRules);
	}

	public boolean meetsThreshold(//
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final UncertainValuesBase<KRatioLabel> calcKratios, //
			final double thresh //
	) {
		double total = 0.0;
		for (final Object lbl : measKratios.getLabels()) {
			assert lbl instanceof KRatioLabel;
			final KRatioLabel measKrl = (KRatioLabel) lbl;
			final double delta = measKratios.getEntry(measKrl) - calcKratios.getEntry(measKrl.asCalculated());
			total += delta * delta;
		}
		System.out.println("Err = " + Math.sqrt(total));
		return total < thresh * thresh;
	}

	public Composition updateEstimate(//
			final Composition est, //
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final UncertainValuesBase<KRatioLabel> calcKratios) throws ArgumentException {
		final Material unkMat = est.getMaterial();
		final Map<Element, Number> res = new HashMap<>();
		for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
			assert measKrl.getUnknown().getMaterial() == unkMat;
			final MassFraction mf = MaterialLabel.buildMassFractionTag(unkMat, measKrl.getElement());
			final double cEst = est.getEntry(mf);
			final double cNext = cEst * measKratios.getEntry(measKrl) //
					/ calcKratios.getEntry(measKrl.asCalculated());
			res.put(measKrl.getElement(), cNext);
		}
		final Composition comp = Composition.massFraction(unkMat, res, mExtraElementRules);
		System.out.println(est.toMassFractionString() + " to " + comp.toMassFractionString());
		return comp;

	}

	final public Composition performIteration(final UncertainValuesBase<KRatioLabel> measKratios) //
			throws ArgumentException {
		final Material unkMat = getModel().getUnknownMaterial();
		Composition est = computeFirstEstimate(unkMat, measKratios);
		for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
			final UncertainValuesBase<EPMALabel> inps = buildInput(est.toMassFraction(), measKratios);
			// Computes the h-model k_meas - k_calc
			final RealVector res = optimized(inps.getValues(getInputLabels()));
			final UncertainValuesBase<KRatioLabel> calcKratios = KRatioLabel.extractKRatios( //
					res, getOutputLabels(), Method.Calculated);
			// Determines how large this difference is...
			if (meetsThreshold(measKratios, calcKratios, THRESH))
				break;
			est = updateEstimate(est, measKratios, calcKratios);
		}
		return est;
	}

	/**
	 * This method performs the iteration and then returns an
	 * UncertainValuesCalculator object with the results of the measurement.
	 *
	 *
	 * @param measKratios
	 * @return {@link UncertainValuesBase}
	 * @throws ArgumentException
	 */
	public final UncertainValuesCalculator<EPMALabel> iterate(
			final UncertainValuesBase<KRatioLabel> measKratios) //
			throws ArgumentException {
		final Composition unkComp = performIteration(measKratios);
		final UncertainValuesBase<EPMALabel> inps = buildInput(unkComp.toMassFraction(), measKratios);
		return new UncertainValuesCalculator<EPMALabel>(this, inps);
	}

}
