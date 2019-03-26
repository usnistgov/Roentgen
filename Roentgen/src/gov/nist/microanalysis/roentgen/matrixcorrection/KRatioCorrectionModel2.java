package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ImplicitMeasurementModel;
import gov.nist.juncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
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
		extends CompositeMeasurementModel<EPMALabel> {

	static private final int MAX_ITERATIONS = 100;
	static private final double THRESH = 1.0e-4;

	private final MatrixCorrectionModel2 mModel;

	private final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> mExtraElementRule;

	private final List<KRatioLabel> mKRatioSet;

	private static class KRatioModel2 extends ImplicitMeasurementModel<EPMALabel, MassFraction> {

		static List<EPMALabel> buildInputs(
				final Set<KRatioLabel> kratios //
		) {
			final ArrayList<EPMALabel> res = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				res.add(krl);
				res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
				res.add(MatrixCorrectionModel2.zafLabel(krl));
			}
			return res;
		}

		static List<MassFraction> buildOutputs(
				final Set<KRatioLabel> kratios //
		) {
			final ArrayList<MassFraction> res = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				res.add(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(), krl.getElement()));
			return new ArrayList<>(res);
		}

		private final List<KRatioLabel> mKRatios;

		public KRatioModel2(
				final Set<KRatioLabel> kratios
		) throws ArgumentException {
			super(buildInputs(kratios), buildOutputs(kratios));
			mKRatios = new ArrayList<>(kratios);
		}

		@Override
		public RealMatrix computeCx(
				final RealVector point
		) {
			final RealMatrix jac = buildEmptyCx();
			for (int oHTag = 0; oHTag < mKRatios.size(); ++oHTag) {
				final KRatioLabel kMeasTag = mKRatios.get(oHTag);
				final Element elm = kMeasTag.getElement();
				final Material unkMat = kMeasTag.getUnknown().getMaterial();
				final Material stdMat = kMeasTag.getStandard().getMaterial();
				final MassFraction mfUnkTag = MaterialLabel.buildMassFractionTag(unkMat, elm);
				final MassFraction mfStdTag = MaterialLabel.buildMassFractionTag(stdMat, elm);
				final ZAFMultiLineLabel zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				// Get inputs
				final double cUnk = getArg(mfUnkTag, point);
				final double cStd = getArg(mfStdTag, point);
				final double zaf = getArg(zafTag, point);
				// assert Math.abs(hi * kMeas) < 1.0e-4 : "h[" + kMeasTag + "] = " + hi;
				setCx(oHTag, kMeasTag, jac, 1.0);
				setCx(oHTag, mfStdTag, jac, (cUnk / (cStd * cStd)) * zaf);
				setCx(oHTag, zafTag, jac, -(cUnk / cStd));
			}
			return jac;
		}

		@Override
		public RealMatrix computeCy(
				final RealVector point
		) {
			final RealMatrix jac = buildEmptyCy();
			for (int oHTag = 0; oHTag < mKRatios.size(); ++oHTag) {
				final KRatioLabel kMeasTag = mKRatios.get(oHTag);
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final MassFraction mfUnkTag = MaterialLabel.buildMassFractionTag(unknown.getMaterial(), elm);
				final MassFraction mfStdTag = MaterialLabel.buildMassFractionTag(standard.getMaterial(), elm);
				final ZAFMultiLineLabel zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				final double cStd = getArg(mfStdTag, point);
				final double zaf = getArg(zafTag, point);
				setCy(oHTag, mfUnkTag, jac, (-1.0 / cStd) * zaf);
			}
			return jac;
		}

	}

	private static List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> buildSteps(
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 mcm, //
			final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> extraElms //
	) throws ArgumentException {
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		final Material unkMat = unk.getMaterial();
		final List<MassFraction> outputs = new ArrayList<>();
		for (final KRatioLabel krl : krs)
			outputs.add(MaterialLabel.buildMassFractionTag(unkMat, krl.getElement()));
		final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends EPMALabel>> res = new ArrayList<>();
		// Perform the matrix correction model to calculate Ci from ki
		res.add(mcm);
		// Perform the implicit model to propagate uncertainty in the Ci
		res.add(new KRatioModel2(krs));
		if (extraElms != null)
			res.add(extraElms);
		// res.add(new Composition.MassFractionToComposition(unkMat, extraElms));
		return res;
	}

	public KRatioCorrectionModel2(
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 model, //
			final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		super("K-Ratio Model[" + model.toString() + "]", buildSteps(krs, model, extraElms), outputLabels);
		mModel = model;
		mKRatioSet = new ArrayList<>(krs);
		mExtraElementRule = extraElms;
	}

	public KRatioCorrectionModel2(
			//
			final Set<KRatioLabel> krs, //
			final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> extraElms, //
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
	public UncertainValuesBase<EPMALabel> buildInput(
			final UncertainValuesBase<KRatioLabel> measKratios //
	) throws ArgumentException {
		Material unkMat = null;
		Set<Composition> stds = new HashSet<>();
		for (KRatioLabel krl : measKratios.getLabels()) {
			final Material mat = krl.getUnknown().getMaterial();
			if (unkMat == null)
				unkMat = mat;
			if (!unkMat.equals(mat))
				throw new ArgumentException("Two different materials [" + mat + " and " + unkMat
						+ " are present in the estimated unknown.");
			stds.add(krl.getStandard().getComposition());
		}
		List<UncertainValuesBase<? extends EPMALabel>> list = new ArrayList<>();
		final UncertainValuesBase<EPMALabel> modInps = mModel.buildInput(unkMat);
		list.add(modInps);
		list.add(measKratios);
		list.add(unkMat.getAtomicWeights());
		for (KRatioLabel krl : measKratios.getLabels()) {
			List<MassFraction> lmf = new ArrayList<>();
			for (MassFraction mf : krl.getStandard().getComposition().massFractionTags())
				if (mModel.outputIndex(mf) == -1)
					lmf.add(mf);
			if (!lmf.isEmpty())
				list.add(krl.getStandard().getComposition().extract(lmf));
		}
		return UncertainValuesBase.<EPMALabel>combine(list, true);
	}

	/**
	 * Builds a {@link ExplicitMeasurementModel} which sequentially combines the
	 * explicit {@link XPPMatrixCorrection2} model with the implicit model
	 * {@link KRatioCorrectionModel2}. It returns a KRatioCorrectionModel2.
	 *
	 * @param keySet
	 * @param extraElms
	 * @param outputLabels
	 * @return {@link KRatioCorrectionModel2}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel(
			//
			final Set<KRatioLabel> keySet, //
			final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		final List<EPMALabel> outputs = new ArrayList<>(outputLabels);
		if (outputs.size() > 0)
			for (final EPMALabel lbl : buildDefaultOutputs(keySet))
				if (!outputs.contains(lbl))
					outputs.add(lbl);
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(keySet, outputs);
		return new KRatioCorrectionModel2(keySet, xpp, extraElms, outputs);
	}

	/**
	 * Builds a {@link ExplicitMeasurementModel} which sequentially combines the
	 * explicit {@link XPPMatrixCorrection2} model with the implicit model
	 * {@link KRatioCorrectionModel2}. It returns a KRatioCorrectionModel2.
	 *
	 * The model returns the values specified by buildDefaultOutputs(...).
	 *
	 * @param keySet
	 * @param extraElms
	 * @return {@link KRatioCorrectionModel2}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel(
			//
			final Set<KRatioLabel> keySet, //
			final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> extraElms //
	) throws ArgumentException {
		return buildXPPModel(keySet, extraElms, buildDefaultOutputs(keySet));
	}

	public static List<EPMALabel> buildDefaultOutputs(
			final Set<KRatioLabel> keySet
	) {
		final FastIndex<EPMALabel> outputs = new FastIndex<>(XPPMatrixCorrection2.buildDefaultOutputs(keySet));
		Material unk = null;
		for (final KRatioLabel krl : keySet) {
			if (unk == null)
				unk = krl.getUnknown().getMaterial();
			assert unk.equals(krl.getUnknown().getMaterial());
			outputs.addIfMissing(krl.asCalculated());
		}
		outputs.addMissing(MaterialLabel.buildMassFractionTags(unk));
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
	 * @param measKratios
	 * @return Composition
	 * @throws ArgumentException
	 */
	public Composition levenbergMarquardtOptimize(
			final UncertainValues<KRatioLabel> measKratios
	) throws ArgumentException {
		final Material unkMat = getUnknownMaterial();
		final UncertainValuesBase<EPMALabel> inps = buildInput(measKratios);
		final RealVector point = inps.getValues(getInputLabels());
		final Map<MassFraction, Double> estComp = computeFirstEstimate(unkMat, measKratios, point);
		// Figure out which elements are in the unknown
		final List<MassFraction> compInp = MaterialLabel.buildMassFractionTags(unkMat);
		final RealVector goal = new ArrayRealVector(measKratios.getDimension());
		final List<KRatioLabel> kCalc = new ArrayList<KRatioLabel>();
		int j = 0;
		for (final KRatioLabel krl : measKratios.getLabels()) {
			kCalc.add(krl.asCalculated());
			goal.setEntry(j, measKratios.getEntry(krl));
			++j;
		}
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction<EPMALabel, EPMALabel> trimmed = //
				new TrimmedNamedMultivariateJacobianFunction<EPMALabel, EPMALabel>(this, compInp, kCalc);
		trimmed.setBaseInputs(buildInput(measKratios));
		final RealVector start = new ArrayRealVector(trimmed.getInputDimension());
		final List<Element> elms = new ArrayList<>();
		for (int i = 0; i < trimmed.getInputDimension(); ++i) {
			final MassFraction mf = (MassFraction) trimmed.getInputLabel(i);
			elms.add(mf.getElement());
			start.setEntry(i, estComp.get(mf));
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
	private Map<MassFraction, Double> computeFirstEstimate(
			final Material unkMat, //
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final RealVector point
	) throws ArgumentException {
		final Map<KRatioLabel, Double> krm = measKratios.getValueMap();
		final Map<MassFraction, Double> est = new HashMap<>();
		for (final Map.Entry<KRatioLabel, Double> me : krm.entrySet()) {
			final KRatioLabel tag = me.getKey();
			assert tag.isMeasured();
			final ElementXRaySet xrs = tag.getXRaySet();
			final Composition std = tag.getStandard().getComposition();
			final MaterialLabel.MassFraction smfTag = MaterialLabel.buildMassFractionTag(//
					std.getMaterial(), xrs.getElement());
			final double cStd = std.getEntry(smfTag);
			final double kR = me.getValue().doubleValue();
			est.put(MaterialLabel.buildMassFractionTag(unkMat, xrs.getElement()), kR * cStd);
		}
		return applyExtraElementRule(unkMat, point, est);
	}

	private Map<MassFraction, Double> applyExtraElementRule(
			Material unkMat, //
			RealVector point, //
			Map<MassFraction, Double> est
	) {
		if (mExtraElementRule != null) {
			RealVector pt = new ArrayRealVector(mExtraElementRule.getInputDimension());
			for (int i = 0; i < pt.getDimension(); ++i) {
				final MaterialLabel ml = mExtraElementRule.getInputLabel(i);
				if (ml instanceof MassFraction)
					pt.setEntry(i, est.get(ml).doubleValue());
				else
					pt.setEntry(i, getArg(ml, point));
			}
			final RealVector res = mExtraElementRule.compute(pt);
			for (int i = 0; i < mExtraElementRule.getOutputDimension(); ++i) {
				final MassFraction mf = mExtraElementRule.getOutputLabel(i);
				est.put(mf, res.getEntry(i));
			}
		}
		return est;
	}

	public boolean meetsThreshold(
			//
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

	public Map<MassFraction, Double> updateEstimate(
			final Material unkMat, //
			final Map<MassFraction, Double> est, //
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final UncertainValuesBase<KRatioLabel> calcKratios, //
			final RealVector point
	) throws ArgumentException {
		final Map<MassFraction, Double> res = new HashMap<>();
		for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
			assert measKrl.getUnknown().getMaterial() == unkMat;
			final MassFraction mf = MaterialLabel.buildMassFractionTag(unkMat, measKrl.getElement());
			final double cEst = est.get(mf);
			final double cNext = cEst * measKratios.getEntry(measKrl) //
					/ calcKratios.getEntry(measKrl.asCalculated());
			res.put(mf, cNext);
		}
		System.out.println(est + " to " + res);
		return applyExtraElementRule(unkMat, point, res);

	}

	final public Map<MassFraction, Double> performIteration(
			final UncertainValuesBase<KRatioLabel> measKratios
	) throws ArgumentException {
		final Material unkMat = getModel().getUnknownMaterial();
		final UncertainValuesBase<EPMALabel> inp = buildInput(measKratios);
		final RealVector point = new ArrayRealVector(getInputDimension());
		for (int i = 0; i < getInputDimension(); ++i) {
			final EPMALabel inpLbl = getInputLabel(i);
			if (inp.hasLabel(inpLbl))
				point.setEntry(i, inp.getEntry(inpLbl));
		}
		Map<MassFraction, Double> est = computeFirstEstimate(unkMat, measKratios, point);
		for (int i = 0; i < getInputDimension(); ++i) {
			final EPMALabel inpLbl = getInputLabel(i);
			if (est.containsKey(inpLbl))
				point.setEntry(i, est.get(inpLbl));
		}
		for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
			// Computes the h-model k_meas - k_calc
			addAdditionalInputs(est);
			final RealVector res = optimized(point);
			final UncertainValuesBase<KRatioLabel> calcKratios = KRatioLabel.extractKRatios( //
					res, getOutputLabels(), Method.Calculated);
			// Determines how large this difference is...
			if (meetsThreshold(measKratios, calcKratios, THRESH))
				break;
			est = updateEstimate(unkMat, est, measKratios, calcKratios, point);
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
			final UncertainValuesBase<KRatioLabel> measKratios
	) throws ArgumentException {
		final Map<MassFraction, Double> bestEst = performIteration(measKratios);
		final UncertainValuesBase<EPMALabel> inps = buildInput(measKratios);
		addAdditionalInputs(bestEst);
		return new UncertainValuesCalculator<EPMALabel>(this, inps);
	}

}
