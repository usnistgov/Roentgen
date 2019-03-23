package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
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

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.CompositeMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ExplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel2;
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
		extends CompositeMeasurementModel<EPMALabel> {

	static private final int MAX_ITERATIONS = 100;
	static private final double THRESH = 1.0e-4;

	private final MatrixCorrectionModel2 mModel;

	private final ExplicitMeasurementModel<? extends MaterialLabel, ? extends MaterialLabel> mExtraElementRules;

	private final List<KRatioLabel> mKRatioSet;

	private static class KRatioModel2 extends ImplicitMeasurementModel2<EPMALabel, MassFraction> {

		static List<EPMALabel> buildInputs(
				final List<KRatioLabel> kratios //
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
				final List<KRatioLabel> kratios //
		) {
			final ArrayList<MassFraction> res = new ArrayList<>();
			for (final KRatioLabel krl : kratios)
				res.add(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(), krl.getElement()));
			return new ArrayList<>(res);
		}

		private final List<KRatioLabel> mKRatios;

		public KRatioModel2(
				final List<KRatioLabel> kratios
		) throws ArgumentException {
			super(buildInputs(kratios), buildOutputs(kratios));
			mKRatios = kratios;
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
				final double cUnk = getOutputValue(mfUnkTag);
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
			final ExplicitMeasurementModel<MaterialLabel, MaterialLabel> extraElms //
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
		final List<KRatioLabel> krl = new ArrayList<>(krs);
		res.add(new KRatioModel2(krl));
		res.add(new Composition.MassFractionToComposition(unkMat, extraElms));
		return res;
	}

	public KRatioCorrectionModel2(
			//
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 model, //
			final ExplicitMeasurementModel<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
			final List<EPMALabel> outputLabels //
	) throws ArgumentException {
		super("K-Ratio Model[" + model.toString() + "]", buildSteps(krs, model, null), outputLabels);
		mModel = model;
		mKRatioSet = new ArrayList<>(krs);
		mExtraElementRules = extraElms;
	}

	public KRatioCorrectionModel2(
			//
			final Set<KRatioLabel> krs, //
			final ExplicitMeasurementModel<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
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
			//
			final UncertainValues<MassFraction> estUnknown, //
			final UncertainValuesBase<KRatioLabel> measKratios //
	) throws ArgumentException {
		final UncertainValuesBase<EPMALabel> inputs = mModel.buildInput(estUnknown);
		final List<UncertainValuesBase<? extends EPMALabel>> list = //
				Arrays.asList(inputs, measKratios);
		return UncertainValuesBase.<EPMALabel>combine(list, false);
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
			final ExplicitMeasurementModel<? extends MaterialLabel, ? extends MaterialLabel> extraElms, //
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
			final ExplicitMeasurementModel<? extends MaterialLabel, ? extends MaterialLabel> extraElms //
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
	public Composition levenbergMarquardtOptimize(
			final UncertainValues<KRatioLabel> kratios
	) //
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
	private Composition computeFirstEstimate(
			//
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

	public Composition updateEstimate(
			//
			final Composition est, //
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final UncertainValuesBase<KRatioLabel> calcKratios
	) throws ArgumentException {
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

	final public Composition performIteration(
			final UncertainValuesBase<KRatioLabel> measKratios
	) //
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
			final UncertainValuesBase<KRatioLabel> measKratios
	) //
			throws ArgumentException {
		final Composition unkComp = performIteration(measKratios);
		final UncertainValuesBase<EPMALabel> inps = buildInput(unkComp.toMassFraction(), measKratios);
		return new UncertainValuesCalculator<EPMALabel>(this, inps);
	}

}
