package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
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
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.XPPMatrixCorrection2;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

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
		extends SerialLabeledMultivariateJacobianFunction {

	private final MatrixCorrectionModel2 mModel;

	private static class KR2HModel //
			extends ImplicitMeasurementModel.HModel //
			implements ILabeledMultivariateFunction {

		private final Set<KRatioLabel> mKRatios;

		static List<MassFraction> buildOutputs(//
				final Set<KRatioLabel> kratios //
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
			final List<MassFraction> res = new ArrayList<>();
			for (final Element elm : elms)
				res.add(MaterialLabel.buildMassFractionTag(unk, elm));
			return res;
		}

		static List<? extends Object> buildInputs(//
				final Set<KRatioLabel> kratios //
		) {
			final Set<Object> res = new HashSet<>();
			for (final KRatioLabel krl : kratios) {
				res.add(krl);
				res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(), krl.getElement()));
				res.add(new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet()));
			}
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
		private KR2HModel(final Set<KRatioLabel> kratios) //
				throws ArgumentException {
			super(buildInputs(kratios), buildOutputs(kratios));
			mKRatios = new HashSet<>(kratios);
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
			for (final KRatioLabel kMeasTag : mKRatios) {
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final Material unkMat = unknown.getMaterial();
				final Material stdMat = standard.getMaterial();
				final Object mfUnkTag = MaterialLabel.buildMassFractionTag(unkMat, elm);
				final Object mfStdTag = MaterialLabel.buildMassFractionTag(stdMat, elm);
				final Object zafTag = new MatrixCorrectionLabel(unknown, standard, kMeasTag.getXRaySet());
				final Object hTag = new ImplicitMeasurementModel.HLabel(mfUnkTag);
				final int iKMeas = inputIndex(kMeasTag);
				final int iMFUnk = inputIndex(mfUnkTag);
				final int iMFStd = inputIndex(mfStdTag);
				final int iZAF = inputIndex(zafTag);
				final int oHTag = outputIndex(hTag);
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
			for (final KRatioLabel kMeasTag : mKRatios) {
				final Element elm = kMeasTag.getElement();
				final Object mfUnkTag = MaterialLabel.buildMassFractionTag(kMeasTag.getUnknown().getMaterial(), elm);
				final Object mfStdTag = MaterialLabel.buildMassFractionTag(kMeasTag.getStandard().getMaterial(), elm);
				final Object zafTag = new MatrixCorrectionLabel(kMeasTag.getUnknown(), kMeasTag.getStandard(),
						kMeasTag.getXRaySet());
				final Object hTag = new ImplicitMeasurementModel.HLabel(mfUnkTag);
				final int iKMeas = inputIndex(kMeasTag);
				final int iMFUnk = inputIndex(mfUnkTag);
				final int iMFStd = inputIndex(mfStdTag);
				final int iZAF = inputIndex(zafTag);
				final int oHTag = outputIndex(hTag);
				final double kMeas = point.getEntry(iKMeas);
				final double cUnk = point.getEntry(iMFUnk);
				final double cStd = point.getEntry(iMFStd);
				final double zaf = point.getEntry(iZAF);
				final double hi = kMeas - (cUnk / cStd) * zaf;
				rv.setEntry(oHTag, hi);
			}
			return rv;
		}

		@Override
		public String toString() {
			return "K-to-C";
		}
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		this(krs, new XPPMatrixCorrection2(krs, variates));
	}


	private static List<LabeledMultivariateJacobianFunction> buildSteps(//
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 mcm //
	) throws ArgumentException {
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		final List<? extends Object> outputs = MaterialLabel.massFractionTags(unk.getMaterial());
		final List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
		// Perform the matrix correction model to calculate Ci from ki
		res.add(mcm);
		// Perform the implicit model to propagate uncertainty in the Ci
		res.add(new ImplicitMeasurementModel(new KR2HModel(krs), outputs));
		return res;
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			final MatrixCorrectionModel2 model //
	) throws ArgumentException {
		super("K-Ratio Model[" + model.toString() + "]", buildSteps(krs, model));
		mModel = model;
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
	public UncertainValues buildInput(final Composition estUnknown, final UncertainValues measKratios) //
			throws ArgumentException {
		assert measKratios.getLabels(KRatioLabel.class).size() == measKratios
				.getDimension() : "Not all the measured k-ratios are k-ratios";
		return UncertainValues.combine(mModel.buildInput(estUnknown), measKratios);
	}

	/**
	 * Builds a {@link LabeledMultivariateJacobianFunction} which sequentially
	 * combines the explicit {@link XPPMatrixCorrection2} model with the implicit
	 * model {@link KRatioCorrectionModel2}. It returns a
	 * Pair&lt;SerialLabeledMultivariateJacobianFunction, UncertainValues&gt;. The
	 * first return value is the model and the second is the input minus the
	 * necessary measured k-ratio values.
	 *
	 *
	 *
	 * @param keySet
	 * @param variates
	 * @param unmeasured A list of UnmeasuredElementRule functions
	 * @return Pair&lt;SerialLabeledMultivariateJacobianFunction,
	 *         UncertainValues&gt;
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel( //
			final Set<KRatioLabel> keySet, //
			final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		return new KRatioCorrectionModel2(keySet, variates);
	}

	public UncertainValues getInputs(//
			final MatrixCorrectionModel2 mcm, //
			final UncertainValues inputs, final UncertainValues krs //
	) throws ArgumentException {
		return UncertainValues.extract(getInputLabels(), UncertainValues.propagate(mcm, inputs), krs, inputs);
	}

	/**
	 * Calculate the first estimate of the composition and build a
	 * {@link MatrixCorrectionDatum} from the data.
	 *
	 * @param kratios
	 * @return {@link Composition}
	 * @throws ArgumentException
	 */
	public static Map<Element, Number> computeFirstEstimate(final UncertainValues kratios) //
			throws ArgumentException {
		final Map<Object, Double> krm = kratios.getValueMap();
		final Map<Element, Number> est = new HashMap<>();
		for (final Map.Entry<Object, Double> me : krm.entrySet()) {
			assert me.getKey() instanceof KRatioLabel;
			if (me.getKey() instanceof KRatioLabel) {
				final KRatioLabel tag = (KRatioLabel) me.getKey();
				assert tag.isMeasured();
				final ElementXRaySet xrs = tag.getXRaySet();
				final Composition std = tag.getStandard().getComposition();
				final MaterialLabel.MassFraction smfTag = MaterialLabel.buildMassFractionTag(std.getMaterial(),
						xrs.getElement());
				final double cStd = std.getEntry(smfTag);
				final double kR = me.getValue().doubleValue();
				est.put(xrs.getElement(), kR * cStd);
			}
		}
		return est;
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
	public Composition optimize(final UncertainValues kratios) //
			throws ArgumentException {
		final Map<Element, Number> unk = computeFirstEstimate(kratios);
		final Material unkMat = mModel.getUnknownMaterial();
		assert unk.keySet().equals(unkMat.getElementSet());
		// Figure out which elements are in the unknown
		final List<MassFraction> compInp = MaterialLabel.buildMassFractionTags(unkMat);
		final RealVector goal = new ArrayRealVector(kratios.getDimension());
		final List<KRatioLabel> kCalc = new ArrayList<KRatioLabel>();
		int j = 0;
		for (final Object krl : kratios.getLabels()) {
			assert krl instanceof KRatioLabel;
			kCalc.add(((KRatioLabel) krl).asCalculated());
			goal.setEntry(j, kratios.getEntry(krl));
			++j;
		}
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction trimmed = //
				new TrimmedNamedMultivariateJacobianFunction(mModel, compInp, kCalc);
		final UncertainValues inputs = mModel.buildInput(Composition.massFraction(unkMat, unk));
		trimmed.initializeConstants(inputs.getValueMap());
		final RealVector start = new ArrayRealVector(trimmed.getInputDimension());
		int i = 0;
		final List<Element> elms = new ArrayList<>();
		for (final Object inLab : trimmed.getInputLabels()) {
			assert inLab instanceof MassFraction;
			final MassFraction mf = (MassFraction) inLab;
			elms.add(mf.getElement());
			start.setEntry(i, unk.get(mf.getElement()).doubleValue());
			++i;
		}
		final ConvergenceChecker<Evaluation> checker = new EvaluationRmsChecker(1.0e-3);
		final LeastSquaresProblem lsm = LeastSquaresFactory.create( //
				trimmed, // The trimmed model
				goal, // Goal of zero
				start, // Starting estimate
				checker, 100, 100);
		final LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		final Optimum res = lmo.optimize(lsm);
		return Composition.massFraction("Quant", new ArrayList<>(elms), res.getPoint(), res.getCovariances(1.0e-10));
	}
}
