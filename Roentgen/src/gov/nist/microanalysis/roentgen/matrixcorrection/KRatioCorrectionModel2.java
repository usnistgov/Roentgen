package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
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
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.Material;

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

		static List<? extends Object> buildOutputs(//
				final Set<KRatioLabel> kratios //
		) throws ArgumentException {
			Material unk = null;
			for (final KRatioLabel krl : kratios) {
				if (unk == null)
					unk = krl.getUnknown().getMaterial();
				if (!unk.equals(krl.getUnknown().getMaterial()))
					throw new ArgumentException("More than one unknown in KRatioCorrectionModel2.HModel");
			}
			return MaterialLabel.massFractionTags(unk);
		}

		static List<? extends Object> buildInputs(//
				final Set<KRatioLabel> kratios //
		) {
			final List<Object> res = new ArrayList<>();
			for (final KRatioLabel krl : kratios) {
				res.add(krl);
				res.add(MaterialLabel.buildMassFractionTag(krl.getStandard().getMaterial(),
						krl.getXRaySet().getElement()));
				res.add(new MatrixCorrectionLabel(krl.getUnknown(), krl.getStandard(), krl.getXRaySet()));
			}
			return res;
		}

		private KR2HModel(//
				final Set<KRatioLabel> kratios //
		) throws ArgumentException {
			super(buildInputs(kratios), buildOutputs(kratios));
			mKRatios = new HashSet<>(kratios);
		}

		@Override
		public Pair<RealVector, RealMatrix> value(final RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
			for (final KRatioLabel kMeasTag : mKRatios) {
				final Element elm = kMeasTag.getXRaySet().getElement();
				final Object mfUnkTag = MaterialLabel.buildMassFractionTag(kMeasTag.getUnknown().getMaterial(), elm);
				final Object mfStdTag = MaterialLabel.buildMassFractionTag(kMeasTag.getStandard().getMaterial(), elm);
				final Object zafTag = new MatrixCorrectionLabel(kMeasTag.getUnknown(), kMeasTag.getStandard(),
						kMeasTag.getXRaySet());
				final Object hTag = new ImplicitMeasurementModel.HLabel(MaterialLabel.buildMassFractionTag(kMeasTag.getUnknown().getMaterial(), elm));
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
				// assert Math.abs(hi) < 1.0e-6 : "h["+kMeasTag + "] = " + hi;
				rv.setEntry(oHTag, hi);
				rm.setEntry(oHTag, iKMeas, 1.0);
				rm.setEntry(oHTag, iMFUnk, (-1.0 / cStd) * zaf);
				rm.setEntry(oHTag, iMFStd, (cUnk / (cStd * cStd)) * zaf);
				rm.setEntry(oHTag, iZAF, -(cUnk / cStd));
			}
			return Pair.create(rv, rm);
		}

		@Override
		public RealVector optimized(RealVector point) {
			final RealVector rv = new ArrayRealVector(getOutputDimension());
			for (final KRatioLabel kMeasTag : mKRatios) {
				final Element elm = kMeasTag.getXRaySet().getElement();
				final Object mfUnkTag = MaterialLabel.buildMassFractionTag(kMeasTag.getUnknown().getMaterial(), elm);
				final Object mfStdTag = MaterialLabel.buildMassFractionTag(kMeasTag.getStandard().getMaterial(), elm);
				final Object zafTag = new MatrixCorrectionLabel(kMeasTag.getUnknown(), kMeasTag.getStandard(),
						kMeasTag.getXRaySet());
				final Object hTag = new ImplicitMeasurementModel.HLabel(
						MaterialLabel.buildMassFractionTag(kMeasTag.getUnknown().getMaterial(), elm));
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
				// assert Math.abs(hi) < 1.0e-6 : "h["+kMeasTag + "] = " + hi;
				rv.setEntry(oHTag, hi);
			}
			return rv;
		}

		@Override
		public String toString() {
			return "H-model";
		}
	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			final Set<MatrixCorrectionModel2.Variate> variates) throws ArgumentException {
		this(krs, new XPPMatrixCorrection2(krs, variates));
	}

	private static List<LabeledMultivariateJacobianFunction> buildSteps(Set<KRatioLabel> krs,
			MatrixCorrectionModel2 mcm) //
			throws ArgumentException {
		final UnknownMatrixCorrectionDatum unk = krs.iterator().next().getUnknown();
		List<? extends Object> outputs = MaterialLabel.massFractionTags(unk.getMaterial());
		List<LabeledMultivariateJacobianFunction> res = new ArrayList<>();
		res.add(mcm);
		res.add(new ImplicitMeasurementModel(new KR2HModel(krs), outputs));
		return res;

	}

	public KRatioCorrectionModel2(//
			final Set<KRatioLabel> krs, //
			MatrixCorrectionModel2 model //
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
	 * Builds the inputs to the {@link MatrixCorrectionModel2} instance. Requires
	 * additional KRatioLabel related to measured items.
	 * 
	 * @return {@link UncertainValues}
	 * @throws ArgumentException
	 */
	public UncertainValues buildInput(Composition estUnknown) //
			throws ArgumentException {
		return mModel.buildInput(estUnknown);
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
	 * @return Pair&lt;SerialLabeledMultivariateJacobianFunction,
	 *         UncertainValues&gt;
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel2 buildXPPModel( //
			final Set<KRatioLabel> keySet, final Set<MatrixCorrectionModel2.Variate> variates //
	) throws ArgumentException {
		return new KRatioCorrectionModel2(keySet, variates);
	}

	public UncertainValues getInputs(final MatrixCorrectionModel2 mcm, final UncertainValues inputs,
			final UncertainValues krs) //
			throws ArgumentException {
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
	public Composition computeEstimate(final UncertainValues kratios) throws ArgumentException {
		final Map<Object, Double> krm = kratios.getValueMap();
		final Map<Element, Number> est = new HashMap<>();
		for (final Map.Entry<Object, Double> me : krm.entrySet()) {
			assert me.getKey() instanceof KRatioLabel;
			if (me.getKey() instanceof KRatioLabel) {
				final KRatioLabel tag = (KRatioLabel) me.getKey();
				assert tag.isMeasured();
				final ElementXRaySet xrs = tag.getXRaySet();
				final Composition std = tag.getStandard().getComposition();
				final MaterialLabel.MassFraction smfTag = MaterialLabel.buildMassFractionTag(std.getMaterial(), xrs.getElement());
				final double cStd = std.getEntry(smfTag);
				final double kR = me.getValue().doubleValue();
				est.put(xrs.getElement(), kR * cStd);
			}
		}
		return Composition.massFraction("Estimate[0]", est);
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
		final Composition unk = computeEstimate(kratios);
		final Set<Element> elms = unk.getElementSet();
		// Figure out which elements are in the unknown
		final List<MaterialLabel.MassFraction> compInp = MaterialLabel.massFractionTags(unk.getMaterial());
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction trimmed = //
				new TrimmedNamedMultivariateJacobianFunction(mModel, compInp, mModel.getOutputLabels());
		final ConvergenceChecker<Evaluation> checker = new EvaluationRmsChecker(1.0e-3);
		final List<? extends Object> trOutTags = trimmed.getOutputLabels();
		final LeastSquaresProblem lsm = LeastSquaresFactory.create( //
				trimmed, // The trimmed model
				new ArrayRealVector(trOutTags.size()), // Goal of zero
				unk.extractValues(trimmed.getInputLabels()), // Starting estimate
				checker, 100, 100);
		final LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		final Optimum res = lmo.optimize(lsm);
		return Composition.massFraction("Quant", new ArrayList<>(elms), res.getPoint(), res.getCovariances(1.0e-10));
	}
}
