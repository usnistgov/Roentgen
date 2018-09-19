package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
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
import org.apache.commons.math3.optim.ConvergenceChecker;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection.Variates;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;

/**
 * Implements the algorithm used to estimate the composition given a set of
 * measured k-ratio values.
 *
 *
 * @author Nicholas
 *
 */
public class KRatioIteration {

	private final Map<ElementXRaySet, MatrixCorrectionDatum> mStandards;
	private final double mBeamEnergy;
	private final double mTakeOffAngle;

	public KRatioIteration(final UncertainValues kratios, final double e0, final double toa,
			final Map<ElementXRaySet, MatrixCorrectionDatum> stds) {
		mStandards = stds;
		mBeamEnergy = e0;
		mTakeOffAngle = toa;

	}

	/**
	 * Calculate the first estimate of the composition and build a
	 * {@link MatrixCorrectionDatum} from the data.
	 *
	 * @param kratios
	 * @return {@link MatrixCorrectionDatum}
	 */
	public MatrixCorrectionDatum computeEstimate(final UncertainValues kratios) {
		final Map<Object, Double> krm = kratios.getValueMap();
		final Map<Element, Number> men = new HashMap<>();
		for (final Map.Entry<Object, Double> me : krm.entrySet()) {
			assert me.getKey() instanceof KRatioTag;
			if (me.getKey() instanceof KRatioTag) {
				final KRatioTag tag = (KRatioTag) me.getKey();
				assert tag.isMeasured();
				final ElementXRaySet xrs = tag.getXRaySet();
				final Composition std = tag.getStandard().getComposition();
				final double cStd = std.getEntry(Composition.buildMassFractionTag(std, xrs.getElement()));
				final double kR = me.getValue().doubleValue();
				men.put(xrs.getElement(), kR * cStd);
			}
		}
		final Composition comp = Composition.massFraction("Estimate", men);
		return new MatrixCorrectionDatum(comp, false, new UncertainValue(mBeamEnergy),
				new UncertainValue(mTakeOffAngle));
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

		final MatrixCorrectionDatum unkMcd = computeEstimate(kratios);

		final Set<Variates> minVariates = XPPMatrixCorrection.minimalVariates();
		final XPPMatrixCorrection xpp = new XPPMatrixCorrection(unkMcd, mStandards, minVariates);
		final KRatioHModel hModel = new KRatioHModel(unkMcd, mStandards);

		// Build the full model with xpp and hModel
		final List<NamedMultivariateJacobianFunction> steps = new ArrayList<>();
		steps.add(xpp);
		steps.add(hModel);
		final SerialNamedMultivariateJacobianFunction model = //
				new SerialNamedMultivariateJacobianFunction("K-ratio iteration", steps);
		// Figure out which elements are in the unknown
		final List<MassFractionTag> compInp = new ArrayList<>();
		final List<Element> elms = new ArrayList<>();
		for (final Element elm : unkMcd.getComposition().getElementSet()) {
			compInp.add(Composition.buildMassFractionTag(unkMcd.getComposition(), elm));
			elms.add(elm);
		}
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction trimmed = //
				new TrimmedNamedMultivariateJacobianFunction(model, compInp, model.getOutputTags());
		// Build up a full set of all required inputs and constants
		final UncertainValues msInp = UncertainValues.build(model.getInputTags(), xpp.buildInput(), kratios);
		{ // Initialize all non-variables as constants...
			final Map<Object, Double> mod = new HashMap<>();
			for (final Object tag : msInp.getTags())
				if (!compInp.contains(tag))
					mod.put(tag, msInp.getEntry(tag));
			trimmed.initializeConstants(mod);
		}
		final ConvergenceChecker<Evaluation> checker = new EvaluationRmsChecker(1.0e-3);
		final List<? extends Object> trOutTags = trimmed.getOutputTags();
		final LeastSquaresProblem lsm = LeastSquaresFactory.create( //
				trimmed, // The trimmed model
				new ArrayRealVector(trOutTags.size()), // Goal of zero
				unkMcd.getComposition().extractValues(trimmed.getInputTags()), // Starting estimate
				checker, 100, 100);
		final LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
		final Optimum res = lmo.optimize(lsm);
		return Composition.massFraction("Quant", elms, res.getPoint(), res.getCovariances(1.0e-10));
	}

}
