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
import org.apache.commons.math3.optim.ConvergenceChecker;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.model.MatrixCorrectionModel2;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;

/**
 * Implements the algorithm used to estimate the composition given a set of
 * measured k-ratio values.
 *
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class KRatioIteration {

	private final Set<KRatioLabel> mKRatios;

	public KRatioIteration(final Set<KRatioLabel> skrl) {
		assert KRatioLabel.areAllSameUnknownElements(skrl);
		mKRatios = new HashSet<>(skrl);

	}

	/**
	 * Calculate the first estimate of the composition and build a
	 * {@link MatrixCorrectionDatum} from the data.
	 *
	 * @param kratios
	 * @return {@link Composition}
	 */
	public Composition computeEstimate(final UncertainValues kratios) {
		final Map<Object, Double> krm = kratios.getValueMap();
		final Map<Element, Number> est = new HashMap<>();
		for (final Map.Entry<Object, Double> me : krm.entrySet()) {
			assert me.getKey() instanceof KRatioLabel;
			if (me.getKey() instanceof KRatioLabel) {
				final KRatioLabel tag = (KRatioLabel) me.getKey();
				assert tag.isMeasured();
				final ElementXRaySet xrs = tag.getXRaySet();
				final Composition std = tag.getStandard().getComposition();
				final MassFractionTag smfTag = Composition.buildMassFractionTag(std, xrs.getElement());
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
		final CompositionFromKRatios2 model = new CompositionFromKRatios2(mKRatios,
				Collections.singleton(MatrixCorrectionModel2.Variate.UnknownComposition));
		// Figure out which elements are in the unknown
		final List<MassFractionTag> compInp = new ArrayList<>();
		for (final Element elm : elms)
			compInp.add(Composition.buildMassFractionTag(unk, elm));
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction trimmed = //
				new TrimmedNamedMultivariateJacobianFunction(model, compInp, model.getOutputLabels());
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
