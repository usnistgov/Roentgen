package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
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
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialLabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.XPPMatrixCorrection2.Variates;
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
	
	public KRatioIteration(Set<KRatioLabel> skrl) {
		assert KRatioLabel.areAllSameUnknown(skrl);
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
		final Map<Element, Number> men = new HashMap<>();
		for (final Map.Entry<Object, Double> me : krm.entrySet()) {
			assert me.getKey() instanceof KRatioLabel;
			if (me.getKey() instanceof KRatioLabel) {
				final KRatioLabel tag = (KRatioLabel) me.getKey();
				assert tag.isMeasured();
				final ElementXRaySet xrs = tag.getXRaySet();
				final Composition std = tag.getStandard().getComposition();
				final double cStd = std.getEntry(Composition.buildMassFractionTag(std, xrs.getElement()));
				final double kR = me.getValue().doubleValue();
				men.put(xrs.getElement(), kR * cStd);
			}
		}
		return Composition.massFraction("Estimate[0]", men);
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
		Set<Element> elms = new HashSet<>();		
		for(KRatioLabel krl : kratios.extractTypeOfLabel(KRatioLabel.class))
			elms.add(krl.getXRaySet().getElement());

		final Set<Variates> minVariates = XPPMatrixCorrection2.minimalVariates();
		final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(mKRatios, minVariates);
		final KRatioHModel2 hModel = new KRatioHModel2(mKRatios);

		// Build the full model with xpp and hModel
		final List<LabeledMultivariateJacobianFunction> steps = new ArrayList<>();
		steps.add(xpp);
		steps.add(hModel);
		final SerialLabeledMultivariateJacobianFunction model = //
				new SerialLabeledMultivariateJacobianFunction("K-ratio iteration", steps);
		// Figure out which elements are in the unknown
		final List<MassFractionTag> compInp = new ArrayList<>();
		for (final Element elm : unk.getElementSet()) {
			compInp.add(Composition.buildMassFractionTag(unk, elm));
		}
		// Trim the inputs down to only the elements in the unknown...
		final TrimmedNamedMultivariateJacobianFunction trimmed = //
				new TrimmedNamedMultivariateJacobianFunction(model, compInp, model.getOutputLabels());
		// Build up a full set of all required inputs and constants
		final UncertainValues msInp = UncertainValues.build(model.getInputLabels(), xpp.buildInput(), kratios);
		{ // Initialize all non-variables as constants...
			final Map<Object, Double> mod = new HashMap<>();
			for (final Object tag : msInp.getLabels())
				if (!compInp.contains(tag))
					mod.put(tag, msInp.getEntry(tag));
			trimmed.initializeConstants(mod);
		}
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
