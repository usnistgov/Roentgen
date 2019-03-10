package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.RealVector;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class Iteration {

	static private final int MAX_ITERATIONS = 100;
	static private final double THRESH = 1.0e-4;

	private final KRatioCorrectionModel2 mModel;

	/**
	 * Calculate the first estimate of the composition and build a
	 * {@link MatrixCorrectionDatum} from the data.
	 *
	 * @param kratios
	 * @return {@link Composition}
	 * @throws ArgumentException
	 */
	public Map<Element, Number> computeFirstEstimate(final UncertainValues kratios) {
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

	public boolean meetsThreshold(final UncertainValues meas, final UncertainValues calc, final double thresh) {
		double total = 0.0;
		for (final Object lbl : meas.getLabels()) {
			assert lbl instanceof KRatioLabel;
			final KRatioLabel measKrl = (KRatioLabel) lbl;
			final double delta = meas.getEntry(measKrl) - calc.getEntry(measKrl.asCalculated());
			total += delta * delta;
		}
		System.out.println("Err = " + Math.sqrt(total));
		return total < thresh * thresh;
	}

	public Map<Element, Number> updateEstimate(//
			final Map<Element, Number> est, //
			final UncertainValues measKratios, //
			final UncertainValues calcKRatios //
	) {
		final Map<Element, Number> res = new HashMap<>();
		for (final Object lbl : measKratios.getLabels()) {
			assert lbl instanceof KRatioLabel;
			final KRatioLabel measKrl = (KRatioLabel) lbl;
			final Element elm = measKrl.getElement();
			final double cEst = est.get(elm).doubleValue();
			final double cNext = cEst * measKratios.getEntry(measKrl) //
					/ calcKRatios.getEntry(measKrl.asCalculated());
			res.put(elm, cNext);
		}
		System.out.println(est + " -> " + res);
		return res;
	}

	/**
	 * @throws ArgumentException
	 *
	 */
	public Iteration(final KRatioCorrectionModel2 krc) throws ArgumentException {
		mModel = krc;
	}

	final Map<Element, Number> iterate(final UncertainValues measKratios) //
			throws ArgumentException {
		final Material unkMat = mModel.getModel().getUnknownMaterial();
		Map<Element, Number> est = computeFirstEstimate(measKratios);
		for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
			final Composition comp = Composition.massFraction(unkMat, est);
			final UncertainValues inps = mModel.buildInput(comp, measKratios);
			final RealVector res = mModel.compute(inps.getValues(mModel.getInputLabels()));
			final UncertainValues calcKratios = KRatioLabel.extractKRatios( //
					res, mModel.getOutputLabels(), Method.Calculated);
			if (meetsThreshold(measKratios, calcKratios, THRESH))
				break;
			est = updateEstimate(est, measKratios, calcKratios);
		}
		return est;
	}

	/**
	 * This method performs the iteration and then computes an UncertainValues
	 * object containing the full results of the measurement.
	 *
	 *
	 * @param measKratios
	 * @return
	 * @throws ArgumentException
	 */
	public final UncertainValues compute(final UncertainValues measKratios) //
			throws ArgumentException {
		final Map<Element, Number> men = iterate(measKratios);
		final Material unkMat = mModel.getModel().getUnknownMaterial();
		final Composition comp = Composition.massFraction(unkMat, men);
		final UncertainValues inps = mModel.buildInput(comp, measKratios);
		return UncertainValues.propagate(mModel, inps);
	}

}
