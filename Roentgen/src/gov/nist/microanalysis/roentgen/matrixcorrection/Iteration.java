package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.RealVector;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesBase;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel.Method;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel.MassFraction;

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
	 * @param measKratios
	 * @return {@link Composition}
	 * @throws ArgumentException
	 */
	public Map<MassFraction, Number> computeFirstEstimate(final Material unkMat,
			final UncertainValuesBase measKratios) {
		final Map<Object, Double> krm = measKratios.getValueMap();
		final Map<MassFraction, Number> est = new HashMap<>();
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
				est.put(MaterialLabel.buildMassFractionTag(unkMat, xrs.getElement()), kR * cStd);
			}
		}
		return est;
	}

	public boolean meetsThreshold(final UncertainValuesBase measKratios, final UncertainValuesBase calcKratios,
			final double thresh) {
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

	public Map<MassFraction, Number> updateEstimate(//
			final Map<MassFraction, Number> est, //
			final UncertainValuesBase measKratios, //
			final UncertainValuesBase calcKratios) {
		final Map<MassFraction, Number> res = new HashMap<>();
		for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
			final Material unkMat = measKrl.getUnknown().getMaterial();
			final MassFraction mf = MaterialLabel.buildMassFractionTag(unkMat, measKrl.getElement());
			final double cEst = est.get(mf).doubleValue();
			final double cNext = cEst * measKratios.getEntry(measKrl) //
					/ calcKratios.getEntry(measKrl.asCalculated());
			res.put(mf, cNext);
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

	final Map<MassFraction, Number> iterate(final UncertainValuesBase measKratios) //
			throws ArgumentException {
		final Material unkMat = mModel.getModel().getUnknownMaterial();
		Map<MassFraction, Number> est = computeFirstEstimate(unkMat, measKratios);
		for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
			final UncertainValuesBase inps = mModel.buildInput(est, measKratios);
			final RealVector res = mModel.compute(inps.getValues(mModel.getInputLabels()));
			final UncertainValuesBase calcKratios = KRatioLabel.extractKRatios( //
					res, mModel.getOutputLabels(), Method.Calculated);
			if (meetsThreshold(measKratios, calcKratios, THRESH))
				break;
			est = updateEstimate(est, measKratios, calcKratios);
		}
		return est;
	}

	/**
	 * This method performs the iteration and then computes an UncertainValuesBase
	 * object containing the full results of the measurement.
	 *
	 *
	 * @param measKratios
	 * @return
	 * @throws ArgumentException
	 */
	public final UncertainValuesBase compute(final UncertainValuesBase measKratios) //
			throws ArgumentException {
		final Map<MassFraction, Number> men = iterate(measKratios);
		final UncertainValuesBase inps = mModel.buildInput(men, measKratios);
		return UncertainValuesBase.propagate(mModel, inps);
	}

}
