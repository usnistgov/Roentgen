package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;

import com.duckandcover.html.Report;
import com.duckandcover.html.Table;
import com.duckandcover.html.Table.Item;

import gov.nist.juncertainty.CompositeMeasurementModel;
import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ImplicitMeasurementModel;
import gov.nist.juncertainty.TrimmedNamedMultivariateJacobianFunction;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.juncertainty.utility.FastIndex;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.math.MathUtilities;
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
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

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
public class KRatioCorrectionModel3 //
		extends CompositeMeasurementModel<EPMALabel> {

	static private final double THRESH = 1.0e-4;
	static private final double C_MIN = 1.0e-5;

	private final MatrixCorrectionModel2 mModel;

	private final ExplicitMeasurementModel<? extends MaterialLabel, MassFraction> mExtraElementRule;

	private final List<KRatioLabel> mKRatioSet;

	private IterationAlgorithm mIteration = new RecordingIterator(new WegsteinIteration());

	private int mMaxIterations = 100;
	private final int mNormIterations = 0;
	private int mConvergedIn = -1;
	private double mBestDeltaK = Double.MAX_VALUE;
	private Map<MassFraction, Double> mBestComposition;
	private Map<EPMALabel, Double> mBestCalculated;

	/**
	 * An interface for iteration algorithms such as SimpleIteration or PAPIteration
	 * to implement.
	 *
	 * @author Nicholas W. M. Ritchie
	 *
	 */
	public interface IterationAlgorithm {

		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calc
		) throws ArgumentException;

		public void reset();
	}

	/**
	 * Based off the description in G. Springer 1976 article "Iterative Procedures
	 * in Electron Probe Analysis Corrections"
	 *
	 * @author Nicholas W. M. Ritchie
	 *
	 */
	public static class WegsteinIteration //
			implements IterationAlgorithm {

		private Map<MassFraction, Double> mPrevResult = null;
		private Map<MassFraction, Double> mPrevFcn = null;

		@Override
		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calc
		) throws ArgumentException {
			final Map<MassFraction, Double> res = new HashMap<>();
			final Map<MassFraction, Double> cn = est;
			final Map<MassFraction, Double> fcn = new HashMap<>();
			if (mPrevResult == null) {
				// One step of simple iteration
				for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
					final MassFraction mfUnk = //
							MaterialLabel.buildMassFractionTag(measKrl.getUnknown().getMaterial(),
									measKrl.getElement());
					final double kCalc = calc.get(measKrl.as(Method.Calculated));
					final double kMeas = measKratios.getEntry(measKrl);
					// Simply gradient iteration
					final double cnu = cn.get(mfUnk);
					// ZAFMultiLineLabel zaf = MatrixCorrectionModel2.zafLabel(measKrl);
					// final double fcnu = calc.get(zaf); // cnu / kCalc;
					final double fcnu = cnu / kCalc;
					fcn.put(mfUnk, fcnu);
					final double cnp1u = kMeas * fcnu;
					res.put(mfUnk, cnp1u > 0.0 ? cnp1u : C_MIN);
				}
			} else {
				for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
					final MassFraction mfUnk = //
							MaterialLabel.buildMassFractionTag(measKrl.getUnknown().getMaterial(),
									measKrl.getElement());
					final double kCalc = calc.get(measKrl.as(Method.Calculated));
					final double cnu = cn.get(mfUnk), cnm1u = mPrevResult.get(mfUnk);
					// ZAFMultiLineLabel zaf = MatrixCorrectionModel2.zafLabel(measKrl);
					// final double fcnu = calc.get(zaf); // cnu / kCalc;
					final double fcnu = cnu / kCalc;
					final double fcnm1u = mPrevFcn.get(mfUnk);
					fcn.put(mfUnk, fcnu);
					final double fp = (fcnu - fcnm1u) / (cnu - cnm1u);
					final double kMeas = measKratios.getEntry(measKrl);
					final double cnp1u = cnu + (kMeas * fcnu - cnu) / (1.0 - Math.min(0.5, kMeas * fp));
					res.put(mfUnk, cnp1u > 0.0 ? cnp1u : C_MIN);
				}
			}
			mPrevResult = est;
			mPrevFcn = fcn;
			return res;
		}

		@Override
		public void reset() {
			mPrevResult = null;
			mPrevFcn = null;
		}

		@Override
		public String toString() {
			return "Wegstein iteration";
		}

	}

	public static class SimpleIteration //
			implements IterationAlgorithm {

		@Override
		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calcKratios //
				) throws ArgumentException {
			final Map<MassFraction, Double> res = new HashMap<>();
			for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
				final MassFraction mfUnk = MaterialLabel.buildMassFractionTag(measKrl.getUnknown().getMaterial(),
						measKrl.getElement());
				final double cEst = est.get(mfUnk);
				final double kCalc = calcKratios.get(measKrl.as(Method.Calculated));
				final double kMeas = measKratios.getEntry(measKrl);
				// Simply gradient iteration
				final double cNew = cEst * kMeas / kCalc;
				res.put(mfUnk, cNew > 0.0 ? cNew : C_MIN);
			}
			return res;
		}

		@Override
		public void reset() {
			// Don't do anything...
		}

		@Override
		public String toString() {
			return "simple iteration";
		}
	}

	public static class RecordingIterator implements IterationAlgorithm {

		private final IterationAlgorithm mBase;

		private UncertainValuesBase<KRatioLabel> mMeasured;
		private ArrayList<Map<MassFraction, Double>> mEstimate;
		private ArrayList<Map<EPMALabel, Double>> mCalculated;

		public RecordingIterator(
				final IterationAlgorithm baseIterator
				) {
			mBase = baseIterator;
			reset();
		}

		@Override
		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calcKratios
				) throws ArgumentException {
			if (mMeasured == null) {
				mMeasured = measKratios;
				mEstimate.add(est);
			}
			final Map<MassFraction, Double> res = mBase.nextEstimate(est, measKratios, calcKratios);
			mCalculated.add(calcKratios);
			mEstimate.add(res);
			return res;
		}

		@Override
		public void reset() {
			mMeasured = null;
			mEstimate = new ArrayList<>();
			mCalculated = new ArrayList<>();
		}

		public void dump(
				final Composition trueComp
				) throws IOException {
			final Report r = new Report("Recording iterator");
			r.addHeader(mBase.toString());
			final Table t = new Table();
			final BasicNumberFormat bnf = new BasicNumberFormat("0.0000");
			final List<MassFraction> mfs = new ArrayList<>();
			final List<KRatioLabel> krls = new ArrayList<>();
			for (final EPMALabel epma : mCalculated.get(0).keySet())
				if (epma instanceof KRatioLabel) {
					final KRatioLabel krl = (KRatioLabel) epma;
					krls.add(krl);
					mfs.add(MaterialLabel.buildMassFractionTag(krl.getUnknown().getMaterial(), krl.getElement()));
				}
			{
				final List<Item> hdr1 = new ArrayList<>();
				hdr1.add(Table.th(" "));
				hdr1.add(Table.th("Estimate", mfs.size()));
				hdr1.add(Table.th("Curr-Prev", mfs.size()));
				hdr1.add(Table.th("True-Prev", mfs.size()));
				hdr1.add(Table.th("Meas/Calc", krls.size()));
				t.addRow(hdr1);
			}
			{
				final List<Item> hdr = new ArrayList<>();
				hdr.add(Table.th("Step"));
				for (final MassFraction mf : mfs)
					hdr.add(Table.th(mf.getElement().getAbbrev()));
				for (final MassFraction mf : mfs)
					hdr.add(Table.th(mf.getElement().getAbbrev()));
				for (final MassFraction mf : mfs)
					hdr.add(Table.th(mf.getElement().getAbbrev()));
				for (final KRatioLabel krl : krls)
					hdr.add(Table.th(krl.getElement().getAbbrev()));
				t.addRow(hdr);
			}
			for (int i = 1; i < mEstimate.size(); ++i) {
				final List<Item> row = new ArrayList<>();
				final Map<MassFraction, Double> prev = mEstimate.get(i - 1);
				final Map<MassFraction, Double> curr = mEstimate.get(i);
				final Map<EPMALabel, Double> calc = mCalculated.get(i - 1);
				row.add(Table.td(i));
				for (final MassFraction mf : mfs)
					row.add(Table.td(bnf.format(prev.get(mf))));

				for (final MassFraction mf : mfs)
					row.add(Table.td(bnf.format(curr.get(mf) - prev.get(mf))));
				for (final MassFraction mf : mfs) {
					final double deltaTrue = trueComp.getMassFraction(mf.getElement()).doubleValue() - prev.get(mf);
					final boolean sameDir = Math.signum(curr.get(mf) - prev.get(mf)) == Math.signum(deltaTrue);
					if (sameDir)
						row.add(Table.td(bnf.format(deltaTrue)));
					else
						row.add(Table.td(bnf.format(deltaTrue), "background-color:pink"));
				}
				for (final KRatioLabel krl : krls)
					row.add(Table
							.td(bnf.format(mMeasured.getEntry(krl.as(KRatioLabel.Method.Measured)) / calc.get(krl))));
				t.addRow(row);
				row.clear();
			}
			r.add(t);
			r.inBrowser(Mode.NORMAL);
		}
	}

	/**
	 * Computes the matrix correction factors for the standard relative to the pure
	 * element.
	 *
	 * @param measKratios
	 * @return
	 * @throws ArgumentException
	 */

	private static Map<Element, Double> computeZAFsp(
			final Collection<KRatioLabel> measKratios
			) throws ArgumentException {
		final Map<Element, Double> res = new HashMap<>();
		for (final KRatioLabel krl : measKratios) {
			final Set<KRatioLabel> stdKrs = new HashSet<>();
			final List<EPMALabel> outputLabels = new ArrayList<>();
			final StandardMatrixCorrectionDatum smcd = krl.getStandard();
			if (smcd.getMaterial().getElementSet().size() > 1) {
				final UnknownMatrixCorrectionDatum umcd = new UnknownMatrixCorrectionDatum(smcd, smcd.getMaterial());
				final Composition pureElm = Composition.pureElement(krl.getElement());
				final StandardMatrixCorrectionDatum pureStd = new StandardMatrixCorrectionDatum(smcd, pureElm);
				final KRatioLabel newKrl = new KRatioLabel(umcd, pureStd, krl.getXRaySet(), Method.Measured);
				stdKrs.add(newKrl);
				final ZAFMultiLineLabel zafLabel = MatrixCorrectionModel2.zafLabel(newKrl);
				outputLabels.add(zafLabel);
				final XPPMatrixCorrection2 xpp = new XPPMatrixCorrection2(stdKrs, outputLabels);
				final UncertainValues<EPMALabel> inputs = xpp.buildInput(umcd.getMaterial());

				xpp.addConstraints(xpp.buildConstraints(inputs));
				xpp.addAdditionalInputs(smcd.getComposition().getValueMap(MassFraction.class));

				final UncertainValuesCalculator<EPMALabel> uvc = new UncertainValuesCalculator<>(xpp, inputs);
				final Map<EPMALabel, Double> valMap = uvc.getValueMap();
				res.put(newKrl.getElement(), valMap.get(zafLabel));
			} else
				res.put(krl.getElement(), 1.0);
		}
		return res;
	}

	public static class PAP1991Iteration //
			implements IterationAlgorithm {

		private Map<Element, Double> mZAFsp = null;

		@Override
		public void reset() {
			mZAFsp = null;
		}

		@Override
		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calc
				) throws ArgumentException {
			if (mZAFsp == null)
				mZAFsp = computeZAFsp(measKratios.getLabels());
			final Map<MassFraction, Double> res = new HashMap<>();
			for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
				final MassFraction mfUnk = MaterialLabel.buildMassFractionTag(measKrl.getUnknown().getMaterial(),
						measKrl.getElement());
				final MassFraction mfStd = MaterialLabel.buildMassFractionTag(measKrl.getStandard().getMaterial(),
						measKrl.getElement());
				final double cEst = est.get(mfUnk);
				final double cStd = measKrl.getStandard().getComposition().getEntry(mfStd);
				final double kCalc = calc.get(measKrl.as(Method.Calculated));
				final double kMeas = measKratios.getEntry(measKrl);
				boolean fallback = true;
				// Pouchou and Pichor iteration from PAP1991 Green book
				// Calculate alpha using the previous step information
				assert mZAFsp != null;
				final double ZAFsp = mZAFsp.get(measKrl.getElement());
				// Compute k for cEst relative to pure element....
				final double kcp = kCalc * cStd * ZAFsp;
				if (Math.abs(1.0 - cEst) > 1.0e-4) {
					final double alpha = (cEst * (1.0 - kcp)) / ((1.0 - cEst) * kcp);
					final double kmp = kMeas * cStd * ZAFsp;
					if (kcp < cEst) {
						// Hyperbolic: C/k = alpha + (1-alpha) C
						final double cNew = (alpha * kmp) / (1.0 + (alpha - 1.0) * kmp);
						if (Double.isFinite(cNew)) {
							res.put(mfUnk, cNew > 0.0 ? cNew : 0.0);
							fallback = false;
						}
					} else {
						// Parabolic: k/C = alpha + (1-alpha) C
						final double[] qs = MathUtilities.quadraticSolver(1.0 - alpha, alpha, -kmp);
						if (qs != null) {
							// assert qs[1] >= 0.0;
							final double cNew = qs[1] >= 0.0 ? qs[1] : qs[0];
							assert cNew >= 0.0 : qs[0] + ", " + qs[1];
							if (cNew >= 0.0) {
								res.put(mfUnk, cNew);
								fallback = false;
							}
						}
					}
				}
				if (fallback) {
					// Simply gradient iteration
					final double cNew = cEst * kMeas / kCalc;
					res.put(mfUnk, cNew > 0.0 ? cNew : 0.0);
				}
			}
			return res;
		}

		@Override
		public String toString() {
			return "PAP1991 iteration";
		}
	}

	public static class AlphaPlusSimpleIteration //
			implements IterationAlgorithm {

		private Map<Element, Double> mZAFsp = null;

		@Override
		public void reset() {
			mZAFsp = null;
		}

		@Override
		public Map<MassFraction, Double> nextEstimate(
				final Map<MassFraction, Double> est, //
				final UncertainValuesBase<KRatioLabel> measKratios, //
				final Map<EPMALabel, Double> calc
				) throws ArgumentException {
			if (mZAFsp == null)
				mZAFsp = computeZAFsp(measKratios.getLabels());
			final Map<MassFraction, Double> res = new HashMap<>();
			final StringBuilder mode = new StringBuilder();
			for (final KRatioLabel measKrl : measKratios.getLabels(KRatioLabel.class)) {
				final MassFraction mfUnk = MaterialLabel.buildMassFractionTag(measKrl.getUnknown().getMaterial(),
						measKrl.getElement());
				final MassFraction mfStd = MaterialLabel.buildMassFractionTag(measKrl.getStandard().getMaterial(),
						measKrl.getElement());
				final double cEst = est.get(mfUnk);
				final double cStd = measKrl.getStandard().getComposition().getEntry(mfStd);
				final double kCalc = calc.get(measKrl.as(Method.Calculated));
				final double kMeas = measKratios.getEntry(measKrl);
				boolean fallback = true;
				// Calculate alpha factor using the previous step information
				if (Math.abs(1.0 - cEst) > 1.0e-4) {
					assert mZAFsp != null;
					final double ZAFsp = mZAFsp.get(measKrl.getElement());
					// Compute k relative to pure element....
					final double kcp = kCalc * cStd * ZAFsp, kmp = kMeas * cStd * ZAFsp;
					final double alpha = (cEst * (1.0 - kcp)) / ((1.0 - cEst) * kcp);
					if (kcp / cEst <= 1.0) {
						// Hyperbolic: C/k = alpha + (1-alpha) C
						final double cNew = (alpha * kmp) / (1.0 + (alpha - 1.0) * kmp);
						if (Double.isFinite(cNew)) {
							res.put(mfUnk, cNew > 0.0 ? cNew : 0.0);
							mode.append(mfUnk.getElement().getAbbrev() + "=H,");
							fallback = false;
						}
					}
				}
				if (fallback) {
					// Simply gradient iteration
					final double cNew = cEst * kMeas / kCalc;
					res.put(mfUnk, cNew > 0.0 ? cNew : 0.0);
					mode.append(mfUnk.getElement().getAbbrev() + "=S,");
				}
			}
			return res;
		}

		@Override
		public String toString() {
			return "Alpha+Simple iteration";
		}
	}

	private static class KRatioModel2 //
			extends ImplicitMeasurementModel<EPMALabel, MassFraction> {

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
		public RealMatrix computeJx(
				final RealVector point
				) {
			final RealMatrix jac = buildEmptyJx();
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
				final double kMeas = getArg(kMeasTag, point);
				// h(x,y) = k_meas*cStd - cUnk*zaf
				setJx(oHTag, kMeasTag, jac, cStd);
				setJx(oHTag, mfStdTag, jac, kMeas);
				setJx(oHTag, zafTag, jac, -cUnk);
			}
			return jac;
		}

		@Override
		public RealVector computeH(
				final RealVector point
				) {
			final RealVector res = buildEmptyH();
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
				final double kMeas = getArg(kMeasTag, point);
				// h(x,y) = k_meas*C_std - Cunk*zaf
				setH(oHTag, res, kMeas * cStd - cUnk * zaf);
			}
			return res;
		}

		@Override
		public RealMatrix computeJy(
				final RealVector point
				) {
			final RealMatrix jac = buildEmptyJy();
			for (int oHTag = 0; oHTag < mKRatios.size(); ++oHTag) {
				final KRatioLabel kMeasTag = mKRatios.get(oHTag);
				final Element elm = kMeasTag.getElement();
				final UnknownMatrixCorrectionDatum unknown = kMeasTag.getUnknown();
				// final StandardMatrixCorrectionDatum standard = kMeasTag.getStandard();
				final MassFraction mfUnkTag = MaterialLabel.buildMassFractionTag(unknown.getMaterial(), elm);
				// final MassFraction mfStdTag =
				// MaterialLabel.buildMassFractionTag(standard.getMaterial(), elm);
				final ZAFMultiLineLabel zafTag = MatrixCorrectionModel2.zafLabel(kMeasTag);
				// final double cStd = getArg(mfStdTag, point);
				final double zaf = getArg(zafTag, point);
				// h(x,y) = k_meas*C_std - Cunk*zaf
				setJy(oHTag, mfUnkTag, jac, -zaf);
			}
			return jac;
		}

		@Override
		public String toString() {
			return "K-ratio H-model";
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

	public KRatioCorrectionModel3(
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

	public KRatioCorrectionModel3(
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
	 * {@link KRatioCorrectionModel3} instance.
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
		final Set<Composition> stds = new HashSet<>();
		for (final KRatioLabel krl : measKratios.getLabels()) {
			final Material mat = krl.getUnknown().getMaterial();
			if (unkMat == null)
				unkMat = mat;
			if (!unkMat.equals(mat))
				throw new ArgumentException("Two different materials [" + mat + " and " + unkMat
						+ " are present in the estimated unknown.");
			stds.add(krl.getStandard().getComposition());
		}
		final List<UncertainValuesBase<? extends EPMALabel>> list = new ArrayList<>();
		final UncertainValuesBase<EPMALabel> modInps = mModel.buildInput(unkMat);
		list.add(modInps);
		list.add(measKratios);
		list.add(unkMat.getAtomicWeights());
		for (final KRatioLabel krl : measKratios.getLabels()) {
			final List<MassFraction> lmf = new ArrayList<>();
			for (final MassFraction mf : krl.getStandard().getComposition().massFractionTags())
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
	 * {@link KRatioCorrectionModel3}. It returns a KRatioCorrectionModel2.
	 *
	 * @param keySet
	 * @param extraElms
	 * @param outputLabels
	 * @return {@link KRatioCorrectionModel3}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel3 buildXPPModel(
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
		return new KRatioCorrectionModel3(keySet, xpp, extraElms, outputs);
	}

	/**
	 * Builds a {@link ExplicitMeasurementModel} which sequentially combines the
	 * explicit {@link XPPMatrixCorrection2} model with the implicit model
	 * {@link KRatioCorrectionModel3}. It returns a KRatioCorrectionModel2.
	 *
	 * The model returns the values specified by buildDefaultOutputs(...).
	 *
	 * @param keySet
	 * @param extraElms
	 * @return {@link KRatioCorrectionModel3}
	 * @throws ArgumentException
	 */
	static public KRatioCorrectionModel3 buildXPPModel(
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
			outputs.addIfMissing(krl.as(Method.Calculated));
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
			kCalc.add(krl.as(Method.Calculated));
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
			final MaterialLabel.MassFraction smfTag = //
					MaterialLabel.buildMassFractionTag(std.getMaterial(), xrs.getElement());
			final double cStd = std.getEntry(smfTag);
			assert cStd >= 0.0;
			final double kR = me.getValue().doubleValue();
			est.put(MaterialLabel.buildMassFractionTag(unkMat, xrs.getElement()), Math.max(0.0, kR * cStd));
		}
		return applyExtraElementRule(unkMat, point, est);
	}

	private Map<MassFraction, Double> applyExtraElementRule(
			final Material unkMat, //
			final RealVector point, //
			final Map<MassFraction, Double> est
			) {
		if (mExtraElementRule != null) {
			final RealVector pt = new ArrayRealVector(mExtraElementRule.getInputDimension());
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

	public double calculateDeltaK(
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final Map<EPMALabel, Double> calc
			) {
		double total = 0.0;
		for (final Object lbl : measKratios.getLabels()) {
			assert lbl instanceof KRatioLabel;
			final KRatioLabel measKrl = (KRatioLabel) lbl;
			final double delta = measKratios.getEntry(measKrl) - calc.get(measKrl.as(Method.Calculated));
			total += delta * delta;
		}
		return Math.sqrt(total);
	}

	/**
	 * Updates the estimated composition by first applying the iteration algorithm
	 * to the k-ratio data and then applying the extra element algorithm.
	 *
	 * @param unkMat      The unknown material
	 * @param est         An estimate of the composition of the unknown material
	 * @param measKratios The measured k-ratios
	 * @param calcKratios The calculated k-ratios associated with the estimated
	 *                    composition
	 * @param point       The point containing input to this from which values for
	 *                    the extra-element function will be extracted.
	 * @return Map&lt;MassFraction, Double&gt;
	 * @throws ArgumentException
	 */
	public Map<MassFraction, Double> updateEstimate(
			final Material unkMat, //
			final Map<MassFraction, Double> est, //
			final UncertainValuesBase<KRatioLabel> measKratios, //
			final Map<EPMALabel, Double> calc, //
			final RealVector point
			) throws ArgumentException {

		final Map<MassFraction, Double> res = mIteration.nextEstimate(est, measKratios, calc);
		return applyExtraElementRule(unkMat, point, res);
	}

	/**
	 * Takes the <code>unknown</code> and computes the k-ratios relative to the
	 * specified set of {@link KRatioLabel}s. Then takes the k-ratio values as input
	 * to the iteration algorithm, finally returning the {@link Composition} from
	 * the iteration algorithm.
	 *
	 * @param unknown Initial composition
	 * @param kratios The set of k-ratios
	 * @return The iterative estimate of <code>unknown</code>
	 * @throws ArgumentException
	 */
	public static Composition roundTripXPP(
			final Composition unknown, //
			final Set<KRatioLabel> kratios
			) throws ArgumentException {
		final UncertainValuesCalculator<EPMALabel> xpp = //
				XPPMatrixCorrection2.buildAnalytical(kratios, unknown.toMassFraction().getValueMap(), false);
		final List<KRatioLabel> asMeas = new ArrayList<>(kratios);
		final List<KRatioLabel> asCalc = KRatioLabel.as(asMeas, KRatioLabel.Method.Calculated);
		final UncertainValuesBase<KRatioLabel> calcKrs = //
				UncertainValues.asUncertainValues(xpp).extract(asCalc);
		final int size = asMeas.size();
		final UncertainValues<KRatioLabel> measKrs = //
				new UncertainValues<>(asMeas, calcKrs.getValues(), MatrixUtils.createRealMatrix(size, size));
		final KRatioCorrectionModel3 iter = KRatioCorrectionModel3.buildXPPModel(new HashSet<>(asMeas), null);
		final UncertainValues<EPMALabel> uvs = UncertainValues.asUncertainValues(iter.iterate(measKrs));
		return Composition.massFraction(iter.getUnknownMaterial(), uvs);
	}
	
	
	/**
	 * Takes the <code>unknown</code> and computes the k-ratios relative to the
	 * specified set of {@link KRatioLabel}s. Then computes the Composition with
	 * the associated uncertainties from the k-ratios.
	 *
	 * @param unknown Initial composition
	 * @param kratios The set of k-ratios
	 * @return The iterative estimate of <code>unknown</code>
	 * @throws ArgumentException
	 */
	public static Composition computeXPP(
			final Composition unknown, //
			final Set<KRatioLabel> kratios
			) throws ArgumentException {
		final UncertainValuesCalculator<EPMALabel> xpp = //
				XPPMatrixCorrection2.buildAnalytical(kratios, unknown.toMassFraction().getValueMap(), false);
		final List<KRatioLabel> asMeas = new ArrayList<>(kratios);
		final List<KRatioLabel> asCalc = KRatioLabel.as(asMeas, KRatioLabel.Method.Calculated);
		final UncertainValuesBase<KRatioLabel> calcKrs = //
				UncertainValues.asUncertainValues(xpp).extract(asCalc);
		final int size = asMeas.size();
		final UncertainValues<KRatioLabel> measKrs = //
				new UncertainValues<>(asMeas, calcKrs.getValues(), MatrixUtils.createRealMatrix(size, size));
		final KRatioCorrectionModel3 krcm = KRatioCorrectionModel3.buildXPPModel(new HashSet<>(asMeas), null);
		final UncertainValuesBase<EPMALabel> inps = krcm.buildInput(measKrs);
		krcm.addAdditionalInputs(unknown.getValueMap());
		UncertainValuesCalculator<EPMALabel> res = new UncertainValuesCalculator<EPMALabel>(krcm, inps);
		return Composition.massFraction(unknown.getMaterial(), res);
	}

	
	

	private Map<EPMALabel, Double> computeOptimized(
			final RealVector point
	) {
		final List<EPMALabel> labels = getOutputLabels();
		final Map<EPMALabel, Double> res = new HashMap<>();
		final RealVector output = computeValue(point.toArray());
		for (int i = 0; i < output.getDimension(); ++i)
			res.put(labels.get(i), output.getEntry(i));
		return res;
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
		mBestDeltaK = Double.MAX_VALUE;
		mBestComposition = est;
		mBestCalculated = null;
		mIteration.reset();
		mConvergedIn = -1;
		for (int iteration = 0; iteration < mMaxIterations; ++iteration) {
			// Computes the h-model k_meas - k_calc
			addAdditionalInputs(est);
			final Map<EPMALabel, Double> res = computeOptimized(point);
			// Determines how large this difference is...
			final double dk = calculateDeltaK(measKratios, res);
			if (dk < mBestDeltaK) {
				mBestDeltaK = dk;
				mBestComposition = est;
				mBestCalculated = res;
			} else {
				System.out.println("Deteriorated in step " + iteration);
				System.out.println("Best = " + mBestDeltaK + " - This = " + dk);
			}
			if (dk < THRESH) {
				System.out.println("Convergence at iteration " + iteration);
				mConvergedIn = iteration;
				break;
			}
			est = updateEstimate(unkMat, est, measKratios, res, point);
			if (iteration < mNormIterations)
				est = MathUtilities.normalize(est);
			if (iteration == mMaxIterations - 1)
				System.out.println("Failed to convergence in " + mMaxIterations + " iterations using " + mIteration);
		}
		return mBestComposition;
	}

	@Override
	public int getStepCount() {
		return mConvergedIn;
	}

	public boolean converged() {
		return mConvergedIn != -1;
	}

	public Map<EPMALabel, Double> getBestCalculated() {
		return mBestCalculated;
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

	public int getMaxIterations() {
		return mMaxIterations;
	}

	public void setMaxIterations(
			final int mmaxiterations
	) {
		mMaxIterations = mmaxiterations;
	}

	public IterationAlgorithm getIterationAlgorithm() {
		return mIteration;
	}

	public void setIterationAlgorithm(
			final IterationAlgorithm iteration
	) {
		mIteration = iteration;
	}
}