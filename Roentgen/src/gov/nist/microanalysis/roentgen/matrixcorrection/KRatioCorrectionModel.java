package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.math.uncertainty.ImplicitMeasurementModel;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.XRaySet.ElementXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * Solves the implicit measurement model uncertainty problem presented by the
 * standard k-ratio protocol for matrix correction:
 * <p>
 * 0 = k<sub>z</sub> -
 * (C<sub>unk,z</sub>/C<sub>std,z</sub>)ZAF<sub>unk,std,z</sub>
 * </p>
 * <p>
 * for z in { Elements }
 * </p>
 * 
 * @author nicholas
 *
 */
public class KRatioCorrectionModel extends ImplicitMeasurementModel {

	private static class hTag extends BaseTag<MatrixCorrectionDatum, ElementXRaySet, Object> {

		protected hTag(MatrixCorrectionDatum unk, ElementXRaySet obj2) {
			super("h", unk, obj2);
		}
	}

	/**
	 * Builds a list of h<sub>z</sub> tags.
	 * 
	 * @param unk
	 * @param stds
	 * @return List&lt;? extends Object&gt;
	 */
	private static List<Object> buildOutputs( //
			MatrixCorrectionDatum unk, //
			Map<ElementXRaySet, MatrixCorrectionDatum> stds //
	) {
		List<Object> res = new ArrayList<>();
		for (ElementXRaySet exrs : stds.keySet())
			res.add(new hTag(unk, exrs));
		return res;
	}

	/**
	 * <p>
	 * Builds a list of C<sub>std,z</sub>, C<sub>unk,z</sub>,
	 * ZAF<sub>unk,std,z</sub> and k<sub>z</sub> for each z in the list of elements
	 * in the unknown as is required to solve
	 * </p>
	 * <p>
	 * 0 = k<sub>z</sub> -
	 * (C<sub>unk,z</sub>/C<sub>std,z</sub>)ZAF<sub>unk,std,z</sub>
	 * </p>
	 * 
	 * @param unkMcd
	 * @param stds
	 * @return List&lt;? extends Object&gt;
	 */
	private static List<? extends Object> buildInputs( //
			MatrixCorrectionDatum unkMcd, //
			Map<ElementXRaySet, MatrixCorrectionDatum> stds //
	) {

		Set<Object> tags = new HashSet<>();
		for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : stds.entrySet()) {
			final ElementXRaySet exrs = me.getKey();
			final MatrixCorrectionDatum stdMcd = me.getValue();
			// Get Cstd,i
			final Composition std = stdMcd.getComposition();
			final Element elm = exrs.getElement();
			final Object smft = Composition.buildMassFractionTag(std, elm);
			tags.add(smft);
			// Get Cunk,i
			final Object umft = Composition.buildMassFractionTag(unkMcd.getComposition(), elm);
			tags.add(umft);
			// Get ZAFi
			final Object mct = new MatrixCorrectionTag(unkMcd, stdMcd, exrs);
			tags.add(mct);
			// Get ki
			final KRatioTag krt = new KRatioTag(unkMcd, stdMcd, exrs);
			tags.add(krt);
		}
		return new ArrayList<>(tags);
	}

	private static class UnknownModel extends NamedMultivariateJacobianFunction {
		private final MatrixCorrectionDatum mUnknown;
		private final Map<ElementXRaySet, MatrixCorrectionDatum> mStandards;

		private UnknownModel(MatrixCorrectionDatum unk, Map<ElementXRaySet, MatrixCorrectionDatum> stds) {
			super(buildInputs(unk, stds), buildOutputs(unk, stds));
			mUnknown = unk;
			mStandards = stds;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
			// for each elmxrayset calculate
			// hi = ki - (Cunk,i/Cstd,i) ZAFunk,std,i
			// dhi/dcunkj = { (-1/Cstdi) ZAFunk,std,i for i=j, 0 otherwise }
			RealVector rv = new ArrayRealVector(getOutputDimension());
			RealMatrix rm = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());

			for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : mStandards.entrySet()) {
				final ElementXRaySet exrs = me.getKey();
				final MatrixCorrectionDatum stdMcd = me.getValue();
				// For all the elements for which this material is a standard
				// Get Cstd,i
				final Composition std = stdMcd.getComposition();
				final Object smft = Composition.buildMassFractionTag(std, exrs.getElement());
				final double c_stdi = point.getEntry(inputIndex(smft));
				// Get Cunk,i
				final Object umft = Composition.buildMassFractionTag(mUnknown.getComposition(), exrs.getElement());
				final double c_unki = point.getEntry(inputIndex(umft));
				// Get ZAFi
				final Object mct = new MatrixCorrectionTag(mUnknown, stdMcd, exrs);
				final double zafi = point.getEntry(inputIndex(mct));
				// Get ki
				final KRatioTag krt = new KRatioTag(mUnknown, stdMcd, exrs);
				final double ki = point.getEntry(inputIndex(krt));
				final int row = outputIndex(new hTag(mUnknown, exrs));
				// Compute output value
				final double hi = ki - (c_unki / c_stdi) * zafi;
				rv.setEntry(row, hi);
				// Compute partials
				final double dhidCunki = (-1.0 / c_stdi) * zafi;
				rm.setEntry(row, inputIndex(umft), dhidCunki);
			}
			return Pair.create(rv, rm);
		}
	}

	private static class StandardsModel extends NamedMultivariateJacobianFunction {

		private final MatrixCorrectionDatum mUnknown;
		private final Map<ElementXRaySet, MatrixCorrectionDatum> mStandards;

		private StandardsModel(MatrixCorrectionDatum unk, Map<ElementXRaySet, MatrixCorrectionDatum> stds)
				throws ArgumentException {
			super(buildInputs(unk, stds), buildOutputs(unk, stds));
			mUnknown = unk;
			mStandards = stds;
		}

		@Override
		public Pair<RealVector, RealMatrix> value(RealVector point) {
			// for each elmxrayset calculate
			// hi = ki - (Cunk,i/Cstd,i) ZAFunk,std,i
			// dhi/dstdj = { (Cunk/(Cstdi^2)) ZAFunk,std,i for i=j, 0 otherwise }
			// dhi/kj = { 1 for i=j, 0 otherwise }
			// dhi/kZAFunk,std,i = { - (Cunk,i/Cstd,i) for i=j, 0 otherwise }
			RealVector rv = new ArrayRealVector(getOutputDimension());
			RealMatrix rm = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());

			for (Map.Entry<ElementXRaySet, MatrixCorrectionDatum> me : mStandards.entrySet()) {
				final ElementXRaySet exrs = me.getKey();
				final MatrixCorrectionDatum stdMcd = me.getValue();
				// For all the elements for which this material is a standard
				// Get Cstd,i
				final Composition std = stdMcd.getComposition();
				final Object smft = Composition.buildMassFractionTag(std, exrs.getElement());
				final double c_stdi = point.getEntry(inputIndex(smft));
				// Get Cunk,i
				final Object umft = Composition.buildMassFractionTag(mUnknown.getComposition(), exrs.getElement());
				final double c_unki = point.getEntry(inputIndex(umft));
				// Get ZAFi
				final Object mct = new MatrixCorrectionTag(mUnknown, stdMcd, exrs);
				final double zafi = point.getEntry(inputIndex(mct));
				// Get ki
				final KRatioTag krt = new KRatioTag(mUnknown, stdMcd, exrs);
				final double ki = point.getEntry(inputIndex(krt));
				final int row = outputIndex(new hTag(mUnknown, exrs));
				// Compute output value
				final double hi = ki - (c_unki / c_stdi) * zafi;
				rv.setEntry(row, hi);
				// Compute partials
				final double dhidCstdi = (c_unki / Math.pow(c_stdi, 2.0)) * zafi;
				final double dhidzafi = (-c_unki / c_stdi);
				final double dhidki = 1.0;
				rm.setEntry(row, inputIndex(smft), dhidCstdi);
				rm.setEntry(row, inputIndex(mct), dhidzafi);
				rm.setEntry(row, inputIndex(krt), dhidki);
			}
			return Pair.create(rv, rm);
		}

	}

	public static List<? extends Object> buildOutputTags(MatrixCorrectionDatum unkMcd) {
		final Composition unk = unkMcd.getComposition();
		ArrayList<Object> res = new ArrayList<>();
		for (Element elm : unk.getElementSet())
			res.add(Composition.buildMassFractionTag(unk, elm));
		return res;
	}

	public KRatioCorrectionModel(MatrixCorrectionDatum unk, Map<ElementXRaySet, MatrixCorrectionDatum> stds)
			throws ArgumentException {
		super(new UnknownModel(unk, stds), new StandardsModel(unk, stds), buildOutputTags(unk));
		final Set<Element> selm1 = new HashSet<Element>();
		final Set<Element> selm2 = new HashSet<Element>(unk.getComposition().getElementSet());
		for (ElementXRaySet exrs : stds.keySet()) {
			final Element elm = exrs.getElement();
			if (selm1.contains(elm))
				throw new ArgumentException("There are duplicate standards for the element " + elm);
			selm1.add(elm);
			if (!selm2.contains(elm))
				throw new ArgumentException("There is a standard for " + elm + " which is not present in the unknown.");
			selm2.remove(elm);
		}
		if (!selm2.isEmpty())
			throw new ArgumentException(
					"The elements " + selm2 + " are present in the unknow but there is no standard.");
	}
}
