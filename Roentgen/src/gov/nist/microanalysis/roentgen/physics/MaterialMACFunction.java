package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.juncertainty.ExplicitMeasurementModel;
import gov.nist.juncertainty.ParallelMeasurementModelBuilder;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.EPMALabel.MaterialMAC;
import gov.nist.microanalysis.roentgen.matrixcorrection.KRatioLabel;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.utility.FastIndex;

/**
 * Compute the material mass absorption coefficient given the elemental mass
 * absorption coefficients.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class MaterialMACFunction //
		extends ExplicitMeasurementModel<EPMALabel, MaterialMAC> {

	public static class DefaultDecorrelationFunction implements DecorrelationFunction {

		private final double mAbove;
		private final double mBelow;

		public DefaultDecorrelationFunction(
				final double above, final double below
		) {
			mAbove = above;
			mBelow = Math.abs(below);
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see gov.nist.microanalysis.roentgen.physics.MaterialMACFunction.
		 * DecorrelationFunction#compute(double, double)
		 */
		@Override
		public double compute(
				final AtomicShell shell, final XRay xray
		) {
			final double MIN_VAL = 0.1;
			double res = 1.0;
			final double delta = xray.getEnergy() - shell.getEdgeEnergy();
			if ((delta > -3.0 * mBelow) && (delta < mAbove)) {
				if (delta > 0.0) // xrayE above edge
					res = MIN_VAL + (1.0 - MIN_VAL) * delta / Math.abs(mAbove);
				else // xrayE below edge
					res = 1.0 - (1.0 - MIN_VAL) * Math.exp(-0.5 * Math.pow(delta / mBelow, 2));
			}
			assert res >= MIN_VAL : res;
			assert res <= 1.0 : res;
			return Math.max(0.0, Math.min(1.0, res * res));
		}
	}

	public static class TotalDecorrelationFunction implements DecorrelationFunction {

		private final double mAbove;
		private final double mBelow;

		public TotalDecorrelationFunction(
				final double above, final double below
		) {
			mAbove = above;
			mBelow = Math.abs(below);
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see gov.nist.microanalysis.roentgen.physics.MaterialMACFunction.
		 * DecorrelationFunction#compute(double, double)
		 */
		@Override
		public double compute(
				final AtomicShell shell, final XRay xray
		) {
			double res = 1.0;
			final double delta = xray.getEnergy() - shell.getEdgeEnergy();
			if ((delta > -2.0 * mBelow) && (delta < mAbove))
				res = 0.0;
			return Math.max(0.0, Math.min(1.0, res * res));
		}
	}

	interface DecorrelationFunction {

		public double compute(
				AtomicShell edgeEnergy, XRay xray
		);

	}

	/**
	 * Constructs a MassAbsorptionCoefficient
	 *
	 * @param comps A list of compositions to compute
	 * @param xrays A list of x-ray energies to compute
	 * @throws ArgumentException
	 */
	public static MaterialMACFunction build(
			final List<Composition> comps, //
			final XRay xray
	) throws ArgumentException {
		final List<Material> mats = new ArrayList<>();
		for (final Composition comp : comps)
			mats.add(comp.getMaterial());
		return new MaterialMACFunction(mats, Collections.emptyList(), xray);
	}

	final static public UncertainValuesCalculator<EPMALabel> compute(
			final List<Composition> comps, //
			final XRay xray
	) throws ArgumentException {
		final List<Material> mats = new ArrayList<>();
		for (final Composition comp : comps)
			mats.add(comp.getMaterial());
		final MaterialMACFunction mmf = new MaterialMACFunction(mats, Collections.emptyList(), xray);
		final UncertainValuesBase<EPMALabel> inputs = mmf.buildInputs(comps, xray);
		return UncertainValuesBase.propagateAnalytical(mmf, inputs);
	}

	/**
	 *
	 * <p>
	 * This function checks the MatrialMAC values to reduce the correlation between
	 * materials that 1) share a common element; and 2) the edge of which is close
	 * in energy to the X-ray energy.
	 * </p>
	 *
	 * <p>
	 * The algorithm considers all the edges associated with all the elements in the
	 * material to determine how close the x-ray energy is to an edge and then
	 * computes a Decorrelation function to determine how much to reduce the
	 * correlation.
	 * </p>
	 *
	 *
	 * @param uvs   An UncertainValues containing more than one MaterialMAC tag
	 * @param xray  X-ray line
	 * @param decor DecorrelationFunction
	 * @return A new UncertainValues object
	 */
	public static UncertainValues<EPMALabel> decorrelate(
			//
			final UncertainValues<EPMALabel> uvs, final DecorrelationFunction decor
	) {
		final UncertainValues<EPMALabel> res = uvs.copy();
		final List<EPMALabel> labels = res.getLabels();
		final int labelCx = labels.size();
		for (int i = 0; i < labelCx; ++i) {
			final Object labeli = labels.get(i);
			if (labeli instanceof MaterialMAC) {
				double mult = 1.0;
				final MaterialMAC mmi = (MaterialMAC) labeli;
				for (final Element elm : mmi.getMaterial().getElementSet()) {
					for (final AtomicShell sh : AtomicShell.forElement(elm))
						mult = Math.min(mult, decor.compute(sh, mmi.getXRay()));
					assert (mult >= 0.0) && (mult <= 1.0) : mult;
					if (mult < 1.0)
						for (int j = i + 1; j < labelCx; ++j)
							if (labels.get(j) instanceof MaterialMAC) {
								final MaterialMAC mmj = (MaterialMAC) labels.get(j);
								if (mmj.getMaterial().getElementSet().contains(elm))
									res.setCovariance(i, j, mult * res.getCovariance(i, j));
							}
				}
			}
		}
		return res;
	}

	/**
	 *
	 * <p>
	 * This function checks the MatrialMAC values to reduce the correlation between
	 * materials that 1) share a common element; and 2) the edge of which is close
	 * in energy to the X-ray energy.
	 * </p>
	 *
	 * <p>
	 * The algorithm considers all the edges associated with all the elements in the
	 * material to determine how close the x-ray energy is to an edge and then
	 * compute a
	 *
	 *
	 * @param uvs   An UncertainValues containing more than one MaterialMAC tag
	 * @param above Range above edge to decorrelate (nominally 400 eV)
	 * @param below Range below edge to decorrelate (nominally 5 eV)
	 * @return
	 */
	static public UncertainValues<EPMALabel> decorrelate(
			//
			final UncertainValues<EPMALabel> uvs, //
			final double above, final double below
	) {
		return decorrelate(uvs, new DefaultDecorrelationFunction(above, below));
	}

	/**
	 *
	 * <p>
	 * This function checks the MatrialMAC values to reduce the correlation between
	 * materials that 1) share a common element; and 2) the edge of which is close
	 * in energy to the X-ray energy.
	 * </p>
	 *
	 * <p>
	 * The algorithm considers all the edges associated with all the elements in the
	 * material to determine how close the x-ray energy is to an edge and then
	 * compute a
	 *
	 *
	 * @param uvs   An UncertainValues containing more than one MaterialMAC tag
	 * @param above Range above edge to decorrelate (nominally 400 eV)
	 * @param below Range below edge to decorrelate (nominally 5 eV)
	 * @return
	 */
	static public UncertainValues<EPMALabel> decorrelate2(
			//
			final UncertainValues<EPMALabel> uvs, //
			final double above, final double below
	) {
		return decorrelate(uvs, new TotalDecorrelationFunction(above, below));
	}

	private static List<EPMALabel> inputTags(
			final List<Material> stds, final List<Material> unk, final XRay xray
	) {
		final List<EPMALabel> res = new ArrayList<>();
		final Set<Element> elms = new HashSet<>();
		for (final Material mat : stds) {
			elms.addAll(mat.getElementSet());
			res.addAll(MaterialLabel.massFractionTags(mat));
		}
		// Don't add MassFraction for unknown
		for (final Material mat : unk)
			elms.addAll(mat.getElementSet());
		for (final Element elm : elms)
			res.add(new ElementalMAC.ElementMAC(elm, xray));
		return res;
	}

	private static List<MaterialMAC> outputTags(
			final List<Material> stds, //
			final List<Material> unks, //
			final XRay xray
	) {
		final List<MaterialMAC> res = new ArrayList<>();
		for (final Material mat : stds)
			res.add(new MaterialMAC(mat, xray));
		for (final Material mat : unks)
			res.add(new MaterialMAC(mat, xray));
		return res;
	}

	/**
	 * Constructs a MassAbsorptionCoefficient
	 *
	 * @param other Materials other than the unknowns
	 * @param unks  Material associated with unknowns
	 * @param xrays A list of x-ray energies to compute
	 * @throws ArgumentException
	 */
	public MaterialMACFunction(
			final List<Material> other, //
			final List<Material> unks, //
			final XRay xray
	) throws ArgumentException {
		super(inputTags(other, unks, xray), outputTags(other, unks, xray));
	}

	public UncertainValuesBase<EPMALabel> buildInputs(
			final List<Composition> comps, //
			final XRay xray //
	) throws ArgumentException {
		final ElementalMAC mac = new ElementalMAC();
		final Map<EPMALabel, Number> men = new HashMap<>();
		for (final EPMALabel lbl : getInputLabels()) {
			if (lbl instanceof ElementalMAC.ElementMAC) {
				final ElementalMAC.ElementMAC eLbl = (ElementalMAC.ElementMAC) lbl;
				final Element elm = eLbl.getElement();
				final XRay xr = eLbl.getXRay();
				men.put(new ElementalMAC.ElementMAC(elm, xr), mac.compute(elm, xr));
			}
		}
		final List<UncertainValuesBase<?>> luvb = new ArrayList<>();
		for (final Composition comp : comps)
			luvb.add(comp.toMassFraction());
		luvb.add(new UncertainValues<EPMALabel>(men));
		return UncertainValues.<EPMALabel>extract(getInputLabels(), luvb);
	}

	@Override
	public RealVector computeValue(
			final double[] point
	) {
		final RealVector vals = buildResult();
		for(MaterialMAC matMacTag : getOutputLabels()) {
			final Material mat = matMacTag.getMaterial();
			double matMac = 0.0;
			for (final Element elm : mat.getElementSet()) {
				final ElementalMAC.ElementMAC elmMacTag = new ElementalMAC.ElementMAC(elm, matMacTag.getXRay());
				final MaterialLabel.MassFraction mfTag = MaterialLabel.buildMassFractionTag(mat, elm);
				matMac += getArg(elmMacTag, point) * getArg(mfTag, point);
			}
			setResult(matMacTag, vals, matMac);
		}
		return vals;
	}

	@Override
	public String toString() {
		return "MaterialMAC" + getOutputLabels();
	}

	/**
	 * Computes the mass absorption coefficient for a set of materials relative to a
	 * specified x-ray. The uncertainties in the mass fractions and in the elemental
	 * MACs are both included.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(
			final RealVector point
	) {
		final RealVector vals = buildResult();
		final RealMatrix jac = buildJacobian();
		for(MaterialMAC matMacTag : getOutputLabels()) {
			final Material mat = matMacTag.getMaterial();
			double matMac = 0.0;
			for (final Element elm : mat.getElementSet()) {
				final ElementalMAC.ElementMAC elmMacTag = new ElementalMAC.ElementMAC(elm, matMacTag.getXRay());
				final MaterialLabel.MassFraction mfTag = MaterialLabel.buildMassFractionTag(mat, elm);
				final double elmMacVal = getArg(elmMacTag, point);
				final double mfVal = getArg(mfTag, point);
				matMac += elmMacVal * mfVal;
				setJacobian(elmMacTag, matMacTag, jac, mfVal);
				setJacobian(mfTag, matMacTag, jac, elmMacVal);
			}
			setResult(matMacTag, vals, matMac);
		}
		return Pair.create(vals, jac);
	}

	public static ExplicitMeasurementModel<EPMALabel, MaterialMAC> buildMaterialMACFunctions(
			final Material estUnknown, //
			final Set<KRatioLabel> kratios
	) throws ArgumentException {
		final Map<Material, CharacteristicXRaySet> comps = new HashMap<>();
		final List<ExplicitMeasurementModel<? extends EPMALabel, ? extends MaterialMAC>> funcs = new ArrayList<>();
		final FastIndex<Material> other = new FastIndex<>();
		final FastIndex<Material> unks = new FastIndex<>();
		{
			final Set<CharacteristicXRay> allCxr = new HashSet<>();
			comps.put(estUnknown, new CharacteristicXRaySet());
			for (final KRatioLabel krl : kratios) {
				other.addIfMissing(krl.getStandard().getMaterial());
				unks.addIfMissing(krl.getUnknown().getMaterial());
				comps.get(estUnknown).addAll(krl.getXRaySet());
				if (krl.getUnknown().hasCoating()) {
					final Material mat = krl.getUnknown().getCoating().getMaterial();
					other.addIfMissing(mat);
					if (!comps.containsKey(mat))
						comps.put(mat, new CharacteristicXRaySet());
					comps.get(mat).addAll(krl.getXRaySet());
				}
				{
					final Material mat = krl.getStandard().getMaterial();
					if (!comps.containsKey(mat))
						comps.put(mat, new CharacteristicXRaySet());
					comps.get(mat).addAll(krl.getXRaySet());
				}
				if (krl.getStandard().hasCoating()) {
					final Material mat = krl.getStandard().getCoating().getMaterial();
					other.addIfMissing(mat);
					if (!comps.containsKey(mat))
						comps.put(mat, new CharacteristicXRaySet());
					comps.get(mat).addAll(krl.getXRaySet());
				}
				allCxr.addAll(krl.getXRaySet().getSetOfCharacteristicXRay());
			}
			for (final CharacteristicXRay cxr : allCxr) {
				final Set<Material> mats = new HashSet<>();
				for (final Map.Entry<Material, CharacteristicXRaySet> me : comps.entrySet())
					if (me.getValue().contains(cxr))
						mats.add(me.getKey());
				funcs.add(new MaterialMACFunction(other, unks, cxr));
			}
		}
		return ParallelMeasurementModelBuilder.join("MaterialMACs", funcs);
	}

}