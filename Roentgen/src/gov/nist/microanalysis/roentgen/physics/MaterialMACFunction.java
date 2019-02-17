package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.NullableRealMatrix;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;
import gov.nist.microanalysis.roentgen.math.uncertainty.ILabeledMultivariateFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;

/**
 * Compute the material mass absorption coefficient given the elemental mass
 * absorption coefficients.
 *
 * @author Nicholas W. M. Ritchie
 *
 */
public class MaterialMACFunction //
		extends LabeledMultivariateJacobianFunction //
		implements ILabeledMultivariateFunction {

	public static class MaterialMAC //
			extends BaseLabel<Material, XRay, Object> {

		public MaterialMAC(final Material mf, final XRay xr) {
			super("[&mu;/&rho;]", mf, xr);
		}

		public Material getMaterial() {
			return getObject1();
		}

		public XRay getXRay() {
			return getObject2();
		}

	}

	private static List<Object> inputTags(final List<? extends Material> comps, final XRay xray) {
		final List<Object> res = new ArrayList<>();
		final Set<Element> elms = new HashSet<>();
		for (final Material comp : comps) {
			elms.addAll(comp.getElementSet());
			res.addAll(MaterialLabel.massFractionTags(comp));
		}
		for (final Element elm : elms)
			res.add(new ElementalMAC.ElementMAC(elm, xray));
		return res;
	}

	/**
	 * Helper for implementing {@link ElementalMAC}.compute(...) and
	 * {@link ElementalMAC}.computeMC(...)
	 *
	 * @param comps List&lt;{@link Composition}&gt;
	 * @param xray  {@link XRay}
	 * @return Pair&lt;{@link UncertainValues},
	 *         {@link ComputeMassAbsorptionCoefficients}&gt;
	 */
	final static public Pair<UncertainValues, MaterialMACFunction> buildCompute(//
			final List<Composition> comps, //
			final XRay xray //
	) {
		final List<Material> mfs = new ArrayList<>();
		for (final Composition comp : comps)
			mfs.add(comp.getMaterial());
		final MaterialMACFunction cmac = new MaterialMACFunction(mfs, xray);
		// Builds an input uncertainty matrix directly
		final List<? extends Object> inp = cmac.getInputLabels();
		final UncertainValues uvs = new UncertainValues(inp);
		final Set<Element> elms = new HashSet<>();
		// Initialize the elemental fractions for each material
		for (final Composition mf : comps) {
			UncertainValues.copy(mf, uvs);
			elms.addAll(mf.getElementSet());
		}
		// Initialize the MACS for all elements relative to xr
		final ElementalMAC mac = new ElementalMAC();
		for (final Element elm : elms)
			uvs.set(new ElementalMAC.ElementMAC(elm, xray), mac.compute(elm, xray));
		return Pair.create(uvs, cmac);
	}

	interface DecorrelationFunction {

		public double compute(AtomicShell edgeEnergy, XRay xray);

	}

	public static class DefaultDecorrelationFunction implements DecorrelationFunction {

		private final double mAbove;
		private final double mBelow;

		public DefaultDecorrelationFunction(final double above, final double below) {
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
		public double compute(final AtomicShell shell, final XRay xray) {
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

		public TotalDecorrelationFunction(final double above, final double below) {
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
		public double compute(final AtomicShell shell, final XRay xray) {
			double res = 1.0;
			final double delta = xray.getEnergy() - shell.getEdgeEnergy();
			if ((delta > -2.0 * mBelow) && (delta < mAbove))
				res = 0.0;
			return Math.max(0.0, Math.min(1.0, res * res));
		}
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
	static public UncertainValues decorrelate( //
			final UncertainValues uvs, //
			final double above, final double below) {
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
	static public UncertainValues decorrelate2( //
			final UncertainValues uvs, //
			final double above, final double below) {
		return decorrelate(uvs, new TotalDecorrelationFunction(above, below));
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
	public static UncertainValues decorrelate( //
			final UncertainValues uvs, final DecorrelationFunction decor) {
		final UncertainValues res = uvs.copy();
		final List<? extends Object> labels = res.getLabels();
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

	final static public UncertainValues compute(final List<Composition> materials, final XRay xray) //
			throws ArgumentException {
		final Pair<UncertainValues, MaterialMACFunction> tmp = buildCompute(materials, xray);
		return UncertainValues.propagate(tmp.getValue(), tmp.getFirst());
	}

	private static List<Object> outputTags(final List<? extends Material> materials, final XRay xray) {
		final List<Object> res = new ArrayList<>();
		for (final Material mf : materials)
			res.add(new MaterialMAC(mf, xray));
		return res;
	}

	/**
	 * Constructs a MassAbsorptionCoefficient
	 *
	 * @param materials A list of materials to compute
	 * @param xrays     A list of x-ray energies to compute
	 */
	public MaterialMACFunction(final List<? extends Material> materials, final XRay xray) {
		super(inputTags(materials, xray), outputTags(materials, xray));
	}

	/**
	 * Computes the mass absorption coefficient for a set of materials relative to a
	 * specified x-ray. The uncertainties in the mass fractions and in the elemental
	 * MACs are both included.
	 *
	 * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
	 */
	@Override
	public Pair<RealVector, RealMatrix> value(final RealVector point) {
		final RealVector res = new ArrayRealVector(getOutputDimension());
		final RealMatrix cov = NullableRealMatrix.build(getOutputDimension(), getInputDimension());
		final List<? extends Object> macTags = getOutputLabels();
		final List<? extends Object> inp = getInputLabels();
		for (int matMacIdx = 0; matMacIdx < macTags.size(); ++matMacIdx) {
			final MaterialMAC matMacTag = (MaterialMAC) macTags.get(matMacIdx);
			final Material comp = matMacTag.getMaterial();
			final XRay xr = matMacTag.getXRay();
			double matMac = 0.0;
			for (final Element elm : comp.getElementSet()) {
				final ElementalMAC.ElementMAC elmMacTag = new ElementalMAC.ElementMAC(elm, xr);
				final MaterialLabel.MassFraction mfTag = MaterialLabel.buildMassFractionTag(comp, elm);
				final int macIdx = inp.indexOf(elmMacTag);
				assert macIdx >= 0 : elmMacTag;
				final int mfIdx = inp.indexOf(mfTag);
				assert mfIdx >= 0 : mfTag + " is missing...";
				final double macVal = point.getEntry(macIdx);
				final double mfVal = point.getEntry(mfIdx);
				matMac += macVal * mfVal;
				cov.setEntry(matMacIdx, macIdx, mfVal);
				cov.setEntry(matMacIdx, mfIdx, macVal);
			}
			res.setEntry(matMacIdx, matMac);
		}
		return Pair.create(res, cov);
	}

	@Override
	public RealVector optimized(final RealVector point) {
		final RealVector res = new ArrayRealVector(getOutputDimension());
		final List<? extends Object> macTags = getOutputLabels();
		final List<? extends Object> inp = getInputLabels();
		for (int i = 0; i < macTags.size(); ++i) {
			final MaterialMAC macTag = (MaterialMAC) macTags.get(i);
			final Material mf = macTag.getMaterial();
			final XRay xr = macTag.getXRay();
			double tmp = 0.0;
			for (final Element elm : mf.getElementSet()) {
				final int p = inp.indexOf(new ElementalMAC.ElementMAC(elm, xr));
				final int q = inp.indexOf(MaterialLabel.buildMassFractionTag(mf, elm));
				assert (p >= 0) && (q >= 0);
				tmp += point.getEntry(p) * point.getEntry(q);
			}
			res.setEntry(i, tmp);
		}
		return res;
	}
}