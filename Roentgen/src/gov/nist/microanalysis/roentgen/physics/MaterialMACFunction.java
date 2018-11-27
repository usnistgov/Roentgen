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
import gov.nist.microanalysis.roentgen.physics.composition.Composition.MassFractionTag;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;

/**
 * Compute the material mass absorption coefficient given the elemental mass
 * absorption coefficients.
 * 
 * @author nicholas
 *
 */
public class MaterialMACFunction extends LabeledMultivariateJacobianFunction implements ILabeledMultivariateFunction {

	public static class MaterialMAC extends BaseLabel<Composition, XRay, Object> {

		public MaterialMAC(final Composition mf, final XRay xr) {
			super("[&mu;/&rho;]", mf, xr);
		}

		public Composition getComposition() {
			return getObject1();
		}

		public XRay getXRay() {
			return getObject2();
		}
	}

	private static List<Object> inputTags(final List<Composition> comps, final XRay xray) {
		final List<Object> res = new ArrayList<>();
		final Set<Element> elms = new HashSet<>();
		for (final Composition comp : comps) {
			assert comp.getNativeRepresentation() == Representation.MassFraction;
			elms.addAll(comp.getElementSet());
		}
		for (final Element elm : elms)
			res.add(new ElementalMAC.ElementMAC(elm, xray));
		// Mass fraction of each element in each material
		for (final Composition comp : comps)
			res.addAll(comp.getLabels());
		return res;
	}

	/**
	 * Helper for implementing {@link ElementalMAC}.compute(...) and
	 * {@link ElementalMAC}.computeMC(...)
	 *
	 * @param comps List&lt;{@link IComposition}&gt;
	 * @param xray  {@link XRay}
	 * @return Pair&lt;{@link UncertainValues},
	 *         {@link ComputeMassAbsorptionCoefficients}&gt;
	 */
	final static public Pair<UncertainValues, MaterialMACFunction> buildCompute(final List<Composition> materials,
			final XRay xray) {
		final List<Composition> mfs = new ArrayList<>();
		for (Composition comp : materials)
			mfs.add(comp.asMassFraction());
		final MaterialMACFunction cmac = new MaterialMACFunction(mfs, xray);
		// Builds an input uncertainty matrix directly
		final List<? extends Object> inp = cmac.getInputLabels();
		final UncertainValues uvs = new UncertainValues(inp);
		final Set<Element> elms = new HashSet<>();
		// Initialize the elemental fractions for each material
		for (final Composition mf : mfs) {
			UncertainValues.copy(mf, uvs);
			elms.addAll(mf.getElementSet());
		}
		// Initialize the MACS for all elements relative to xr
		final ElementalMAC mac = new ElementalMAC();
		for (final Element elm : elms)
			uvs.set(new ElementalMAC.ElementMAC(elm, xray), mac.compute(elm, xray));
		return Pair.create(uvs, cmac);
	}

	/**
	 * 
	 * This function checks the MatrialMAC values to reduce the correlation between
	 * materials that 1) share a common element; and 2) the edge of which is close
	 * in energy to the X-ray energy.
	 * 
	 * 
	 * @param uvs   An UncertainValues containing more than one MaterialMAC tag
	 * @param xray  X-ray line
	 * @param above Range above edge to decorrelate (nominally 400 eV)
	 * @param below Range below edge to decorrelate (nominally 5 eV)
	 * @return
	 */
	static public UncertainValues decorrelate( //
			UncertainValues uvs, //
			final XRay xray, //
			final double above, //
			final double below) {
		final List<? extends Object> labels = uvs.getLabels();
		final double xrayE = xray.getEnergy();
		for (Object label : labels) {
			double mult = 1.0;
			if (label instanceof MaterialMAC) {
				final MaterialMAC mm1 = (MaterialMAC) label;
				for (Element elm : mm1.getComposition().getElementSet()) {
					for (final AtomicShell sh : AtomicShell.forElement(elm)) {
						final double ee = sh.getEdgeEnergy();
						final double delta = xrayE - ee;
						if (delta > 0.0) {
							// xrayE above edge
							final double tmp = delta > Math.abs(above) ? 1.0 : 0.1 + 0.9 * delta / Math.abs(above);
							if (tmp < mult)
								mult = tmp;
						} else {
							// xrayE below edge
							final double tmp = delta < -3.0 * Math.abs(below) ? 1.0
									: 1.0 - 0.9 * Math.exp(-0.5 * Math.pow(delta / below, 2));
							if (tmp < mult)
								mult = tmp;
						}
					}
					if (mult < 1.0) {
						assert (mult >= 0.0) && (mult < 1.0);
						boolean past = false;
						for (Object label2 : labels) {
							if (label2 == label) {
								past = true;
								continue;
							}
							if (past && (label2 instanceof MaterialMAC)) {
								MaterialMAC mm2 = (MaterialMAC) label2;
								if (mm2.getComposition().getElementSet().contains(elm)) {
									final double cov = mult * mult * uvs.getCovariance(label, label2);
									uvs.setCovariance(label2, label, cov);
								}
							}
						}
					}
				}

			}
		}
		return uvs;
	}

	final static public UncertainValues compute(List<Composition> materials, final XRay xray) //
			throws ArgumentException {
		Pair<UncertainValues, MaterialMACFunction> tmp = buildCompute(materials, xray);
		return UncertainValues.propagate(tmp.getValue(), tmp.getFirst());
	}

	private static List<Object> outputTags(final List<Composition> materials, final XRay xray) {
		final List<Object> res = new ArrayList<>();
		for (final Composition mf : materials)
			res.add(new MaterialMAC(mf, xray));
		return res;
	}

	/**
	 * Constructs a MassAbsorptionCoefficient
	 *
	 * @param materials A list of materials to compute
	 * @param xrays     A list of x-ray energies to compute
	 */
	public MaterialMACFunction(final List<Composition> materials, final XRay xray) {
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
		final RealMatrix cov = new NullableRealMatrix(getOutputDimension(), getInputDimension());
		final List<? extends Object> macTags = getOutputLabels();
		final List<? extends Object> inp = getInputLabels();
		for (int matMacIdx = 0; matMacIdx < macTags.size(); ++matMacIdx) {
			final MaterialMAC matMacTag = (MaterialMAC) macTags.get(matMacIdx);
			final Composition comp = matMacTag.getComposition();
			final XRay xr = matMacTag.getXRay();
			double matMac = 0.0;
			for (final Element elm : comp.getElementSet()) {
				final ElementalMAC.ElementMAC elmMacTag = new ElementalMAC.ElementMAC(elm, xr);
				final MassFractionTag mfTag = Composition.buildMassFractionTag(comp, elm);
				final int macIdx = inp.indexOf(elmMacTag);
				assert macIdx >= 0 : elmMacTag;
				final int mfIdx = inp.indexOf(mfTag);
				assert mfIdx >= 0 : mfTag;
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
	public RealVector optimized(RealVector point) {
		final RealVector res = new ArrayRealVector(getOutputDimension());
		final List<? extends Object> macTags = getOutputLabels();
		final List<? extends Object> inp = getInputLabels();
		for (int i = 0; i < macTags.size(); ++i) {
			final MaterialMAC macTag = (MaterialMAC) macTags.get(i);
			final Composition mf = macTag.getComposition();
			final XRay xr = macTag.getXRay();
			double tmp = 0.0;
			for (final Element elm : mf.getElementSet()) {
				final int p = inp.indexOf(new ElementalMAC.ElementMAC(elm, xr));
				final int q = inp.indexOf(Composition.buildMassFractionTag(mf, elm));
				assert (p >= 0) && (q >= 0);
				tmp += point.getEntry(p) * point.getEntry(q);
			}
			res.setEntry(i, tmp);
		}
		return res;
	}
}