package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.MassAbsorptionCoefficient.MaterialMAC;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;

public class ComputeMACs
   extends
   NamedMultivariateJacobianFunction {

   private static List<Object> inputTags(final List<Composition> comps, final XRay xray) {
      final List<Object> res = new ArrayList<>();
      final Set<Element> elms = new HashSet<>();
      for(final Composition comp : comps) {
         assert comp.getNativeRepresentation() == Representation.MassFraction;
         elms.addAll(comp.getElementSet());
      }
      for(final Element elm : elms)
         res.add(new MassAbsorptionCoefficient.ElementMAC(elm, xray));
      // Mass fraction of each element in each material
      for(final Composition comp : comps)
         res.addAll(comp.getTags());
      return res;
   }

   /**
    * Helper for implementing {@link MassAbsorptionCoefficient}.compute(...) and
    * {@link MassAbsorptionCoefficient}.computeMC(...)
    *
    * @param comps List&lt;{@link IComposition}&gt;
    * @param xray {@link XRay}
    * @return Pair&lt;{@link UncertainValues},
    *         {@link ComputeMassAbsorptionCoefficients}&gt;
    */
   final static public Pair<UncertainValues, ComputeMACs> buildCompute(final List<Composition> materials, final XRay xray) {
      final ComputeMACs cmac = new ComputeMACs(materials, xray);
      // Builds an input uncertainty matrix directly
      final List<? extends Object> inp = cmac.getInputTags();
      final UncertainValues uvs = new UncertainValues(inp);
      final Set<Element> elms = new HashSet<>();
      // Initialize the elemental fractions for each material
      for(final Composition mf : materials) {
         UncertainValues.copy(mf, uvs);
         elms.addAll(mf.getElementSet());
      }
      // Initialize the MACS for all elements relative to xr
      final MassAbsorptionCoefficient mac = MassAbsorptionCoefficient.instance();
      for(final Element elm : elms)
         uvs.set(new MassAbsorptionCoefficient.ElementMAC(elm, xray), mac.compute(elm, xray));
      return Pair.create(uvs, cmac);
   }

   private static List<Object> outputTags(final List<Composition> materials, final XRay xray) {
      final List<Object> res = new ArrayList<>();
      for(final Composition mf : materials)
         res.add(new MaterialMAC(mf, xray));
      return res;
   }

   public static UncertainValues convert(final List<Composition> materials, final XRay xray) {
      final List<Composition> mfs = new ArrayList<>();
      for(final Composition c : materials)
         mfs.add(c.asMassFraction());
      final Pair<UncertainValues, ComputeMACs> pr = buildCompute(mfs, xray);
      return UncertainValues.propagate(pr.getSecond(), pr.getFirst());
   }

   /**
    * Constructs a MassAbsorptionCoefficient
    *
    * @param materials A list of materials to compute
    * @param xrays A list of x-ray energies to compute
    */
   public ComputeMACs(final List<Composition> materials, final XRay xray) {
      super(inputTags(materials, xray), outputTags(materials, xray));
   }

   /**
    * Computes the mass absorption coefficient for a set of materials relative
    * to a specified x-ray. The uncertainties in the mass fractions and in the
    * elemental MACs are both included.
    *
    * @see org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction#value(org.apache.commons.math3.linear.RealVector)
    */
   @Override
   public Pair<RealVector, RealMatrix> value(final RealVector point) {
      final RealVector res = new ArrayRealVector(getOutputDimension());
      final RealMatrix cov = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());
      final List<? extends Object> macTags = getOutputTags();
      final List<? extends Object> inp = getInputTags();
      for(int i = 0; i < macTags.size(); ++i) {
         final MaterialMAC macTag = (MaterialMAC) macTags.get(i);
         final Composition mf = macTag.getComposition();
         final XRay xr = macTag.getXRay();
         double tmp = 0.0;
         for(final Element elm : mf.getElementSet()) {
            final int p = inp.indexOf(new MassAbsorptionCoefficient.ElementMAC(elm, xr));
            final int q = inp.indexOf(Composition.buildMassFractionTag(mf, elm));
            assert (p >= 0) && (q >= 0);
            tmp += point.getEntry(p) * point.getEntry(q);
            cov.setEntry(i, p, point.getEntry(q));
            cov.setEntry(i, q, point.getEntry(p));
         }
         res.setEntry(i, tmp);
      }
      return Pair.create(res, cov);
   }
}