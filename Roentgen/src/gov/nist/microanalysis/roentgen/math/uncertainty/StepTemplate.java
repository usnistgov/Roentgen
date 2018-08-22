package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * This class is intended to serve as a template for implementing
 * {@link NamedMultivariateJacobian}. It shows an effective way to implement the
 * input and output tags and the <code>value(...)</code> function
 *
 * @author Nicholas
 */
public class StepTemplate
   extends
   NamedMultivariateJacobianFunction {

   private final Composition mComposition;
   private final AtomicShell mShell;

   public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
      final List<Object> res = new ArrayList<>();
      res.add(new BaseTag<Composition, Object, Object>("P", comp));
      return res;
   }

   public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
      final List<Object> res = new ArrayList<>();
      res.add(new BaseTag<Composition, Object, Object>("F", comp));
      res.add(new BaseTag<Composition, AtomicShell, Object>("G", comp, shell));
      return res;
   }

   public StepTemplate(final Composition comp, final AtomicShell shell) {
      super(buildInputs(comp, shell), buildOutputs(comp, shell));
      mComposition = comp;
      mShell = shell;
   }

   @Override
   public Pair<RealVector, RealMatrix> value(final RealVector point) {
      final int iF = inputIndex(new BaseTag<Composition, Object, Object>("F", mComposition));
      final int iG = inputIndex(new BaseTag<Composition, AtomicShell, Object>("G", mComposition, mShell));

      final double F = point.getEntry(iF);
      final double G = point.getEntry(iG);

      final int oP = outputIndex(new BaseTag<Composition, Object, Object>("P", mComposition));

      final RealVector rv = new ArrayRealVector(getOutputDimension());
      final RealMatrix rm = new Array2DRowRealMatrix(getOutputDimension(), getInputDimension());

      final double P = F * G;
      // Calculate and assign values
      rv.setEntry(oP, P);
      // Calculate and assign Jacobian partials
      rm.setEntry(oP, iF, G);
      rm.setEntry(oP, iG, F);

      return Pair.create(rv, rm);
   }
}
