package gov.nist.microanalysis.roentgen.math.uncertainty;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;

public class ImplicitMeasurementModel {

   private final NamedMultivariateJacobianFunction mCy;
   private final NamedMultivariateJacobianFunction mCx;

   private static final boolean NAIVE = true;

   private static List<? extends Object> combine(List<? extends Object> a, List<? extends Object> b) {
      ArrayList<Object> res = new ArrayList<Object>();
      res.addAll(a);
      for(Object bo : b)
         if(!res.contains(bo))
            res.add(bo);
      return res;
   }

   public ImplicitMeasurementModel(final NamedMultivariateJacobianFunction cy, NamedMultivariateJacobianFunction cx) {
      mCy = cy;
      mCx = cx;
   }

   public List<? extends Object> getInputTags() {
      return combine(mCy.getInputTags(), mCx.getInputTags());
   }

   public List<? extends Object> getOutputTags() {
      return mCy.getOutputTags();
   }

   /**
    * Solve Cy Uy CyT = Cx Ux CxT for Uy
    * 
    * @param ux
    * @param y
    * @return UncertainValues
    */
   public UncertainValues computeUy(UncertainValues ux, RealVector y) {
      if(NAIVE) {

      }
      return null;
   }

}
