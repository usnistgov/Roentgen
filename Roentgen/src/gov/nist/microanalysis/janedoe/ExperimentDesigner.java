package gov.nist.microanalysis.janedoe;

import com.duckandcover.html.IToHTML;

public class ExperimentDesigner
   implements
   IToHTML {

   public enum DesignMode {
      /**
       * Determine the best value for a single factor
       */
      Comparative,
      /**
       * Determine most important factors
       */
      Sensitivity,
      /**
       * Maximize the output
       */
      Optimization,
      /**
       * Model the output
       */
      Regression,
      /**
       * Determine a consensus value for Y
       */
      Consensus,
      /**
       * Pass-fail test
       */
      PassFail
   };

   private final DesignMode mMode;

   public ExperimentDesigner(final DesignMode mode) {
      mMode = mode;
   }

   /**
    * Assign the number of levels associated with the continuous factor labeled
    * with tag
    *
    * @param tag
    * @return ExperimentDesigner
    */
   public ExperimentDesigner assignLevels(final Object tag, final double[] levels) {

      return this;
   }

   /**
    * Assign the indices of the discrete levels.
    *
    * @param tag
    * @param indices
    * @return ExperimentDesigner
    */
   public ExperimentDesigner assignLevels(final Object tag, final int[] indices) {

      return this;
   }

   /**
    * Specify the maximum number of measurements in budget.
    *
    * @param n
    * @return
    */
   public ExperimentDesigner assignMaxEvaluations(final int n) {

      return this;
   }

   public DesignMode getMode() {
      return mMode;
   }

   @Override
   public String toHTML(final Mode mode) {
      // TODO Auto-generated method stub
      return null;
   }

}
