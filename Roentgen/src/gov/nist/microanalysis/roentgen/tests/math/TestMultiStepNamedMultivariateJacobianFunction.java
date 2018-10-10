package gov.nist.microanalysis.roentgen.tests.math;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.SerialNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobian;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import junit.framework.TestCase;

/**
 * @author Nicholas
 */
public class TestMultiStepNamedMultivariateJacobianFunction
   extends
   TestCase {

   private final List<? extends Object> mInputs = Arrays.asList(new String[] {
      "X0",
      "X1",
      "X2"
   });
   private final List<? extends Object> mOut1 = Arrays.asList(new String[] {
      "F1_0",
      "F1_1",
      "F1_2"
   });
   private final List<? extends Object> mInputs2 = Arrays.asList(new String[] {
      "X0",
      "X1",
      "X2",
      "F1_0",
      "F1_1",
      "F1_2"
   });
   private final List<? extends Object> mOut2 = Arrays.asList(new String[] {
      "F2_0",
      "F2_1"
   });

   private final LabeledMultivariateJacobianFunction mStep1 = new LabeledMultivariateJacobianFunction(mInputs, mOut1) {

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final RealVector fv = new ArrayRealVector(getOutputDimension());
         final RealMatrix jac = MatrixUtils.createRealMatrix(getOutputDimension(), getOutputDimension());
         final int ix0 = inputIndex("X0");
         final int ix1 = inputIndex("X1");
         final int ix2 = inputIndex("X2");
         final int iF1_0 = outputIndex("F1_0");
         fv.setEntry(iF1_0, 13.0 + 5.0 * Math.pow(point.getEntry(ix0), 2.7) + 7.0 * Math.pow(point.getEntry(ix1), 1.3)
               + 11.0 * Math.pow(point.getEntry(ix2), -3.2));
         jac.setEntry(iF1_0, ix0, 5.0 * 2.7 * Math.pow(point.getEntry(ix0), 1.7));
         jac.setEntry(iF1_0, ix1, 7.0 * 1.3 * Math.pow(point.getEntry(ix1), 0.3));
         jac.setEntry(iF1_0, ix2, 11.0 * (-3.2) * Math.pow(point.getEntry(ix2), -4.2));

         final int iF1_1 = outputIndex("F1_1");
         fv.setEntry(iF1_1, 2.0 + 4.0 * Math.pow(point.getEntry(ix0), 2.0) + 6.0 * Math.pow(point.getEntry(ix1), 2.0)
               + 8.0 * Math.pow(point.getEntry(ix2), 2.0));
         jac.setEntry(iF1_1, ix0, 2.0 * 4.0 * point.getEntry(ix0));
         jac.setEntry(iF1_1, ix1, 2.0 * 6.0 * point.getEntry(ix1));
         jac.setEntry(iF1_1, ix2, 2.0 * 8.0 * point.getEntry(ix2));

         final int iF1_2 = outputIndex("F1_2");
         fv.setEntry(iF1_2, 1.0 + 0.5 * Math.pow(point.getEntry(ix0), 2.0) + 0.25 * Math.pow(point.getEntry(ix1), 4.0)
               + 0.125 * Math.pow(point.getEntry(ix2), 8.0));
         jac.setEntry(iF1_2, ix0, point.getEntry(ix0));
         jac.setEntry(iF1_2, ix1, Math.pow(point.getEntry(ix1), 3.0));
         jac.setEntry(iF1_2, ix2, Math.pow(point.getEntry(ix2), 7.0));
         return Pair.create(fv, jac);
      }
   };

   private final LabeledMultivariateJacobianFunction mStep2 = new LabeledMultivariateJacobianFunction(mInputs2, mOut2) {

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final double x0 = point.getEntry(0);
         final double x1 = point.getEntry(1);
         final double x2 = point.getEntry(2);
         final double f1_0 = point.getEntry(3);
         final double f1_1 = point.getEntry(4);
         final double f1_2 = point.getEntry(5);
         final RealVector fv = new ArrayRealVector(2);
         final RealMatrix jac = MatrixUtils.createRealMatrix(2, 6);
         fv.setEntry(0, Math.log(7.0 * f1_1 + f1_2) + 3.0 * x1 + 11.0 * x2);
         jac.setEntry(0, 1, 3.0);
         jac.setEntry(0, 2, 11.0);
         jac.setEntry(0, 4, (1.0 / (7.0 * f1_1 + f1_2)) * 7.0);
         jac.setEntry(0, 5, (1.0 / (7.0 * f1_1 + f1_2)) * 1.0);

         fv.setEntry(1, Math.log(f1_2 + 3.0 * f1_0) + 7.0 * x0 + 2.0 * x2);
         jac.setEntry(1, 0, 7.0);
         jac.setEntry(1, 2, 2.0);
         jac.setEntry(1, 3, (1.0 / (f1_2 + 3.0 * f1_0)) * 3.0);
         jac.setEntry(1, 5, (1.0 / (f1_2 + 3.0 * f1_0)) * 1.0);

         return Pair.create(fv, jac);
      }
   };

   private final RealVector mValues0 = new ArrayRealVector(new double[] {
      2.3,
      3.5,
      1.3
   });
   private final RealMatrix mCov0 = new Array2DRowRealMatrix(new double[][] {
      {
         0.1,
         0.03,
         -0.013
      },
      {
         0.03,
         0.7,
         0.002
      },
      {
         -0.013,
         0.002,
         0.3
      }
   });

   private final UncertainValues mInput0 = new UncertainValues(mInputs, mValues0, mCov0);

   private final RealVector mValues1 = new ArrayRealVector(new double[] {
      2.3,
      3.5,
      1.3
   });
   private final RealMatrix mCov1 = new Array2DRowRealMatrix(new double[][] {
      {
         0.01,
         0.003,
         -0.0013
      },
      {
         0.003,
         0.07,
         0.0002
      },
      {
         -0.0013,
         0.0002,
         0.03
      }
   });

   private final UncertainValues mInput1 = new UncertainValues(mInputs, mValues1, mCov1);

   public void test1()
         throws IOException {
      final UncertainValues uv = UncertainValues.propagate(mStep1, mInput0);
      final Report rep = new Report("Step 1");
      rep.addHeader("Inputs");
      rep.add(mInput0);
      rep.addHeader("Jacobian");
      rep.add(new LabeledMultivariateJacobian(mStep1, mValues0));
      rep.addHeader("Outputs 1");
      rep.add(uv);
      rep.inBrowser(Mode.NORMAL);
      // Calculated with Mathematica notebook MultiStepNMJF test1.nb
      final double[] rv = new double[] {
         100.812,
         110.18,
         42.1803
      };
      final double[][] rc = new double[][] {
         {
            533.887,
            483.696,
            455.926
         },
         {
            483.696,
            1438.36,
            1330.66
         },
         {
            455.926,
            1330.66,
            1305.74
         }
      };
      final UncertainValues res = new UncertainValues(mOut1, new ArrayRealVector(rv), MatrixUtils.createRealMatrix(rc));
      assertTrue(uv.equals(res, 0.01));
   }

   public void test2()
         throws IOException, ArgumentException {
      final LabeledMultivariateJacobianFunction[] steps = new LabeledMultivariateJacobianFunction[] {
         mStep1,
         mStep2
      };
      final SerialNamedMultivariateJacobianFunction msnmjf = new SerialNamedMultivariateJacobianFunction("Test1", Arrays.asList(steps));
      final UncertainValues uv = UncertainValues.propagate(msnmjf, mInput1);
      final Report rep = new Report("Step 1 and 2");
      rep.addHeader("Inputs");
      rep.add(mInput1);
      rep.addHeader("Jacobian");
      rep.add(LabeledMultivariateJacobian.compute(msnmjf, mValues0));
      rep.addHeader("Outputs 1 and 2");
      rep.add(uv);
      rep.addHeader("MC Outputs 1 and 2");
      rep.add(UncertainValues.propagateMC(msnmjf, mInput1, 100000));
      rep.inBrowser(Mode.NORMAL);

      // Calculated with Mathematica notebook MultiStepNMJF test1.nb
      // double[] rv = new double[] { 100.812, 110.18, 42.1803, 31.501, 24.5424
      // };
      // double[][] rc = new double[][] { { 533.887, 483.696, 455.926 }, {
      // 483.696,
      // 1438.36, 1330.66 },
      // { 455.926, 1330.66, 1305.74 } };
      // UncertainValues res = new UncertainValues(mOut1, new
      // ArrayRealVector(rv),
      // MatrixUtils.createRealMatrix(rc));
      // assertTrue(uv.equals(res, 0.01));
   }

}
