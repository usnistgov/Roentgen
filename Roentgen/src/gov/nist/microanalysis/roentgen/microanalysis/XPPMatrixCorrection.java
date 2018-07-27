package gov.nist.microanalysis.roentgen.microanalysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;
import gov.nist.microanalysis.roentgen.math.uncertainty.MultiStepNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.NamedMultivariateJacobianFunctionBuilder;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.AtomicShell;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.MassAbsorptionCoefficient;
import gov.nist.microanalysis.roentgen.physics.XRaySet.CharacteristicXRaySet;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>This class implements the XPP matrix correction algorithm as described in the
 * "Green Book" (Heinrich & Newbury 1991) by Pouchou and Pichoir. It also
 * implements all the Jacobians necessary to implement an uncertainty
 * calculation to propagate uncertainties in the input parameters into output
 * parameters. It breaks the calculation into a series of much simpler steps for
 * which the Jacobian can be readily calculated.<p>
 * 
 * <p>Since it is derived from {@link NamedMultivariateJacobianFunction}, it computes
 * not only the value (ie. the k-ratio) but also sensitivity matrix (Jacobian)
 * that maps uncertainty in the input parameters into uncertainty in the output
 * parameters.<p>
 * 
 * <p>The following parameters are assumed to have associated uncertainties.
 * <ul>
 *    <li>The composition of the standards (user provided)</li>
 *    <li>The composition of the unknown (user provided)</li>
 *    <li>The beam energy (user provided)</li>
 *    <li>The take-off angle (user provided)</li>
 *    <li>The mass absorption coefficient (computed using FFAST)</li>
 *    <li>The mean ionization coefficient (computed using the equation in PAP1991)</li>
 * </ul>
 * 
 * @author Nicholas
 */
public class XPPMatrixCorrection
   extends
   MultiStepNamedMultivariateJacobianFunction
   implements
   IToHTML {

   private final Map<Composition, CharacteristicXRaySet> mStandards;
   private final Composition mUnknown;

   public static final String CHI = "&chi;";
   public static final String F_CHI = "F(&chi;)";
   public static final String F_CHI_F = "F(&chi;)/F";
   public static final String EPS = "&epsilon;";
   public static final String PHI0 = "&phi;<sub>0</sub>";
   public static final String ONE_OVER_S = "<sup>1</sup>/<sub>S</sub>";
   public static final String QLA = "<html>Q<sub>l</sub><sup>a</sup>";
   public static final String ZBARB = "Z<sub>barb</sub>";
   public static final String RBAR = "R<sub>bar</sub>";
   public static final String E_0 = "E<sub>0</sub>";

   private static class ElementTag
      extends
      BaseTag {

      private ElementTag(final String name, final Element obj) {
         super(name, obj);
      }
   }

   public static class CompositionTag
      extends
      BaseTag {
      public CompositionTag(final String name, final Composition comp) {
         super(name, "<html>" + comp.toHTML(Mode.TERSE));
      }
   }

   public static class ResultTag
      extends
      BaseTag {

      private ResultTag(final Composition unk, final Composition std, final CharacteristicXRay cxr) {
         super("ZA", "<html>" + unk.toHTML(Mode.TERSE), "<html>" + std.toHTML(Mode.TERSE), cxr);
      }
   }

   private static class CompositionTag2<H>
      extends
      BaseTag {

      private CompositionTag2(final String name, final Composition obj1, final H obj2) {
         super(name, "<html>" + obj1.toHTML(Mode.TERSE), obj2);
      }
   }

   private static final double eVtokeV(final double eV) {
      return 0.001 * eV;
   }

   /**
    * <p>
    * Do some very basic checks on the input and output indices.
    * <ul>
    * <li>Not duplicated..</li>
    * <li>Not missing</li>
    * </ul>
    * </p>
    *
    * @param idxs
    */
   private static void checkIndices(final int... idxs) {
      final int len = idxs.length;
      for(int i = 0; i < len; ++i) {
         assert idxs[i] >= 0;
         for(int j = i + 1; j < len; ++j)
            assert idxs[i] != idxs[j];
      }
   }

   private static class StepMJZBarb // C1
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;

      public StepMJZBarb(final Composition comp) {
         super(buildInputs(comp), buildOutputs(comp));
         mComposition = comp;
      }

      public static List<? extends Object> buildOutputs(final Composition comp) {
         final List<Object> res = new ArrayList<>();
         res.add(new CompositionTag("M", comp));
         res.add(new CompositionTag("J", comp));
         res.add(new CompositionTag(ZBARB, comp));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp) {
         final List<Object> res = new ArrayList<>();
         for(final Element elm : comp.getElementSet()) {
            res.add(Composition.buildMassFractionTag(comp, elm));
            res.add(meanIonizationTag(elm));
         }
         return res;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {

         final List<Element> elms = new ArrayList<>(mComposition.getElementSet());
         final int[] iCi = new int[elms.size()];
         final int[] iJi = new int[elms.size()];
         final double[] Ci = new double[elms.size()];
         final double[] ZoA = new double[elms.size()];
         final double[] Z = new double[elms.size()];
         final double[] Ji = new double[elms.size()];
         for(int i = 0; i < Ci.length; ++i) {
            final Element elm = elms.get(i);
            iCi[i] = inputIndex(Composition.buildMassFractionTag(mComposition, elm));
            iJi[i] = inputIndex(meanIonizationTag(elm));
            Ci[i] = point.getEntry(iCi[i]);
            Z[i] = elm.getAtomicNumber();
            Ji[i] = point.getEntry(iJi[i]);
            ZoA[i] = Z[i] / elm.getAtomicWeight();
         }
         checkIndices(iCi);
         checkIndices(iJi);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oM = outputIndex(new CompositionTag("M", mComposition));
         final int oJ = outputIndex(new CompositionTag("J", mComposition));
         final int oZbarb = outputIndex(new CompositionTag(ZBARB, mComposition));
         checkIndices(oM, oJ, oZbarb);

         // Calculate M & partials
         double M = 0.0;
         for(int i = 0; i < Ci.length; ++i)
            M += Ci[i] * ZoA[i];
         rv.setEntry(oM, M); // c1
         for(int i = 0; i < iCi.length; ++i)
            rm.setEntry(oM, iCi[i], ZoA[i]); // c1

         // Calculate J and partials
         double lnJ = 0.0;
         for(int i = 0; i < elms.size(); ++i)
            lnJ += Math.log(Ji[i]) * Ci[i] * ZoA[i];
         lnJ /= M;
         final double J = Math.exp(lnJ); // keV
         rv.setEntry(oJ, J); // C2
         for(int i = 0; i < iCi.length; ++i) {
            rm.setEntry(oJ, iJi[i], Ci[i] * J * ZoA[i] / (M * Ji[i])); // C3
            double k1 = 0.0, k2 = 0.0;
            for(int j = 0; j < iCi.length; ++j) {
               if(j != i) {
                  k1 += Ci[j] * ZoA[j];
                  k2 -= Ci[j] * ZoA[j] * Math.log(Ji[j]);
               }
            }
            rm.setEntry(oJ, iCi[i], J * ZoA[i] * (k1 * Math.log(Ji[i]) + k2) / (M * M)); // C3
         }
         // Calculate Zbarb and partials
         double Zbt = 0.0;
         for(int i = 0; i < elms.size(); ++i)
            Zbt += Ci[i] * Math.sqrt(Z[i]); // C2
         final double Zbarb = Math.pow(Zbt, 2.0);
         rv.setEntry(oZbarb, Zbarb);
         for(int i = 0; i < elms.size(); ++i)
            rm.setEntry(oZbarb, iCi[i], 2.0 * Math.sqrt(Z[i]) * Zbt); // Ci
         return Pair.create(rv, rm);
      }
   }

   private static class StepQlaE0OneOverS // C1
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(ONE_OVER_S, comp, shell));
         res.add(tag2(QLA, comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp) {
         final List<Object> res = new ArrayList<>();
         res.add(new CompositionTag("M", comp));
         res.add(new CompositionTag("J", comp));
         res.add(beamEnergyTag(comp));
         return res;
      }

      public StepQlaE0OneOverS(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iE0 = inputIndex(beamEnergyTag(mComposition));
         final int iJ = inputIndex(new CompositionTag("J", mComposition));
         final int iM = inputIndex(new CompositionTag("M", mComposition));
         checkIndices(iE0, iJ, iM);

         final double e0 = point.getEntry(iE0);
         final double J = point.getEntry(iJ);
         final double M = point.getEntry(iM);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final double Ea = eVtokeV(mShell.getEdgeEnergy());
         final double u0 = e0 / Ea;
         assert u0 > 1.0;

         double m = 1.0;
         final double Za = mShell.getElement().getAtomicNumber();
         switch(mShell.getFamily().getPrincipleQuantumNumber()) {
            case 1: // K
               m = 0.86 + 0.12 * Math.exp(-1.0 * Math.pow(Za / 5, 2.0));
               break;
            case 2: // L
               m = 0.82;
               break;
            case 3: // M
               m = 0.78;
               break;
         }

         final int oOneOverS = outputIndex(tag2(ONE_OVER_S, mComposition, mShell));
         final int oQlaE0 = outputIndex(tag2(QLA, mComposition, mShell));
         checkIndices(oOneOverS, oQlaE0);

         final double logU0 = Math.log(u0);
         // QlaE0 and partials
         final double QlaE0 = logU0 / (Math.pow(u0, m) * Math.pow(Ea, 2.0));
         rv.setEntry(oQlaE0, QlaE0); // C1
         rm.setEntry(oQlaE0, iE0, Math.pow(u0, -1.0 - m) * (1.0 - m * logU0) / Math.pow(Ea, 3.0)); // C1

         // oneOverS and partials
         final double D1 = 6.6e-6, D2 = 1.12e-5 * (1.35 - 0.45 * J), D3 = 2.2e-6 / J;
         final double P1 = 0.78, P2 = 0.1, P3 = 0.25 * J - 0.5;
         final double T1 = 1.0 + P1 - m, T2 = 1.0 + P2 - m, T3 = 1.0 + P3 - m;
         final double dD2dJ = 1.12e-5 * (-0.45), dD3dJ = -1.0 * 2.2e-6 / (J * J), dP3dJ = 0.25, dT3dJ = 0.25;
         final double v0ou0 = Ea / J; // c1

         final double T1_2 = Math.pow(T1, 2), T2_2 = Math.pow(T2, 2), T3_2 = Math.pow(T3, 2), T3_3 = Math.pow(T3, 3);
         final double v0ou0_P1 = Math.pow(v0ou0, P1), v0ou0_P2 = Math.pow(v0ou0, P2), v0ou0_P3 = Math.pow(v0ou0, P3);
         final double u0_T1 = Math.pow(u0, T1), u0_T2 = Math.pow(u0, T2), u0_T3 = Math.pow(u0, T3);
         final double logV0oU0 = Math.log(v0ou0);
         final double dOoSdJ = (T3
               * (T3_2
                     * (-(D1 * v0ou0_P1 * (-1.0 + P1) * T2_2 * (1 - u0_T1 + u0_T1 * T1 * logU0))
                           - v0ou0_P2 * (-1.0 + P2) * T1_2 * D2 * (1 - u0_T2 + u0_T2 * T2 * logU0)
                           + v0ou0_P2 * J * T1_2 * (1 - u0_T2 + u0_T2 * T2 * logU0) * dD2dJ)
                     - (-1.0 + u0_T3) * v0ou0_P3 * J * T1_2 * T2_2 * dD3dJ
                     + u0_T3 * v0ou0_P3 * J * T1_2 * T2_2 * logU0 * T3 * dD3dJ)
               - v0ou0_P3 * T1_2 * T2_2 * D3 * (-2.0 * (-1.0 + u0_T3) * J * dT3dJ
                     - u0_T3 * logU0 * T3_2 * (1 - P3 + J * logV0oU0 * dP3dJ + J * logU0 * dT3dJ) + T3 * (-1.0 + u0_T3
                           - (-1.0 + u0_T3) * P3 + (-1.0 + u0_T3) * J * logV0oU0 * dP3dJ + 2.0 * u0_T3 * J * logU0 * dT3dJ)))
               / (Ea * M * T1_2 * T2_2 * T3_3);

         final double dOoSdE0 = (J * (D1 * u0_T1 * v0ou0_P1 + u0_T2 * v0ou0_P2 * D2 + u0_T3 * v0ou0_P3 * D3) * logU0)
               / (e0 * Ea * M);

         final double dOoSdM = -((J * ((D1 * v0ou0_P1 * (1 - u0_T1 + u0_T1 * T1 * logU0)) / T1_2
               + (v0ou0_P2 * D2 * (1.0 - u0_T2 + u0_T2 * T2 * logU0)) / T2_2
               + (v0ou0_P3 * D3 * (1.0 - u0_T3 + u0_T3 * logU0 * T3)) / T3_2)) / (Ea * Math.pow(M, 2)));

         final double oneOverS = ((J / (Ea * M)) //
               * ((D1 * v0ou0_P1 * (1.0 - u0_T1 + T1 * u0_T1 * logU0)) / T1_2
                     + (D2 * v0ou0_P2 * (1.0 - u0_T2 + T2 * u0_T2 * logU0)) / T2_2 //
                     + (D3 * v0ou0_P3 * (1.0 - u0_T3 + T3 * u0_T3 * logU0)) / T3_2));
         rv.setEntry(oOneOverS, oneOverS); // C2
         rm.setEntry(oOneOverS, iM, dOoSdM); // C2
         rm.setEntry(oOneOverS, iJ, dOoSdJ); // C2
         rm.setEntry(oOneOverS, iE0, dOoSdE0); // C2

         return Pair.create(rv, rm);
      }
   }

   private static class StepRphi0
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2("R", comp, shell));
         res.add(tag2(PHI0, comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp) {
         final List<Object> res = new ArrayList<>();
         res.add(new CompositionTag(ZBARB, comp));
         res.add(new CompositionTag(E_0, comp));
         return res;
      }

      public StepRphi0(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iZbarb = inputIndex(new CompositionTag(ZBARB, mComposition));
         final int iE0 = inputIndex(new CompositionTag(E_0, mComposition));
         checkIndices(iZbarb, iE0);

         final double Ea = eVtokeV(mShell.getEdgeEnergy());
         final double e0 = point.getEntry(iE0);
         final double u0 = e0 / Ea;
         final double Zbarb = point.getEntry(iZbarb);

         final double etabar = 1.75e-3 * Zbarb + 0.37 * (1.0 - Math.exp(-0.015 * Math.pow(Zbarb, 1.3)));
         final double Wbar = 0.595 + etabar / 3.7 + Math.pow(etabar, 4.55);
         final double q = (2.0 * Wbar - 1.0) / (1.0 - Wbar);
         final double Ju0 = 1.0 + u0 * (Math.log(u0) - 1.0);
         final double Gu0 = (u0 - 1.0 - (1.0 - 1.0 / Math.pow(u0, 1.0 + q)) / (1.0 + q)) / ((2.0 + q) * Ju0);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oR = outputIndex(tag2("R", mComposition, mShell));
         final int oPhi0 = outputIndex(tag2(PHI0, mComposition, mShell));
         checkIndices(oR, oPhi0);

         final double R = 1.0 - etabar * Wbar * (1.0 - Gu0);
         rv.setEntry(oR, R);

         final double detabardZbarb = 0.00175 + (0.007215 * Math.pow(Zbarb, 0.3)) / Math.exp(0.015 * Math.pow(Zbarb, 1.3));

         final double dWbardZbarb = (0.27027027027027023 + 4.55 * Math.pow(etabar, 3.55)) * detabardZbarb;

         final double dqdZbarb = dWbardZbarb / Math.pow(-1.0 + Wbar, 2.0);

         final double dGu0dZbarb = -((Math.pow(u0, -1.0 - q) * (3.0 - 4.0 * Math.pow(u0, 1.0 + q) + Math.pow(u0, 2.0 + q)
               + 2.0 * Math.log(u0) + (2.0 - 4.0 * Math.pow(u0, 1.0 + q) + 2.0 * Math.pow(u0, 2.0 + q) + 3.0 * Math.log(u0)) * q
               + ((-1.0 + u0) * Math.pow(u0, 1.0 + q) + Math.log(u0)) * Math.pow(q, 2.0)) * dqdZbarb)
               / (Ju0 * Math.pow(1.0 + q, 2.0) * Math.pow(2.0 + q, 2.0)));

         final double dRdZbarb = -((1 - Gu0) * Wbar * detabardZbarb) + etabar * Wbar * dGu0dZbarb
               - etabar * (1.0 - Gu0) * dWbardZbarb;

         final double dJu0de0 = Math.log(u0) / Ea;

         final double dGu0de0 = ((1.0 / Ea - Ea / (Math.pow(e0, 2.0) * Math.pow(u0, q))) * Ju0
               + (1.0 - e0 / Ea + 1.0 / (1.0 + q) - Ea / (Math.pow(u0, q) * (e0 + e0 * q))) * dJu0de0)
               / ((2.0 + q) * Math.pow(Ju0, 2.0));

         final double dRde0 = etabar * Wbar * dGu0de0;

         rm.setEntry(oR, iZbarb, dRdZbarb);
         rm.setEntry(oR, iE0, dRde0);

         final double r = 2.0 - 2.3 * etabar;
         final double phi0 = 1.0 + 3.3 * (1.0 - 1.0 / Math.pow(u0, r)) * Math.pow(etabar, 1.2);
         rv.setEntry(oPhi0, phi0);
         rm.setEntry(oPhi0, iZbarb, detabardZbarb * (Math.pow(etabar, 0.2) * (3.96 - 3.96 * Math.pow(u0, -2 + 2.3 * etabar))
               - 7.59 * Math.pow(etabar, 1.2) * Math.pow(u0, -2 + 2.3 * etabar) * Math.log(u0)));
         rm.setEntry(oPhi0, iE0, (3.3 * Math.pow(etabar, 1.2) * r * Math.pow(u0, -1.0 - r)) / Ea);
         return Pair.create(rv, rm);
      }
   }

   private static class StepFRbar
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(RBAR, comp, shell));
         res.add(tag2("F", comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(new CompositionTag(ZBARB, comp));
         res.add(new CompositionTag(E_0, comp));
         res.add(tag2(ONE_OVER_S, comp, shell));
         res.add(tag2(QLA, comp, shell));
         res.add(tag2(PHI0, comp, shell));
         res.add(tag2("R", comp, shell));
         return res;
      }

      public StepFRbar(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp, shell), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iE0 = inputIndex(new CompositionTag(E_0, mComposition));
         final int iZbarb = inputIndex(new CompositionTag(ZBARB, mComposition));
         final int iOneOverS = inputIndex(tag2(ONE_OVER_S, mComposition, mShell));
         final int iQlaE0 = inputIndex(tag2(QLA, mComposition, mShell));
         final int iPhi0 = inputIndex(tag2(PHI0, mComposition, mShell));
         final int iR = inputIndex(tag2("R", mComposition, mShell));
         checkIndices(iE0, iZbarb, iOneOverS, iQlaE0, iPhi0, iR);

         final double Ea = eVtokeV(mShell.getEdgeEnergy());
         final double u0 = point.getEntry(iE0) / Ea;
         final double Zbarb = point.getEntry(iZbarb);
         final double oneOverS = point.getEntry(iOneOverS);
         final double QlaE0 = point.getEntry(iQlaE0);
         final double phi0 = point.getEntry(iPhi0);
         final double R = point.getEntry(iR);

         final double logZbarb = Math.log(Zbarb);
         final double X = 1.0 + 1.3 * logZbarb;
         final double Y = 0.2 + 0.005 * Zbarb;
         final double pu42 = Math.pow(u0, 0.42);
         final double FoRbar = 1.0 + (X * Math.log(1.0 + Y * (1.0 - 1.0 / pu42)) / Math.log(1.0 + Y));

         final double Zbarb240 = 240. + Zbarb;
         final double kk = Math.log(Zbarb240 / 200.0);
         final double dFoRbardZbarb = (((-200.0 + 200.0 * pu42 + (-260.0 + 260.0 * pu42) * logZbarb) * kk)
               / (-40.0 - Zbarb + pu42 * Zbarb240)
               + Math.log(1.2 + (-0.2 - Zbarb / 200.0) / pu42 + Zbarb / 200.0)
                     * ((-200.0 - 260.0 * logZbarb) / Zbarb240 + (260.0 * kk) / Zbarb))
               / (200.0 * kk * kk);
         final double dFoRbardu0 = (0.546 * (40.0 + Zbarb) * (0.7692307692307692 + logZbarb))
               / (u0 * (-40.0 - Zbarb + pu42 * Zbarb240) * kk);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oRbar = outputIndex(tag2(RBAR, mComposition, mShell));
         final int oF = outputIndex(tag2("F", mComposition, mShell));
         checkIndices(oRbar, oF);

         // F and partials
         final double F = R * oneOverS / QlaE0; // C1
         rv.setEntry(oF, F);
         final double dFdR = oneOverS / QlaE0; // C1
         rm.setEntry(oF, iR, dFdR);
         final double dFdOneOverS = R / QlaE0; // C1
         rm.setEntry(oF, iOneOverS, dFdOneOverS);
         final double dFdQlaE0 = -1.0 * R * oneOverS / Math.pow(QlaE0, 2.0); // C1
         rm.setEntry(oF, iQlaE0, dFdQlaE0);

         // Rbar and partials
         double Rbar = Double.NaN;
         double dRbardF = Double.NaN;
         if(FoRbar >= phi0) {
            Rbar = F / FoRbar; // C1
            dRbardF = 1.0 / FoRbar;
            final double dRbardFoRbar = -F / Math.pow(FoRbar, 2);
            rm.setEntry(oRbar, iE0, dRbardFoRbar * dFoRbardu0 / Ea);
            rm.setEntry(oRbar, iZbarb, dRbardFoRbar * dFoRbardZbarb);
         } else {
            Rbar = F / phi0;
            dRbardF = 1.0 / phi0;
            rm.setEntry(oRbar, iPhi0, -F / Math.pow(phi0, 2));
         }
         rv.setEntry(oRbar, Rbar);
         rm.setEntry(oRbar, iR, dRbardF * dFdR);
         rm.setEntry(oRbar, iOneOverS, dRbardF * dFdOneOverS);
         rm.setEntry(oRbar, iQlaE0, dRbardF * dFdQlaE0);

         return Pair.create(rv, rm);
      }
   }

   private static class StepP
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2("P", comp, shell));
         res.add(tag2("b", comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(RBAR, comp, shell));
         res.add(tag2("F", comp, shell));
         res.add(new CompositionTag(ZBARB, comp));
         res.add(new CompositionTag(E_0, comp));
         res.add(tag2(PHI0, comp, shell));
         return res;
      }

      StepP(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp, shell), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iZbarb = inputIndex(new CompositionTag(ZBARB, mComposition));
         final int iE0 = inputIndex(new CompositionTag(E_0, mComposition));
         final int iRbar = inputIndex(tag2(RBAR, mComposition, mShell));
         final int iF = inputIndex(tag2("F", mComposition, mShell));
         final int iphi0 = inputIndex(tag2(PHI0, mComposition, mShell));
         checkIndices(iZbarb, iE0, iRbar, iF, iphi0);

         final double Zbarb = point.getEntry(iZbarb);
         final double e0 = point.getEntry(iE0);
         final double Ea = eVtokeV(mShell.getEdgeEnergy());
         final double u0 = e0 / Ea;
         final double Rbar = point.getEntry(iRbar);
         final double F = point.getEntry(iF);
         final double phi0 = point.getEntry(iphi0);

         // P and partials
         final double kg = Math.exp((-1.0 / 15.0) * Zbarb * (u0 - 1.0));
         final double g = 0.22 * Math.log(4.0 * Zbarb) * (1.0 - 2.0 * kg);
         final double h = 1.0 - 10.0 * (1.0 - 1.0 / (1.0 + 0.1 * u0)) / Math.pow(Zbarb, 2.0);
         final double b = Math.sqrt(2.0) * (1.0 + Math.sqrt(1.0 - Rbar * phi0 / F)) / Rbar;
         final double sqrt2 = Math.sqrt(2.0), Rbar2 = Math.pow(Rbar, 2.0);

         // Two different ways to compute gh4
         final double gh4_1 = g * Math.pow(h, 4.0);
         final double gh4_2 = 0.9 * b * Rbar2 * (b - 2.0 * phi0 / F);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oP = outputIndex(tag2("P", mComposition, mShell));
         final int ob = outputIndex(tag2("b", mComposition, mShell));
         checkIndices(oP, ob);

         if(gh4_1 < gh4_2) {
            // Depends on Zbarb, u0, F, Rbar
            final double P = gh4_1 * F / Rbar2;
            rv.setEntry(oP, P);
            final double k1 = Math.exp(((-1.0 + u0) * Zbarb) / 15.0);
            final double Zbarb2 = Math.pow(Zbarb, 2), Zbarb3 = Math.pow(Zbarb, 3.0);
            final double k2 = 10.0 * Zbarb2 + u0 * (-10.0 + Zbarb2);
            final double u02 = Math.pow(u0, 2.0);
            rm.setEntry(oP, iZbarb, (11.0 * F * Math.pow(k2, 3.0)
                  * (15.0 * (-2.0 + k1) * k2 + 2.0 * (-10.0 * Zbarb3 + u02 * Zbarb * (-10.0 + Zbarb2)
                        + u0 * (-1200.0 + 600.0 * k1 + 10.0 * Zbarb + 9.0 * Zbarb3)) * Math.log(4.0 * Zbarb)))
                  / (750.0 * k1 * Rbar2 * Math.pow(10.0 + u0, 4.0) * Math.pow(Zbarb, 9.0)));
            rm.setEntry(oP, iE0, ((11.0 * F * Math.pow(k2, 3.0)
                  * (-3000.0 * k1 + u02 * Zbarb * (-10.0 + Zbarb2) + 20.0 * u0 * Zbarb * (-5.0 + Zbarb2)
                        + 100.0 * (60.0 + Zbarb3))
                  * Math.log(4.0 * Zbarb)) / (375.0 * k1 * Rbar2 * Math.pow(10.0 + u0, 5.0) * Math.pow(Zbarb, 8.0))) / Ea);
            rm.setEntry(oP, iF, P / F);
            rm.setEntry(oP, iRbar, -2.0 * P / Rbar);
         } else {
            // Depends on F, Rbar, b, phi0
            final double P = gh4_2 * F / Rbar2;
            rv.setEntry(oP, P);
            final double k2 = Math.sqrt(1.0 - (phi0 * Rbar) / F);
            rm.setEntry(oP, iphi0, (-0.9 * sqrt2
                  * (-3.0 * phi0 * Math.pow(Rbar, 3.0) + F * Rbar2 * (2.0 + 2.0 * k2) + Math.pow(F, 2.0) * (1. + k2) * sqrt2))
                  / (F * Math.pow(Rbar, 3.0) * k2));
            rm.setEntry(oP, iF, (-0.9 * sqrt2
                  * (Math.pow(phi0, 2.0) * Math.pow(Rbar, 4.0) + Math.pow(F, 3.0) * (-4.0 - 4.0 * k2) * sqrt2
                        + Math.pow(F, 2.0) * phi0 * Rbar * (3.0 + k2) * sqrt2))
                  / (Math.pow(F, 2.0) * Math.pow(Rbar, 4.0) * k2));
            rm.setEntry(oP, iRbar, (-7.2 * sqrt2
                  * (0.125 * Math.pow(phi0, 2.0) * Math.pow(Rbar, 4.0) + F * phi0 * Math.pow(Rbar, 3.0) * (-0.25 - 0.25 * k2)
                        + Math.pow(F, 2.0) * phi0 * Rbar * (-0.875 - 0.375 * k2) * sqrt2
                        + Math.pow(F, 3.0) * (1.0 + 1.0 * k2) * sqrt2))
                  / (F * Math.pow(Rbar, 5.0) * k2));
         }

         final double k3 = Math.sqrt(1. - (phi0 * Rbar) / F);

         rv.setEntry(ob, b);
         rm.setEntry(ob, iRbar, (-1. * (-0.5 * phi0 * Rbar + F * (1. + k3)) * sqrt2) / (F * Rbar2 * k3));
         rm.setEntry(ob, iphi0, -sqrt2 / (2. * F * k3));
         rm.setEntry(ob, iF, (phi0 * sqrt2) / (2. * Math.pow(F, 2) * k3));

         return Pair.create(rv, rm);
      }

   }

   private static class Stepa // C1
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2("a", comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(PHI0, comp, shell));
         res.add(tag2("P", comp, shell));
         res.add(tag2(RBAR, comp, shell));
         res.add(tag2("F", comp, shell));
         res.add(tag2("b", comp, shell));
         return res;
      }

      public Stepa(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp, shell), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iphi0 = inputIndex(tag2(PHI0, mComposition, mShell));
         final int iF = inputIndex(tag2("F", mComposition, mShell));
         final int iP = inputIndex(tag2("P", mComposition, mShell));
         final int iRbar = inputIndex(tag2(RBAR, mComposition, mShell));
         final int ib = inputIndex(tag2("b", mComposition, mShell));
         checkIndices(iphi0, iF, iP, iRbar, ib);

         final double phi0 = point.getEntry(iphi0);
         final double P = point.getEntry(iP);
         final double Rbar = point.getEntry(iRbar);
         final double F = point.getEntry(iF);
         final double b = point.getEntry(ib);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oa = outputIndex(tag2("a", mComposition, mShell));
         checkIndices(oa);

         final double b2 = Math.pow(b, 2.0);
         final double den = Math.pow(phi0 + b2 * F * Rbar - 2.0 * b * F, 2.0);
         final double a = (P + b * (2.0 * phi0 - b * F)) / (b * F * (2.0 - b * Rbar) - phi0); // C1

         rv.setEntry(oa, a); // C1
         rm.setEntry(oa, iP, 1.0 / (b * F * (2.0 - b * Rbar) - phi0)); // C1
         rm.setEntry(oa, ib, (-2.0 * (F * P + phi0 * phi0 - b * F * (phi0 + P * Rbar) + b2 * F * (F - phi0 * Rbar))) / den); // C1
         rm.setEntry(oa, iphi0, (3.0 * b2 * F + P - 2.0 * b2 * b * F * Rbar) / den); // C1
         rm.setEntry(oa, iF, (b * (P * (-2.0 + b * Rbar) + b * phi0 * (2.0 * b * Rbar - 3.0))) / den); // C1
         rm.setEntry(oa, iRbar, (b2 * F * (P + 2.0 * b * phi0 - b2 * F)) / den); // C1

         return Pair.create(rv, rm);
      }
   };

   static public class StepEps // C1
      extends
      NamedMultivariateJacobianFunction {

      private static final double MIN_EPS = 1.0e-6;

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(EPS, comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2("a", comp, shell));
         res.add(tag2("b", comp, shell));
         return res;
      }

      public StepEps(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp, shell), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int ia = inputIndex(tag2("a", mComposition, mShell));
         final int ib = inputIndex(tag2("b", mComposition, mShell));
         checkIndices(ia, ib);

         final double a = point.getEntry(ia);
         final double b = point.getEntry(ib);

         final int oeps = outputIndex(tag2(EPS, mComposition, mShell));
         checkIndices(oeps);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         if(Math.abs((a - b) / b) >= MIN_EPS) {
            final double eps = (a - b) / b;
            rv.setEntry(oeps, eps); // C1
            rm.setEntry(oeps, ia, 1.0 / b); // C1
            rm.setEntry(oeps, ib, -a / (b * b)); // C1
         } else {
            rv.setEntry(oeps, MIN_EPS); // C1
            rm.setEntry(oeps, ia, 0.0); // C1
            rm.setEntry(oeps, ib, 0.0); // C1
         }
         return Pair.create(rv, rm);
      }
   }

   private static class StepAB // C1
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mComposition;
      private final AtomicShell mShell;

      public static List<? extends Object> buildOutputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2("A", comp, shell));
         res.add(tag2("B", comp, shell));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final AtomicShell shell) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(PHI0, comp, shell));
         res.add(tag2("F", comp, shell));
         res.add(tag2("P", comp, shell));
         res.add(tag2("b", comp, shell));
         res.add(tag2(EPS, comp, shell));
         return res;
      }

      public StepAB(final Composition comp, final AtomicShell shell) {
         super(buildInputs(comp, shell), buildOutputs(comp, shell));
         mComposition = comp;
         mShell = shell;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iphi0 = inputIndex(tag2(PHI0, mComposition, mShell));
         final int iF = inputIndex(tag2("F", mComposition, mShell));
         final int iP = inputIndex(tag2("P", mComposition, mShell));
         final int ib = inputIndex(tag2("b", mComposition, mShell));
         final int ieps = inputIndex(tag2(EPS, mComposition, mShell));
         checkIndices(iphi0, iF, iP, ib, ieps);

         final double phi0 = point.getEntry(iphi0);
         final double P = point.getEntry(iP);
         final double F = point.getEntry(iF);
         final double b = point.getEntry(ib);
         final double eps = point.getEntry(ieps);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());
         final int oA = outputIndex(tag2("A", mComposition, mShell));
         final int oB = outputIndex(tag2("B", mComposition, mShell));
         checkIndices(oA, oB);

         final double B = (b * b * F * (1.0 + eps) - P - phi0 * b * (2.0 + eps)) / eps; // C2

         rv.setEntry(oB, B);
         rm.setEntry(oB, iphi0, (-b * (2.0 + eps)) / eps); // C2
         rm.setEntry(oB, iP, -1.0 / eps); // C2
         rm.setEntry(oB, iF, (b * b * (1.0 + eps)) / eps); // C2
         rm.setEntry(oB, ib, (2.0 * b * F * (1.0 + eps) - phi0 * (2.0 + eps)) / eps); // C2
         rm.setEntry(oB, ieps, (P + 2 * b * phi0 - b * b * F) / (eps * eps)); // C2

         final double k1 = (1.0 + eps) / (eps * eps); // C1

         rv.setEntry(oA, k1 * (b * (b * F - 2.0 * phi0) - P) / b); // C1
         rm.setEntry(oA, iphi0, -2.0 * k1); // C1
         rm.setEntry(oA, iP, -k1 / b); // C1
         rm.setEntry(oA, iF, k1 * b); // C1
         rm.setEntry(oA, ib, k1 * (b * b * F + P) / (b * b)); // C1
         rm.setEntry(oA, ieps, (((2.0 + eps) * (P + 2.0 * b * phi0 - b * b * F)) / (b * Math.pow(eps, 3.0)))); // C1
         return Pair.create(rv, rm);
      }
   }

   private static class StepAaBb
      extends
      MultiStepNamedMultivariateJacobianFunction {

      private static List<NamedMultivariateJacobianFunction> buildSteps(final Composition comp, final AtomicShell shell) {
         final List<NamedMultivariateJacobianFunction> res = new ArrayList<>();
         res.add(new StepRphi0(comp, shell));
         res.add(new StepFRbar(comp, shell));
         res.add(new StepP(comp, shell));
         res.add(new Stepa(comp, shell));
         res.add(new StepEps(comp, shell));
         res.add(new StepAB(comp, shell));
         return res;
      }

      public StepAaBb(final Composition comp, final AtomicShell shell)
            throws ArgumentException {
         super("StepAaBb", buildSteps(comp, shell));
      }

   }

   private static class StepChi
      extends
      NamedMultivariateJacobianFunction { // C1

      private final Composition mComposition;
      private final CharacteristicXRay mXRay;

      public StepChi(final Composition comp, final CharacteristicXRay cxr) {
         super(buildInputs(comp, cxr), buildOutputs(comp, cxr));
         mComposition = comp;
         mXRay = cxr;
      }

      public static List<? extends Object> buildOutputs(final Composition comp, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(CHI, comp, cxr));
         return res;
      }

      public static List<? extends Object> buildInputs(final Composition comp, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         for(final Element elm : comp.getElementSet()) {
            res.add(macTag(elm, cxr));
            res.add(Composition.buildMassFractionTag(comp, elm));
         }
         res.add(new CompositionTag("TOA", comp));
         return res;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int itoa = inputIndex(new CompositionTag("TOA", mComposition));
         checkIndices(itoa);
         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         double chi = 0.0;
         final int ochi = outputIndex(tag2(CHI, mComposition, mXRay));
         checkIndices(ochi);

         final double toa = point.getEntry(itoa);
         assert (toa > 0.0) && (toa < 0.5 * Math.PI);
         final double csc = 1.0 / Math.sin(toa);

         for(final Element elm : mComposition.getElementSet()) {
            final int imac = inputIndex(macTag(elm, mXRay));
            final int imf = inputIndex(Composition.buildMassFractionTag(mComposition, elm));
            final double mac = point.getEntry(imac);
            final double mf = point.getEntry(imf);
            final double tmp = mac * mf * csc; // C1
            rm.setEntry(ochi, imac, mf * csc); // C1
            rm.setEntry(ochi, imf, mac * csc); // C1
            chi += tmp;
         }
         rm.setEntry(ochi, itoa, -1.0 * chi / Math.tan(toa)); // C1
         rv.setEntry(ochi, chi);

         return Pair.create(rv, rm);
      }
   }

   private static class StepFx
      extends
      NamedMultivariateJacobianFunction {

      public static List<? extends Object> buildInputs(final Composition comp, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         final AtomicShell shell = cxr.getInner();
         res.add(tag2("A", comp, shell));
         // res.add(tag2("a", comp, shell));
         res.add(tag2("B", comp, shell));
         res.add(tag2("b", comp, shell));
         res.add(tag2(PHI0, comp, shell));
         res.add(tag2(EPS, comp, shell));
         res.add(tag2(CHI, comp, cxr));
         return res;
      }

      public static List<? extends Object> buildOutputs(final Composition comp, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(F_CHI, comp, cxr));
         return res;
      }

      private final Composition mMaterial;
      private final CharacteristicXRay mXRay;

      public StepFx(final Composition unk, final CharacteristicXRay cxr) {
         super(buildInputs(unk, cxr), buildOutputs(unk, cxr));
         mMaterial = unk;
         mXRay = cxr;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iA = inputIndex(tag2("A", mMaterial, mXRay.getInner()));
         final int iB = inputIndex(tag2("B", mMaterial, mXRay.getInner()));
         final int ib = inputIndex(tag2("b", mMaterial, mXRay.getInner()));
         final int iPhi0 = inputIndex(tag2(PHI0, mMaterial, mXRay.getInner()));
         final int iChi = inputIndex(tag2(CHI, mMaterial, mXRay));
         final int ieps = inputIndex(tag2(EPS, mMaterial, mXRay.getInner()));
         checkIndices(iA, iB, ib, iPhi0, iChi, ieps);

         final double A = point.getEntry(iA);
         final double B = point.getEntry(iB);
         final double b = point.getEntry(ib);
         final double phi0 = point.getEntry(iPhi0);
         final double chi = point.getEntry(iChi);
         final double eps = point.getEntry(ieps);

         final int oFx = outputIndex(tag2(F_CHI, mMaterial, mXRay));
         checkIndices(oFx);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final double k0 = b + chi;
         final double k1 = chi + b * (1 + eps);

         rv.setEntry(oFx, (B / k0 - (A * b * eps) / k1 + phi0) / k0); // C2
         rm.setEntry(oFx, iA, -(1.0 / k0) + 1.0 / k1); // C2
         rm.setEntry(oFx, iB, Math.pow(k0, -2.0)); // C2
         rm.setEntry(oFx, ib, (-(B / k0) + (A * b * eps) / k1
               + k0 * (-(B / Math.pow(k0, 2.0)) - (A * chi * eps) / Math.pow(k1, 2.0)) - phi0) / Math.pow(k0, 2.0)); // C2
         rm.setEntry(oFx, iPhi0, 1.0 / k0); // C2
         rm.setEntry(oFx, iChi, (-2.0 * B + k0 * ((A * b * eps * (2.0 * k0 + b * eps)) / Math.pow(k1, 2.0) - phi0))
               / Math.pow(k0, 3.0)); // C2
         rm.setEntry(oFx, ieps, -((A * b) / Math.pow(k1, 2.0))); // C2
         return Pair.create(rv, rm);
      }
   }

   private static class StepZA
      extends
      NamedMultivariateJacobianFunction {

      private final Composition mUnknown;
      private final Composition mStandard;
      private final CharacteristicXRay mXRay;

      public static List<? extends Object> buildInputs(final Composition unk, final Composition std, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(F_CHI, unk, cxr));
         res.add(tag2("F", unk, cxr.getInner()));
         res.add(tag2(F_CHI, std, cxr));
         res.add(tag2("F", std, cxr.getInner()));
         return res;
      }

      public static List<? extends Object> buildOutputs(final Composition unk, final Composition std, final CharacteristicXRay cxr) {
         final List<Object> res = new ArrayList<>();
         res.add(tag2(F_CHI_F, unk, cxr));
         res.add(tag2(F_CHI_F, std, cxr));
         res.add(zaTag(unk, std, cxr));
         return res;
      }

      public StepZA(final Composition unk, final Composition std, final CharacteristicXRay cxr) {
         super(buildInputs(unk, std, cxr), buildOutputs(unk, std, cxr));
         mUnknown = unk;
         mStandard = std;
         mXRay = cxr;
      }

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
         final int iFxu = inputIndex(tag2(F_CHI, mUnknown, mXRay));
         final int iFxs = inputIndex(tag2(F_CHI, mStandard, mXRay));
         final int iFu = inputIndex(tag2("F", mUnknown, mXRay.getInner()));
         final int iFs = inputIndex(tag2("F", mStandard, mXRay.getInner()));
         checkIndices(iFxu, iFxs, iFu, iFs);

         final double Fxu = point.getEntry(iFxu);
         final double Fxs = point.getEntry(iFxs);
         final double Fu = point.getEntry(iFu);
         final double Fs = point.getEntry(iFs);

         final RealVector rv = new ArrayRealVector(getOutputDimension());
         final RealMatrix rm = MatrixUtils.createRealMatrix(getOutputDimension(), getInputDimension());

         final int oFxFu = outputIndex(tag2(F_CHI_F, mUnknown, mXRay));
         final int oFxFs = outputIndex(tag2(F_CHI_F, mStandard, mXRay));
         final int oZA = outputIndex(zaTag(mUnknown, mStandard, mXRay));
         checkIndices(oFxFu, oFxFs, oZA);

         rv.setEntry(oZA, Fxu / Fxs); // C2
         rm.setEntry(oZA, iFxu, 1.0 / Fxs); // C2
         rm.setEntry(oZA, iFxs, -1.0 * Fxu / (Fxs * Fxs)); // C2

         rv.setEntry(oFxFu, Fxu / Fu); // C2
         rm.setEntry(oFxFu, iFu, -Fxu / (Fu * Fu)); // C2
         rm.setEntry(oFxFu, iFxu, 1.0 / Fu); // C2

         rv.setEntry(oFxFs, Fxs / Fs); // C2
         rm.setEntry(oFxFs, iFs, -Fxs / (Fs * Fs)); // C2
         rm.setEntry(oFxFs, iFxs, 1.0 / Fs); // C2

         return Pair.create(rv, rm);
      }
   }

   private static CharacteristicXRaySet extractXRays(final Map<Composition, CharacteristicXRaySet> stds) {
      final CharacteristicXRaySet res = new CharacteristicXRaySet();
      for(final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet())
         res.addAll(me.getValue());
      return res;
   }

   /**
    * This class performs the full XPP calculation for a single composition and
    * the associated set of characteristic x-rays. It breaks the calculation
    * into a series of steps which involve 1) only the composition; 2) the
    * composition and the atomic shells; and 3) the composition and the
    * characteristic x-rays. It is designed to break the calculation down into
    * small independent blocks each of which is relatively fast to calculate.
    * This minimizes and postpones the need to build and copy large Jacobian
    * matrices and is thus slightly more computationally efficient.
    *
    * @author nritchie
    */
   private static final class StepXPP
      extends
      MultiStepNamedMultivariateJacobianFunction {

      private static List<NamedMultivariateJacobianFunction> buildSteps(final Composition comp, final CharacteristicXRaySet exrs)
            throws ArgumentException {
         final List<NamedMultivariateJacobianFunction> res = new ArrayList<>();
         res.add(new StepMJZBarb(comp));
         final Set<AtomicShell> shells = new HashSet<>();
         for(final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
            shells.add(cxr.getInner());
         {
            final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
            for(final AtomicShell shell : shells)
               step.add(new StepQlaE0OneOverS(comp, shell));
            res.add(NamedMultivariateJacobianFunctionBuilder.join("QlaOoS", step));
         }
         {
            final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
            for(final AtomicShell shell : shells)
               step.add(new StepAaBb(comp, shell));
            res.add(NamedMultivariateJacobianFunctionBuilder.join("AaBb", step));
         }
         {
            final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
            for(final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
               step.add(new StepChi(comp, cxr));
            res.add(NamedMultivariateJacobianFunctionBuilder.join("Chi", step));
         }
         {
            final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
            for(final CharacteristicXRay cxr : exrs.getSetOfCharacteristicXRay())
               step.add(new StepFx(comp, cxr));
            res.add(NamedMultivariateJacobianFunctionBuilder.join("Fx", step));
         }
         return res;
      }

      public StepXPP(final Composition comp, final CharacteristicXRaySet cxrs)
            throws ArgumentException {
         super("XPP[" + comp + "]", buildSteps(comp, cxrs));
      }

   };

   /**
    * This method joins together the individual Composition computations into a
    * single final step that outputs the final matrix corrections.
    *
    * @param unk
    * @param stds
    * @return List&lt;NamedMultivariateJacobianFunction&gt;
    * @throws ArgumentException
    */
   private static List<NamedMultivariateJacobianFunction> buildSteps(final Composition unk, final Map<Composition, CharacteristicXRaySet> stds)
         throws ArgumentException {
      final List<NamedMultivariateJacobianFunction> res = new ArrayList<>();
      {
         final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
         // final Set<CharacteristicXRay> scxr = new HashSet<>();
         final CharacteristicXRaySet cxrs = new CharacteristicXRaySet();
         for(final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet()) {
            step.add(new StepXPP(me.getKey(), me.getValue()));
            // scxr.addAll(me.getValue().getSetOfCharacteristicXRay());
            cxrs.addAll(me.getValue());
         }
         step.add(new StepXPP(unk, cxrs));
         res.add(NamedMultivariateJacobianFunctionBuilder.join("XPP", step));
      }
      {
         final List<NamedMultivariateJacobianFunction> step = new ArrayList<>();
         for(final Map.Entry<Composition, CharacteristicXRaySet> me : stds.entrySet())
            for(final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay())
               step.add(new StepZA(unk, me.getKey(), cxr));
         res.add(NamedMultivariateJacobianFunctionBuilder.join("ZA", step));
      }
      return res;
   }
   
   // public Composition 
   
   
   

   static public Object beamEnergyTag(final Composition comp) {
      return new CompositionTag(E_0, comp);
   }

   static public Object takeOffAngleTag(final Composition comp) {
      return new CompositionTag("TOA", comp);
   }

   static public Object CompositionTag(final Composition comp, final Element elm) {
      return Composition.buildMassFractionTag(comp, elm);
   }

   static public Object meanIonizationTag(final Element elm) {
      return new ElementTag("J", elm);
   }

   static public Object macTag(final Element elm, final CharacteristicXRay cxr) {
      return new MassAbsorptionCoefficient.ElementMAC(elm, cxr);
   }

   static public <H> Object tag2(final String name, final Composition comp, final H other) {
      return new CompositionTag2<>(name, comp, other);
   }

   static public Object zaTag(final Composition unk, final Composition std, final CharacteristicXRay cxr) {
      return new ResultTag(unk, std, cxr);
   }

   /**
    * @param unk Composition
    * @param std Composition
    * @param scxr A set of {@link CharacteristicXRay} objects. It is important
    *           that the edge energies must all be below the beam energy.
    * @throws ArgumentException
    */
   public XPPMatrixCorrection(final Composition unk, final Map<Composition, CharacteristicXRaySet> stds)
         throws ArgumentException {
      super("XPP Matrix Correction", buildSteps(unk, stds));
      mStandards = stds;
      mUnknown = unk;
   }
   
   /**
    * @param unk Composition
    * @param std Composition
    * @param scxr A {@link CharacteristicXRay} object. It is important that the
    *           edge energies must be below the beam energy.
    * @throws ArgumentException
    */
   public XPPMatrixCorrection(final Composition unk, final Composition std, final CharacteristicXRaySet exrs)
         throws ArgumentException {
      this(unk, Collections.singletonMap(std, exrs));
   }

   /**
    * @param unk Composition
    * @param std Composition
    * @param scxr A {@link CharacteristicXRay} object. It is important that the
    *           edge energies must be below the beam energy.
    * @throws ArgumentException
    */
   public XPPMatrixCorrection(final Composition unk, final Composition std, final CharacteristicXRay cxr)
         throws ArgumentException {
      this(unk, std, CharacteristicXRaySet.build(cxr));
   }

   public UncertainValue computeJi(final Element elm) {
      final double z = elm.getAtomicNumber();
      final double j = 1.0e-3 * z * (10.04 + 8.25 * Math.exp(-z / 11.22));
      return new UncertainValue(j, "J[" + elm.getAbbrev() + "]", 0.03 * j);
   }

   /**
    * Many of the input parameters are computed or tabulated. Some are input as
    * experimental conditions like beam energy.
    *
    * @param conditions
    * @return
    * @throws ArgumentException
    */
   public UncertainValues buildInputs(final UncertainValues conditions)
         throws ArgumentException {
      final Set<Element> elms = new TreeSet<>();
      elms.addAll(mUnknown.getElementSet());
      for(final Composition std : mStandards.keySet())
         elms.addAll(std.getElementSet());
      final UncertainValues mip = buildMeanIonizationPotentials(elms);
      final UncertainValues macs = buildMassAbsorptionCoefficients(elms, extractXRays(mStandards));
      final UncertainValues[] parts = new UncertainValues[4 + mStandards.size()];
      parts[0] = conditions;
      parts[1] = mip;
      parts[2] = macs;
      parts[3] = mUnknown;
      int i = 0;
      for(final Composition std : mStandards.keySet()) {
         parts[4 + i] = std.asMassFraction();
         ++i;
      }
      return UncertainValues.build(getInputTags(), parts);
   }

   private UncertainValues buildMeanIonizationPotentials(final Set<Element> elms) {
      final RealVector vals = new ArrayRealVector(elms.size());
      final RealVector var = new ArrayRealVector(vals.getDimension());
      final List<Object> tags = new ArrayList<>();
      int i = 0;
      for(final Element elm : elms) {
         tags.add(meanIonizationTag(elm));
         final UncertainValue j = computeJi(elm);
         vals.setEntry(i, j.doubleValue());
         var.setEntry(i, j.variance());
         ++i;
      }
      assert tags.size() == vals.getDimension();
      final UncertainValues mip = new UncertainValues(tags, vals, var);
      return mip;
   }

   private UncertainValues buildMassAbsorptionCoefficients(final Set<Element> elms, final CharacteristicXRaySet scxr) {
      final RealVector vals = new ArrayRealVector(elms.size() * scxr.size());
      final RealVector var = new ArrayRealVector(vals.getDimension());
      final List<Object> tags = new ArrayList<>();
      final MassAbsorptionCoefficient macInstance = MassAbsorptionCoefficient.instance();
      int i = 0;
      for(final Element elm : elms) {
         for(final CharacteristicXRay cxr : scxr.getSetOfCharacteristicXRay()) {
            tags.add(macTag(elm, cxr));
            final UncertainValue j = macInstance.compute(elm, cxr);
            vals.setEntry(i, j.doubleValue());
            var.setEntry(i, 0.01 * j.variance());
            ++i;
         }
      }
      assert tags.size() == vals.getDimension();
      return new UncertainValues(tags, vals, var);
   }

   @Override
   public String toHTML(final Mode mode) {
      final Report r = new Report("XPP");
      r.addHeader("XPP Matrix Correction");
      final Table tbl = new Table();
      tbl.addRow(Table.td("Unknown"), Table.td(mUnknown), Table.td());
      for(final Map.Entry<Composition, CharacteristicXRaySet> me : mStandards.entrySet()) {
         tbl.addRow(Table.td("Standard"), Table.td(me.getKey()), Table.td());
         for(final CharacteristicXRay cxr : me.getValue().getSetOfCharacteristicXRay())
            tbl.addRow(Table.td(), Table.td("X-Ray Line"), Table.td(cxr.toHTML(Mode.NORMAL)));
      }
      r.add(tbl);
      r.addSubHeader("Inputs / Outputs");
      r.addHTML(super.toHTML(mode));
      return r.toHTML(mode);
   }
}
