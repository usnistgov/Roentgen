package gov.nist.microanalysis.roentgen.physics;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>
 * Calculates the mass absorption coefficient by performing a log-log
 * interpolation on the Chantler-2005 database.
 * </p>
 * <p>
 * Chantler, C.T., Olsen, K., Dragoset, R.A., Chang, J., Kishore, A.R.,
 * Kotochigova, S.A., and Zucker, D.S. (2005), <i>X-Ray Form Factor, Attenuation
 * and Scattering Tables (version 2.1).</i> [Online] Available:
 * http://physics.nist.gov/ffast [Wednesday, 12-Oct-2016 08:02:04 EDT]. National
 * Institute of Standards and Technology, Gaithersburg, MD.
 * </p>
 * <p>
 * Originally published as Chantler, C.T., J. Phys. Chem. Ref. Data
 * <b>29</b>(4), 597-1048 (2000); and Chantler, C.T., J. Phys. Chem. Ref. Data
 * <b>24</b>, 71-643 (1995).
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author nritchie
 * @version $Rev: 280 $
 */
public class MassAbsorptionCoefficient {

   private static MassAbsorptionCoefficient mInstance = new MassAbsorptionCoefficient();

   public static MassAbsorptionCoefficient instance() {
      return mInstance;
   }

   private final Map<Element, PolynomialSplineFunction> mData;

   public static class ElementMAC
      implements
      IToHTML {

      private final Element mElement;
      private final XRay mXRay;

      public ElementMAC(final Element elm, final XRay xr) {
         mElement = elm;
         mXRay = xr;
      }

      @Override
      public String toString() {
         return "\u03BC[" + mElement.getAbbrev() + "," + mXRay.toString() + "]";
      }

      /**
       * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
       */
      @Override
      public String toHTML(final Mode mode) {
         return "[&mu;/&rho;][" + mElement.getAbbrev() + "," + mXRay.toHTML(Mode.TERSE) + "]";
      }

      @Override
      public int hashCode() {
         return Objects.hash(mElement, mXRay);
      }

      @Override
      public boolean equals(final Object obj) {
         if(this == obj)
            return true;
         if(obj == null)
            return false;
         if(getClass() != obj.getClass())
            return false;
         final ElementMAC other = (ElementMAC) obj;
         return mElement.equals(other.mElement) && (mXRay.equals(other.mXRay));
      }
   }

   public static class MaterialMAC
      implements
      IToHTML {
      private final Composition mComp;
      private final XRay mXRay;

      public MaterialMAC(final Composition mf, final XRay xr) {
         mComp = mf;
         mXRay = xr;
      }

      public Composition getComposition() {
         return mComp;
      }

      public XRay getXRay() {
         return mXRay;
      }

      @Override
      public int hashCode() {
         return Objects.hash(mComp, mXRay);
      }

      @Override
      public boolean equals(final Object obj) {
         if(this == obj)
            return true;
         if(obj == null)
            return false;
         if(getClass() != obj.getClass())
            return false;
         final MaterialMAC other = (MaterialMAC) obj;
         return (mComp == other.mComp) && mXRay.equals(other.mXRay);
      }

      /**
       * @param mode
       * @return
       * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
       */
      @Override
      public String toHTML(final Mode mode) {
         return "&mu;[" + mComp.toHTML(Mode.TERSE) + ", " + mXRay.toHTML(Mode.TERSE) + "]";
      }

      @Override
      public String toString() {
         return "MAC[" + HTML.stripTags(mComp.toHTML(Mode.TERSE)) + ", " + mXRay.toString() + "]";
      }

   }

   public static final double MIN_E = 0.0;
   public static final double MAX_E = 432.9451e3;

   private static final Map<Element, PolynomialSplineFunction> loadData() {
      final InputStream is = MassAbsorptionCoefficient.class.getResourceAsStream("Chantler2005MAC.csv");
      final BufferedReader br = new BufferedReader(new InputStreamReader(is));
      final Map<Element, PolynomialSplineFunction> res = new HashMap<>();
      try {
         String keVline = br.readLine();
         String macLine = keVline != null ? br.readLine() : null;
         Element elm = Element.Hydrogen;
         final LinearInterpolator uvi = new LinearInterpolator();
         while((keVline != null) && (macLine != null)) {
            final String[] eVItems = keVline.split(",");
            final String[] macItems = macLine.split(",");
            assert eVItems.length == macItems.length;
            final double[] eVres = new double[eVItems.length];
            final double[] macRes = new double[eVItems.length];
            for(int i = 0; i < eVItems.length; ++i) {
               eVres[i] = Math.log(1000.0 * Double.parseDouble(eVItems[i].trim()));
               macRes[i] = Math.log(Double.parseDouble(macItems[i].trim()));
            }
            // The tables are linearly interpolated in log-log space.
            final PolynomialSplineFunction psf = uvi.interpolate(eVres, macRes);
            res.put(elm, psf);
            keVline = br.readLine();
            macLine = keVline != null ? br.readLine() : null;
            elm = elm.next();
         }
      }
      catch(final Exception e) {
         final String msg = "Unable to read Chantler2005 MAC data file.";
         Globals.getLogger().error(msg, e);
         throw new Error(msg, e);
      }
      return res;
   }

   /**
    * Constructs a MassAbsorptionCoefficient2
    */
   private MassAbsorptionCoefficient() {
      mData = loadData();
   }

   /**
    * @param el The Element
    * @param eV The energy in eV
    * @return double
    */
   private double fractionalUncertainty(final Element el, final double eV) {
      double err = 0.0;
      if(eV < 200.0)
         err = 1.50; // 100-200%
      else if(eV < 500)
         err = 0.5 + (((1.0 - 0.5) * (eV - 200)) / (500.0 - 200.0)); // 50-100%
      else if(eV < 1000)
         err = 0.05 + (((0.20 - 0.05) * (eV - 500)) / (1000.0 - 500.0));
      for(final AtomicShell sh : AtomicShell.forElement(el)) {
         final double ee = sh.getEdgeEnergy();
         if(ee > eV)
            continue;
         final double delta = (eV - ee) / eV;
         assert delta >= 0.0;
         if(delta < 0.001)
            err = Math.max(err, 0.5);
         else
            switch(sh.getShell()) {
               case K:
                  if(delta < 0.1)
                     err = Math.max(err, 0.15); // 10-20%
                  else if(eV < (1.1 * ee))
                     err = Math.max(err, 0.03); // 3%
                  else
                     err = Math.max(err, 0.01); // 1%
                  break;
               case L1:
               case M1:
               case M2:
               case M3:
               case N1: // Special case for N shells added by
                  // NWMR
               case N3:
               case N4:
               case N5:
                  if(delta < 0.15)
                     err = Math.max(err, 0.225); // 15-30%
                  else if(delta < 0.4)
                     err = Math.max(err, 0.04); // 4%
                  else
                     err = Math.max(err, 0.01); // 1%
                  break;
               case L2:
               case L3:
               case M4:
               case M5:
               case N6: // Special case for N shells added
                  // by NWMR
               case N7:
                  if(delta < 0.15)
                     err = Math.max(err, 0.3); // 20-40%
                  else if(delta < 0.40)
                     err = Math.max(err, 0.04); // 4%
                  else
                     err = Math.max(err, 0.01); // 1%
                  break;
               default:
                  break;
            }
      }
      return err;
   }

   /**
    * Computes the mass absorption coefficient for the specified element at the
    * specified x-ray energy.
    *
    * @param el Element
    * @param energy in eV
    * @return UncertainValue In cm^2/g
    */
   public UncertainValue compute(final Element el, final double energy) {
      assert mData.get(el) != null;
      assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
      assert energy < MAX_E : "Energy too high.";
      final double res = Math.exp(mData.get(el).value(Math.log(energy)));
      return new UncertainValue(res, new ElementMAC(el, new XRay(energy)), res * fractionalUncertainty(el, energy));
   }

   /**
    * Returns a PolynomialSplineFunction that implements a cubic spline
    * interpolation function over the mass absorption coefficent values for the
    * specified element. The values returned by this function are likely to
    * agree well but not perfectly with those computed using the other compute
    * functions in this class which use a log-log interpolations.
    *
    * @param elm
    * @return UnivariateFunction
    */
   public PolynomialSplineFunction getMACFunction(final Element elm) {
      final PolynomialSplineFunction uvf = mData.get(elm);
      final double[] x = uvf.getKnots();
      final double[] y = new double[x.length];
      for(int i = 0; i < x.length; ++i)
         y[i] = uvf.value(x[i]);
      return (new SplineInterpolator()).interpolate(x, y);
   }

   /**
    * Returns an array containing the energies at which values of the MACs were
    * provided that served as the basis for the interpolation.
    *
    * @param elm
    * @return double[]
    */
   public double[] getKnots(final Element elm) {
      final double[] res = mData.get(elm).getKnots();
      for(int i = 0; i < res.length; ++i)
         res[i] = Math.exp(res[i]);
      return res;
   }

   /**
    * @param el Element
    * @param energy in eV
    * @return true if a MAC is available for this element
    */
   public boolean isAvailable(final Element el, final XRay xr) {
      assert xr != null;
      final double energy = xr.getEnergy();
      assert el != null;
      if(!mData.containsKey(el))
         return false;
      return (energy > MIN_E) && (energy <= MAX_E);
   }

   /**
    * Computes the mass absorption coefficient for the specified element at the
    * specified x-ray energy.
    *
    * @param el Element
    * @param XRayPhoton
    * @return UncertainValue In cm^2/g
    */
   public UncertainValue compute(final Element el, final XRay xr) {
      assert mData.get(el) != null;
      final double energy = xr.getEnergy();
      assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
      assert energy < MAX_E : "Energy too high.";
      final double res = Math.exp(mData.get(el).value(Math.log(energy)));
      return new UncertainValue(res, new ElementMAC(el, xr), res * fractionalUncertainty(el, energy));
   }

   /**
    * Returns the elemental MAC associated with the specified material and xray
    * energy in eV.
    *
    * @param element
    * @param energy (eV)
    * @return MAC in cm^2/g
    */
   public double computeQ(final Element el, final double energy) {
      assert mData.get(el) != null;
      assert energy > MIN_E : "Can't compute a MAC for a negative or zero x-ray energy";
      assert energy < MAX_E : "Energy too high.";
      return Math.exp(mData.get(el).value(Math.log(energy)));
   }

   /**
    * Computes the MAC associated with the specified material and xray energy in
    * eV.
    *
    * @param material
    * @param energy (eV)
    * @return MAC in cm^2/g
    */
   public double computeQ(final Composition comp, final double energy) {
      Preconditions.checkArgument(energy > MIN_E, "Can't compute a MAC for a negative or zero x-ray energy");
      Preconditions.checkArgument(energy < MAX_E, "Energy too high.");
      final Set<Element> elms = comp.getElementSet();
      double res = 0.0;
      for(final Element elm : elms)
         res += computeQ(elm, energy) * comp.getValue(elm).doubleValue();
      return res;
   }

   /**
    * Computes the MAC associated with the specified material and xray.
    *
    * @param material
    * @param xr
    * @return MAC in cm^2/g
    */
   public double computeQ(final Composition material, final XRay xr) {
      return computeQ(material, xr.getEnergy());
   }
}
