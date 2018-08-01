package gov.nist.microanalysis.roentgen.tests.spectrum;

import java.awt.Color;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.math.MathUtilities;
import gov.nist.microanalysis.roentgen.math.MathUtilities.IDoubleAsColor;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.spectrum.AdaptiveGaussianFilter;
import gov.nist.microanalysis.roentgen.spectrum.AdaptiveTophatFilter;
import gov.nist.microanalysis.roentgen.spectrum.EDSSpectrum;
import gov.nist.microanalysis.roentgen.spectrum.EnergyCalibration;
import gov.nist.microanalysis.roentgen.spectrum.LineshapeCalibration;

/**
 * <p>
 * Tests classes derived from
 * {@link gov.nist.microanalysis.roentgen.spectrum.EDSFittingFilter}
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestEDSFittingFilter {

   /**
    * Constructs a TestEDSFittingFilter
 * @throws Exception 
    */
   @Test
   public void testGaussianExtent() throws Exception {
      final EnergyCalibration ec = EnergyCalibration.Linear(0.0, 10.0, 2048);
      final LineshapeCalibration ls = new LineshapeCalibration.Gaussian(130.0, LineshapeCalibration.Gaussian.SDD_EV_PER_EH);
      final AdaptiveGaussianFilter agf = new AdaptiveGaussianFilter(2048, ec, ls);
      System.out.println("Gaussian");
      System.out.println(agf.extent(200, 1.0e-4));
      System.out.println(agf.extent(2000, 1.0e-4));
      System.out.println(agf.extent(10000, 1.0e-4));

      System.out.println("Mg = " + agf.extents(Element.Magnesium, 20.0e3, 1.0e-4).toString());
      System.out.println("Cu = " + agf.extents(Element.Copper, 20.0e3, 1.0e-4).toString());
      System.out.println("Ca = " + agf.extents(Element.Calcium, 20.0e3, 1.0e-4).toString());
      System.out.println("Ag = " + agf.extents(Element.Silver, 20.0e3, 1.0e-4).toString());
      System.out.println("Au = " + agf.extents(Element.Gold, 20.0e3, 1.0e-4).toString());
      System.out.println("Pb = " + agf.extents(Element.Lead, 20.0e3, 1.0e-4).toString());
      
      final EDSSpectrum spec = ReadSpectrum.fromResource("Fe_ref1.msa");
      final AdaptiveTophatFilter agf2 = new AdaptiveTophatFilter(spec.size(), spec.getEnergyCalibration(), ls);
      final Pair<RealVector, RealMatrix> pr = agf2.evaluate(spec.getData());
      final Report r = new Report("Fe_ref1.msa - Tophat");
      final RealMatrix jac = pr.getSecond();
      final double ext = Math.max(Math.abs(MathUtilities.min(jac)), Math.abs(MathUtilities.max(jac)));
      final IDoubleAsColor colorize = new MathUtilities.PositiveNegativeAsColor(-ext, Color.RED, ext, Color.GREEN, Color.YELLOW, 1.0);
      r.addImage(MathUtilities.RealMatrixAsBitmap(pr.getValue(), colorize), "Jacobian matrix - Gaussian");
      r.inBrowser(Mode.VERBOSE);
   }

   @Test
   public void testTophatExtent()
         throws Exception {
      final EnergyCalibration ec = EnergyCalibration.Linear(0.0, 10.0, 2048);
      final LineshapeCalibration ls = new LineshapeCalibration.Gaussian(130.0, LineshapeCalibration.Gaussian.SDD_EV_PER_EH);
      final AdaptiveTophatFilter agf = new AdaptiveTophatFilter(2048, ec, ls);
      System.out.println("Tophat");
      System.out.println(agf.extent(200, 1.0e-2));
      System.out.println(agf.extent(200, 1.0e-3));
      System.out.println(agf.extent(200, 1.0e-4));
      System.out.println(agf.extent(2000, 1.0e-4));
      System.out.println(agf.extent(10000, 1.0e-4));

      agf.setUseFastExtent(true);
      System.out.println("Fast = True");
      System.out.println("Mg = " + agf.extents(Element.Magnesium, 20.0e3, 1.0e-4).toString());
      System.out.println("Cu = " + agf.extents(Element.Copper, 20.0e3, 1.0e-4).toString());
      System.out.println("Ca = " + agf.extents(Element.Calcium, 20.0e3, 1.0e-4).toString());
      System.out.println("Ag = " + agf.extents(Element.Silver, 20.0e3, 1.0e-4).toString());
      System.out.println("Au = " + agf.extents(Element.Gold, 20.0e3, 1.0e-4).toString());
      System.out.println("Pb = " + agf.extents(Element.Lead, 20.0e3, 1.0e-4).toString());
      agf.setUseFastExtent(false);
      System.out.println("Fast = False");
      System.out.println("Mg = " + agf.extents(Element.Magnesium, 20.0e3, 1.0e-4).toString());
      System.out.println("Cu = " + agf.extents(Element.Copper, 20.0e3, 1.0e-4).toString());
      System.out.println("Ca = " + agf.extents(Element.Calcium, 20.0e3, 1.0e-4).toString());
      System.out.println("Ag = " + agf.extents(Element.Silver, 20.0e3, 1.0e-4).toString());
      System.out.println("Au = " + agf.extents(Element.Gold, 20.0e3, 1.0e-4).toString());
      System.out.println("Pb = " + agf.extents(Element.Lead, 20.0e3, 1.0e-4).toString());

      final EDSSpectrum spec = ReadSpectrum.fromResource("Fe_ref1.msa");
      final AdaptiveTophatFilter agf2 = new AdaptiveTophatFilter(spec.size(), spec.getEnergyCalibration(), ls);
      final Pair<RealVector, RealMatrix> pr = agf2.evaluate(spec.getData());
      final Report r = new Report("Fe_ref1.msa - Tophat");
      final RealMatrix jac = pr.getSecond();
      final double ext = Math.max(Math.abs(MathUtilities.min(jac)), Math.abs(MathUtilities.max(jac)));
      final IDoubleAsColor colorize = new MathUtilities.PositiveNegativeAsColor(-ext, Color.RED, ext, Color.GREEN, Color.YELLOW, 1.0);
      r.addImage(MathUtilities.RealMatrixAsBitmap(pr.getValue(), colorize), "Jacobian matrix - Tophat");
      r.inBrowser(Mode.VERBOSE);
   }

}
