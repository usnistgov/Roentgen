package gov.nist.microanalysis.roentgen.tests.physics;

import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

public class CompositionTest2 {

   private static final boolean mHTML = true;

   @Test
   public void testStoich1()
         throws ParseException {
      final Composition sc = Composition.parse("NaAlSi3O8");
      final Set<Element> elementSet = sc.getElementSet();
      Assert.assertTrue(elementSet.contains(Element.Sodium));
      Assert.assertTrue(elementSet.contains(Element.Aluminum));
      Assert.assertTrue(elementSet.contains(Element.Silicon));
      Assert.assertTrue(elementSet.contains(Element.Oxygen));

      Assert.assertEquals(1.0, sc.getValue(Element.Sodium).doubleValue(), 1e-12);
      Assert.assertEquals(1.0, sc.getValue(Element.Aluminum).doubleValue(), 1e-12);
      Assert.assertEquals(3.0, sc.getValue(Element.Silicon).doubleValue(), 1e-12);
      Assert.assertEquals(8.0, sc.getValue(Element.Oxygen).doubleValue(), 1e-12);

      final UncertainValues vals = sc;
      Assert.assertEquals(1.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Sodium)), 1e-12);
      Assert.assertEquals(1.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Aluminum)), 1e-12);
      Assert.assertEquals(3.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Silicon)), 1e-12);
      Assert.assertEquals(8.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Oxygen)), 1e-12);

      Assert.assertTrue(sc.toHTML(Mode.NORMAL).equals("NaAlSi<sub>3</sub>O<sub>8</sub>"));

      final Composition af = sc.asAtomFraction();
      Assert.assertEquals(1.0 / 13.0, af.getValue(Element.Sodium).doubleValue(), 1.0e-9);
      Assert.assertEquals(1.0 / 13.0, af.getValue(Element.Aluminum).doubleValue(), 1.0e-9);
      Assert.assertEquals(3.0 / 13.0, af.getValue(Element.Silicon).doubleValue(), 1.0e-9);
      Assert.assertEquals(8.0 / 13.0, af.getValue(Element.Oxygen).doubleValue(), 1.0e-9);

      Assert.assertEquals(0.0, af.getValue(Element.Sodium).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Aluminum).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Silicon).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Oxygen).uncertainty(), 1.0e-9);

      final Composition mf = af.asMassFraction();
      Assert.assertEquals(mf.getValue(Element.Oxygen).doubleValue(), 0.488116, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Sodium).doubleValue(), 0.0876726, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Aluminum).doubleValue(), 0.102895, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Silicon).doubleValue(), 0.321316, 0.000001);

      if(mHTML) {
         try {
            final File f = File.createTempFile("AtomicFraction", ".html");
            try (FileWriter fr = new FileWriter(f)) {
               fr.append(HTML.createPageHeader("Atomic Fraction Report", null));
               fr.append(HTML.header("Stoichiometric Compound"));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(sc.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(sc.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(sc.toHTML(Mode.VERBOSE));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(af.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(af.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(af.toHTML(Mode.VERBOSE));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(mf.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(mf.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(mf.toHTML(Mode.VERBOSE));
               fr.append(HTML.createPageEnd());
            }
            Desktop.getDesktop().browse(f.toURI());
         }
         catch(final IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }

   @Test
   public void testStoich2()
         throws ParseException {
      final Composition sc = Composition.parse("Ca5(PO4)3(OH)");
      final Set<Element> elementSet = sc.getElementSet();
      Assert.assertTrue(elementSet.contains(Element.Calcium));
      Assert.assertTrue(elementSet.contains(Element.Phosphorus));
      Assert.assertTrue(elementSet.contains(Element.Oxygen));
      Assert.assertTrue(elementSet.contains(Element.Hydrogen));

      Assert.assertEquals(5.0, sc.getValue(Element.Calcium).doubleValue(), 1e-12);
      Assert.assertEquals(3.0, sc.getValue(Element.Phosphorus).doubleValue(), 1e-12);
      Assert.assertEquals(13.0, sc.getValue(Element.Oxygen).doubleValue(), 1e-12);
      Assert.assertEquals(1.0, sc.getValue(Element.Hydrogen).doubleValue(), 1e-12);

      final UncertainValues vals = sc;
      Assert.assertEquals(5.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Calcium)), 1e-12);
      Assert.assertEquals(3.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Phosphorus)), 1e-12);
      Assert.assertEquals(13.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Oxygen)), 1e-12);
      Assert.assertEquals(1.0, vals.getEntry(Composition.buildStoichiometryTag(sc, Element.Hydrogen)), 1e-12);

      Assert.assertTrue(sc.toHTML(Mode.TERSE).equals("Ca<sub>5</sub>(PO<sub>4</sub>)<sub>3</sub>(OH)"));
      Assert.assertTrue(sc.toHTML(Mode.NORMAL).equals("Ca<sub>5</sub>(PO<sub>4</sub>)<sub>3</sub>(OH)"));

      final Composition af = sc.asAtomFraction();

      Assert.assertEquals(5.0 / 22.0, af.getValue(Element.Calcium).doubleValue(), 1.0e-9);
      Assert.assertEquals(3.0 / 22.0, af.getValue(Element.Phosphorus).doubleValue(), 1.0e-9);
      Assert.assertEquals(13.0 / 22.0, af.getValue(Element.Oxygen).doubleValue(), 1.0e-9);
      Assert.assertEquals(1.0 / 22.0, af.getValue(Element.Hydrogen).doubleValue(), 1.0e-9);

      Assert.assertEquals(0.0, af.getValue(Element.Calcium).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Phosphorus).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Oxygen).uncertainty(), 1.0e-9);
      Assert.assertEquals(0.0, af.getValue(Element.Hydrogen).uncertainty(), 1.0e-9);

      final Composition mf = af.asMassFraction();
      Assert.assertEquals(mf.getValue(Element.Calcium).doubleValue(), 0.398936, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Phosphorus).doubleValue(), 0.184987, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Oxygen).doubleValue(), 0.41407, 0.000001);
      Assert.assertEquals(mf.getValue(Element.Hydrogen).doubleValue(), 0.0020066, 0.000001);

      if(mHTML) {
         try {
            final File f = File.createTempFile("AtomicFraction", ".html");
            try (FileWriter fr = new FileWriter(f)) {
               fr.append(HTML.createPageHeader("Atomic Fraction Report", null));
               fr.append(HTML.header("Stoichiometric Compound"));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(sc.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(sc.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(sc.toHTML(Mode.VERBOSE));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(af.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(af.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(af.toHTML(Mode.VERBOSE));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(mf.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(mf.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(mf.toHTML(Mode.VERBOSE));
               fr.append(HTML.createPageEnd());
            }
            Desktop.getDesktop().browse(f.toURI());
         }
         catch(final IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }

   @Test
   public void testMassFraction() {
      final Map<Element, Number> massFracs = new HashMap<>();
      final double total = 0.98;
      massFracs.put(Element.Calcium, new UncertainValue(0.398936 / total, "dCa", 0.012));
      massFracs.put(Element.Phosphorus, new UncertainValue(0.184987 / total, "dP", 0.008));
      massFracs.put(Element.Oxygen, new UncertainValue(0.41407 / total, "dO", 0.015));
      massFracs.put(Element.Hydrogen, new UncertainValue(0.0020066 / total, "dH", 0.001));
      final Composition mf = Composition.massFraction("Apatite", massFracs);

      final Composition af = mf.asAtomFraction();
      Assert.assertEquals(5.0 / 22.0, af.getValue(Element.Calcium).doubleValue(), 1.0e-6);
      Assert.assertEquals(3.0 / 22.0, af.getValue(Element.Phosphorus).doubleValue(), 1.0e-6);
      Assert.assertEquals(13.0 / 22.0, af.getValue(Element.Oxygen).doubleValue(), 1.0e-6);
      Assert.assertEquals(1.0 / 22.0, af.getValue(Element.Hydrogen).doubleValue(), 1.0e-6);

      Assert.assertEquals(0.0000767132, af.getCovariance(Element.Calcium, Element.Calcium), 1.0e-7);
      Assert.assertEquals(0.0000430936, af.getCovariance(Element.Phosphorus, Element.Phosphorus), 1.0e-7);
      Assert.assertEquals(0.000273062, af.getCovariance(Element.Oxygen, Element.Oxygen), 1.0e-6);
      Assert.assertEquals(0.000450103, af.getCovariance(Element.Hydrogen, Element.Hydrogen), 1.0e-6);
      Assert.assertEquals(-0.0000624525, af.getCovariance(Element.Hydrogen, Element.Phosphorus), 1.0e-6);
      Assert.assertEquals(0.0000176265, af.getCovariance(Element.Calcium, Element.Phosphorus), 1.0e-6);
      Assert.assertEquals(-0.000284053, af.getCovariance(Element.Oxygen, Element.Hydrogen), 1.0e-6);

      if(mHTML) {
         try {
            final File f = File.createTempFile("MassFraction", ".html");
            try (FileWriter fr = new FileWriter(f)) {
               fr.append(HTML.createPageHeader("Mass Fraction Report", null));
               fr.append(HTML.header("Mass Fraction"));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(mf.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(mf.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(mf.toHTML(Mode.VERBOSE));
               fr.append(HTML.subHeader("TERSE"));
               fr.append(af.toHTML(Mode.TERSE));
               fr.append(HTML.subHeader("NORMAL"));
               fr.append(af.toHTML(Mode.NORMAL));
               fr.append(HTML.subHeader("VERBOSE"));
               fr.append(af.toHTML(Mode.VERBOSE));
               final Table t = new Table();
               t.addRow(Table.th("Element"), Table.thc("Verbose"), Table.thc("Normal"));
               for(final Element elm : af.getElementSet())
                  t.addRow(Table.td(elm.toHTML(Mode.TERSE)), Table.tdc(af.getValue(elm).toHTML(Mode.VERBOSE)), Table.tdc(af.getValue(elm).toHTML(Mode.NORMAL)));
               fr.append(t.toHTML(Mode.NORMAL));
               fr.append(HTML.createPageEnd());
            }
            Desktop.getDesktop().browse(f.toURI());
         }
         catch(final IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }

   @Test
   public void testAtomicFraction() {
      final Map<Element, Number> atFracs = new HashMap<>();
      atFracs.put(Element.Oxygen, new UncertainValue(0.593981, "dO", 0.05));
      atFracs.put(Element.Magnesium, new UncertainValue(0.106595, "dMg", 0.005));
      atFracs.put(Element.Aluminum, new UncertainValue(0.0404141, "dAl", 0.001));
      atFracs.put(Element.Silicon, new UncertainValue(0.167755, "dSi", 0.002));
      atFracs.put(Element.Calcium, new UncertainValue(0.0604423, "dCa", 0.0021));
      atFracs.put(Element.Iron, new UncertainValue(0.0308124, "dFe", 0.0012));
      final Composition af = Composition.atomFraction("K412", atFracs);

      final Composition mf = af.asMassFraction();
      Assert.assertEquals(0.431202, mf.getValue(Element.Oxygen).doubleValue(), 1.0e-5);
      Assert.assertEquals(0.117554, mf.getValue(Element.Magnesium).doubleValue(), 1.0e-5);
      Assert.assertEquals(0.049477, mf.getValue(Element.Aluminum).doubleValue(), 1.0e-5);
      Assert.assertEquals(0.21377, mf.getValue(Element.Silicon).doubleValue(), 1.0e-5);
      Assert.assertEquals(0.109914, mf.getValue(Element.Calcium).doubleValue(), 1.0e-5);
      Assert.assertEquals(0.0780755, mf.getValue(Element.Iron).doubleValue(), 1.0e-5);

      Assert.assertEquals(0.0209244, mf.getValue(Element.Oxygen).uncertainty(), 1.0e-5);
      Assert.assertEquals(0.00650561, mf.getValue(Element.Magnesium).uncertainty(), 1.0e-5);
      Assert.assertEquals(0.00217441, mf.getValue(Element.Aluminum).uncertainty(), 1.0e-5);
      Assert.assertEquals(0.00817155, mf.getValue(Element.Silicon).uncertainty(), 1.0e-5);
      Assert.assertEquals(0.00529588, mf.getValue(Element.Calcium).uncertainty(), 1.0e-5);
      Assert.assertEquals(0.00402649, mf.getValue(Element.Iron).uncertainty(), 1.0e-5);

      Assert.assertEquals(-0.0000857086, mf.getCovariance(Element.Oxygen, Element.Calcium), 1.0e-7);
      Assert.assertEquals(9.2027e-6, mf.getCovariance(Element.Magnesium, Element.Iron), 1.0e-7);
      Assert.assertEquals(9.68555e-6, mf.getCovariance(Element.Iron, Element.Calcium), 1.0e-7);
      Assert.assertEquals(0.0000274102, mf.getCovariance(Element.Silicon, Element.Magnesium), 1.0e-7);
      Assert.assertEquals(0.0000139519, mf.getCovariance(Element.Aluminum, Element.Silicon), 1.0e-7);
      Assert.assertEquals(0.0, mf.getCovariance(Element.Oxygen, Element.Antimony), 1.e-12);
      Assert.assertEquals(0.0, mf.getCovariance(Element.Antimony, Element.Silver), 1.e-12);

      if(mHTML) {
         final Report rep = new Report("Atomic Fraction");
         try {
            rep.addHeader("Atom Fraction");
            rep.addSubHeader("TERSE");
            rep.addHTML(mf.toHTML(Mode.TERSE));
            rep.addSubHeader("NORMAL");
            rep.addHTML(mf.toHTML(Mode.NORMAL));
            rep.addSubHeader("VERBOSE");
            rep.addHTML(mf.toHTML(Mode.VERBOSE));
            rep.addSubHeader("TERSE");
            rep.addHTML(af.toHTML(Mode.TERSE));
            rep.addSubHeader("NORMAL");
            rep.addHTML(af.toHTML(Mode.NORMAL));
            rep.addSubHeader("VERBOSE");
            rep.addHTML(af.toHTML(Mode.VERBOSE));
            final Table t = new Table();
            t.addRow(Table.th("Element"), Table.thc("Verbose"), Table.thc("Normal"));
            for(final Element elm : mf.getElementSet())
               t.addRow(Table.td(elm.toHTML(Mode.TERSE)), Table.tdc(mf.getValue(elm).toHTML(Mode.VERBOSE)), Table.tdc(mf.getValue(elm).toHTML(Mode.NORMAL)));
            rep.addHTML(t.toHTML(Mode.NORMAL));
         }
         catch(final Exception e) {
            rep.addHTML(HTML.error(e.getMessage()));
         }
         try {
            rep.inBrowser(Mode.NORMAL);
         }
         catch(final IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
      }
   }

   @Test
   public void testAtomicFractionMC() throws ArgumentException {
      final Map<Element, Number> atFracs = new HashMap<>();
      atFracs.put(Element.Oxygen, new UncertainValue(0.593981, "dO", 0.05));
      atFracs.put(Element.Magnesium, new UncertainValue(0.106595, "dMg", 0.005));
      atFracs.put(Element.Aluminum, new UncertainValue(0.0404141, "dAl", 0.001));
      atFracs.put(Element.Silicon, new UncertainValue(0.167755, "dSi", 0.002));
      atFracs.put(Element.Calcium, new UncertainValue(0.0604423, "dCa", 0.0021));
      atFracs.put(Element.Iron, new UncertainValue(0.0308124, "dFe", 0.0012));
      final Composition af = Composition.atomFraction("K412", atFracs);

      final UncertainValues uv = af;
      final UncertainValues mup = UncertainValues.propagateMC(new Composition.AtomFractionToMassFraction(af), uv, 160000);
      final Composition mf = af.asMassFraction();
      // for(Element elm )

      if(mHTML) {
         final Report rep = new Report("Atomic Fraction MC");
         try {
            rep.addHeader("MC: Atom Fraction to Mass Fraction");
            rep.addSubHeader("TERSE");
            rep.addHTML(mup.toHTML(Mode.TERSE));
            rep.addSubHeader("NORMAL");
            rep.addHTML(mup.toHTML(Mode.NORMAL));
            rep.addSubHeader("VERBOSE - Analytic");
            rep.addHTML(mf.toHTML(Mode.VERBOSE));
            rep.addSubHeader("VERBOSE - Monte Carlo");
            rep.addHTML(mup.toHTML(Mode.VERBOSE, Composition.MASS_FRACTION_FORMAT));
         }

         catch(final Exception e) {
            rep.addHTML(HTML.error(e.getMessage()));
         }
         try {
            rep.inBrowser(Mode.NORMAL);
         }
         catch(final IOException e) {
            e.printStackTrace();
         }
      }
   }
}
