package gov.nist.microanalysis.roentgen.tests.physics;

import static org.junit.Assert.assertEquals;

import java.awt.Color;
import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;
import com.duckandcover.html.Table;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.LabeledMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition.Representation;
import gov.nist.microanalysis.roentgen.physics.composition.Material;
import gov.nist.microanalysis.roentgen.physics.composition.MaterialLabel;
import gov.nist.microanalysis.roentgen.swing.LinearToColor;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class CompositionTest2 {

	private static final boolean mHTML = true;

	@Test
	public void testStoich1() throws ParseException, ArgumentException {
		final Composition sc = Composition.parse("NaAlSi3O8");
		final Set<Element> elementSet = sc.getElementSet();
		Assert.assertTrue(elementSet.contains(Element.Sodium));
		Assert.assertTrue(elementSet.contains(Element.Aluminum));
		Assert.assertTrue(elementSet.contains(Element.Silicon));
		Assert.assertTrue(elementSet.contains(Element.Oxygen));

		Assert.assertEquals(1.0, sc.getStoichiometry(Element.Sodium).doubleValue(), 1e-12);
		Assert.assertEquals(1.0, sc.getStoichiometry(Element.Aluminum).doubleValue(), 1e-12);
		Assert.assertEquals(3.0, sc.getStoichiometry(Element.Silicon).doubleValue(), 1e-12);
		Assert.assertEquals(8.0, sc.getStoichiometry(Element.Oxygen).doubleValue(), 1e-12);

		final Material mat = sc.getMaterial();
		Assert.assertEquals(1.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Sodium)), 1e-12);
		Assert.assertEquals(1.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Aluminum)), 1e-12);
		Assert.assertEquals(3.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Silicon)), 1e-12);
		Assert.assertEquals(8.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Oxygen)), 1e-12);

		Assert.assertTrue(sc.toHTML(Mode.NORMAL).equals("NaAlSi<sub>3</sub>O<sub>8</sub>"));

		Assert.assertEquals(1.0 / 13.0, sc.getAtomFraction(Element.Sodium).doubleValue(), 1.0e-9);
		Assert.assertEquals(1.0 / 13.0, sc.getAtomFraction(Element.Aluminum).doubleValue(), 1.0e-9);
		Assert.assertEquals(3.0 / 13.0, sc.getAtomFraction(Element.Silicon).doubleValue(), 1.0e-9);
		Assert.assertEquals(8.0 / 13.0, sc.getAtomFraction(Element.Oxygen).doubleValue(), 1.0e-9);

		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Sodium).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Aluminum).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Silicon).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Oxygen).uncertainty(), 1.0e-9);

		Assert.assertEquals(sc.getMassFraction(Element.Oxygen).doubleValue(), 0.488116, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Sodium).doubleValue(), 0.0876726, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Aluminum).doubleValue(), 0.102895, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Silicon).doubleValue(), 0.321316, 0.000001);

		if (mHTML) {
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
					fr.append(HTML.createPageEnd());
				}
				Desktop.getDesktop().browse(f.toURI());
			} catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Test
	public void testStoich2() throws ParseException, ArgumentException {
		final Composition sc = Composition.parse("Ca5(PO4)3(OH)");
		final Set<Element> elementSet = sc.getElementSet();
		Assert.assertTrue(elementSet.contains(Element.Calcium));
		Assert.assertTrue(elementSet.contains(Element.Phosphorus));
		Assert.assertTrue(elementSet.contains(Element.Oxygen));
		Assert.assertTrue(elementSet.contains(Element.Hydrogen));

		Assert.assertEquals(5.0, sc.getStoichiometry(Element.Calcium).doubleValue(), 1e-12);
		Assert.assertEquals(3.0, sc.getStoichiometry(Element.Phosphorus).doubleValue(), 1e-12);
		Assert.assertEquals(13.0, sc.getStoichiometry(Element.Oxygen).doubleValue(), 1e-12);
		Assert.assertEquals(1.0, sc.getStoichiometry(Element.Hydrogen).doubleValue(), 1e-12);

		final Material mat = sc.getMaterial();
		Assert.assertEquals(5.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Calcium)), 1e-12);
		Assert.assertEquals(3.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Phosphorus)), 1e-12);
		Assert.assertEquals(13.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Oxygen)), 1e-12);
		Assert.assertEquals(1.0, sc.getEntry(MaterialLabel.buildStoichiometryTag(mat, Element.Hydrogen)), 1e-12);

		Assert.assertTrue(sc.toHTML(Mode.TERSE).equals("Ca<sub>5</sub>(PO<sub>4</sub>)<sub>3</sub>(OH)"));
		Assert.assertTrue(sc.toHTML(Mode.NORMAL).equals("Ca<sub>5</sub>(PO<sub>4</sub>)<sub>3</sub>(OH)"));

		Assert.assertEquals(5.0 / 22.0, sc.getAtomFraction(Element.Calcium).doubleValue(), 1.0e-9);
		Assert.assertEquals(3.0 / 22.0, sc.getAtomFraction(Element.Phosphorus).doubleValue(), 1.0e-9);
		Assert.assertEquals(13.0 / 22.0, sc.getAtomFraction(Element.Oxygen).doubleValue(), 1.0e-9);
		Assert.assertEquals(1.0 / 22.0, sc.getAtomFraction(Element.Hydrogen).doubleValue(), 1.0e-9);

		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Calcium).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Phosphorus).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Oxygen).uncertainty(), 1.0e-9);
		Assert.assertEquals(0.0, sc.getAtomFraction(Element.Hydrogen).uncertainty(), 1.0e-9);

		Assert.assertEquals(sc.getMassFraction(Element.Calcium).doubleValue(), 0.398936, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Phosphorus).doubleValue(), 0.184987, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Oxygen).doubleValue(), 0.41407, 0.000001);
		Assert.assertEquals(sc.getMassFraction(Element.Hydrogen).doubleValue(), 0.0020066, 0.000001);

		if (mHTML) {
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
					fr.append(HTML.createPageEnd());
				}
				Desktop.getDesktop().browse(f.toURI());
			} catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Test
	public void testMassFraction() throws ArgumentException {
		LabeledMultivariateJacobianFunction.sDump = System.out;
		try {
		final Map<Element, Number> massFracs = new HashMap<>();
		final double total = 0.98;
		massFracs.put(Element.Calcium, new UncertainValue(0.398936 / total, "dCa", 0.012));
		massFracs.put(Element.Phosphorus, new UncertainValue(0.184987 / total, "dP", 0.008));
		massFracs.put(Element.Oxygen, new UncertainValue(0.41407 / total, "dO", 0.015));
		massFracs.put(Element.Hydrogen, new UncertainValue(0.0020066 / total, "dH", 0.001));
		final Composition mf = Composition.massFraction("Apatite", massFracs);

		Assert.assertEquals(5.0 / 22.0, mf.getAtomFraction(Element.Calcium).doubleValue(), 1.0e-6);
		Assert.assertEquals(3.0 / 22.0, mf.getAtomFraction(Element.Phosphorus).doubleValue(), 1.0e-6);
		Assert.assertEquals(13.0 / 22.0, mf.getAtomFraction(Element.Oxygen).doubleValue(), 1.0e-6);
		Assert.assertEquals(1.0 / 22.0, mf.getAtomFraction(Element.Hydrogen).doubleValue(), 1.0e-6);

		final Material mat = mf.getMaterial();
		final MaterialLabel.AtomFraction aft_ca = MaterialLabel.buildAtomFractionTag(mat, Element.Calcium);
		final MaterialLabel.AtomFraction aft_p = MaterialLabel.buildAtomFractionTag(mat, Element.Phosphorus);
		final MaterialLabel.AtomFraction aft_o = MaterialLabel.buildAtomFractionTag(mat, Element.Oxygen);
		final MaterialLabel.AtomFraction aft_h = MaterialLabel.buildAtomFractionTag(mat, Element.Hydrogen);
		Assert.assertEquals(0.0000767132, mf.getCovariance(aft_ca, aft_ca), 1.0e-7);
		Assert.assertEquals(0.0000430936, mf.getCovariance(aft_p, aft_p), 1.0e-7);
		Assert.assertEquals(0.000273062, mf.getCovariance(aft_o, aft_o), 1.0e-6);
		Assert.assertEquals(0.000450103, mf.getCovariance(aft_h, aft_h), 1.0e-6);
		Assert.assertEquals(-0.0000624525, mf.getCovariance(aft_h, aft_p), 1.0e-6);
		Assert.assertEquals(0.0000176265, mf.getCovariance(aft_ca, aft_p), 1.0e-6);
		Assert.assertEquals(-0.000284053, mf.getCovariance(aft_o, aft_h), 1.0e-6);

		if (mHTML) {
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
					final Table t = new Table();
					t.addRow(Table.th("Element"), Table.thc("Verbose"), Table.thc("Normal"));
					for (final Element elm : mf.getElementSet())
						t.addRow(Table.td(elm.toHTML(Mode.TERSE)),
								Table.tdc(mf.getAtomFraction(elm).toHTML(Mode.VERBOSE)),
								Table.tdc(mf.getAtomFraction(elm).toHTML(Mode.NORMAL)));
					fr.append(t.toHTML(Mode.NORMAL));
					fr.append(HTML.createPageEnd());
				}
				Desktop.getDesktop().browse(f.toURI());
			} catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		}
		finally {
			LabeledMultivariateJacobianFunction.sDump=null;
		}
		
	}

	@Test
	public void testAtomicFraction() throws ArgumentException {
		final Map<Element, Number> atFracs = new HashMap<>();
		atFracs.put(Element.Oxygen, new UncertainValue(0.593981, "dO", 0.05));
		atFracs.put(Element.Magnesium, new UncertainValue(0.106595, "dMg", 0.005));
		atFracs.put(Element.Aluminum, new UncertainValue(0.0404141, "dAl", 0.001));
		atFracs.put(Element.Silicon, new UncertainValue(0.167755, "dSi", 0.002));
		atFracs.put(Element.Calcium, new UncertainValue(0.0604423, "dCa", 0.0021));
		atFracs.put(Element.Iron, new UncertainValue(0.0308124, "dFe", 0.0012));
		final Composition af = Composition.atomFraction("K412", atFracs);

		Assert.assertEquals(0.431202, af.getMassFraction(Element.Oxygen).doubleValue(), 1.0e-5);
		Assert.assertEquals(0.117554, af.getMassFraction(Element.Magnesium).doubleValue(), 1.0e-5);
		Assert.assertEquals(0.049477, af.getMassFraction(Element.Aluminum).doubleValue(), 1.0e-5);
		Assert.assertEquals(0.21377, af.getMassFraction(Element.Silicon).doubleValue(), 1.0e-5);
		Assert.assertEquals(0.109914, af.getMassFraction(Element.Calcium).doubleValue(), 1.0e-5);
		Assert.assertEquals(0.0780755, af.getMassFraction(Element.Iron).doubleValue(), 1.0e-5);

		Assert.assertEquals(0.0209244, af.getMassFraction(Element.Oxygen).uncertainty(), 1.0e-5);
		Assert.assertEquals(0.00650561, af.getMassFraction(Element.Magnesium).uncertainty(), 1.0e-5);
		Assert.assertEquals(0.00217441, af.getMassFraction(Element.Aluminum).uncertainty(), 1.0e-5);
		Assert.assertEquals(0.00817155, af.getMassFraction(Element.Silicon).uncertainty(), 1.0e-5);
		Assert.assertEquals(0.00529588, af.getMassFraction(Element.Calcium).uncertainty(), 1.0e-5);
		Assert.assertEquals(0.00402649, af.getMassFraction(Element.Iron).uncertainty(), 1.0e-5);

		final Material mat = af.getMaterial();
		final MaterialLabel.MassFraction mft_o = MaterialLabel.buildMassFractionTag(mat, Element.Oxygen);
		final MaterialLabel.MassFraction mft_ca = MaterialLabel.buildMassFractionTag(mat, Element.Calcium);
		final MaterialLabel.MassFraction mft_mg = MaterialLabel.buildMassFractionTag(mat, Element.Magnesium);
		final MaterialLabel.MassFraction mft_fe = MaterialLabel.buildMassFractionTag(mat, Element.Iron);
		final MaterialLabel.MassFraction mft_si = MaterialLabel.buildMassFractionTag(mat, Element.Silicon);
		final MaterialLabel.MassFraction mft_al = MaterialLabel.buildMassFractionTag(mat, Element.Aluminum);
		final MaterialLabel.MassFraction mft_ag = MaterialLabel.buildMassFractionTag(mat, Element.Silver);
		final MaterialLabel.MassFraction mft_sb = MaterialLabel.buildMassFractionTag(mat, Element.Antimony);
		Assert.assertEquals(-0.0000857086, af.getCovariance(mft_o, mft_ca), 1.0e-7);
		Assert.assertEquals(9.2027e-6, af.getCovariance(mft_mg, mft_fe), 1.0e-7);
		Assert.assertEquals(9.68555e-6, af.getCovariance(mft_fe, mft_ca), 1.0e-7);
		Assert.assertEquals(0.0000274102, af.getCovariance(mft_si, mft_mg), 1.0e-7);
		Assert.assertEquals(0.0000139519, af.getCovariance(mft_al, mft_si), 1.0e-7);
		Assert.assertEquals(0.0, af.getCovariance(mft_o, mft_sb), 1.e-12);
		Assert.assertEquals(0.0, af.getCovariance(mft_sb, mft_ag), 1.e-12);

		if (mHTML) {
			final Report rep = new Report("Atomic Fraction");
			try {
				rep.addHeader("Atom Fraction");
				rep.addSubHeader("TERSE");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.TERSE));
				rep.addSubHeader("NORMAL");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.NORMAL));
				rep.addSubHeader("VERBOSE");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.VERBOSE));
				rep.addSubHeader("TERSE");
				rep.addHTML(af.toHTML(Mode.TERSE));
				rep.addSubHeader("NORMAL");
				rep.addHTML(af.toHTML(Mode.NORMAL));
				rep.addSubHeader("VERBOSE");
				rep.addHTML(af.toHTML(Mode.VERBOSE));
				final Table t = new Table();
				t.addRow(Table.th("Element"), Table.thc("Verbose"), Table.thc("Normal"));
				for (final Element elm : af.getElementSet())
					t.addRow(Table.td(elm.toHTML(Mode.TERSE)), Table.tdc(af.getMassFraction(elm).toHTML(Mode.VERBOSE)),
							Table.tdc(af.getMassFraction(elm).toHTML(Mode.NORMAL)));
				rep.addHTML(t.toHTML(Mode.NORMAL));
			} catch (final Exception e) {
				rep.addHTML(HTML.error(e.getMessage()));
			}
			try {
				rep.inBrowser(Mode.NORMAL);
			} catch (final IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Test
	public void testAtomicFractionMC() //
			throws ArgumentException {
		final Map<Element, UncertainValue> atFracs = new HashMap<>();
		atFracs.put(Element.Oxygen, new UncertainValue(0.593981, "dO", 0.05));
		atFracs.put(Element.Magnesium, new UncertainValue(0.106595, "dMg", 0.005));
		atFracs.put(Element.Aluminum, new UncertainValue(0.0404141, "dAl", 0.001));
		atFracs.put(Element.Silicon, new UncertainValue(0.167755, "dSi", 0.002));
		atFracs.put(Element.Calcium, new UncertainValue(0.0604423, "dCa", 0.0021));
		atFracs.put(Element.Iron, new UncertainValue(0.0308124, "dFe", 0.0012));
		final Composition af = Composition.atomFraction("K412", atFracs);

		// From DTSA-II
		assertEquals(0.4312, af.getMassFraction(Element.Oxygen).doubleValue(), 0.0001);
		assertEquals(0.1176, af.getMassFraction(Element.Magnesium).doubleValue(), 0.0001);
		assertEquals(0.0495, af.getMassFraction(Element.Aluminum).doubleValue(), 0.0001);
		assertEquals(0.2138, af.getMassFraction(Element.Silicon).doubleValue(), 0.0001);
		assertEquals(0.1099, af.getMassFraction(Element.Calcium).doubleValue(), 0.0001);
		assertEquals(0.0781, af.getMassFraction(Element.Iron).doubleValue(), 0.0001);

		// Test MC against analytic
		final UncertainValues mup = af.propagateMC(160000);
		for (final Object label : af.massFractionTags())
			assertEquals(af.getValue(label).fractionalUncertainty(), mup.getUncertainValue(label).fractionalUncertainty(), 0.001);

		if (mHTML) {
			final Report rep = new Report("Atomic Fraction MC");
			try {
				rep.addHeader("MC: Atom Fraction to Mass Fraction");
				rep.addHTML(af.toHTML(Mode.VERBOSE));
				rep.addSubHeader("TERSE");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.TERSE));
				rep.addSubHeader("NORMAL");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.NORMAL));
				rep.addSubHeader("VERBOSE - Analytic");
				rep.addHTML(af.toHTML(Representation.MassFraction, Mode.VERBOSE));
				rep.addSubHeader("VERBOSE - Monte Carlo");
				rep.addHTML(mup.toHTML(Mode.VERBOSE, Composition.MASS_FRACTION_FORMAT));
				
				final ValueToLog3 v2l = new ValueToLog3(1.0);
				final LinearToColor l2c = new LinearToColor(1.0, Color.blue, Color.red);
				rep.addImage(af.reorder(mup.getLabels()).asCovarianceBitmap(8, v2l, l2c), "Jacobian");
				rep.addImage(mup.asCovarianceBitmap(8, v2l, l2c), "Monte Carlo");
				
			}

			catch (final Exception e) {
				rep.addHTML(HTML.error(e.getMessage()));
			}
			try {
				rep.inBrowser(Mode.NORMAL);
			} catch (final IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Test
	public void testMixture() throws ArgumentException, ParseException, IOException {

		final String name = "K412";

		final double fSiO2 = 0.4535;
		final double fFeO = 0.0996;
		final double fMgO = 0.1933;
		final double fCaO = 0.1525;
		final double fAl2O3 = 0.0927;
		final double unc = 0.0020;
		final Map<Composition, Number> mcn = new HashMap<>();
		mcn.put(Composition.parse("SiO2"), new UncertainValue(fSiO2, unc));
		mcn.put(Composition.parse("FeO"), new UncertainValue(fFeO, unc));
		mcn.put(Composition.parse("MgO"), new UncertainValue(fMgO, unc));
		mcn.put(Composition.parse("CaO"), new UncertainValue(fCaO, unc));
		mcn.put(Composition.parse("Al2O3"), new UncertainValue(fAl2O3, unc));
		final Composition mix = Composition.combine(name, mcn, false);

		final UncertainValues jres = UncertainValues.force(mix);

		final RealVector dinp = mix.getValues().mapMultiply(0.0001);
		mix.setCalculator(new UncertainValuesCalculator.FiniteDifference(dinp));
		final UncertainValues dres = UncertainValues.force(mix);

		final UncertainValues mcres = mix.propagateMC(160000);

		for (final Object row : dres.getLabels()) {
			final double mm = mix.getVariance(row);
			final double dm = dres.getVariance(row);
			if (mix.getVariance(row) > 1.0e-10) {
				if (Math.abs(mm - dm) > 0.1 * Math.abs(mm))
					System.out.println(row + " " + mm + " != " + dm);
				for (final Object col : dres.getLabels())
					if (mix.getVariance(col) > 1.0e-10) {
						final double mc = mix.getCorrelationCoefficient(row, col);
						final double dc = dres.getCorrelationCoefficient(row, col);
						if (Math.abs(mc - dc) > 0.1 * Math.abs(mc))
							System.out.println("Row: " + row + "   Col: " + col);
					}
			}
		}

		final Map<Element, Number> res = new HashMap<Element, Number>();
		res.put(Element.Silicon, new UncertainValue(fSiO2 * 0.4674350, unc * 0.4674350));
		res.put(Element.Iron, new UncertainValue(fFeO * 0.777305, unc * 0.777305));
		res.put(Element.Magnesium, new UncertainValue(fMgO * 0.603036, unc * 0.603036));
		res.put(Element.Calcium, new UncertainValue(fCaO * 0.714691, unc * 0.714691));
		res.put(Element.Aluminum, new UncertainValue(fAl2O3 * 0.529251, unc * 0.529251));
		res.put(Element.Oxygen, new UncertainValue(fSiO2 * 0.532565 + // Si
				fFeO * 0.222695 + // Fe
				fMgO * 0.396964 + // Mg
				fCaO * 0.285309 + // Ca
				fAl2O3 * 0.470749, // Al
				unc * 0.532565 + // Si
						unc * 0.222695 + // Fe
						unc * 0.396964 + // Mg
						unc * 0.285309 + // Ca
						unc * 0.470749)); // Al
		final Composition uv = Composition.massFraction(name, res);

		final Report r = new Report("Mixture");
		r.addHeader("Jacobian");
		r.add(jres, Mode.NORMAL);
		r.addHeader("Delta");
		r.add(dres, Mode.NORMAL);
		r.addHeader("Monte Carlo");
		r.add(mcres, Mode.NORMAL);
		r.addHTML(mix.toSimpleHTML(new BasicNumberFormat("0.00E0")));
		r.addHTML(mix.getAnalyticalTotal().toHTML(Mode.TERSE));
		r.addHeader("Mass Fractions");
		r.add(uv);
		r.addHTML(uv.toSimpleHTML(new BasicNumberFormat("0.00E0")));
		r.addHTML(uv.getAnalyticalTotal().toHTML(Mode.TERSE));
		r.addHeader("Compare");
		final ValueToLog3 v2l = new ValueToLog3(1.0);
		final LinearToColor l2c = new LinearToColor(1.0, Color.blue, Color.red);
		r.addImage(dres.asCovarianceBitmap(8, v2l, l2c), "Delta");
		r.addImage(jres.asCovarianceBitmap(8, v2l, l2c), "Mixture");
		r.addImage(UncertainValues.compareAsBitmap(dres, jres, v2l, 8), "Mass Fractions");

		r.inBrowser(Mode.VERBOSE);
	}

	@Test
	public void runAnorthoclase() //
			throws ArgumentException, ParseException, IOException {
		final String name = "Anorthoclase";
		final Map<Composition, Number> mcn = new HashMap<>();
		final Composition sanidine = Composition.parse("K(AlSi3O8)");
		mcn.put(sanidine, new UncertainValue(0.40, 0.01));
		final Composition albite = Composition.parse("Na(AlSi3O8)");
		mcn.put(albite, new UncertainValue(0.60, 0.01));
		final Composition combiner = Composition.combine(name, mcn, true);
		
		final Report r = new Report("Anorthoclase");
		r.addHeader("Mixture");
		r.addSubHeader("Jacobian");
		r.add(albite);
		r.add(sanidine);
		r.add(combiner);
		UncertainValues mix = UncertainValues.force(combiner);
		combiner.setCalculator(new UncertainValuesCalculator.FiniteDifference(combiner.getInputValues().mapMultiply(0.001)));
		UncertainValues dmix = UncertainValues.force(combiner);
		r.addHTML(mix.toSimpleHTML(new BasicNumberFormat("0.00E0")));
		r.addSubHeader("Delta");
		r.add(dmix);
		r.addImage(
				UncertainValues.compareAsBitmap(UncertainValues.extract(dmix.getLabels(), mix), dmix,
						new LinearToColor(1.0, Color.blue, Color.red), 8),
				"Comparing analytical with finite difference.");

		r.inBrowser(Mode.VERBOSE);
	}

}
