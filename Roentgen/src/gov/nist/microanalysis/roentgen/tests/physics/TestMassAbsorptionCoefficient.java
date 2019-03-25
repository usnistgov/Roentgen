package gov.nist.microanalysis.roentgen.tests.physics;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Test;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.juncertainty.UncertainValue;
import gov.nist.juncertainty.UncertainValues;
import gov.nist.juncertainty.UncertainValuesBase;
import gov.nist.juncertainty.UncertainValuesCalculator;
import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.EPMALabel;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction;
import gov.nist.microanalysis.roentgen.physics.XRay;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;
import gov.nist.microanalysis.roentgen.swing.ValueToLog3;

/**
 * <p>
 * Tests for the {@link ElementalMAC} class
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 244 $
 */
public class TestMassAbsorptionCoefficient {

	private static final boolean REPORT = true;

	@Test
	public void testCompute() throws ParseException, ArgumentException {
		final Composition fe = Composition.parse("Fe");
		testMac(fe, Element.Oxygen, XRayTransition.KL3, 3625.375742042211, 208.35034389516583);
		testMac(Element.Iron, Element.Oxygen, XRayTransition.KL3, 3625.375742042211, 208.35034389516583);

		final Composition sio2 = Composition.parse("SiO2");
		testMac(sio2, Element.Uranium, XRayTransition.M5N7, 483.741, 4.0156);
		testMac(sio2, Element.Uranium, XRayTransition.M3N5, 349.892, 2.91998);

		final Composition u3o8 = Composition.parse("U3O8");
		testMac(u3o8, Element.Oxygen, XRayTransition.KA1, 6528.85, 365.556);
		testMac(u3o8, Element.Uranium, XRayTransition.LA1, 63.5378, 0.632347);
	}

	@Test
	public void testIsAvailable() {
		final ElementalMAC alg = new ElementalMAC();
		for (final Element elm : Element.range(Element.Hydrogen, Element.Uranium)) {
			assertEquals(alg.isAvailable(elm, new XRay(100.0)), true);
			assertEquals(alg.isAvailable(elm, new XRay(100.0e3)), true);
			assertEquals(alg.isAvailable(elm, new XRay(5.0e5)), false);
		}
	}

	@Test
	public void testK412() throws ParseException, IOException, ArgumentException {

		final Composition k412 = Composition.combine("K412", false, //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.4541, 0.0077)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.0994, 0.0018)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1966, 0.0025)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1544, 0.0015)), //
				Pair.create(Composition.parse("Al2O3"), new UncertainValue(0.0934, 0.0029)));
		final Composition mf1 = Composition.parse("MgO");
		final Composition mf2 = Composition.parse("Al2O3");
		final Composition mf3 = Composition.parse("SiO2");
		final Composition mf4 = Composition.parse("CaF2");
		final Composition mf5 = Composition.parse("Fe");
		final Composition mf6 = Composition.combine("K411", false, //
				Pair.create(Composition.parse("SiO2"), new UncertainValue(0.5389, 0.0096)), //
				Pair.create(Composition.parse("FeO"), new UncertainValue(0.1448, 0.0027)), //
				Pair.create(Composition.parse("MgO"), new UncertainValue(0.1512, 0.0020)), //
				Pair.create(Composition.parse("CaO"), new UncertainValue(0.1549, 0.0015)));
		final Composition[] mfs = new Composition[] { k412, mf1, mf2, mf3, mf4, mf5, mf6 };
		final CharacteristicXRay[] cxrs = { CharacteristicXRay.create(Element.Oxygen, XRayTransition.KA1),
				CharacteristicXRay.create(Element.Magnesium, XRayTransition.KA1),
				CharacteristicXRay.create(Element.Aluminum, XRayTransition.KA1),
				CharacteristicXRay.create(Element.Silicon, XRayTransition.KA1),
				CharacteristicXRay.create(Element.Calcium, XRayTransition.KA1),
				CharacteristicXRay.create(Element.Iron, XRayTransition.KA1) };
		final Report r = new Report("Test K412");
		try {
			for (final CharacteristicXRay cxr : cxrs) {
				final List<Composition> comps = Arrays.asList(mfs);
				final MaterialMACFunction mmf = MaterialMACFunction.build(comps, cxr);
				final UncertainValuesBase<EPMALabel> inps = mmf.buildInputs(comps, cxr);
				final UncertainValuesCalculator<EPMALabel> uv = UncertainValuesBase.propagateAnalytical(mmf, inps);
				// final UncertainValuesBase mc =
				// UncertainValues.propagateDeltaOrdered(pr.getSecond(), pr.getFirst(),
				// pr.getFirst().getValues().mapMultiply(0.001));
				r.addHeader(cxr);
				r.addSubHeader("Taylor");
				r.add(uv);
				r.addHTML(uv.toHTML(Mode.NORMAL));
				// r.addSubHeader("Delta");
				// r.add(mc);
				final ValueToLog3 V2L3 = new ValueToLog3(1.0);
				r.addImage(uv.asCovarianceBitmap(20, V2L3, V2L3), "Taylor");
				// r.addImage(mc.asCovarianceBitmap(8, V2L3, L2C), "Delta");
			}
		} finally {
			r.inBrowser(Mode.VERBOSE);
		}

	}

	public void testMac(final Composition comp, final Element elm, final XRayTransition tr, final double mac,
			final double uMac) {
		final XRay xr = CharacteristicXRay.find(CharacteristicXRay.forElement(elm), tr);
		assertEquals(xr != null, true);
		// assertEquals(UncertainValue.mean(alg.compute(comp, xr)), mac, 0.1);
		// assertEquals(UncertainValue.uncertainty(alg.compute(comp, xr)), uMac,
		// 0.1);
	}

	public void testMac(final Element mat, final Element elm, final XRayTransition tr, final double mac,
			final double uMac) {
		final ElementalMAC alg = new ElementalMAC();
		final XRay xr = CharacteristicXRay.find(CharacteristicXRay.forElement(elm), tr);
		assertEquals(xr != null, true);
		assertEquals(UncertainValue.mean(alg.compute(mat, xr)), mac, 0.1);
		assertEquals(UncertainValue.uncertainty(alg.compute(mat, xr)), uMac, 0.1);
	}

	@Test
	public void testMultiMaterial() throws ParseException, IOException, ArgumentException {
		final Composition mf1 = Composition.parse("Fe2O3");
		final Composition mf2 = Composition.parse("FeO2");
		final Composition mf3 = Composition.parse("FeO");
		final Composition mf4 = Composition.parse("Fe3Al2(SiO4)3");
		final Composition mf5 = Composition.parse("Al2O3");
		final Composition mf6 = Composition.parse("CaF2");
		final Composition[] mfs = new Composition[] { mf1, mf2, mf3, mf4, mf5, mf6 };
		final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Iron, XRayTransition.LA1);
		final UncertainValuesCalculator<EPMALabel> uvc = MaterialMACFunction.compute(Arrays.asList(mfs), cxr);

		final UncertainValues<EPMALabel> uv = UncertainValues.asUncertainValues(uvc);
		uvc.setCalculator(uvc.new MonteCarlo(16000));
		final UncertainValues<EPMALabel> mc = UncertainValues.asUncertainValues(uvc);

		if (REPORT) {
			final Report r = new Report("TestMultiMAC");
			for (final Composition mf : mfs)
				r.add(mf);
			r.addHeader("Taylor");
			r.add(uv);
			r.addHeader("Monte-Carlo");
			r.add(mc);
			r.inBrowser(Mode.VERBOSE);
		}
		for (int p = 0; p < uv.getDimension(); ++p) {
			Assert.assertEquals(uv.getEntry(p), mc.getEntry(p), 0.02 * mc.getEntry(p));
			for (int q = p + 1; q < uv.getDimension(); ++q) {
				// System.out.println(uv.getTag(p) + ", " + uv.getTag(q) + ","
				// + Math.abs(uv.getCovariance(p, q) - mc.getCovariance(p, q)) + ","
				// + 0.03 * Math.sqrt(uv.getVariance(p) * uv.getVariance(q)));
				Assert.assertEquals(uv.getCovariance(p, q), mc.getCovariance(p, q),
						0.05 * Math.sqrt(uv.getVariance(p) * uv.getVariance(q) + 1.0e8));
			}
		}

	}

}
