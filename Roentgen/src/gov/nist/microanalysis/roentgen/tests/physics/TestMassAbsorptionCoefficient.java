package gov.nist.microanalysis.roentgen.tests.physics;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Test;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.ArgumentException;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.physics.CharacteristicXRay;
import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.ElementalMAC;
import gov.nist.microanalysis.roentgen.physics.MaterialMACFunction;
import gov.nist.microanalysis.roentgen.physics.XRay;
import gov.nist.microanalysis.roentgen.physics.XRayTransition;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

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
	public void testCompute() throws ParseException {
		final Composition fe = Composition.parse("Fe").asMassFraction();
		testMac(fe, Element.Oxygen, XRayTransition.KL3, 3625.375742042211, 208.35034389516583);
		testMac(Element.Iron, Element.Oxygen, XRayTransition.KL3, 3625.375742042211, 208.35034389516583);

		final Composition sio2 = Composition.parse("SiO2").asMassFraction();
		testMac(sio2, Element.Uranium, XRayTransition.M5N7, 483.741, 4.0156);
		testMac(sio2, Element.Uranium, XRayTransition.M3N5, 349.892, 2.91998);

		final Composition u3o8 = Composition.parse("U3O8").asMassFraction();
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
	public void testMultiMaterial() throws ParseException, IOException, ArgumentException {
		final Composition mf1 = Composition.parse("Fe2O3").asMassFraction();
		final Composition mf2 = Composition.parse("FeO2").asMassFraction();
		final Composition mf3 = Composition.parse("FeO").asMassFraction();
		final Composition mf4 = Composition.parse("Fe3Al2(SiO4)3").asMassFraction();
		final Composition mf5 = Composition.parse("Al2O3").asMassFraction();
		final Composition mf6 = Composition.parse("CaF2").asMassFraction();
		final Composition[] mfs = new Composition[] { mf1, mf2, mf3, mf4, mf5, mf6 };
		final CharacteristicXRay cxr = CharacteristicXRay.create(Element.Iron, XRayTransition.LA1);
		final Pair<UncertainValues, MaterialMACFunction> pr = MaterialMACFunction.buildCompute(Arrays.asList(mfs), cxr);
		final UncertainValues uv = UncertainValues.propagate(pr.getSecond(), pr.getFirst());
		final RealVector minCov = new ArrayRealVector(pr.getFirst().getDimension());
		minCov.set(1.0);
		final UncertainValues mc = UncertainValues.propagateMC(pr.getSecond(),
				UncertainValues.forceMinCovariance(pr.getFirst(), minCov), 16000);
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
						0.05 * Math.sqrt(uv.getVariance(p) * uv.getVariance(q)));
			}
		}

	}

}
