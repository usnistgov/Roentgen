package gov.nist.microanalysis.roentgen.tests.math;

import org.junit.Assert;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValue;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;
import junit.framework.TestCase;

/**
 * <p>
 * These test UncertainValue against UncertainValue and against the output of
 * the program
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 199 $
 */
public class TestUncertainValue extends TestCase {

	private final UncertainValue mA = new UncertainValue(1.24, "A", 0.3);
	private final UncertainValue mA2a = makeA2a();

	private final UncertainValue mB = new UncertainValue(8.82, "B", 1.2);
	private final UncertainValue mB2a = makeB2a();

	private final UncertainValue mC = new UncertainValue(-9.3, "C", 2.1);
	private final UncertainValue mC2a = makeC2a();

	private static final UncertainValue makeA2a() {
		final UncertainValue res = new UncertainValue(1.24);
		res.assignComponent("V1", Math.sqrt(0.05));
		res.assignComponent("V2", Math.sqrt(0.04));
		return res;
	}

	private static final UncertainValue makeB2a() {
		final UncertainValue res = new UncertainValue(8.82);
		res.assignComponent("V1", Math.sqrt(0.84));
		res.assignComponent("V2", Math.sqrt(0.6));
		return res;
	}

	private static final UncertainValue makeC2a() {
		final UncertainValue res = new UncertainValue(-9.3);
		// 2.1 * 2.1 = 4.2 + 0.21 = 4.41
		res.assignComponent("V1", Math.sqrt(3.0));
		res.assignComponent("V3", Math.sqrt(1.41));
		return res;
	}

	private static void assertEquals(final Number n1, final Number n2, final double delta) {
		Assert.assertEquals(n1.doubleValue(), n2.doubleValue(), delta);
		if (n1.getClass().equals(n2.getClass())) {
			if (n1 instanceof UncertainValue)
				Assert.assertEquals(((UncertainValue) n1).uncertainty(), ((UncertainValue) n2).uncertainty(), delta);
		} else {
			final double unc1 = n1 instanceof UncertainValue ? ((UncertainValue) n1).doubleValue() : 0.0;
			final double unc2 = n2 instanceof UncertainValue ? ((UncertainValue) n2).doubleValue() : 0.0;
			Assert.assertEquals(unc1, unc2, delta);
		}
	}

	public void testA() {
		assertEquals(mA, mA2a, 1.0e-10);
		assertEquals(mA.uncertainty(), mA2a.uncertainty(), 1.0e-8);
		assertEquals(mA, mA2a, 1.0e-8);
		assertEquals(mA2a.format(new BasicNumberFormat("0.000")), "1.240\u00B10.300");
		assertEquals(mA2a.formatLong(new BasicNumberFormat("0.000")), "1.240\u00B10.224(V1)\u00B10.200(V2)");
	}

	public void testB() {
		assertEquals(mB, mB2a, 1.0e-10);
		assertEquals(mB.uncertainty(), mB2a.uncertainty(), 1.0e-8);
		assertEquals(mB, mB2a, 1.0e-8);
		assertEquals(mB2a.format(new BasicNumberFormat("0.000")), "8.820\u00B11.200");
		assertEquals(mB2a.formatLong(new BasicNumberFormat("0.000")), "8.820\u00B10.917(V1)\u00B10.775(V2)");
	}

	public void testC() {
		assertEquals(mC, mC2a, 1.0e-10);
		assertEquals(mC.uncertainty(), mC2a.uncertainty(), 1.0e-8);
		assertEquals(mC, mC2a, 1.0e-8);
		assertEquals(mC2a.format(new BasicNumberFormat("0.000")), "-9.300\u00B12.100");
		assertEquals(mC2a.formatLong(new BasicNumberFormat("0.000")), "-9.300\u00B11.732(V1)\u00B11.187(V3)");
	}
};