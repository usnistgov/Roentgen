package gov.nist.microanalysis.roentgen.tests.math;

import org.junit.Assert;

import gov.nist.juncertainty.UncertainValue;
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

	private final UncertainValue mA = new UncertainValue(1.24, 0.3);
	private final UncertainValue mA2 = new UncertainValue(mA);

	private final UncertainValue mB = new UncertainValue(8.82, 0.0);
	private final UncertainValue mB2 = new UncertainValue(8.82);

	private final UncertainValue mC = new UncertainValue(1.24, 2.1);

	public void testA() {
		Assert.assertTrue(mA.equals(mA2));
		Assert.assertTrue(mB.equals(mB2));
		Assert.assertEquals(mA.uncertainty(), 0.3, 1.0e-10);
		Assert.assertEquals(mB.uncertainty(), 0.0, 1.0e-10);
		Assert.assertFalse(mA.equals(mB));
		Assert.assertEquals(mA.fractionalUncertainty(), 0.3 / 1.24, 1.0e-10);
		Assert.assertEquals(mB.fractionalUncertainty(), 0.0, 1.0e-10);

		Assert.assertEquals(mA.compareTo(mB), -1);
		Assert.assertEquals(mB.compareTo(mA), 1);
		Assert.assertEquals(mA.compareTo(mC), -1);
		Assert.assertEquals(mC.compareTo(mA), 1);
		Assert.assertTrue(mA.equals(mA2));

		Assert.assertEquals(mA.uncertainty(), mA2.uncertainty(), 1.0e-8);
		
		Assert.assertEquals(Math.pow(mA.uncertainty(),2.0), mA2.variance(), 1.0e-10);

		Assert.assertTrue(mA.isUncertain());
		Assert.assertFalse(mB.isUncertain());
		Assert.assertFalse(UncertainValue.isUncertain(Double.valueOf(100.0)));
		Assert.assertTrue(UncertainValue.isUncertain(mA));
		Assert.assertFalse(UncertainValue.isUncertain(mB));
		Assert.assertFalse(UncertainValue.isSpecialNumber(Double.valueOf(100.0)));
		Assert.assertTrue(UncertainValue.isSpecialNumber(mA));
		
		Assert.assertEquals(3.3*1.24, mA.multiply(3.3).doubleValue(), 1.0e-10);
		Assert.assertEquals(2.3*0.3, mA.multiply(2.3).uncertainty(), 1.0e-10);
		
		Assert.assertEquals(mA2.format(new BasicNumberFormat("0.000")), "1.240\u00B10.300");
		Assert.assertEquals(mA2.formatLong(new BasicNumberFormat("0.000")), "1.240\u00B10.300");
		Assert.assertEquals(mB2.formatLong(new BasicNumberFormat("0.000")), "8.820\u00B10.000");
		Assert.assertEquals(mC.format(new BasicNumberFormat("0.00")), "1.24\u00B12.10");
	}

};