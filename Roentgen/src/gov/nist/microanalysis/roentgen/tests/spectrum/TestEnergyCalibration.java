package gov.nist.microanalysis.roentgen.tests.spectrum;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import gov.nist.microanalysis.roentgen.spectrum.EnergyCalibration;
import gov.nist.microanalysis.roentgen.spectrum.EnergyCalibration.Polynomial;

/**
 * <p>
 * Description
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 199 $
 */
public class TestEnergyCalibration {

	@Test
	public void testLinear() {
		final Polynomial l = EnergyCalibration.Linear(0.0, 10.0, 2048);
		assertEquals(l.averageEnergyForChannel(0), 5.0, 1.0e-6);
		assertEquals(l.averageEnergyForChannel(1000), 10005.0, 1.0e-6);
		assertEquals(l.channelIndex(10005.0), 1000);
		assertEquals(l.channelIndex(0.0), 0);
		assertEquals(l.channelIndex(10.0), 1);
		assertEquals(l.channelWidth(100), 10.0, 1.0e-6);
		assertEquals(l.maxEnergyForChannel(100), 1010.0, 1.e-6);
		assertEquals(l.minEnergyForChannel(100), 1000.0, 1.e-6);
		final RealVector er = l.range(0, 2048);
		assertEquals(er.getDimension(), 2048);
		for (int i = 0; i < 2048; i++)
			assertEquals(er.getEntry(i), i * 10.0, 1.0e6);
		assertEquals(l.equals(EnergyCalibration.Linear(0.0, 10.0, 2048)), true);
		assertEquals(l.equals(EnergyCalibration.Linear(0.0, 10.0, 1024)), false);
		assertEquals(l.equals(EnergyCalibration.Linear(1.0, 10.0, 2048)), false);
		assertEquals(l.equals(EnergyCalibration.Linear(0.0, 9.900, 2048)), false);
		assertEquals(l.equals(null), false);
		assertEquals(l.equals(EnergyCalibration.Quadratic(0.0, 10.0, 0.0, 2048)), false);
	}

	@Test
	public void testQuadratic() {
		final Polynomial q = EnergyCalibration.Quadratic(-100.0, 10.0, 1.0e-6, 2048);
		assertEquals(q.averageEnergyForChannel(0), -94.9999999, 1.0e-5);
		assertEquals(q.averageEnergyForChannel(1000), 9906.0010005, 1.0e-4);
		assertEquals(q.channelIndex(10005.0), 1010);
		assertEquals(q.channelIndex(0.0), 9);
		assertEquals(q.channelIndex(10.0), 10);
		assertEquals(q.channelWidth(100), 10.000201, 1.0e-5);
		assertEquals(q.maxEnergyForChannel(100), 910.010201, 1.e-6);
		assertEquals(q.minEnergyForChannel(100), 900.01, 1.e-6);
		final RealVector er = q.range(0, 2048);
		assertEquals(er.getDimension(), 2048);
		for (int i = 0; i < 2048; i++)
			assertEquals(er.getEntry(i), -100.0 + (i * (10.0 + (i * 1.0e-6))), 1.0e-6);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.0, 1.0e-6 }, 2048)), true);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-101.0, 10.0, 1.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-100.0, 10.01, 1.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-100.0, 10.0, 2.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Linear(-100.0, 10.0, 2048)), false);
		assertEquals(q.equals(null), false);
	}

	@Test
	public void testQuadratic2() {
		final Polynomial q = EnergyCalibration.Quadratic(-100.0, 10.0, -1.0e-6, 2048);
		assertEquals(q.averageEnergyForChannel(0), -94.9999999, 1.0e-5);
		assertEquals(q.averageEnergyForChannel(1000), 9903.9989995, 1.0e-4);
		assertEquals(q.channelIndex(10005.0), 1010);
		assertEquals(q.channelIndex(0.0), 10);
		assertEquals(q.channelIndex(10.0), 11);
		assertEquals(q.channelWidth(100), 9.9997989999, 1.0e-5);
		assertEquals(q.maxEnergyForChannel(100), 909.989799, 1.e-6);
		assertEquals(q.minEnergyForChannel(100), 899.99, 1.e-6);
		final RealVector er = q.range(0, 2048);
		assertEquals(er.getDimension(), 2048);
		for (int i = 0; i < 2048; i++)
			assertEquals(er.getEntry(i), -100.0 + (i * (10.0 + (i * -1.0e-6))), 1.0e-6);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.0, -1.0e-6 }, 2048)), true);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-101.0, 10.0, -1.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-100.0, 10.01, -1.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Quadratic(-100.0, 10.0, -2.0e-6, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Linear(-100.0, 10.0, 2048)), false);
		assertEquals(q.equals(null), false);
	}

	@Test
	public void testPolynomial() {
		final Polynomial q = new Polynomial(new double[] { -100.0, 10.0, 1.0e-6, -1.0e-9 }, 2048);
		assertEquals(q.averageEnergyForChannel(0), -94.9999999, 1.0e-5);
		assertEquals(q.averageEnergyForChannel(1000), 9904.9995, 1.0e-4);
		assertEquals(q.channelIndex(10005.0), 1010);
		assertEquals(q.channelIndex(0.0), 9);
		assertEquals(q.channelIndex(10.0), 10);
		assertEquals(q.channelWidth(100), 10.0001706, 1.0e-5);
		assertEquals(q.maxEnergyForChannel(100), 910.0091707, 1.e-6);
		assertEquals(q.minEnergyForChannel(100), 900.009, 1.e-6);
		assertEquals(q.invert(5005.0), 510.487, 0.001);
		final RealVector er = q.range(0, 2048);
		assertEquals(er.getDimension(), 2048);
		for (int i = 0; i < 2048; i++)
			assertEquals(er.getEntry(i), -100.0 + (i * (10.0 + (i * (1.0e-6 - (i * 1.0e-9))))), 1.0e-3);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.0, 1.0e-6, -1.0e-9 }, 2048)), true);
		assertEquals(q.equals(new Polynomial(new double[] { -101.0, 10.0, 1.0e-6 - 1.0e-9 }, 2048)), false);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.01, 1.0e-6 - 1.0e-9 }, 2048)), false);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.0, 2.0e-6 - 1.0e-9 }, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Linear(-100.0, 10.0, 2048)), false);
		assertEquals(q.equals(null), false);
	}

	public double e(final double a, final double b, final double c, final double ch) {
		return (a * ch) + (Math.sqrt(ch) * b) + c;
	}

	@Test
	public void testQuadraticInSqrt() {
		final EnergyCalibration.QuadraticInSqrt q = new EnergyCalibration.QuadraticInSqrt(-100.0, 0.01, 10.0, 2048);
		assertEquals(q.averageEnergyForChannel(0), 0.5 * (e(10.0, 0.01, -100.0, 0) + e(10.0, 0.01, -100.0, 1.0)),
				1.0e-5);
		assertEquals(q.averageEnergyForChannel(1000),
				0.5 * (e(10.0, 0.01, -100.0, 1000) + e(10.0, 0.01, -100.0, 1001.0)), 1.0e-4);
		assertEquals(q.channelIndex(10005.0), 1010);
		assertEquals(q.channelIndex(0.0), 9);
		assertEquals(q.channelIndex(10.0), 10);
		assertEquals(q.channelWidth(100), 910.100498756 - 900.100, 1.0e-5);
		assertEquals(q.maxEnergyForChannel(100), 910.100498756, 1.e-6);
		assertEquals(q.minEnergyForChannel(100), 900.1, 1.e-6);
		assertEquals(q.invert(5005.0), 510.4774, 0.001);
		final RealVector er = q.range(0, 2048);
		assertEquals(er.getDimension(), 2048);
		for (int i = 0; i < 2048; i++)
			assertEquals(er.getEntry(i), -100.0 + (Math.sqrt(i) * 0.01) + (10.0 * i), 1.0e-3);
		assertEquals(q.equals(new EnergyCalibration.QuadraticInSqrt(-100.0, 0.01, 10.0, 2048)), true);
		assertEquals(q.equals(new EnergyCalibration.QuadraticInSqrt(-100.0, 0.02, 10.0, 2048)), false);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.01, 1.0e-6 - 1.0e-9 }, 2048)), false);
		assertEquals(q.equals(new Polynomial(new double[] { -100.0, 10.0, 2.0e-6 - 1.0e-9 }, 2048)), false);
		assertEquals(q.equals(EnergyCalibration.Linear(-100.0, 10.0, 2048)), false);
		assertEquals(q.equals(null), false);
	}
}
