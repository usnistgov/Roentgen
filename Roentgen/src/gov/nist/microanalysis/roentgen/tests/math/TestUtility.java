package gov.nist.microanalysis.roentgen.tests.math;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import gov.nist.microanalysis.roentgen.math.Histogram;
import gov.nist.microanalysis.roentgen.math.Utility;
import junit.framework.TestCase;

public class TestUtility extends TestCase {

	public void testExpRand() {
		final SummaryStatistics ss = new SummaryStatistics();
		final Histogram h = new Histogram(0.0, 5.0, 50);
		for (int i = 0; i < 10000000; ++i) {
			final double val = Utility.expRand();
			ss.addValue(val);
			h.add(val);
		}
		// System.out.println("Mean=" + ss.getMean());
		// System.out.println("Variance=" + ss.getVariance());
		assertEquals(1.0, ss.getMean(), 0.01);
		assertEquals(1.0, ss.getVariance(), 0.01);
		double cdf = 0.0;
		final double k = 1.0 / h.totalCounts();
		for (int i = 0; i < h.binCount(); ++i) {
			cdf += k * h.counts(i);
			// System.out.println(cdf + " <->" + Double.toString(1.0 -
			// Math.exp(-h.maxValue(i))));
			assertEquals(cdf, 1.0 - Math.exp(-h.maxValue(i)), 0.001);
		}
	}

	public void testQuadratic() {
		// (2x+4)(x-5)=0
		// 2x^2-6x-20=0
		// x-> -2, 5
		final double[] qs1 = Utility.quadraticSolver(2.0, -6.0, -20.0);
		assertEquals(-2.0, qs1[0], 1.0e-10);
		assertEquals(5.0, qs1[1], 1.0e-10);
		final double[] qs2 = Utility.quadraticSolver(2.0, -6.0, 20.0);
		assertNull(qs2);
		// (x-2)(x-2) = (x^2-4x+4)
		// x-> 2, 2
		final double[] qs3 = Utility.quadraticSolver(1.0, -4.0, 4.0);
		assertEquals(2.0, qs3[0], 1.0e-10);
		assertEquals(2.0, qs3[1], 1.0e-10);
	}
}
