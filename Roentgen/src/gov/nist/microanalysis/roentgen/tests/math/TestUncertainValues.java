package gov.nist.microanalysis.roentgen.tests.math;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;

/**
 * @author Nicholas W. M. Ritchie
 *
 */
public class TestUncertainValues {

	@Test
	public void testUVs1() {
		final List<String> labels = Arrays.asList("B", "A", "C", "F", "E", "D", "G");
		final RealVector vals = new ArrayRealVector(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 });
		final RealVector vars = new ArrayRealVector(
				new double[] { 0.3 * 0.3, 12.0 * 12.0, 0.4 * 0.4, 1.1 * 1.1, 3.3 * 3.3, 0.2 * 0.2, 0.27 * 0.27 });
		final RealMatrix corrCoeff = MatrixUtils.createRealMatrix(new double[][] { //
				{ 0.0, 0.1, 0.02, 0.0, -0.2, 0.0, -0.8 }, //
				{ 0.0, 0.0, -0.3, 0.0, 0.3, 0.0, 0.6 }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3 }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33 }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, //
		});
		final double EPS = 1.0e-8;
		final UncertainValues uvs = new UncertainValues(labels, vals, vars, corrCoeff);
		assertEquals(0.1, uvs.getCorrelationCoefficient(0, 1), EPS);
		assertEquals(0.1, uvs.getCorrelationCoefficient(1, 0), EPS);
		assertEquals(0.33, uvs.getCorrelationCoefficient(6, 5), EPS);
		assertEquals(0.33, uvs.getCorrelationCoefficient(5, 6), EPS);

		assertEquals(2.0, uvs.getEntry("A"), EPS);
		assertEquals(2.0, uvs.getEntry(1), EPS);
		assertEquals(4.0, uvs.getEntry("F"), EPS);
		assertEquals(4.0, uvs.getEntry(3), EPS);
		assertEquals(7.0, uvs.getEntry("G"), EPS);
		assertEquals(7.0, uvs.getEntry(6), EPS);
		assertEquals(1.0, uvs.getEntry("B"), EPS);
		assertEquals(1.0, uvs.getEntry(0), EPS);
		assertEquals(12.0 * 12.0, uvs.getVariance("A"), EPS);
		assertEquals(12.0, uvs.getUncertainty("A"), EPS);
		assertEquals(3.3 * 3.3, uvs.getVariance("E"), EPS);
		assertEquals(3.3, uvs.getUncertainty("E"), EPS);
		final double covCF = uvs.getCovariance("C", "F");
		final double covCG = uvs.getCovariance("C", "G");
		final double covBG = uvs.getCovariance("B", "G");
		assertEquals(0.4 * 1.0 * 0.0, covCF, EPS);
		assertEquals(0.4 * 0.27 * -0.3, covCG, EPS);
		assertEquals(0.3 * 0.27 * -0.8, covBG, EPS);

		final UncertainValues uvss = uvs.sort();

		assertEquals(2.0, uvss.getEntry("A"), EPS);
		assertEquals(4.0, uvss.getEntry("F"), EPS);
		assertEquals(7.0, uvss.getEntry("G"), EPS);
		assertEquals(1.0, uvss.getEntry("B"), EPS);
		assertEquals(12.0 * 12.0, uvss.getVariance("A"), EPS);
		assertEquals(12.0, uvss.getUncertainty("A"), EPS);
		assertEquals(3.3 * 3.3, uvss.getVariance("E"), EPS);
		assertEquals(3.3, uvss.getUncertainty("E"), EPS);
		assertEquals(covCF, uvss.getCovariance("C", "F"), EPS);
		assertEquals(covCG, uvss.getCovariance("C", "G"), EPS);
		assertEquals(covBG, uvss.getCovariance("B", "G"), EPS);

		assertEquals(0.33, uvss.getCorrelationCoefficient("D", "G"), EPS);
		assertEquals(0.33, uvs.getCorrelationCoefficient("D", "G"), EPS);
		assertEquals(-0.3, uvss.getCorrelationCoefficient("A", "C"), EPS);
		assertEquals(-0.3, uvs.getCorrelationCoefficient("A", "C"), EPS);
		assertEquals(-0.3, uvss.getCorrelationCoefficient("C", "A"), EPS);
		assertEquals(-0.3, uvs.getCorrelationCoefficient("C", "A"), EPS);

		final UncertainValues uvsd = uvs.blockDiagnonalize();
		assertEquals(2.0, uvsd.getEntry("A"), EPS);
		assertEquals(4.0, uvsd.getEntry("F"), EPS);
		assertEquals(7.0, uvsd.getEntry("G"), EPS);
		assertEquals(1.0, uvsd.getEntry("B"), EPS);
		assertEquals(12.0 * 12.0, uvsd.getVariance("A"), EPS);
		assertEquals(12.0, uvsd.getUncertainty("A"), EPS);
		assertEquals(3.3 * 3.3, uvsd.getVariance("E"), EPS);
		assertEquals(3.3, uvsd.getUncertainty("E"), EPS);
		assertEquals(covCF, uvsd.getCovariance("C", "F"), EPS);
		assertEquals(covCG, uvsd.getCovariance("C", "G"), EPS);
		assertEquals(covBG, uvsd.getCovariance("B", "G"), EPS);

	}

	@Test
	public void testUVs2() {
		final List<String> labels = Arrays.asList("B", "A", "C", "F", "E", "D", "G");
		final RealVector vals = new ArrayRealVector(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 });
		final RealVector vars = new ArrayRealVector(
				new double[] { 0.3 * 0.3, 12.0 * 12.0, 0.4 * 0.4, 1.1 * 1.1, 3.3 * 3.3, 0.2 * 0.2, 0.27 * 0.27 });
		final RealMatrix corrCoeff = MatrixUtils.createRealMatrix(new double[][] { //
				{ 0.0, 0.1, 0.02, 0.0, 0.0, 0.0, 0.0 }, //
				{ 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0 }, //
				{ 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0 }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, -0.11, 0.0, }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33 }, //
				{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, //
		});
		final double EPS = 1.0e-8;
		final UncertainValues uvs = new UncertainValues(labels, vals, vars, corrCoeff);
		final UncertainValues uvss = uvs.sort();

		final UncertainValues uvsd = uvss.blockDiagnonalize();
		assertEquals(uvs.getEntry("A"), uvsd.getEntry("A"), EPS);
		assertEquals(uvs.getEntry("F"), uvsd.getEntry("F"), EPS);
		assertEquals(uvs.getEntry("G"), uvsd.getEntry("G"), EPS);
		assertEquals(uvs.getEntry("B"), uvsd.getEntry("B"), EPS);
		assertEquals(uvs.getVariance("A"), uvsd.getVariance("A"), EPS);
		assertEquals(uvs.getUncertainty("A"), uvsd.getUncertainty("A"), EPS);
		assertEquals(uvs.getVariance("E"), uvsd.getVariance("E"), EPS);
		assertEquals(uvs.getUncertainty("E"), uvsd.getUncertainty("E"), EPS);

		System.out.println(uvsd.getLabels());
		System.out.println(uvsd.getCovariances());

	}

}
