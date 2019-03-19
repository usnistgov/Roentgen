package gov.nist.microanalysis.roentgen.tests.math;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

import com.duckandcover.html.IToHTML.Mode;
import com.duckandcover.html.Report;

import gov.nist.microanalysis.roentgen.math.EstimateUncertainValues;
import gov.nist.microanalysis.roentgen.math.SafeMultivariateNormalDistribution;
import gov.nist.microanalysis.roentgen.math.uncertainty.UncertainValues;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

public class TestSafeMultivariateNormalDistribution {

	private final boolean REPORT = true;

	@Test
	public void test1() throws IOException {
		final RealVector vals = MatrixUtils.createRealVector(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 });
		final int tagCount = vals.getDimension();
		final List<String> tags = new ArrayList<>();
		for (int i = 0; i < tagCount; ++i)
			tags.add("Tag" + i);
		final RealMatrix cov = MatrixUtils.createRealMatrix(tagCount, tagCount);
		cov.setEntry(0, 0, 0.5);
		cov.setEntry(1, 1, 0.0);
		cov.setEntry(2, 2, 1.0);
		cov.setEntry(3, 3, 2.0);
		cov.setEntry(4, 4, 3.0);

		cov.setEntry(2, 3, 0.3);
		cov.setEntry(3, 2, 0.3);
		cov.setEntry(4, 2, 0.8);
		cov.setEntry(2, 4, 0.8);
		cov.setEntry(3, 4, -1.2);
		cov.setEntry(4, 3, -1.2);
		final UncertainValues<String> uvs = new UncertainValues<String>(tags, vals, cov);
		final SafeMultivariateNormalDistribution smnd = new SafeMultivariateNormalDistribution(vals, cov);
		final EstimateUncertainValues<String> euvs = new EstimateUncertainValues<String>(tags);
		euvs.add(smnd, 1600000);
		if (REPORT) {
			final Report r = new Report("Test1 - SafeMultivariateNormalDistribution");
			final BasicNumberFormat nf = new BasicNumberFormat("0.00");
			r.addHTML(uvs.toHTML(Mode.NORMAL, nf));
			r.addHTML(euvs.toHTML(Mode.NORMAL, nf));
			r.inBrowser(Mode.NORMAL);
		}
		for (int r = 0; r < tagCount; ++r) {
			Assert.assertEquals(uvs.getEntry(r), euvs.getEntry(r), 0.01);
			for (int c = r; c < tagCount; ++c)
				Assert.assertEquals(uvs.getCovariance(r, c), euvs.getCovariance(r, c), 0.05);
		}
	}

	@Test
	public void test2() throws IOException {
		final RealVector vals = MatrixUtils.createRealVector(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 });
		final int tagCount = vals.getDimension();
		final List<String> tags = new ArrayList<>();
		for (int i = 0; i < tagCount; ++i)
			tags.add("Tag" + i);
		final RealMatrix cov = MatrixUtils.createRealMatrix(tagCount, tagCount);
		cov.setEntry(0, 0, 0.5);
		cov.setEntry(1, 1, 0.0);
		cov.setEntry(2, 2, 1.0);
		cov.setEntry(3, 3, 2.0);
		cov.setEntry(4, 4, 3.0);
		cov.setEntry(5, 5, 1.0);
		cov.setEntry(6, 6, 1.0);

		cov.setEntry(2, 3, 0.3);
		cov.setEntry(3, 2, 0.3);
		cov.setEntry(4, 2, 0.8);
		cov.setEntry(2, 4, 0.8);
		cov.setEntry(3, 4, -1.2);
		cov.setEntry(4, 3, -1.2);

		cov.setEntry(5, 6, -0.3);
		cov.setEntry(6, 5, -0.3);
		final UncertainValues<String> uvs = new UncertainValues<String>(tags, vals, cov);
		final SafeMultivariateNormalDistribution smnd = new SafeMultivariateNormalDistribution(vals, cov);
		final EstimateUncertainValues<String> euvs = new EstimateUncertainValues<String>(tags);
		euvs.add(smnd, 1600000);
		if (REPORT) {
			final Report r = new Report("Test1 - SafeMultivariateNormalDistribution");
			final BasicNumberFormat nf = new BasicNumberFormat("0.00");
			r.addHTML(uvs.toHTML(Mode.NORMAL, nf));
			r.addHTML(euvs.toHTML(Mode.NORMAL, nf));
			r.inBrowser(Mode.NORMAL);
		}
		for (int r = 0; r < tagCount; ++r) {
			Assert.assertEquals(uvs.getEntry(r), euvs.getEntry(r), 0.01);
			for (int c = r; c < tagCount; ++c)
				Assert.assertEquals(uvs.getCovariance(r, c), euvs.getCovariance(r, c), 0.05);
		}
	}

}
