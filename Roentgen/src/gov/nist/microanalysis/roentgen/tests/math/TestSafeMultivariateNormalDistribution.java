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

	private final boolean REPORT = false;

	@Test
	public void test1() throws IOException {
		List<String> tags = new ArrayList<>();
		for (int i = 0; i < 5; ++i)
			tags.add("Tag" + i);
		RealVector vals = MatrixUtils.createRealVector(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 });
		RealMatrix cov = MatrixUtils.createRealMatrix(vals.getDimension(), vals.getDimension());
		cov.setEntry(0, 0, 0.5);
		cov.setEntry(1, 1, 0.0);
		cov.setEntry(2, 2, 1.0);
		cov.setEntry(3, 3, 2.0);
		cov.setEntry(4, 4, 3.0);

		cov.setEntry(2, 3, 0.3);
		cov.setEntry(3, 2, 0.3);
		cov.setEntry(4, 2, 0.8);
		cov.setEntry(2, 4, 0.8);
		cov.setEntry(3, 4, 1.2);
		cov.setEntry(4, 3, 1.2);
		UncertainValues uvs = new UncertainValues(tags, vals, cov);
		SafeMultivariateNormalDistribution smnd = new SafeMultivariateNormalDistribution(vals, cov);
		EstimateUncertainValues euv = new EstimateUncertainValues(tags);
		euv.add(smnd, 1600000);
		UncertainValues euvs = euv.estimateDistribution();
		for(int r=0;r<5;++r) {
			Assert.assertEquals(uvs.getEntry(r), euvs.getEntry(r),0.01);
			for (int c=r;c<5;++c)
				Assert.assertEquals(uvs.getCovariance(r, c), euvs.getCovariance(r, c),0.05);
		}
		if (REPORT) {
			Report r = new Report("Test1 - SafeMultivariateNormalDistribution");
			BasicNumberFormat nf = new BasicNumberFormat("0.00");
			r.addHTML(uvs.toHTML(Mode.NORMAL, nf));
			r.addHTML(euv.estimateDistribution().toHTML(Mode.NORMAL, nf));
			r.inBrowser(Mode.NORMAL);
		}
	}

}
