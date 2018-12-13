package gov.nist.microanalysis.roentgen.tests.spectrum;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.spectrum.AdaptiveGaussianFilter;
import gov.nist.microanalysis.roentgen.spectrum.AdaptiveTophatFilter;
import gov.nist.microanalysis.roentgen.spectrum.EDSSpectrum;
import gov.nist.microanalysis.roentgen.spectrum.LineshapeCalibration;
import gov.nist.microanalysis.roentgen.spectrum.SpectrumFileReader;

/**
 * <p>
 * Tests the SpectrumFileReader class
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 312 $
 */
public class TestSpectrumFileReader {

	private File makeFile(final String resName) throws IOException {
		try (final InputStream is = getClass().getResourceAsStream(resName)) {
			final File f = File.createTempFile("spec", ".msa");
			try (final FileOutputStream fos = new FileOutputStream(f)) {
				final byte[] tmp = new byte[4096];
				while (is.available() > 0) {
					final int cx = is.read(tmp);
					if (cx > 0)
						fos.write(tmp, 0, cx);
				}
			}
			return f;
		}
	}

	@Test
	public void testEMSAFileReader() throws Exception {
		final SpectrumFileReader sfr = Globals.getSpectrumFileReader();
		final File f = makeFile("Al_ref1.msa");
		try {
			final EDSSpectrum spec = sfr.read(f, 0);
			assertEquals(spec.size(), 2048);
			assertEquals(spec.getBeamEnergy(), 20.0e3, 1.0);
			assertEquals(spec.getEnergyCalibration().channelWidth(0), 10.0, 0.001);
			assertEquals(spec.getEnergyCalibration().channelWidth(2047), 10.0, 0.001);
			assertEquals(spec.getEnergyCalibration().averageEnergyForChannel(0), 5.0, 0.001);
			assertEquals(spec.getEnergyCalibration().averageEnergyForChannel(1024), 10245.0, 0.01);
			assertEquals(!Double.isNaN(spec.getLiveTime()), true);
			assertEquals(spec.getLiveTime(), 60.0, 1.0e-5);
			final RealVector data = spec.getData();
			assertEquals(data.getEntry(63), 80.58932308980361, 0.01);
			assertEquals(data.getEntry((232 * 4) - 1), 7.854691297950559, 0.01);
			if (false) {
				final LineshapeCalibration lsc = new LineshapeCalibration.Gaussian(130.0,
						LineshapeCalibration.Gaussian.SDD_EV_PER_EH);
				final AdaptiveTophatFilter atf = new AdaptiveTophatFilter(spec.size(), spec.getEnergyCalibration(),
						lsc);
				final Pair<RealVector, RealMatrix> res = atf.evaluate(spec.getData());
				final RealVector vals = res.getFirst();
				final AdaptiveGaussianFilter agf = new AdaptiveGaussianFilter(spec.size(), spec.getEnergyCalibration(),
						lsc);
				final Pair<RealVector, RealMatrix> res2 = agf.evaluate(spec.getData());
				final RealVector vals2 = res2.getFirst();
				for (int i = 0; i < vals.getDimension(); ++i) {
					final StringBuffer sb = new StringBuffer();
					sb.append(i);
					sb.append("\t");
					sb.append(data.getEntry(i));
					sb.append("\t");
					sb.append(vals.getEntry(i));
					sb.append("\t");
					sb.append(vals2.getEntry(i));
					System.out.println(sb.toString());
				}
			}
		} finally {
			f.delete();
		}
	}

	@Test
	public void testEMSAFileReader2() throws Exception {
		final SpectrumFileReader sfr = Globals.getSpectrumFileReader();
		final File f = makeFile("MgO_ref1.msa");
		try {
			final EDSSpectrum spec = sfr.read(f, 0);
			assertEquals(spec.size(), 2048);
			assertEquals(spec.getBeamEnergy(), 20.0e3, 1.0);
			assertEquals(spec.getEnergyCalibration().channelWidth(0), 10.0, 0.001);
			assertEquals(spec.getEnergyCalibration().channelWidth(2047), 10.0, 0.001);
			assertEquals(spec.getEnergyCalibration().averageEnergyForChannel(0), 5.0, 0.001);
			assertEquals(spec.getEnergyCalibration().averageEnergyForChannel(1024), 10245.0, 0.01);
			assertEquals(!Double.isNaN(spec.getLiveTime()), true);
			assertEquals(spec.getLiveTime(), 60.0, 1.0e-5);
			final RealVector data = spec.getData();
			// assertEquals(data.getEntry(63), 80.58932308980361, 0.01);
			// assertEquals(data.getEntry((232 * 4) - 1), 7.854691297950559, 0.01);

			if (true) {
				final LineshapeCalibration lsc = new LineshapeCalibration.Gaussian(130.0,
						LineshapeCalibration.Gaussian.SDD_EV_PER_EH);
				final AdaptiveTophatFilter atf = new AdaptiveTophatFilter(spec.size(), spec.getEnergyCalibration(),
						lsc);
				final Pair<RealVector, RealMatrix> res = atf.evaluate(spec.getData());
				final RealVector vals = res.getFirst();
				final AdaptiveGaussianFilter agf = new AdaptiveGaussianFilter(spec.size(), spec.getEnergyCalibration(),
						lsc);
				final Pair<RealVector, RealMatrix> res2 = agf.evaluate(spec.getData());
				final RealVector vals2 = res2.getFirst();

				final double sc = vals.getMaxValue() / vals2.getMaxValue();
				for (int i = 0; i < vals.getDimension(); ++i) {
					final StringBuffer sb = new StringBuffer();
					sb.append(i);
					sb.append("\t");
					sb.append(data.getEntry(i));
					sb.append("\t");
					sb.append(vals.getEntry(i));
					sb.append("\t");
					sb.append(vals2.getEntry(i) * sc);
					System.out.println(sb.toString());
				}
			}

		} finally {
			f.delete();
		}
	}

}
