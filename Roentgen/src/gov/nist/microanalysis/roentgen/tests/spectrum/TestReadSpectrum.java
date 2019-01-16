package gov.nist.microanalysis.roentgen.tests.spectrum;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.spectrum.EDSSpectrum;
import gov.nist.microanalysis.roentgen.spectrum.SpectrumFileReader;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestReadSpectrum {

	private static File makeFile(final String resName) throws IOException {
		try (final InputStream is = TestReadSpectrum.class.getResourceAsStream(resName)) {
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

	public static EDSSpectrum fromResource(final String resName) throws Exception {
		final SpectrumFileReader sfr = Globals.getSpectrumFileReader();
		final File f = makeFile(resName);
		return sfr.read(f, 0);
	}

}
