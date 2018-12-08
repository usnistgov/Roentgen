package gov.nist.microanalysis.roentgen.spectrum;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.swing.filechooser.FileFilter;

/**
 * <p>
 * Attempts to determine the format of the spectrum file and use the appropriate
 * algorithm to read the spectrum out of the file.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 312 $
 */
public class SpectrumFileReader {

	private final ArrayList<ISpectrumReaderFactory> mReaderFactories = new ArrayList<>();

	public SpectrumFileReader() {
	}

	public void register(final ISpectrumReaderFactory reader) {
		for (final ISpectrumReaderFactory fact : mReaderFactories)
			if (fact.getClass().equals(reader.getClass()))
				return;
		mReaderFactories.add(reader);
	}

	public List<FileFilter> getSpectrumFileFilters() {
		final ArrayList<FileFilter> res = new ArrayList<>();
		for (final ISpectrumReaderFactory reader : mReaderFactories)
			res.add(reader.getFileFilter());
		return res;
	}

	public EDSSpectrum read(final InputStream fis, final int i) throws Exception {
		final File tmp = File.createTempFile("tmp", "sfr");
		try (FileOutputStream fos = new FileOutputStream(tmp)) {
			final byte[] buffer = new byte[64 * 1024];
			for (int len = fis.read(buffer); len > 0; len = fis.read(buffer))
				fos.write(buffer, 0, len);
		}
		final List<EDSSpectrum> res = readInternal(tmp);
		tmp.delete();
		return res.get(i);
	}

	private List<EDSSpectrum> readInternal(final File file) throws Exception {
		final ArrayList<EDSSpectrum> res = new ArrayList<>();
		for (final ISpectrumReaderFactory fact : mReaderFactories) {
			final ISpectrumReader reader = fact.get();
			int count;
			try (FileInputStream fis = new FileInputStream(file)) {
				count = reader.getSpectrumCount(fis);
			}
			for (int i = 0; i < count; ++i)
				try (FileInputStream fis = new FileInputStream(file)) {
					final EDSSpectrum spec = reader.read(fis, i);
					// assert spec.getDetector() != null;
					res.add(spec);
				}
		}
		return res;
	}

	public List<EDSSpectrum> read(final File f) throws Exception {
		return readInternal(f);
	}

	public EDSSpectrum read(final File f, final int i) throws Exception {
		try (FileInputStream fis = new FileInputStream(f)) {
			final EDSSpectrum res = read(fis, i);
			// res.setSource(f, i);
			return res;
		}
	}

	public void unregister(final ISpectrumReaderFactory reader) {
		for (final ISpectrumReaderFactory fact : mReaderFactories)
			if (fact.getClass().equals(reader.getClass())) {
				mReaderFactories.remove(fact);
				break;
			}
	}
}