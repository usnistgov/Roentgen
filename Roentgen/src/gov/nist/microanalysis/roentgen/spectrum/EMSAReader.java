package gov.nist.microanalysis.roentgen.spectrum;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.sql.Date;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>
 * A basic file reader for the EMSA/ISO standard spectrum file.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 312 $
 */
public class EMSAReader implements ISpectrumReader {

	private final NumberFormat mDefaultFormat = NumberFormat.getInstance(Locale.US);

	private double parseDouble(String value) {
		try {
			if (value.indexOf('+') != -1)
				value = value.replaceAll("\\+", "");
			if (value.indexOf('e') != -1)
				value = value.replaceAll("\\e", "E");
			return mDefaultFormat.parse(value).doubleValue();
		} catch (final ParseException e) {
			return Double.NaN;
		}
	}

	private int findMonth(final String month) {
		final String[] months = { "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec" };
		for (int i = 0; i < months.length; ++i)
			if (month.equalsIgnoreCase(months[i]))
				return i + 1;
		return 1;
	}

	/**
	 * @param str      = "hh:mm" or "hh:mm:ss"
	 * @param acquired
	 * @return Date
	 */
	private Date parseTime(final String str, final Date acquired) {
		int day, month, year, hour = 9, min = 0, sec = 0;
		final Calendar c = Calendar.getInstance();
		if (acquired != null)
			c.setTime(acquired);
		day = c.get(Calendar.DAY_OF_MONTH);
		month = c.get(Calendar.MONTH);
		year = c.get(Calendar.YEAR);
		try {
			final String[] items = str.split(":");
			if (items.length >= 2) {
				hour = Integer.parseInt(items[0].trim());
				min = Integer.parseInt(items[1].trim());
				if (items.length > 2)
					sec = Integer.parseInt(items[2].trim()); // seconds
			}
		} catch (final Exception e) {
		}
		c.set(year, month, day, hour, min, sec);
		return new Date(c.getTimeInMillis());
	}

	/**
	 * @param          value="dd-mmm-yyyy"
	 * @param acquired
	 * @return Date
	 */
	private Date parseDate(final String value, final Date acquired) {
		final Calendar c = Calendar.getInstance();
		// Default to today
		int day = c.get(Calendar.DAY_OF_MONTH);
		int month = c.get(Calendar.MONTH);
		int year = c.get(Calendar.YEAR);
		int hour = 0, min = 0, sec = 0;
		if (acquired != null) {
			// Get previously set time
			c.setTime(acquired);
			hour = c.get(Calendar.HOUR_OF_DAY);
			min = c.get(Calendar.MINUTE);
			sec = c.get(Calendar.SECOND);
		}
		try {
			// According to the EMSA standard (dd-mmm-yyyy)
			final String[] items = value.split("-");
			if (items.length < 3)
				throw new Exception("Misformatted date in EMSA file: " + value);
			day = Integer.parseInt(items[0].trim());
			month = findMonth(items[1].trim()) - 1;
			year = Integer.parseInt(items[2].trim());
			if (year < 100)
				year += 2000;
			if (year < 200)
				year += 1900;
		} catch (final Exception ex1) {
			// According to the locale
			try {
				final java.util.Date dt = DateFormat.getInstance().parse(value);
				c.setTimeInMillis(dt.getTime());
				day = c.get(Calendar.DAY_OF_MONTH);
				month = c.get(Calendar.MONTH);
				year = c.get(Calendar.YEAR);
			} catch (final Exception e) {
				System.err.println("Unable to parse date: " + value);
			}
		}
		c.set(year, month, day, hour, min, sec);
		return new Date(c.getTimeInMillis());
	}

	private String[] split(final String line) {
		final int p1 = line.indexOf(':');
		final String tag = p1 != -1 ? line.substring(0, p1).toUpperCase().trim() : "";
		final int p2 = tag.indexOf('-', 2);
		final String[] res = new String[p2 > 2 ? 3 : 2];
		res[0] = p2 > 2 ? tag.substring(0, p2).trim() : tag;
		if (p2 > 2)
			res[2] = tag.substring(p2 + 1).trim();
		res[1] = line.substring(p1 + 1).trim();
		return res;
	}

	private Pair<Element, Double> parsePair(final String elStr) {
		if (elStr.startsWith("(") && elStr.endsWith(")")) {
			final int p = elStr.indexOf(":");
			if (p > 0) {
				final Element elm = Element.parse(elStr.substring(1, p));
				final Double qty = Double.parseDouble(elStr.substring(p + 1, elStr.length() - 1));
				return new Pair<>(elm, qty);
			}
		}
		return null;
	}

	private Composition parseComposition(final String compStr) {
		// Name,(el1:qty1),(el2:qty2),...,density
		try {
			final String[] items = compStr.split(",");
			final Map<Element, Number> men = new TreeMap<>();
			@SuppressWarnings("unused")
			Double den = Double.NaN;
			for (int i = 1; i < items.length; ++i) {
				final Pair<Element, Double> pr = parsePair(items[i].trim());
				if (pr != null)
					men.put(pr.getFirst(), 0.01 * pr.getSecond());
				else if (i == (items.length - 1))
					den = Double.parseDouble(items[i].trim());
			}
			final Composition mf = Composition.massFraction(items[0], men);
			// if(!Double.isNaN(den)) mf.setDensity(den);
			return mf;
		} catch (final NumberFormatException e) {
			return null;
		}
	}

	@Override
	public EDSSpectrum read(final InputStream fis, final int i) throws IOException {
		if (i > 1)
			return null;
		final Reader rd = new InputStreamReader(fis, Charset.forName("US-ASCII"));
		// Number always use '.' as decimal separator
		try (final BufferedReader br = new BufferedReader(rd)) {
			String line = br.readLine().trim();
			// It seems that DTSA expects a blank line up front!?! and
			// LISPIX obliges
			if (line.length() == 0)
				line = br.readLine().trim();
			String[] items = split(line);
			if (items[0].equals("#FORMAT")) {
				if (!(items[1].equalsIgnoreCase("EMSA/MAS SPECTRAL DATA FILE")
						|| items[1].equalsIgnoreCase("EMSA/MAS SPECTRAL DATA STANDARD"))) {
					System.err.println("The format header in this EMSA file is spurious: " + items[1]);
					return null;
				}
			} else
				return null;
			items = split(br.readLine());
			if (items[0].equals("#VERSION")) {
				if ((parseDouble(items[1]) != 1.0) && (!items[1].equalsIgnoreCase("TC202v1.0")))
					System.err.println("The EMSA file version number was not 1.0. It was " + items[1]);
			} else
				return null;
			line = br.readLine().trim();
			items = split(line);
			final HashMap<String, String> props = new HashMap<>();
			while (!items[0].equals("#SPECTRUM")) {
				if (items[0].length() > 0)
					props.put(items[0], items[1]);
				line = br.readLine().trim();
				items = split(line);
			}
			final boolean isXY = "XY".equals(props.get("#DATATYPE"));
			int nCh;
			{
				String value = props.get("#NPOINTS");
				if (value.startsWith("+"))
					value = value.substring(1);
				nCh = Integer.parseInt(value);
			}
			final double zeroOffset = parseDouble(props.get("#OFFSET"));
			final double chWidth = parseDouble(props.get("#XPERCHAN"));
			double baseXUnit = 1.0;
			if (props.containsKey("#XUNITS") && props.get("#XUNITS").equalsIgnoreCase("kev"))
				baseXUnit = 1000.0;
			final double[] data = new double[nCh];
			line = br.readLine();
			if ((line != null) && (!line.startsWith("#ENDOFDATA"))) {
				int start = 0, dataCounter = 0;
				int xyCounter = 0;
				while (br.ready() && (dataCounter < data.length)) {
					String item;
					final int end = line.indexOf(',', start);
					if (end == -1) { // no ","
						item = line.substring(start).trim();
						start = line.length();
					} else {
						item = line.substring(start, end).trim();
						start = end + 1;
					}
					if (item.length() == 0) {
						line = br.readLine();
						start = 0;
						if ((line == null) || line.startsWith("#ENDOFDATA"))
							break;
						continue;
					}
					if ((!isXY) || ((xyCounter % 2) == 1)) {
						data[dataCounter] = parseDouble(item);
						dataCounter++;
					}
					++xyCounter;
				}
				if (dataCounter != data.length)
					System.err.println("The number of data points was fewer than the reported number of channels.");
			}
			final String title = props.containsKey("#TITLE") ? props.get("#TITLE")
					: "Spectrum " + Arrays.hashCode(data);
			final EDSSpectrum spec = new EDSSpectrum(title, chWidth * baseXUnit, zeroOffset * baseXUnit, data, true);
			Date acquired = parseDate(props.get("#DATE"), null);
			acquired = parseTime(props.get("#TIME"), acquired);
			if (acquired != null)
				spec.setAcquired(acquired);
			if (props.containsKey("#BEAMKV"))
				spec.setBeamEnergy(parseDouble(props.get("#BEAMKV")) * 1000.0);
			if (props.containsKey("#PROBECUR"))
				spec.setProbeCurrent(parseDouble(props.get("#PROBECUR")));
			final HashMap<String, Double> stage = new HashMap<>();
			if (props.containsKey("#XTILTSTGE"))
				stage.put(EDSSpectrum.X_TILT, parseDouble(props.get("#XTILTSTGE")));
			if (props.containsKey("#YTILTSTGE"))
				stage.put(EDSSpectrum.Y_TILT, parseDouble(props.get("#YTILTSTGE")));
			if (props.containsKey("#XPOSITION"))
				stage.put(EDSSpectrum.X_AXIS, parseDouble(props.get("#XPOSITION")));
			if (props.containsKey("#YPOSITION"))
				stage.put(EDSSpectrum.Y_AXIS, parseDouble(props.get("#YPOSITION")));
			if (props.containsKey("#ZPOSITION"))
				stage.put(EDSSpectrum.Z_AXIS, parseDouble(props.get("#ZPOSITION")));
			if (stage.size() > 0)
				spec.setPosition(stage);
			if (props.containsKey("#LIVETIME"))
				spec.setLiveTime(parseDouble(props.get("#LIVETIME")));
			if (props.containsKey("#REALTIME"))
				spec.setRealTime(parseDouble(props.get("#REALTIME")));
			// if(props.containsKey("##DET_HASH"))
			// spec.setOtherProperty("DET_HASH",
			// UUID.fromString(props.get("##DET_HASH")));
			if (props.containsKey("##D2STDCMP")) // DTSA-II custom tag
				spec.setComposition(parseComposition(props.get("##D2STDCMP")));
			// Cases: nextDatum then "," or nextDatum then EOL or EOL
			return spec;
		}
	}

	public boolean matches(final InputStream fis) {
		try (final InputStreamReader isr = new InputStreamReader(fis, Charset.forName("US-ASCII"))) {
			try (final BufferedReader br = new BufferedReader(isr)) {
				String line = br.readLine();
				if (line != null)
					line = line.trim();
				while ((line != null) && (line.length() == 0)) {
					line = br.readLine();
					if (line != null)
						line = line.trim();
				}

				final boolean res = (line == null) || line.substring(0, 7).equalsIgnoreCase("#FORMAT");
				return res;
			}

		} catch (final Exception e) {
			return false;
		}

	}

	@Override
	public int getSpectrumCount(final InputStream fis) {
		return matches(fis) ? 1 : 0;
	}

	public static class EMSAFileReaderFactory implements ISpectrumReaderFactory {

		/**
		 * @see com.duckandcover.roentgen.inputOutput.ISpectrumReaderFactory#get()
		 */
		@Override
		public ISpectrumReader get() {
			return new EMSAReader();
		}

		/**
		 * @see com.duckandcover.roentgen.inputOutput.ISpectrumReaderFactory#getFileFilter()
		 */
		@Override
		public FileFilter getFileFilter() {
			return new FileNameExtensionFilter("EMSA spectrum file", "msa", "emsa", "txt");
		}

	};

}