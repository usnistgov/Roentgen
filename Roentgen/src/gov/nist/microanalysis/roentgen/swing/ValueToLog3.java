package gov.nist.microanalysis.roentgen.swing;

import java.awt.Color;

public class ValueToLog3 implements IValueToColor {

	final double mMax;

	public ValueToLog3(final double max) {
		mMax = max;
	}

	@Override
	public Color map(final double value) {
		final int i = Math.min(Math.max(0, (int) (-255.0 * Math.log10(value / mMax) / 3.0)), 255);
		if (i <= 86)
			return new Color((200 * i) / 86, (200 * i) / 86, 60 + ((195 * i) / 86));
		if (i <= 170) {
			final int j = i - 87;
			return new Color((200 * j) / (170 - 87), (200 * j) / (170 - 87), 60 + ((195 * j) / (170 - 87)));
		}
		if (i <= 254) {
			final int j = i - 171;
			new Color(60 + ((195 * j) / (254 - 171)), (200 * j) / (254 - 171), (200 * j) / (254 - 171));
		}
		return new Color(255, 255, 255);
	}
}
