package gov.nist.microanalysis.roentgen.swing;

import java.awt.Color;

public class LinearToColor implements IValueToColor {

	final int[] mPosColor;
	final int[] mNegColor;
	final double mScale;

	public LinearToColor(final double scale, final Color positiveColor, final Color negativeColor) {
		mPosColor = new int[] { positiveColor.getRed(), positiveColor.getGreen(), positiveColor.getBlue() };
		mNegColor = new int[] { negativeColor.getRed(), negativeColor.getGreen(), negativeColor.getBlue() };
		mScale = scale;
	}

	@Override
	public Color map(final double value) {
		final double v = Math.min(1.0, Math.max(-1.0, value / mScale));
		final int[] rgb = v > 0.0 ? mPosColor : mNegColor;
		return new Color(rgb[0], rgb[1], rgb[2], (int) (255 * Math.abs(v)));
	}

}
