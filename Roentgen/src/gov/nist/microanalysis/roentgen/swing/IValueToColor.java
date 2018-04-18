package gov.nist.microanalysis.roentgen.swing;

import java.awt.Color;

/**
 * Maps a numeric value to a color. Useful for plotting values as colors.
 *
 * @author Nicholas
 */
public interface IValueToColor {

   public Color map(double value);

}
