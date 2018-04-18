package com.duckandcover.swing;

import java.awt.Color;
import java.awt.SystemColor;
import java.awt.image.IndexColorModel;

import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;

/**
 * <p>
 * Some functions to make for a more consistent user interface.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class SwingUtils {

   private final static Color sLineColor = createLineColor();

   private static final int balance(final int i1, final int i2) {
      return (i1 + 3 * i2) / 4;
   }

   private static final Color createLineColor() {
      final Color c1 = SystemColor.controlShadow;
      final Color c2 = SystemColor.control;
      return new Color(balance(c1.getRed(), c2.getRed()), balance(c1.getGreen(), c2.getGreen()), balance(c1.getBlue(), c2.getBlue()), balance(c1.getAlpha(), c2.getAlpha()));
   }

   public static TitledBorder createTitledBorder(final String name) {
      return BorderFactory.createTitledBorder(createDefaultBorder(), name);
   }

   public static Border createDefaultBorder() {
      return BorderFactory.createMatteBorder(1, 1, 1, 1, sLineColor);
   }

   public static Border createEmptyBorder() {
      return BorderFactory.createEmptyBorder(5, 5, 5, 5);
   }

   public static IndexColorModel createLog3BandColorModel() {
      final byte[] red = new byte[256];
      final byte[] green = new byte[256];
      final byte[] blue = new byte[256];
      for(int i = 1; i <= 86; ++i) {
         red[i] = green[i] = (byte) ((200 * i) / 86);
         blue[i] = (byte) (60 + ((195 * i) / 86));
      }
      for(int i = 0; i <= (170 - 87); ++i) {
         red[i + 87] = blue[i + 87] = (byte) ((200 * i) / (170 - 87));
         green[i + 87] = (byte) (60 + ((195 * i) / (170 - 87)));
      }
      for(int i = 0; i <= (254 - 171); ++i) {
         red[i + 171] = (byte) (60 + ((195 * i) / (254 - 171)));
         green[i + 171] = blue[i + 171] = (byte) ((200 * i) / (254 - 171));
      }
      red[255] = green[255] = (byte) 255;
      blue[255] = 0;
      return new IndexColorModel(8, 256, red, green, blue);
   }
}
