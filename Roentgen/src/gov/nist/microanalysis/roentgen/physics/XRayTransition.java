package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.Arrays;

import com.google.common.base.Preconditions;

/**
 * <p>
 * This abstraction represents the two Shell objects that participate in an
 * x-ray emission event.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 228 $
 */
public enum XRayTransition implements Comparable<XRayTransition> {
   KL3(Shell.K, Shell.L3, "K&alpha;<sub>1</sub>"), // Ka1
   KL2(Shell.K, Shell.L2, "K&alpha;<sub>2</sub>"), // Ka2
   KL1(Shell.K, Shell.L1), // For Li, Be
   KM3(Shell.K, Shell.M3, "K&beta;<sub>1</sub>"), // KB1,
   KN3(Shell.K, Shell.N3, "K<super>i</super>&beta;<sub>2</sub>"), // KiB2,
   KN2(Shell.K, Shell.N2, "K<super>ii</super>&beta;<sub>2</sub>"), // KiiB2
   KM2(Shell.K, Shell.M2, "K&beta;<sub>3</sub>"), // KB3
   KN5(Shell.K, Shell.N5, "K<super>i</super>&beta;<sub>4</sub>"), // KiB4
   KN4(Shell.K, Shell.N4, "K<super>ii</super>&beta;<sub>4</sub>"), // KiiB4
   KM5(Shell.K, Shell.M5, "K<super>i</super>&beta;<sub>5</sub>"), // KiB5
   KM4(Shell.K, Shell.M4, "K<super>ii</super>&beta;<sub>5</sub>"), // KiiB5
   KO2(Shell.K, Shell.O2), //
   KO3(Shell.K, Shell.O3), //
   KO4(Shell.K, Shell.O4), //
   KO5(Shell.K, Shell.O5), //
   KP2(Shell.K, Shell.P2), //
   KP3(Shell.K, Shell.P3), //
   L3N2(Shell.L3, Shell.N2), // L3-N2,
   L3N3(Shell.L3, Shell.N3), // L3-N3,
   L3O2(Shell.L3, Shell.O2), // L3-O2,
   L3O3(Shell.L3, Shell.O3), // L3-O3,
   L3P1(Shell.L3, Shell.P1), // L3-P1,
   L3M5(Shell.L3, Shell.M5, "L&alpha;<sub>1</sub>"), // LA1-M5,
   L3M4(Shell.L3, Shell.M4, "L&alpha;<sub>2</sub>"), // LA2-M4,
   L3N4(Shell.L3, Shell.N4, "L&beta;<sub>15</sub>"), // LB15-N4,
   L3N5(Shell.L3, Shell.N5, "L&beta;<sub>2</sub>"), // LB2-N5,
   L3O4(Shell.L3, Shell.O4, "L&beta;<sub>5</sub>"), // LB5-O4_O5,
   L3O5(Shell.L3, Shell.O5), //
   L3N1(Shell.L3, Shell.N1, "L&beta;<sub>6</sub>"), // LB6-N1,
   L3O1(Shell.L3, Shell.O1, "L&beta;<sub>7</sub>"), // LB7-O1,
   L3M1(Shell.L3, Shell.M1, "L<italic>l</italic>"), // Ll-M1,
   L3M3(Shell.L3, Shell.M3, "L<italic>s</italic>"), // Ls-M3,
   L3M2(Shell.L3, Shell.M2, "L<italic>t</italic>"), // Lt-M2,
   L3N6(Shell.L3, Shell.N6, "L<italic>u</italic>"), // Lu-N4_N7,
   L3N7(Shell.L3, Shell.N7), //
   L3P3(Shell.L3, Shell.P3), //
   L2M2(Shell.L2, Shell.M2), // L2-M2,
   L2M5(Shell.L2, Shell.M5), // L2-M5,
   L2N2(Shell.L2, Shell.N2), // L2-N2,
   L2N3(Shell.L2, Shell.N3), // L2-N3,
   L2N5(Shell.L2, Shell.N5), // L2-N5,
   L2O2(Shell.L2, Shell.O2), // L2-O2,
   L2O3(Shell.L2, Shell.O3), // L2-O3,
   L2P1(Shell.L2, Shell.P1), // L2-P1,
   L2P2(Shell.L2, Shell.P2), // L2-P2,
   L2M4(Shell.L2, Shell.M4, "L&beta;<sub>1</sub>"), // LB1-M4,
   L2M3(Shell.L2, Shell.M3, "L&beta;<sub>17</sub>"), // LB17-M3,
   L2N4(Shell.L2, Shell.N4, "L&gamma;<sub>1</sub>"), // LG1-N4,
   L2N1(Shell.L2, Shell.N1, "L&gamma;<sub>5</sub>"), // LG5-N1,
   L2O4(Shell.L2, Shell.O4, "L&gamma;<sub>6</sub>"), // LG6-O4,
   L2O1(Shell.L2, Shell.O1, "L&gamma;<sub>8</sub>"), // LG8-O1,
   L2M1(Shell.L2, Shell.M1, "L&eta;"), // Ln,
   L2N6(Shell.L2, Shell.N6, "L&nu;"), // Lv-N6,
   L2N7(Shell.L2, Shell.N7), //
   L2P3(Shell.L2, Shell.P3), //
   L3P2(Shell.L3, Shell.P2), //
   L3Q1(Shell.L3, Shell.Q1), //
   L1M1(Shell.L1, Shell.M1), // L1-M1,
   L1N1(Shell.L1, Shell.N1), // L1-N1,
   L1N4(Shell.L1, Shell.N4), // L1-N4,
   L1O1(Shell.L1, Shell.O1), // L1-O1,
   L1O4(Shell.L1, Shell.O4), // L1-O4,
   L1M4(Shell.L1, Shell.M4, "L&beta;<sub>10</sub>"), // LB10-M4,
   L1M3(Shell.L1, Shell.M3, "L&beta;<sub>3</sub>"), // LB3-M3,
   L1M2(Shell.L1, Shell.M2, "L&beta;<sub>4</sub>"), // LB4-M2,
   L1M5(Shell.L1, Shell.M5, "L&beta;<sub>9</sub>"), // LB9-M5,
   L1N2(Shell.L1, Shell.N2, "L&gamma;<sub>2</sub>"), // LG2-N2,
   L1N5(Shell.L1, Shell.N5, "L&gamma;<sub>11</sub>"), // LG11-N5,
   L1N3(Shell.L1, Shell.N3, "L&gamma;<sub>3</sub>"), // LG3-N3,
   L1O5(Shell.L1, Shell.O5), //
   L1O3(Shell.L1, Shell.O3, "L&gamma;<sub>4</sub>"), // LG4-O3,
   L1O2(Shell.L1, Shell.O2, "L&gamma;<sub>4'</sub>"), // LG4`-O2,
   L1P2(Shell.L1, Shell.P2), //
   L1P3(Shell.L1, Shell.P3), //
   L2Q1(Shell.L2, Shell.Q1), //
   L2O6(Shell.L2, Shell.O6), //
   L3O6(Shell.L3, Shell.O6), //
   L3O7(Shell.L3, Shell.O7), //
   M1N2(Shell.M1, Shell.N2), // M1-N2,
   M1N3(Shell.M1, Shell.N3), // M1-N3,
   M1O2(Shell.M1, Shell.O2), // M1-O2,
   M1O3(Shell.M1, Shell.O3), // M1-O3,
   M1P2(Shell.M1, Shell.P2), //
   M1P3(Shell.M1, Shell.P3), //
   M2M4(Shell.M2, Shell.M4), // M2-M4,
   M2N1(Shell.M2, Shell.N1), // M2-N1,
   M2N4(Shell.M2, Shell.N4), // M2-N4,
   M2O1(Shell.M2, Shell.O1), // M2-O1
   M2O4(Shell.M2, Shell.O4), // M2-O4,
   M2P1(Shell.M2, Shell.P1), // M2-P1
   M2Q1(Shell.M2, Shell.Q1), // M2-Q1
   M3M4(Shell.M3, Shell.M4), // M3-M4,
   M3M5(Shell.M3, Shell.M5), // M3-M5,
   M3N1(Shell.M3, Shell.N1), // M3-N1,
   M3N4(Shell.M3, Shell.N4), // M3-N4,
   M3O1(Shell.M3, Shell.O1), // M3-O1,
   M3O4(Shell.M3, Shell.O4), // M3-O4,
   M3O5(Shell.M3, Shell.O5), // M3-O5,
   M3N5(Shell.M3, Shell.N5, "M&gamma;"), // MG-N5,
   M3P1(Shell.M3, Shell.P1), // M3P1
   M3Q1(Shell.M3, Shell.Q1), //
   M4N3(Shell.M4, Shell.N3), // M4-N3,
   M4O2(Shell.M4, Shell.O2), // M4-O2,
   M4O3(Shell.M4, Shell.O3), // M4-O3,
   M4O6(Shell.M4, Shell.O6), //
   M4N6(Shell.M4, Shell.N6), // MB-N6,
   M4N2(Shell.M4, Shell.N2, "M&zeta;<sub>2</sub>"), // MZ2-N2,
   M4P2(Shell.M4, Shell.P2), //
   M4P3(Shell.M4, Shell.P3), //
   M5O3(Shell.M5, Shell.O3), // M5-O3,
   M5P3(Shell.M5, Shell.P3), //
   M5N7(Shell.M5, Shell.N7, "M&alpha;<sub>1</sub>"), // MA1-N7,
   M5N6(Shell.M5, Shell.N6, "M&alpha;<sub>2</sub>"), // MA2-N6,
   M5N3(Shell.M5, Shell.N3, "M&zeta;<sub>1</sub>"), // MZ1-N3,
   M5O6(Shell.M5, Shell.O6), //
   M5O7(Shell.M5, Shell.O7), //
   N4N6(Shell.N4, Shell.N6), // N4-N6,
   N5N6(Shell.N5, Shell.N6), // N5-N6
   N1O2(Shell.N1, Shell.O2), //
   N1O3(Shell.N1, Shell.O3), //
   N1P2(Shell.N1, Shell.P2), //
   N1P3(Shell.N1, Shell.P3), //
   N2O1(Shell.N2, Shell.O1), //
   N2O4(Shell.N2, Shell.O4), //
   N2P1(Shell.N2, Shell.P1), //
   N2Q1(Shell.N2, Shell.Q1), //
   N3O1(Shell.N3, Shell.O1), //
   N3O4(Shell.N3, Shell.O4), //
   N3O5(Shell.N3, Shell.O5), //
   N3P1(Shell.N3, Shell.P1), //
   N3Q1(Shell.N3, Shell.Q1), //
   N4O2(Shell.N4, Shell.O2), //
   N4O3(Shell.N4, Shell.O3), //
   N4P2(Shell.N4, Shell.P2), //
   N4P3(Shell.N4, Shell.P3), //
   N4O6(Shell.N4, Shell.O6), //
   N5O3(Shell.N5, Shell.O3), //
   N5P3(Shell.N5, Shell.P3), //
   N5O6(Shell.N5, Shell.O6), //
   N5O7(Shell.N5, Shell.O7), //
   N6O4(Shell.N6, Shell.O4), //
   N6O5(Shell.N6, Shell.O5), //
   N7O5(Shell.N7, Shell.O5), //
   O1P2(Shell.O1, Shell.P2), //
   O1P3(Shell.O1, Shell.P3), //
   O2P1(Shell.O2, Shell.P1), //
   O2Q1(Shell.O2, Shell.Q1), //
   O3P1(Shell.O3, Shell.P1), //
   O3Q1(Shell.O3, Shell.Q1);

   // Aliases for Seigbahn names
   public static final XRayTransition KA1 = KL3;
   public static final XRayTransition KA2 = KL2;
   public static final XRayTransition KB1 = KM3;
   public static final XRayTransition KiB2 = KN3;
   public static final XRayTransition KiiB2 = KN2;
   public static final XRayTransition KB3 = KM2;
   public static final XRayTransition KiB4 = KN5;
   public static final XRayTransition KiiB4 = KN4;
   public static final XRayTransition KiB5 = KM5;
   public static final XRayTransition KiiB5 = KM4;
   public static final XRayTransition LA1 = L3M5;
   public static final XRayTransition LA2 = L3M4;
   public static final XRayTransition LB15 = L3N4;
   public static final XRayTransition LB2 = L3N5;
   public static final XRayTransition LB5 = L3O4;
   public static final XRayTransition LB6 = L3N1;
   public static final XRayTransition LB7 = L3O1;
   public static final XRayTransition Ll = L3M1;
   public static final XRayTransition Ls = L3M3;
   public static final XRayTransition Lt = L3M2;
   public static final XRayTransition Lu = L3N6;
   public static final XRayTransition LB1 = L2M4;
   public static final XRayTransition LB17 = L2M3;
   public static final XRayTransition LG1 = L2N4;
   public static final XRayTransition LG5 = L2N1;
   public static final XRayTransition LG6 = L2O4;
   public static final XRayTransition LG8 = L2O1;
   public static final XRayTransition Ln = L2M1;
   public static final XRayTransition Lv = L2N6;
   public static final XRayTransition LB10 = L1M4;
   public static final XRayTransition LB3 = L1M3;
   public static final XRayTransition LB4 = L1M2;
   public static final XRayTransition LB9 = L1M5;
   public static final XRayTransition LG2 = L1N2;
   public static final XRayTransition LG11 = L1N5;
   public static final XRayTransition LG3 = L1N3;
   public static final XRayTransition LG4 = L1O3;
   public static final XRayTransition LG4p = L1O2;
   public static final XRayTransition MG = M3N5;
   public static final XRayTransition MZ2 = M4N2;
   public static final XRayTransition MA1 = M5N7;
   public static final XRayTransition MA2 = M5N6;
   public static final XRayTransition MZ1 = M5N3;

   private final Shell mInner;
   private final Shell mOuter;
   private final String mSeigbahn;

   public final static XRayTransition[] K_ALPHA = new XRayTransition[] {
      KL3,
      KL2,
      KL1
   };

   public final static XRayTransition[] K_BETA = new XRayTransition[] {
      KM3,
      KN3,
      KN2,
      KM2,
      KN5,
      KN4,
      KM5,
      KM4
   };

   public final static XRayTransition[] K_OTHER = new XRayTransition[] {
      KO2,
      KO3,
      KO4,
      KO5,
      KP2,
      KP3
   };

   public final static XRayTransition[] K_ALL = combine(K_ALPHA, K_BETA, K_OTHER);

   public final static XRayTransition[] L_ALPHA = new XRayTransition[] {
      L3M5,
      L3M4
   };

   public final static XRayTransition[] L_BETA = new XRayTransition[] {
      L3N4,
      L3N5,
      L3O4,
      L3N1,
      L3O1,
      L2M4,
      L2M3,
      L1M4,
      L1M3,
      L1M2,
      L1M5
   };

   public final static XRayTransition[] L_GAMMA = new XRayTransition[] {
      L2N4,
      L2N1,
      L2O4,
      L2O1,
      L1N2,
      L1N5,
      L1N3,
      L1O3,
      L1O5,
      L1O2
   };

   public final static XRayTransition[] L_OTHER = new XRayTransition[] {
      L3N2,
      L3N3,
      L3O2,
      L3O3,
      L3P1,
      L3O5,
      L3M1,
      L3M3,
      L3M2,
      L3N6,
      L3N7,
      L3O6,
      L3O7,
      L2M2,
      L2M5,
      L2N2,
      L2N3,
      L2N5,
      L2O2,
      L2O3,
      L2P2,
      L2M1,
      L2N6,
      L2N7,
      L1M1,
      L1N1,
      L1N4,
      L1O1,
      L1O4,
      L2O6,
      L2P1,
      L2Q1,
      L3Q1,
      L1P2,
      L1P3,
      L2P3,
      L3P2,
      L3P3
   };

   public static final XRayTransition[] L_ALL = combine(L_ALPHA, L_BETA, L_GAMMA, L_OTHER);

   public static final XRayTransition[] M_ALPHA = new XRayTransition[] {
      M5N7,
      M5N6
   };

   public static final XRayTransition[] M_GAMMA = new XRayTransition[] {
      M3N5
   };

   public static final XRayTransition[] M_ZETA = new XRayTransition[] {
      M4N2,
      M5N3
   };

   public static final XRayTransition[] M_OTHER = new XRayTransition[] {
      M1N2,
      M1N3,
      M1O2,
      M1O3,
      M2M4,
      M2N1,
      M2N4,
      M2O1,
      M2O4,
      M2Q1,
      M3M4,
      M3M5,
      M3N1,
      M3N4,
      M3O1,
      M3O4,
      M3O5,
      M3Q1,
      M4N3,
      M4O2,
      M4O3,
      M4N6,
      M4O6,
      M4P2,
      M4P3,
      M5O3,
      M2P1,
      M3P1,
      M1P2,
      M1P3,
      M5P3,
      M5O6,
      M5O7
   };

   public static final XRayTransition[] M_ALL = combine(M_ALPHA, M_GAMMA, M_ZETA, M_OTHER);

   public static final XRayTransition[] N_ALL = new XRayTransition[] {
      N1O2,
      N1O3,
      N2O1,
      N2O4,
      N1P2,
      N1P3,
      N2P1,
      N2Q1,
      N3O1,
      N3O4,
      N3O5,
      N3P1,
      N3Q1,
      N4O2,
      N4O3,
      N4N6,
      N4P2,
      N4P3,
      N4O6,
      N5N6,
      N5O3,
      N5O7,
      N5P3,
      N5O6,
      N6O4,
      N6O5,
      N7O5
   };

   public static final XRayTransition[] O_ALL = new XRayTransition[] {
      O1P2,
      O1P3,
      O2P1,
      O2Q1,
      O3P1,
      O3Q1
   };

   public static XRayTransition[] combine(final XRayTransition[]... xrtas) {
      int len = 0;
      for(final XRayTransition[] xrta : xrtas)
         len += xrta.length;
      final XRayTransition[] res = new XRayTransition[len];
      int i = 0;
      for(final XRayTransition[] xrta : xrtas)
         for(final XRayTransition xrt : xrta) {
            res[i] = xrt;
            ++i;
         }
      return res;
   }

   private XRayTransition(final Shell inner, final Shell outer) {
      this(inner, outer, null);
   }

   private XRayTransition(final Shell inner, final Shell outer, final String seigbahn) {
      Preconditions.checkNotNull(inner, "Inner shell must be non-null");
      Preconditions.checkNotNull(outer, "Outer shell must be non-null");
      Preconditions.checkArgument(inner.compareTo(outer) < 0, "The inner shell must be strictly inside the outer shell.");
      mInner = inner;
      mOuter = outer;
      mSeigbahn = seigbahn;
   }

   public static XRayTransition parse(final String str) {
      final String[] parts = str.split("-");
      if(parts.length == 2) {
         final Shell inner = Shell.parse(parts[0]);
         final Shell outer = Shell.parse(parts[1]);
         if((inner != null) && (outer != null) && (inner.ordinal() < outer.ordinal()))
            return XRayTransition.find(outer, inner);
      }
      return null;
   }

   @Override
   public String toString() {
      return mInner.toString() + "-" + mOuter.toString();
   }

   /**
    * Returns the transition name in HTML format.
    *
    * @return String in HTML
    */
   public String toHTML() {
      return mInner.toHTML() + "&#8209;" + mOuter.toHTML();
   }

   /**
    * Converts the Transition to the standard Seigbahn notation.
    *
    * @return String in HTML format
    */
   public String toSeigbahn() {
      return mSeigbahn != null ? mSeigbahn : toHTML();
   }

   /**
    * The inner shell. (The initially ionized shell)
    *
    * @return Shell
    */
   public Shell getInner() {
      return mInner;
   }

   /**
    * The outer shell (The donor shell)
    *
    * @return Shell
    */
   public Shell getOuter() {
      return mOuter;
   }

   private boolean isMember(final XRayTransition[] tra) {
      for(final XRayTransition trx : tra)
         if(equals(trx))
            return true;
      return false;

   }

   static public XRayTransition find(final Shell outer, final Shell inner) {
      for(final XRayTransition tr : values())
         if(outer.equals(tr.mOuter) && inner.equals(tr.mInner))
            return tr;
      return null;

   }

   public boolean isKAlpha() {
      return isMember(K_ALPHA);
   }

   public boolean isKBeta() {
      return isMember(K_BETA);
   }

   public boolean isK() {
      return isMember(K_ALL);
   }

   public boolean isL() {
      return isMember(L_ALL);
   }

   public boolean isLAlpha() {
      return isMember(L_ALPHA);
   }

   public boolean isLBeta() {
      return isMember(L_BETA);
   }

   public boolean isLGamma() {
      return isMember(L_GAMMA);
   }

   public boolean isLOther() {
      return isMember(L_ALL) && !(isLAlpha() || isLBeta() || isLGamma());
   }

   public boolean isMAlpha() {
      return isMember(M_ALPHA);
   }

   public boolean isMGamma() {
      return isMember(M_GAMMA);
   }

   public boolean isMZeta() {
      return isMember(M_ZETA);
   }

   public boolean isMOther() {
      return isMember(M_ALL) && !(isMAlpha() || isMGamma() || isMZeta());
   }

   public boolean isM() {
      return isMember(M_ALL);
   }

   public static void main(final String[] args) {
      for(final XRayTransition xrt : values()) {
         final String tmp = xrt.mInner.toString() + xrt.mOuter.toString();
         if(!tmp.equals(xrt.name()))
            System.out.println(xrt);
      }
      {
         final ArrayList<XRayTransition> al = new ArrayList<>(Arrays.asList(combine(K_ALL, L_ALL, M_ALL, N_ALL, O_ALL)));
         for(final XRayTransition xrt : values())
            al.remove(xrt);
         System.out.println(al);
      }
      {
         final ArrayList<XRayTransition> al = new ArrayList<>(Arrays.asList(values()));
         for(final XRayTransition xrt : combine(K_ALL, L_ALL, M_ALL, N_ALL, O_ALL))
            al.remove(xrt);
         System.out.println(al);
      }

   }

}
