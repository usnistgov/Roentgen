package gov.nist.microanalysis.roentgen.physics;

import java.util.Arrays;
import java.util.List;

import com.duckandcover.html.IToHTML;

/**
 * <p>
 * An enumeration representing the elements from Hydrogen(Z=1) up to Lawrencium
 * (Z=103).
 * </p>
 * <a>The atomic weight data comes from ... The weight data is nominal based on
 * an assumed isotopic mix which is typical of terrestrial matter.
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author nritchie
 * @version $Rev: 285 $
 */
public enum Element implements IToHTML {

   Hydrogen(1, "H", 1.00794), //
   Helium(2, "He", 4.0026), //
   Lithium(3, "Li", 6.941), //
   Beryllium(4, "Be", 9.01218), //
   Boron(5, "B", 10.811), //
   Carbon(6, "C", 12.0107), //
   Nitrogen(7, "N", 14.0067), //
   Oxygen(8, "O", 15.9994), //
   Fluorine(9, "F", 18.998), //
   Neon(10, "Ne", 20.1797), //
   Sodium(11, "Na", 22.9898), //
   Magnesium(12, "Mg", 24.305), //
   Aluminum(13, "Al", 26.9815), //
   Silicon(14, "Si", 28.0855), //
   Phosphorus(15, "P", 30.9738), //
   Sulfur(16, "S", 32.065), //
   Chlorine(17, "Cl", 35.453), //
   Argon(18, "Ar", 39.948), //
   Potassium(19, "K", 39.0983), //
   Calcium(20, "Ca", 40.078), //
   Scandium(21, "Sc", 44.9559), //
   Titanium(22, "Ti", 47.867), //
   Vanadium(23, "V", 50.9415), //
   Chromium(24, "Cr", 51.9961), //
   Manganese(25, "Mn", 54.938), //
   Iron(26, "Fe", 55.845), //
   Cobalt(27, "Co", 58.933), //
   Nickel(28, "Ni", 58.6934), //
   Copper(29, "Cu", 63.546), //
   Zinc(30, "Zn", 65.409), //
   Gallium(31, "Ga", 69.723), //
   Germanium(32, "Ge", 72.64), //
   Arsenic(33, "As", 74.9216), //
   Selenium(34, "Se", 78.96), //
   Bromine(35, "Br", 79.904), //
   Krypton(36, "Kr", 83.798), //
   Rubidium(37, "Rb", 85.4678), //
   Strontium(38, "Sr", 87.62), //
   Yttrium(39, "Y", 88.905), //
   Zirconium(40, "Zr", 91.224), //
   Niobium(41, "Nb", 92.906), //
   Molybdenum(42, "Mo", 95.94), //
   Technetium(43, "Tc", 98), //
   Ruthenium(44, "Ru", 101.07), //
   Rhodium(45, "Rh", 102.906), //
   Palladium(46, "Pd", 106.42), //
   Silver(47, "Ag", 107.868), //
   Cadmium(48, "Cd", 112.411), //
   Indium(49, "In", 114.818), //
   Tin(50, "Sn", 118.71), //
   Antimony(51, "Sb", 121.76), //
   Tellurium(52, "Te", 127.6), //
   Iodine(53, "I", 126.904), //
   Xenon(54, "Xe", 131.293), //
   Cesium(55, "Cs", 132.904), //
   Barium(56, "Ba", 137.327), //
   Lanthanum(57, "La", 138.906), //
   Cerium(58, "Ce", 140.116), //
   Praseodymium(59, "Pr", 140.907), //
   Neodymium(60, "Nd", 144.24), //
   Promethium(61, "Pm", 145), //
   Samarium(62, "Sm", 150.36), //
   Europium(63, "Eu", 151.964), //
   Gadolinium(64, "Gd", 157.25), //
   Terbium(65, "Tb", 158.925), //
   Dysprosium(66, "Dy", 160.5), //
   Holmium(67, "Ho", 164.93), //
   Erbium(68, "Er", 167.259), //
   Thulium(69, "Tm", 168.934), //
   Ytterbium(70, "Yb", 173.04), //
   Lutetium(71, "Lu", 174.967), //
   Hafnium(72, "Hf", 178.49), //
   Tantalum(73, "Ta", 180.948), //
   Tungsten(74, "W", 183.84), //
   Rhenium(75, "Re", 186.207), //
   Osmium(76, "Os", 190.23), //
   Iridium(77, "Ir", 192.217), //
   Platinum(78, "Pt", 195.078), //
   Gold(79, "Au", 196.967), //
   Mercury(80, "Hg", 200.59), //
   Thallium(81, "Tl", 204.383), //
   Lead(82, "Pb", 207.2), //
   Bismuth(83, "Bi", 208.98), //
   Polonium(84, "Po", 209), //
   Astatine(85, "At", 210), //
   Radon(86, "Rn", 222), //
   Francium(87, "Fr", 223), //
   Radium(88, "Ra", 226), //
   Actinium(89, "Ac", 227), //
   Thorium(90, "Th", 232.038), //
   Protactinium(91, "Pa", 231.036), //
   Uranium(92, "U", 238.029), //
   Neptunium(93, "Np", 237), //
   Plutonium(94, "Pu", 244), //
   Americium(95, "Am", 243), //
   Curium(96, "Cm", 247), //
   Berkelium(97, "Bk", 247), //
   Californium(98, "Cf", 251), //
   Einsteinium(99, "Es", 252), //
   Fermium(100, "Fm", 257), //
   Mendelevium(101, "Md", 258), //
   Nobelium(102, "No", 259), //
   Lawrencium(103, "Lr", 262);

   private final int mAtomicNumber;
   private final String mAbbrev;
   private final double mWeight;
   static final public int MIN_Z = 1;
   static final public int MAX_Z = 103;

   private Element(final int z, final String abbrev, final double weight) {
      mAtomicNumber = z;
      mAbbrev = abbrev;
      mWeight = weight;
   }

   /**
    * Finds an Element by atomic number.
    *
    * @param z The atomic number
    * @return Element or null if z corresponds to no element
    */
   public static Element byAtomicNumber(final int z) {
      final Element[] vals = Element.values();
      if((z > 0) && (z <= vals.length)) {
         final Element res = Element.values()[z - 1];
         assert res.mAtomicNumber == z;
         return res;
      } else
         return null;
   }

   /**
    * @param str
    * @return An Element or null if no matches
    */
   public static Element parse(final String str) {
      try {
         return byAtomicNumber(Integer.parseInt(str));
      }
      catch(final NumberFormatException ex) {
         for(final Element elm : Element.values())
            if(elm.toString().equalsIgnoreCase(str) || elm.mAbbrev.equals(str))
               return elm;
      }
      return null;
   }

   /**
    * @param str
    * @return An Element or null if no matches
    */
   public static Element parseAbbrev(final String str) {
      for(final Element elm : Element.values())
         if(elm.mAbbrev.equals(str))
            return elm;
      return null;
   }

   /**
    * Returns the atomic number associated with the element.
    *
    * @return Z the atomic number
    */
   public int getAtomicNumber() {
      return mAtomicNumber;
   }

   /**
    * Returns the atomic weight associated with the element.
    *
    * @return A the atomic weight in AMU
    */
   public double getAtomicWeight() {
      return mWeight;
   }

   /**
    * Returns the standard one or two letter abbreviation for the element.
    *
    * @return
    */
   public String getAbbrev() {
      return mAbbrev;
   }

   /**
    * Creates a list of elements from first to last in atomic number order.
    *
    * @param first Inclusive
    * @param last Inclusive
    * @return List&lt;Element&gt;
    */
   static public List<Element> range(final Element first, final Element last) {
      final Element[] res = new Element[last.getAtomicNumber() - first.getAtomicNumber() + 1];
      Element curr = first;
      for(int i = 0; i < res.length; ++i) {
         res[i] = curr;
         curr = curr.next();
      }
      return Arrays.asList(res);
   }

   static public List<Element> allElements() {
      return Arrays.asList(Element.values());
   }

   public boolean hasNext() {
      return byAtomicNumber(mAtomicNumber + 1) != null;
   }

   public Element next() {
      return byAtomicNumber(mAtomicNumber + 1);
   }

   public boolean hasPrevious() {
      return byAtomicNumber(mAtomicNumber - 1) != null;
   }

   public Element previous() {
      return byAtomicNumber(mAtomicNumber - 1);
   }

   @Override
   public String toHTML(final Mode mode) {
      switch(mode) {
         default:
         case TERSE:
            return getAbbrev();
         case NORMAL:
            return toString();
         case VERBOSE: {

            final StringBuffer sb = new StringBuffer();
            sb.append("<table><tr><th>Property</th><th>Value</th></tr>");
            sb.append("<tr><td>Name</td><td>" + toString() + "</td></tr>");
            sb.append("<tr><td>Abbreviation</td><td>" + getAbbrev() + "</td></tr>");
            sb.append("<tr><td>Atomic Number</td><td>" + Integer.toString(getAtomicNumber()) + "</td></tr>");
            sb.append("<tr><td>Atomic Weight</td><td>" + Double.toString(getAtomicWeight()) + " AMU</td></tr>");
            final AtomicShell[] shells = AtomicShell.forElement(this);
            if(shells != null) {
               sb.append("<tr><th colspan=\"2\">Atomic shells</th></tr>");
               for(final AtomicShell sh : shells)
                  sb.append("<tr><td>" + sh.toHTML(Mode.NORMAL) + "</td><td>" + Double.toString(sh.getEdgeEnergy())
                        + " eV</td></tr>");

            }
            final CharacteristicXRay[] cxrs = CharacteristicXRay.forElement(this);
            if(cxrs != null) {
               sb.append("<tr><th colspan=\"2\">Characteristic x-ray transitions</th></tr>");
               for(final CharacteristicXRay cxr : cxrs)
                  sb.append("<tr><td>" + cxr.toHTML(Mode.NORMAL) + "</td><td>" + cxr.toHTML(Mode.VERBOSE) + "</td></tr>");
            }
            sb.append("</table>");
            return sb.toString();
         }
      }
   }
}