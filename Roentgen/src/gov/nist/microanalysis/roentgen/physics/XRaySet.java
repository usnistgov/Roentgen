package gov.nist.microanalysis.roentgen.physics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.geometry.euclidean.oned.Interval;
import org.apache.commons.math3.geometry.partitioning.Region.Location;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.duckandcover.html.Table;
import com.duckandcover.lazy.SimplyLazy;
import com.google.common.base.Objects;
import com.google.common.base.Preconditions;

import gov.nist.microanalysis.roentgen.physics.Shell.Principle;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p> XRaySet implements a collection of XRay objects. All the XRay objects
 * must conform to XRaySet.IMemberTest to be a member of the XRaySet. This
 * typically means that the XRay object must be within a certain energy range of
 * the other XRay objects within the XRaySet. Implement IMemberTest to ensure
 * that all x-rays are within a certain number of FWHM of the other x-rays in
 * the set or to ensure that all x-rays are members of the same family or some
 * other useful criteria. </p> <p> Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 307 $
 */
public class XRaySet {

   /**
    * <p> Implement this interface to determine whether an XRay should be a
    * member of this TransitionSet. </p> <p> Copyright Nicholas W. M. Ritchie
    * 2014-2016 </p>
    *
    * @author nritchie
    * @version $Rev: 307 $
    */
   public interface IMemberTest {
      /**
       * Should this XRay be a member of this TransitionSet?
       *
       * @param ts
       * @param xr
       * @return true if the XRay should be a member, false otherwise.
       */
      public boolean belongs(XRaySet ts, XRay xr);
   }

   public static IMemberTest TRUE = new IMemberTest() {

      @Override
      public boolean belongs(final XRaySet ts, final XRay xr) {
         return true;
      }

   };

   public static class SimpleRangeMemberTest
      implements
      IMemberTest {
      private final double mdELow;
      private final double mdEHigh;

      public SimpleRangeMemberTest(final double dELow, final double dEHigh) {
         mdELow = dELow;
         mdEHigh = dEHigh;
      }

      @Override
      public boolean belongs(final XRaySet xrs, final XRay xr) {
         final double e = xr.getEnergy();
         final Interval iv = xrs.getEnergyInterval();
         return (iv.checkPoint(e - mdELow, 0.0) == Location.INSIDE) && (iv.checkPoint(e + mdEHigh, 0.0) == Location.INSIDE);
      }
   }

   /**
    * <p> Use this member test to return all the CharacteristicXRay objects
    * which result from ionizations in the specified Shell objects. </p> <p>
    * Copyright Nicholas W. M. Ritchie 2014-2016 </p>
    *
    * @author nritchie
    * @version $Rev: 307 $
    */
   public static class FamilyMemberTest
      implements
      IMemberTest {

      public Shell[] mShells;

      public FamilyMemberTest(final Shell[] shells) {
         mShells = shells.clone();
      }

      @Override
      public boolean belongs(final XRaySet xrs, final XRay xr) {
         if(xr instanceof CharacteristicXRay) {
            final Shell sh = ((CharacteristicXRay) xr).getTransition().getInner();
            for(final Shell sh2 : mShells)
               if(sh2.equals(sh))
                  return true;
         }
         return false;
      }
   }

   /**
    * <p> Is the edge energy for the specified CharacteristicXRay between eMin
    * and eMax? </p> <p> Copyright Nicholas W. M. Ritchie 2014-2016 </p>
    *
    * @author nritchie
    * @version $Rev: 307 $
    */
   public static class EdgeMemberTest
      implements
      IMemberTest {

      private final double mEMin, mEMax;

      public EdgeMemberTest(final double eMin, final double eMax) {
         Preconditions.checkArgument(eMin >= 0.0);
         Preconditions.checkArgument(eMin < eMax);
         mEMin = eMin;
         mEMax = eMax;
      }

      public EdgeMemberTest(final double eMax) {
         Preconditions.checkArgument(eMax > 0);
         mEMin = 0.0;
         mEMax = eMax;
      }

      @Override
      public boolean belongs(final XRaySet xrs, final XRay xr) {
         if(xr instanceof CharacteristicXRay) {
            final double ee = ((CharacteristicXRay) xr).getEdgeEnergy();
            return (ee >= mEMin) && (ee <= mEMax);
         }
         return false;
      }
   }

   protected HashSet<XRay> mSet = new HashSet<>();
   protected final ArrayList<IMemberTest> mTests;
   protected final SimplyLazy<Interval> mInterval = new SimplyLazy<Interval>() {

      @Override
      protected Interval initialize() {
         Interval res = null;
         for(final XRay tr : mSet) {
            final double energy = tr.getEnergy();
            if(res.checkPoint(energy, 0.0) != Location.INSIDE)
               res = new Interval(Math.min(res.getInf(), energy), Math.max(res.getSup(), energy));
         }
         return res;
      }

   };
   protected final SimplyLazy<String> mName = new SimplyLazy<String>() {

      @Override
      protected String initialize() {
         return buildName();
      }
   };

   public XRaySet(final IMemberTest mt) {
      Preconditions.checkArgument(mt != null);
      mTests = new ArrayList<>();
      add(mt);
   }

   public XRaySet(final Collection<IMemberTest> mtc) {
      mTests = new ArrayList<>();
      for(final IMemberTest mt : mtc)
         add(mt);
   }

   protected void add(final IMemberTest mt) {
      Preconditions.checkArgument(mSet.size() == 0);
      Preconditions.checkArgument(mt != null);
      mTests.add(mt);
      mName.reset();
   }

   public Interval getEnergyInterval() {
      return mInterval.get();
   }

   final public boolean add(final XRay xr) {
      final boolean res = belongs(xr);
      if(res) {
         mSet.add(xr);
         mInterval.reset();
         mName.reset();
      }
      return res;
   }

   final public void addAll(XRaySet xrs) {
      for(XRay xr : xrs.getXRaySet())
         add(xr);
   }

   protected boolean belongs(final XRay xr) {
      for(final IMemberTest mt : mTests)
         if(!mt.belongs(this, xr))
            return false;
      return true;
   }

   public Set<? extends XRay> getXRaySet() {
      return new TreeSet<>(mSet);
   }

   public void remove(final Collection<? extends XRay> xrc) {
      mSet.removeAll(xrc);
      mName.reset();
      mInterval.reset();
   }

   public int size() {
      return mSet.size();
   }

   @Override
   public XRaySet clone() {
      final XRaySet xrs = new XRaySet(mTests);
      xrs.mSet.addAll(mSet);
      xrs.mInterval.reset();
      xrs.mName.reset();
      return xrs;
   }

   public int compareTo(final XRaySet xrs) {
      final Set<XRay> all = new TreeSet<>();
      all.addAll(getXRaySet());
      all.addAll(xrs.getXRaySet());
      for(final XRay xr : all) {
         if(!mSet.contains(xr))
            return -1;
         if(!xrs.mSet.contains(xr))
            return 1;
      }
      return 0;
   }

   public boolean contains(final XRay cxr) {
      return mSet.contains(cxr);
   }

   /**
    * Merges all the XRay objects in the collection that meet the IMemberTest
    * objects associated with this set. The method repeated attempts to add the
    * XRay objects until no more can be added so that contiguous ranges of XRay
    * objects can be constructed.
    *
    * @param xrc
    * @return int The number of XRay objects added
    */
   final public int merge(final Collection<? extends XRay> xrc) {
      final HashSet<XRay> added = new HashSet<>();
      int cx;
      do {
         cx = 0;
         for(final XRay xr : xrc)
            if((!added.contains(xr)) && add(xr)) {
               added.add(xr);
               ++cx;
            }
      } while(cx > 0);
      return added.size();
   }

   @Override
   public int hashCode() {
      return Objects.hashCode(mSet);
   }

   @Override
   public boolean equals(final Object obj) {
      if(this == obj)
         return true;
      if(obj == null)
         return false;
      if(getClass() != obj.getClass())
         return false;
      return Objects.equal(mSet, ((XRaySet) obj).mSet);
   }

   protected String buildName() {
      return "X-rays[N=" + Integer.toString(size()) + ", E=" + mInterval.toString() + "]";
   }

   @Override
   public String toString() {
      return mName.get();
   }

   /**
    * <p> A set of CharacteristicXRay objects. </p> <p> Copyright Nicholas W. M.
    * Ritchie 2014-2016 </p>
    *
    * @author nritchie
    * @version $Rev: 307 $
    */
   public static class CharacteristicXRaySet
      extends
      XRaySet {
      private static class IsCharacteristic
         implements
         IMemberTest {

         @Override
         public boolean belongs(final XRaySet ts, final XRay xr) {
            return xr instanceof CharacteristicXRay;
         }
      }

      public CharacteristicXRaySet(final Collection<IMemberTest> mtc) {
         super(mtc);
         add(new IsCharacteristic());
      }

      /**
       * Constructs a CharacteristicXRaySet that can contain any type of
       * CharacteristicXRay.
       */
      public CharacteristicXRaySet() {
         super(new IsCharacteristic());
      }

      public static CharacteristicXRaySet build(CharacteristicXRay cxr) {
         final CharacteristicXRaySet res = new CharacteristicXRaySet();
         res.add(cxr);
         return res;
      }

      /**
       * Returns a set containing all Element objects for which a
       * CharacteristicXRay in in this set.
       *
       * @return Set&lt;Element&gt;
       */
      public Set<Element> getElementSet() {
         final Set<Element> res = new HashSet<>();
         for(final XRay xr : mSet)
            res.add(((CharacteristicXRay) xr).getElement());
         return res;
      }

      /**
       * Return all the SingleElementTransitionSets that can be associated with
       * the specified element.
       *
       * @param elm
       * @return List&lt;SingleElementTransitionSet&gt;
       */
      public List<ElementXRaySet> getSingleElementTransitionSets(final Element elm) {
         final HashSet<XRay> tmp = new HashSet<>();
         for(final XRay xr : mSet) {
            assert xr instanceof CharacteristicXRay;
            final CharacteristicXRay cxr = (CharacteristicXRay) xr;
            if(cxr.getElement().equals(elm))
               tmp.add(cxr);
         }
         final ArrayList<ElementXRaySet> res = new ArrayList<>();
         while(tmp.size() > 0) {
            final ElementXRaySet exrs = new ElementXRaySet(elm, mTests);
            exrs.merge(tmp);
            assert exrs.size() > 0;
            tmp.removeAll(exrs.getXRaySet());
            res.add(exrs);
         }
         return res;
      }

      @Override
      public CharacteristicXRaySet clone() {
         final CharacteristicXRaySet xrs = new CharacteristicXRaySet(mTests);
         xrs.mSet.addAll(mSet);
         xrs.mInterval.reset();
         xrs.mName.reset();
         return xrs;
      }

      @Override
      protected String buildName() {
         final StringBuffer elms = new StringBuffer();
         for(final Element elm : getElementSet())
            elms.append(elm.getAbbrev() + ",");
         final BasicNumberFormat df = new BasicNumberFormat("0");
         final Interval ii = mInterval.get();
         return "Characteristic x-rays[" + elms + "N=" + Integer.toString(size()) + ", E=[" + df.format(ii.getInf()) + ","
               + df.format(ii.getSup()) + "]]";
      }

      public Set<CharacteristicXRay> getSetOfCharacteristicXRay() {
         final TreeSet<CharacteristicXRay> res = new TreeSet<>();
         for(final XRay xr : mSet)
            res.add((CharacteristicXRay) xr);
         return res;
      }

      public CharacteristicXRay getBrightest() {
         CharacteristicXRay brightest = null;
         double weight = -1.0;
         for(final XRay xr : mSet) {
            final CharacteristicXRay cxr = (CharacteristicXRay) xr;
            if(cxr.getWeight() > weight) {
               brightest = cxr;
               weight = cxr.getWeight();
            }
         }
         return brightest;
      }

      /**
       * Returns the sum of the weights of the transitions in this
       * {@link CharacteristicXRaySet}.
       *
       * @return double (nominally between 0 and 1 assuming that the set comes
       *         from a single line family (K, L or M).)
       */
      public double getSumOfWeights() {
         double weight = 0.0;
         for(final XRay xr : mSet)
            weight += ((CharacteristicXRay) xr).getWeight();
         return weight;
      }

   }

   /**
    * <p> ElementXRaySet represents a set of CharacteristicXRay objects
    * associated with one and only one Element. </p> <p> Copyright Nicholas W.
    * M. Ritchie 2014-2016 </p>
    *
    * @author nritchie
    * @version $Rev: 307 $
    */
   public static class ElementXRaySet
      extends
      CharacteristicXRaySet
      implements
      IToHTML {

      private final Element mElement;

      private static final class ElementMemberTest
         implements
         IMemberTest {

         private final Element mElement;

         ElementMemberTest(final Element elm) {
            mElement = elm;
         }

         @Override
         public boolean belongs(final XRaySet xrs, final XRay xr) {
            return (xr instanceof CharacteristicXRay) && ((CharacteristicXRay) xr).getElement().equals(mElement);
         }
      };

      /**
       * Constructs a ElementXRaySet associated with the specified Element that
       * fulfills the IMemberTest objects in mtc.
       *
       * @param elm
       * @param mtc
       */
      public ElementXRaySet(final Element elm, final Collection<IMemberTest> mtc) {
         super(mtc);
         add(new ElementMemberTest(elm));
         mElement = elm;
      }

      public ElementXRaySet(final Collection<CharacteristicXRay> ccxr) {
         this(ccxr.iterator().next().getElement());
         for(final CharacteristicXRay cxr : ccxr)
            add(cxr);
      }

      public ElementXRaySet(final Element elm, final CharacteristicXRaySet cxrs) {
         this(elm);
         for(final CharacteristicXRay cxr : cxrs.getSetOfCharacteristicXRay())
            if(cxr.getElement().equals(elm))
               add(cxr);
      }

      /**
       * Constructs a ElementXRaySet associated with the specified Element.
       *
       * @param elm
       */
      public ElementXRaySet(final Element elm) {
         super();
         add(new ElementMemberTest(elm));
         mElement = elm;
      }

      /**
       * Constructs a ElementXRaySet associated with the specified Element.
       *
       * @param elm
       */
      public ElementXRaySet(final CharacteristicXRay cxr) {
         super();
         add(new ElementMemberTest(cxr.getElement()));
         mElement = cxr.getElement();
         add(cxr);
      }

      /**
       * A set containing all transitions below maxE for the specified element.
       *
       * @param elm
       * @param maxE in eV
       * @return ElementXRaySet
       */
      public static ElementXRaySet buildElementXRaySet(final Element elm, final double maxE) {
         final ElementXRaySet res = new ElementXRaySet(elm);
         final CharacteristicXRay[] cxrs = CharacteristicXRay.forElement(elm);
         for(final XRay xr : cxrs)
            if(xr.getEnergy() < maxE)
               res.add(xr);
         return res;
      }

      /**
       * A set containing all transitions below maxE for the specified element.
       *
       * @param elm
       * @param maxE in eV
       * @return ElementXRaySet
       */
      public static ElementXRaySet buildElementXRaySet(final AtomicShell ashell) {
         final ElementXRaySet res = new ElementXRaySet(ashell.getElement());
         final CharacteristicXRay[] cxrs = CharacteristicXRay.forAtomicShell(ashell);
         for(final XRay xr : cxrs)
            res.add(xr);
         return res;
      }

      /**
       * The Element with which this ElementXRaySet is associated.
       *
       * @return Element
       */
      public Element getElement() {
         return mElement;
      }

      /**
       * Create a List of all the ElementXRaySet objects associated with the
       * specified element such that each ElementXRaySet meets the specified
       * collection of IMemberTest objects.
       *
       * @param elm
       * @param mt
       * @return List&lt;SingleElementTransitionSet&gt;
       */
      public static List<ElementXRaySet> create(final Element elm, final Collection<IMemberTest> mtc) {
         final List<ElementXRaySet> res = new ArrayList<>();
         for(final CharacteristicXRay cxr : CharacteristicXRay.forElement(elm)) {
            boolean added = false;
            for(final ElementXRaySet sets : res)
               if(sets.add(cxr)) {
                  added = true;
                  break;
               }
            if(!added) {
               // Check whether it could be added to a new, empty SETS
               final ElementXRaySet sets = new ElementXRaySet(elm, mtc);
               if(sets.add(cxr))
                  res.add(sets);
            }
         }
         boolean merged;
         do {
            merged = false;
            for(int i = res.size() - 1; i >= 0; --i) {
               for(int j = i - 1; j >= 0; --j) {
                  final Set<? extends XRay> tmp = res.get(j).getXRaySet();
                  if(res.get(i).merge(tmp) > 0) {
                     merged = true;
                     res.get(j).remove(res.get(i).getXRaySet());
                     if(res.get(j).size() == 0)
                        res.remove(j);
                     break;
                  }
               }
               if(merged)
                  break;
            }
         } while(merged);
         return res;
      }

      @Override
      public ElementXRaySet clone() {
         final ElementXRaySet xrs = new ElementXRaySet(mElement, mTests);
         xrs.mSet.addAll(mSet);
         xrs.mInterval.reset();
         xrs.mName.reset();
         return xrs;
      }

      @Override
      protected String buildName() {
         final StringBuffer res = new StringBuffer();
         res.append("[");
         // Present in this XRaySet
         boolean k = false, ka = false, kb = false, l = false, la = false, lb = false, lo = false, m = false, ma = false,
               mz = false, mo = false;
         // Present in the element
         boolean cka = false, ckb = false, cla = false, clb = false, clo = false, cma = false, cmz = false, cmo = false;
         final Set<CharacteristicXRay> allChar = ElementXRaySet.buildElementXRaySet(mElement, Double.MAX_VALUE).getSetOfCharacteristicXRay();
         for(final CharacteristicXRay cx : allChar) {
            final XRayTransition tr = cx.getTransition();
            cka |= tr.isKAlpha();
            ckb |= tr.isKBeta();
            cla |= tr.isLAlpha();
            clb |= tr.isLBeta();
            tr.isLGamma();
            clo |= tr.isLOther();
            cma |= tr.isMAlpha();
            cmz |= tr.isMZeta();
            cmo |= tr.isMOther();
            if(contains(cx)) {
               ka |= tr.isKAlpha();
               kb |= tr.isKBeta();
               k |= tr.isK();
               la |= tr.isLAlpha();
               lb |= tr.isLBeta();
               tr.isLGamma();
               lo |= tr.isLOther();
               l |= tr.isL();
               ma |= tr.isMAlpha();
               mz |= tr.isMZeta();
               mo |= tr.isMOther();
               m |= tr.isM();
            }
         }
         if((!cka || ka) && (!ckb || kb) && (!cla || la) && (!clb || lb) && (!cma || ma) && (!cmz || mz))
            res.append(mElement.getAbbrev() + " all");
         else {
            final int len = res.length();
            if(k) {
               if((ka || kb) && ((!cka || ka) && (!ckb || kb)))
                  res.append(mElement.getAbbrev() + " K");
               else {
                  if(ka)
                     res.append(mElement.getAbbrev() + " K\u03B1");
                  if(kb)
                     res.append(mElement.getAbbrev() + " K\u03B2");
                  if(!(ka || kb))
                     res.append(mElement.getAbbrev() + " K misc");
               }
            }
            if(l)
               if((!cla || la) && (!clb || lb) && (!clo || lo))
                  res.append(mElement.getAbbrev() + " L");
               else {
                  if(la)
                     res.append(mElement.getAbbrev() + " L\u03B1");
                  if(lb)
                     res.append(mElement.getAbbrev() + " L\u03B2");
                  if(lo)
                     res.append(mElement.getAbbrev() + " L other");
                  if(!(la || lb || lo))
                     res.append(mElement.getAbbrev() + " L misc");
               }
            if(m) {
               if((!cma || ma) && (!cmz || mz) && (!cmo || mo))
                  res.append(mElement.getAbbrev() + " M");
               if(ma)
                  res.append(mElement.getAbbrev() + " M&alpha;");
               if(mz)
                  res.append(mElement.getAbbrev() + " M&zeta;");
               if(mo)
                  res.append(mElement.getAbbrev() + " M other");
               if(!(ma || mz || mo))
                  res.append(mElement.getAbbrev() + " M misc");
            }
            if(res.length() == len) {
               res.append(getBrightest().toString());
               if(mSet.size() > 1) {
                  res.append(" + ");
                  res.append(mSet.size() - 1);
                  res.append(" others");
               }
            }
         }
         res.append("]");
         return res.toString();
      }

      public Set<AtomicShell> getShells() {
         final HashSet<AtomicShell> res = new HashSet<>();
         for(final XRay xr : mSet) {
            assert xr instanceof CharacteristicXRay;
            final CharacteristicXRay cxr = (CharacteristicXRay) xr;
            res.add(cxr.getInner());
         }
         return res;
      }

      @Override
      public int hashCode() {
         return Objects.hashCode(mSet, mElement);
      }

      @Override
      public boolean equals(final Object obj) {
         return super.equals(obj) && (((ElementXRaySet) obj).mElement == this.mElement);
      }

      @Override
      public String toHTML(final Mode mode) {
         if(mode == Mode.TERSE) {
            final CharacteristicXRay br = this.getBrightest();
            return br.toHTML(Mode.TERSE) + " + " + Integer.toString(mSet.size() - 1) + " others";
         } else if(mode == Mode.NORMAL)
            return HTML.escape(buildName());
         else {
            final Table t = new Table();
            t.addRow(Table.th("Transition"), Table.th("Details"));
            for(final CharacteristicXRay cxr : getSetOfCharacteristicXRay()) {
               t.addRow(Table.td(cxr.toHTML(Mode.TERSE)), Table.td(cxr.toHTML(Mode.VERBOSE)));

            }
            return t.toHTML(Mode.VERBOSE);
         }
      }
   }

   /**
    * Build the ElementXRaySet associated with the specified Element.
    * 
    * @param elm The element
    * @return ElementXRaySet
    */
   static public ElementXRaySet build(Element elm) {
      final ElementXRaySet exrs = new ElementXRaySet(elm);
      for(CharacteristicXRay cxr : CharacteristicXRay.forElement(elm))
         exrs.add(cxr);
      return exrs;
   }

   /**
    * Build the ElementXRaySet associated with the specified Element and energy
    * interval.
    * 
    * @param elm The element
    * @param iv The interval of edge energies which to include
    * @return ElementXRaySet
    */
   static public ElementXRaySet build(Element elm, Interval iv) {
      ElementXRaySet exrs = new ElementXRaySet(elm);
      for(CharacteristicXRay cxr : CharacteristicXRay.forElement(elm))
         if(iv.checkPoint(cxr.getEdgeEnergy(), 0.0) == Location.INSIDE)
            exrs.add(cxr);
      return exrs;
   }

   /**
    * Build the ElementXRaySet associated with the specified Element and
    * Shell(s).
    * 
    * @param elm
    * @param shells
    * @return ElementXRaySet
    */
   static public ElementXRaySet build(Element elm, Shell[] shells) {
      ElementXRaySet exrs = new ElementXRaySet(elm);
      for(CharacteristicXRay cxr : CharacteristicXRay.forElement(elm))
         for(Shell sh : shells)
            if(cxr.getInner().getShell() == sh)
               exrs.add(cxr);
      return exrs;
   }

   /**
    * Build the ElementXRaySet associated with the specified Element and
    * Shell(s).
    * 
    * @param elm
    * @param shells
    * @return ElementXRaySet
    */
   static public ElementXRaySet build(Element elm, Principle prin) {
      ElementXRaySet exrs = new ElementXRaySet(elm);
      for(CharacteristicXRay cxr : CharacteristicXRay.forElement(elm))
         if(cxr.getInner().getShell().getFamily() == prin)
            exrs.add(cxr);
      return exrs;
   }

   /**
    * Merge an array of {@link CharacteristicXRaySet} into a single one.
    * 
    * @param cxrss
    * @return A new CharacteristicXRaySet instance.
    */
   static public CharacteristicXRaySet merge(CharacteristicXRaySet... cxrss) {
      CharacteristicXRaySet res = new CharacteristicXRaySet();
      for(CharacteristicXRaySet cxrs : cxrss)
         res.addAll(cxrs);
      return res;
   }

}
