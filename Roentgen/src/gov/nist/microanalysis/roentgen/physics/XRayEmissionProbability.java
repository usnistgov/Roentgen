package gov.nist.microanalysis.roentgen.physics;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.lang3.concurrent.ConcurrentException;
import org.apache.commons.lang3.concurrent.LazyInitializer;

import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.Globals;
import gov.nist.microanalysis.roentgen.math.Utility;
import gov.nist.microanalysis.roentgen.utility.BasicNumberFormat;

/**
 * <p>
 * When an ionization occurs, it can lead to a cascade of x-rays, Auger and
 * Coster-Kronig transitions. We care primarily about x-rays so I've extracted
 * the x-ray data from the ENDF-VII database. The database only contains the
 * events occuring as an immediate result of the ionization (not the subsequent
 * cascade of ionizations as the ionization percolates towards the continuum).
 * So the data set was processed to account for the cascade and all x-rays that
 * could occur as the result of a single ionization are included in the file
 * "relax_endf6.csv" along with the fractional likelyhood of occurance. The
 * fractional likelyhood is diminished by the likelyhood of being ionized and
 * the fluorescence yield (radiative/non-radiative). The likelyhood accounts for
 * Auger and Coster-Kronig transitions.
 * </p>
 * <p>
 * Loads the x-ray probability data as a function of ionized atomic shell once
 * and then feeds it out on a per atomic shell basis as needed.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class XRayEmissionProbability
   implements
   IToHTML {

   private static final LazyInitializer<Map<AtomicShell, Map<CharacteristicXRay, Double>>> mData = new LazyInitializer<Map<AtomicShell, Map<CharacteristicXRay, Double>>>() {

      @Override
      protected Map<AtomicShell, Map<CharacteristicXRay, Double>> initialize() {
         final Map<AtomicShell, Map<CharacteristicXRay, Double>> res = new HashMap<>();
         try {
            final InputStream is = getClass().getResourceAsStream("relax_endf6.csv");
            final InputStreamReader isr = new InputStreamReader(is, Charset.forName("US-ASCII"));
            try (final BufferedReader br = new BufferedReader(isr)) {
               AtomicShell prev = null;
               Map<CharacteristicXRay, Double> data = new HashMap<>();
               for(String str = br.readLine(); str != null; str = br.readLine()) {
                  if(str.startsWith("//"))
                     continue;
                  final String[] items = str.split(",");
                  final Element elm = Element.byAtomicNumber(Integer.parseInt(items[0]));
                  if(elm.equals(Element.Neptunium))
                     break;
                  final AtomicShell src = AtomicShell.find(elm, Shell.values()[Integer.parseInt(items[1])]);
                  if(!src.equals(prev)) {
                     if(prev != null)
                        res.put(prev, Collections.unmodifiableMap(data));
                     data = new HashMap<>();
                     prev = src;
                  }
                  final Shell inner = Shell.values()[Integer.parseInt(items[2])];
                  final Shell outer = Shell.values()[Integer.parseInt(items[3])];
                  final XRayTransition xrt = XRayTransition.find(outer, inner);

                  final CharacteristicXRay cxr = CharacteristicXRay.create(elm, xrt);
                  if(cxr != null) {
                     final double f = Double.parseDouble(items[4]);
                     data.put(cxr, f);
                  }
               }
            }
         }
         catch(final Exception e) {
            System.out.println(e.getMessage());
            Globals.getLogger().error(e.getMessage());
         }
         return Collections.unmodifiableMap(res);
      }

   };

   private final AtomicShell mIonized;
   private final Map<CharacteristicXRay, Double> mProbabilities;

   public static boolean isAvailable(final AtomicShell sh) {
      try {
         return mData.get().containsKey(sh);
      }
      catch(final ConcurrentException e) {
         Globals.getLogger().catching(e);
      }
      return false;
   }

   public static Set<AtomicShell> getAvailableShells(final Element elm) {
      final Set<AtomicShell> res = new TreeSet<>();
      for(final AtomicShell sh : AtomicShell.forElement(elm))
         if(isAvailable(sh))
            res.add(sh);
      return res;
   }

   public XRayEmissionProbability(final AtomicShell ionized) {
      mIonized = ionized;
      Map<CharacteristicXRay, Double> prob = null;
      try {
         prob = mData.get().get(ionized);
         assert prob != null;
      }
      catch(final ConcurrentException e) {
         Globals.getLogger().catching(e);
         System.exit(1);
      }
      mProbabilities = prob;
   }

   public AtomicShell getIonized() {
      return mIonized;
   }

   public Map<CharacteristicXRay, Double> getProbabilties() {
      return mProbabilities;
   }

   public double getFluorescenceYield() {
      double sum = 0.0;
      for(final Map.Entry<CharacteristicXRay, Double> me : mProbabilities.entrySet())
         if(me.getKey().getInner().equals(mIonized))
            sum += me.getValue().doubleValue();
      return sum;
   }

   /**
    * The relative weights of the various CharacteristicXRay emission lines
    * produced when this shell is ionized. The weights are normalized to the
    * most probable being a weight of unity.
    *
    * @param direct true to only include those transitions directly from this
    *           shell
    * @return Map&lt;CharacteristicXRay, Double&gt;
    */
   public Map<CharacteristicXRay, Double> getWeights(final boolean direct) {
      final Map<CharacteristicXRay, Double> prob = getProbabilties();
      final Map<CharacteristicXRay, Double> res = new TreeMap<>();
      final double max = Utility.max(prob.values());
      for(final CharacteristicXRay cxr : prob.keySet())
         if(!direct || (cxr.getInner().equals(mIonized)))
            res.put(cxr, prob.get(cxr).doubleValue() / max);
      return res;
   }

   /**
    * Calculates weights of lines based on the specified beam energy.
    *
    * @param elm Element
    * @param e0 Beam energy (eV)
    * @return Map&lt;CharacteristicXRay, Double&gt;
    */
   public static Map<CharacteristicXRay, Double> getWeights(final Element elm, final double e0) {
      final Map<CharacteristicXRay, Double> res = new TreeMap<>();
      final Set<AtomicShell> shells = XRayEmissionProbability.getAvailableShells(elm);
      double maxP = 0.0;
      for(final AtomicShell shell : shells) {
         if(IonizationCrossSection.isSupported(shell)) {
            final double icx = IonizationCrossSection.value(shell, e0);
            final XRayEmissionProbability xep = new XRayEmissionProbability(shell);
            for(final Map.Entry<CharacteristicXRay, Double> me : xep.getProbabilties().entrySet()) {
               final CharacteristicXRay cxr = me.getKey();
               final double p = me.getValue() * icx + res.getOrDefault(cxr, 0.0);
               res.put(cxr, p);
               maxP = Math.max(maxP, p);
            }
         }
      }
      for(final Map.Entry<CharacteristicXRay, Double> me : res.entrySet())
         me.setValue(me.getValue() / maxP);
      return res;
   }

   @Override
   public String toHTML(final Mode mode) {
      final StringBuffer sb = new StringBuffer();
      if(mode == Mode.TERSE) {
         sb.append("<p>X-ray emission probabilities for " + mIonized.toHTML(Mode.TERSE) + "</p>");
      } else if((mode == Mode.NORMAL) || (mode == Mode.VERBOSE)) {
         final BasicNumberFormat nf = new BasicNumberFormat("#,##0");
         final BasicNumberFormat nf2 = new BasicNumberFormat("0.000E0");
         final Map<CharacteristicXRay, Double> wgts = getWeights(false);
         sb.append("<table>");
         sb.append("<tr><th>Transition</th>");
         if(mode == Mode.VERBOSE)
            sb.append("<th>Inner</th><th>Outer</th><th>Energy<br/>(eV)</th>");
         sb.append("<th>Normalized<br/>Probability<th>Probability</th></tr>");
         for(final Map.Entry<CharacteristicXRay, Double> me : mProbabilities.entrySet()) {
            sb.append("<tr><td>");
            final CharacteristicXRay cxr = me.getKey();
            sb.append(cxr.toHTML(Mode.TERSE));
            sb.append("</td><td>");
            if(mode == Mode.VERBOSE) {
               sb.append(cxr.getInner().toHTML(Mode.TERSE));
               sb.append("</td><td>");
               sb.append(cxr.getOuter().toHTML(Mode.TERSE));
               sb.append("</td><td>");
               sb.append(nf.format(cxr.getEnergy()));
               sb.append("</td><td>");
            }
            sb.append(nf2.format(wgts.get(cxr)));
            sb.append("</td><td>");
            sb.append(nf2.format(me.getValue()));
            sb.append("</td></tr>");
         }
         sb.append("</table>");
      }
      return sb.toString();
   }

}
