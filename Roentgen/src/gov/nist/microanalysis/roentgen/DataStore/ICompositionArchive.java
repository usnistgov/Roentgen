package gov.nist.microanalysis.roentgen.DataStore;

import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.Pair;

import gov.nist.microanalysis.roentgen.physics.Element;
import gov.nist.microanalysis.roentgen.physics.composition.Composition;

/**
 * <p>
 * This interface provides a generic mechanism to look up the definition of a
 * named Composition in a database or some such.
 * </p>
 *
 * @author Nicholas
 * @version $Rev: 312 $
 */
public interface ICompositionArchive {

   /**
    * This interface provides a mechanism to look up the definition of a named
    * Composition in a database or some such. @param name @return A Composition
    * or null if none defined for this name.
    *
    * @throws Exception
    */
   public Curated<Composition> findComposition(String name)
         throws Exception;

   /**
    * Searches the database for the Compositions for which mass fractions fall
    * on the specified intervals.
    *
    * @param criteria A map of Element to mass faction Interval
    * @return All Compositions matching the criteria
    */
   public Set<Curated<Composition>> findMatchingComposition(Map<Element, Pair<Double, Double>> criteria)
         throws Exception;

   /**
    * Contents of database by name
    *
    * @throws Exception
    */
   public Set<String> getCompositionNames()
         throws Exception;

   /**
    * Add the specified {@link Composition} to the database with the specified
    * name. If name already has been defined, then the old definition is
    * overwritten.
    *
    * @param name
    * @param comp
    * @return Curated&lt;Composition&gt;
    */
   public Curated<Composition> add(String name, Composition comp);

}