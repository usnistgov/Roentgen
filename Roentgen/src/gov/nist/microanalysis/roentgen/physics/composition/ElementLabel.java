package gov.nist.microanalysis.roentgen.physics.composition;

import java.util.Objects;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * <p>
 * Acts as a mechanism to tag variables associated with elements.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class ElementLabel
   implements
   IToHTML {

   private final Element mElement;
   private final String mPrefix;
   private Composition mSource;

   /**
    * Constructs a ElementTag
    */
   protected ElementLabel(final Element elm, final String prefix) {
      this(elm, prefix, null);
   }

   /**
    * Constructs a ElementTag
    */
   protected ElementLabel(final Element elm, final String prefix, final Composition src) {
      mElement = elm;
      mPrefix = prefix;
      mSource = src;
   }

   public Element getElement() {
      return mElement;
   }

   public Composition getSource() {
      return mSource;
   }

   public void setSource(final Composition src) {
      mSource = src;
   }

   public int compareTo(final ElementLabel et) {
      int r = mPrefix.compareTo(et.mPrefix);
      if((r == 0) && (mSource != et.mSource))
         r = mSource.toString().compareTo(et.mSource.toString());
      return r == 0 ? mElement.compareTo(et.mElement) : r;
   }

   @Override
   public int hashCode() {
      return Objects.hash(mElement, mPrefix, mSource);
   }

   @Override
   public boolean equals(final Object obj) {
      if(this == obj)
         return true;
      if(obj == null)
         return false;
      if(getClass() != obj.getClass())
         return false;
      final ElementLabel other = (ElementLabel) obj;
      if(!mElement.equals(other.mElement))
         return false;
      if(!mPrefix.equals(other.mPrefix))
         return false;
      if(mSource != other.mSource)
         return false;
      return true;
   }

   @Override
   public String toString() {
      return HTML.stripTags(toHTML(Mode.TERSE));
   }

   /**
    * @see gov.nist.microanalysis.roentgen.html.IToHTML#toHTML(gov.nist.microanalysis.roentgen.Representation.IToHTML.Mode)
    */
   @Override
   public String toHTML(final Mode mode) {
      if(mSource == null)
         return mPrefix + "[" + mElement.getAbbrev() + "]";
      else
         return mPrefix + "[" + mSource.toHTML(Mode.TERSE) + ", " + mElement.getAbbrev() + "]";
   }
}
