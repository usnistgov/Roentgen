package gov.nist.microanalysis.roentgen.math.uncertainty;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;
import com.google.common.base.Objects;

/**
 * <p>
 * The BaseTag class makes it easy to construct tags for use with
 * {@link NamedMultivariateJacobianFunction} and {@link UncertainValues} class
 * instances.
 * </p>
 * <p>
 * Tags need to be derived from a class that implements hashCode() and equals().
 * This class meets those requirements so long as the Object arguments also
 * implement hashCode() and equals().
 * </p>
 *
 * @author Nicholas
 */
public class BaseTag
   implements
   IToHTML {

   final String mName;
   final Object mObject1;
   final Object mObject2;
   final Object mObject3;

   /**
    * Creates a BaseTag object with an name from three objects. The objects must
    * implement hashCode() and equals().
    *
    * @param name
    * @param obj1
    * @param obj2
    * @param obj3
    */
   protected BaseTag(final String name, final Object obj1, final Object obj2, final Object obj3) {
      mName = name;
      mObject1 = obj1;
      mObject2 = obj2;
      mObject3 = obj3;
   }

   /**
    * Creates a BaseTag object with an name from two objects. The objects must
    * implement hashCode() and equals().
    *
    * @param name
    * @param obj1
    * @param obj2
    */
   protected BaseTag(final String name, final Object obj1, final Object obj2) {
      this(name, obj1, obj2, null);
   }

   /**
    * Creates a BaseTag object with an name from one object. The object must
    * implement hashCode() and equals().
    *
    * @param name
    * @param obj
    */
   protected BaseTag(final String name, final Object obj) {
      this(name, obj, null, null);
   }

   protected Object getObject1() {
      return mObject1;
   }

   protected Object getObject2() {
      return mObject2;
   }

   protected Object getObject3() {
      return mObject3;
   }

   @Override
   public int hashCode() {
      return Objects.hashCode(mName, mObject1, mObject2, mObject3);
   }

   @Override
   public boolean equals(final Object obj) {
      if(this == obj)
         return true;
      if(obj == null)
         return false;
      if(getClass() != obj.getClass())
         return false;
      final BaseTag other = (BaseTag) obj;
      return Objects.equal(mName, other.mName) && Objects.equal(mObject1, other.mObject1)
            && Objects.equal(mObject2, other.mObject2) && Objects.equal(mObject3, other.mObject3);
   }

   @Override
   public String toString() {
      return HTML.stripTags(toHTML(Mode.TERSE));
   }

   @Override
   public String toHTML(final Mode mode) {
      final StringBuffer sb = new StringBuffer();
      sb.append(mName);
      sb.append("[");
      if(mObject1 != null)
         sb.append(HTML.toHTML(mObject1, mode));
      if(mObject2 != null) {
         sb.append(",");
         sb.append(HTML.toHTML(mObject2, mode));
      }
      if(mObject3 != null) {
         sb.append(",");
         sb.append(HTML.toHTML(mObject3, mode));
      }
      sb.append("]");
      return sb.toString();
   }

}
