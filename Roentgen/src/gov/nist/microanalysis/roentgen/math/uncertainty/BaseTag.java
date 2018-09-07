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
public class BaseTag<H, I, J> //
		implements IToHTML {

	final String mName;
	final H mObject1;
	final I mObject2;
	final J mObject3;

	final transient int mHashCode;

	/**
	 * Creates a BaseTag object with an name from three objects. The objects must
	 * implement hashCode() and equals().
	 *
	 * @param name
	 * @param obj1
	 * @param obj2
	 * @param obj3
	 */
	protected BaseTag(final String name, final H obj1, final I obj2, final J obj3) {
		mName = name;
		mObject1 = obj1;
		mObject2 = obj2;
		mObject3 = obj3;
		// Calculate this once to speed access
		mHashCode = Objects.hashCode(mName, mObject1, mObject2, mObject3);

	}

	/**
	 * Creates a BaseTag object with an name from two objects. The objects must
	 * implement hashCode() and equals().
	 *
	 * @param name
	 * @param obj1
	 * @param obj2
	 */
	protected BaseTag(final String name, final H obj1, final I obj2) {
		this(name, obj1, obj2, null);
	}

	/**
	 * Creates a BaseTag object with an name from one object. The object must
	 * implement hashCode() and equals().
	 *
	 * @param name
	 * @param obj
	 */
	protected BaseTag(final String name, final H obj) {
		this(name, obj, null, null);
	}

	protected H getObject1() {
		return mObject1;
	}

	protected I getObject2() {
		return mObject2;
	}

	protected J getObject3() {
		return mObject3;
	}

	@Override
	public int hashCode() {
		return mHashCode;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		final BaseTag<?, ?, ?> other = (BaseTag<?, ?, ?>) obj;
		return Objects.equal(mName, other.mName) && //
				Objects.equal(mObject1, other.mObject1) && //
				Objects.equal(mObject2, other.mObject2) && //
				Objects.equal(mObject3, other.mObject3);
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
		if (mObject1 != null)
			sb.append(HTML.toHTML(mObject1, Mode.TERSE));
		if (mObject2 != null) {
			sb.append(",");
			sb.append(HTML.toHTML(mObject2, Mode.TERSE));
		}
		if (mObject3 != null) {
			sb.append(",");
			sb.append(HTML.toHTML(mObject3, Mode.TERSE));
		}
		sb.append("]");
		return sb.toString();
	}

}
