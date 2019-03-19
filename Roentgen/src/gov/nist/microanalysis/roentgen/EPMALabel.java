/**
 *
 */
package gov.nist.microanalysis.roentgen;

import java.util.Objects;

import com.duckandcover.html.HTML;
import com.duckandcover.html.IToHTML;

import gov.nist.microanalysis.roentgen.physics.XRay;
import gov.nist.microanalysis.roentgen.physics.composition.Material;

/**
 * @author nicho
 *
 */
public class EPMALabel implements IToHTML {

	final private String mName;

	public EPMALabel(String name) {
		mName = name;
	}

	@Override
	public int hashCode() {
		return Objects.hash(mName);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		EPMALabel other = (EPMALabel) obj;
		return Objects.equals(mName, other.mName);
	}

	@Override
	public String toHTML(Mode mode) {
		return mName;
	}

	public String toString() {
		return mName;
	}

	static public class BaseLabel<H, I, J> //
			extends EPMALabel {

		final H mObject1;
		final I mObject2;
		final J mObject3;

		final int mHashCode;

		/**
		 * Creates a BaseLabel object with an name from one object. The object must
		 * implement hashCode() and equals().
		 *
		 * @param name
		 * @param obj
		 */
		public BaseLabel(final String name, final H obj) {
			this(name, obj, null, null);
		}

		/**
		 * Creates a BaseLabel object with an name from two objects. The objects must
		 * implement hashCode() and equals().
		 *
		 * @param name
		 * @param obj1
		 * @param obj2
		 */
		public BaseLabel(final String name, final H obj1, final I obj2) {
			this(name, obj1, obj2, null);
		}

		/**
		 * Creates a BaseLabel object with an name from three objects. The objects must
		 * implement hashCode() and equals().
		 *
		 * @param name
		 * @param obj1
		 * @param obj2
		 * @param obj3
		 */
		public BaseLabel(final String name, final H obj1, final I obj2, final J obj3) {
			super(name);
			assert obj1 != null;
			mObject1 = obj1;
			mObject2 = obj2;
			assert !((obj2 == null) && (obj3 != null));
			mObject3 = obj3;
			// Calculate this once to speed access
			mHashCode = super.hashCode() ^ Objects.hash(mObject1, mObject2, mObject3);

		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			final BaseLabel<?, ?, ?> other = (BaseLabel<?, ?, ?>) obj;
			return super.equals(obj) && Objects.equals(mObject1, other.mObject1) && //
					Objects.equals(mObject2, other.mObject2) && //
					Objects.equals(mObject3, other.mObject3);
		}

		@Override
		public int hashCode() {
			return mHashCode;
		}

		@Override
		public String toHTML(final Mode mode) {
			final StringBuffer sb = new StringBuffer();
			sb.append(super.toString());
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

		@Override
		public String toString() {
			return HTML.stripTags(toHTML(Mode.TERSE));
		}

		protected H getObject1() {
			assert mObject1 != null;
			return mObject1;
		}

		protected I getObject2() {
			assert mObject1 != null;
			assert mObject2 != null;
			return mObject2;
		}

		protected J getObject3() {
			assert mObject1 != null;
			assert mObject2 != null;
			assert mObject3 != null;
			return mObject3;
		}

	}

	public static class MaterialMAC //
			extends EPMALabel.BaseLabel<Material, XRay, Object> {

		public MaterialMAC(final Material mf, final XRay xr) {
			super("[&mu;/&rho;]", mf, xr);
		}

		public Material getMaterial() {
			return getObject1();
		}

		public XRay getXRay() {
			return getObject2();
		}

	}

}
