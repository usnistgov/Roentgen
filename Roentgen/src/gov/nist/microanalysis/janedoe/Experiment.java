package gov.nist.microanalysis.janedoe;

import org.apache.commons.math3.geometry.euclidean.oned.Interval;

import com.duckandcover.html.IToHTML;

/**
 * @author Nicholas
 */
public class Experiment implements IToHTML {

	public Experiment(final String name, final String measuredVariable) {

	}

	public Experiment addDiscreteFactor(final Object tag, final String[] levels) {

		return this;
	}

	public Experiment addContinuousFactor(final Object tag, final Interval iv) {

		return this;
	}

	@Override
	public String toHTML(final Mode mode) {
		// TODO Auto-generated method stub
		return null;
	}

}
