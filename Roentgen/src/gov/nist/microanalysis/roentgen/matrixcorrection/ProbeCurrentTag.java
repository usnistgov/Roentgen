package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Date;
import java.util.Optional;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseTag;

public class ProbeCurrentTag extends BaseTag<Integer, Object, Object> {

	private final Optional<Date> mTimestamp;

	protected ProbeCurrentTag(int index, Date timestamp) {
		super("PC<sub>" + index + "</sub>", Integer.valueOf(index));
		mTimestamp = Optional.of(timestamp);
	}

	protected ProbeCurrentTag(int index) {
		super("PC<sub>" + index + "</sub>", Integer.valueOf(index));
		mTimestamp = Optional.empty();
	}

	public Optional<Date> getTimestamp() {
		return mTimestamp;
	}

}
