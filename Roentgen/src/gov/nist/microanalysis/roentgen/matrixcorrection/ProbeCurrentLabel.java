package gov.nist.microanalysis.roentgen.matrixcorrection;

import java.util.Date;
import java.util.Optional;

import gov.nist.microanalysis.roentgen.math.uncertainty.BaseLabel;

public class ProbeCurrentLabel extends BaseLabel<Integer, Object, Object> {

	private final Optional<Date> mTimestamp;

	protected ProbeCurrentLabel(final int index, final Date timestamp) {
		super("PC<sub>" + index + "</sub>", Integer.valueOf(index));
		mTimestamp = Optional.of(timestamp);
	}

	protected ProbeCurrentLabel(final int index) {
		super("PC<sub>" + index + "</sub>", Integer.valueOf(index));
		mTimestamp = Optional.empty();
	}

	public Optional<Date> getTimestamp() {
		return mTimestamp;
	}

}
