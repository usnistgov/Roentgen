package com.duckandcover.lazy;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * <p>
 * Used to delay the computation of a computationally expensive object until
 * necessary. Derived classes compute the object in the function compute(). Use
 * this instead of the similar Apache
 * {@link org.apache.commons.lang3.concurrent.LazyInitializer} to avoid the
 * potential for throwing ConcurrentExceptions.
 * </p>
 * <p>
 * Avoid using SimplyLazy for computing anything that may throw an uncaught
 * exception as this may lead to SimplyLazy repeatedly attempting to recompute
 * the quantity and failing again and again.
 * </p>
 * <p>
 * <h4>Uses</h4>
 * <ul>
 * <li>SimplyLazy can make programs more responsive particularly at startup if
 * calculations can be delayed as long as possible.</li>
 * <li>SimplyLazy can also facilitate resources shared between threads as the
 * first thread will create the object and subsequent threads will wait until
 * the first thread completes. This will temporarily block threads leading to a
 * short-term inefficiency.</li>
 * <li>SimplyLazy can help with serialization as it can facilitate initializing
 * computed objects after the main object has been loaded from disk.</li>
 * </ul>
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2017
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 309 $
 */
public abstract class SimplyLazy<H> {

	@XStreamOmitField // Make certain this is not saved as it can be computed
	private transient volatile H mValue;

	/**
	 * Constructs a LazyEvaluate which will evaluate compute when an instance of
	 * mValue is required.
	 */
	public SimplyLazy() {
		mValue = null;
	}

	/**
	 * Constructs a SimplyLazy in which the value of mValue is already assigned.
	 *
	 * @param value
	 */
	public SimplyLazy(
			final H value
	) {
		mValue = value;
	}

	/**
	 * Clears the cached value. If get() is called again, the cached value will be
	 * recomputed via a call to compute().
	 */
	public void reset() {
		synchronized (this) {
			mValue = null;
		}
	}

	/**
	 * get() will return the cached value when one exists or call compute to assign
	 * the cached value and then return it.
	 *
	 * @return H An instance of the object
	 */
	public H get() {
		H res = mValue;
		if (res == null)
			synchronized (this) {
				res = mValue;
				if (res == null) {
					res = initialize();
					mValue = res;
				}
			}
		return res;
	}

	public boolean initialized() {
		synchronized (this) {
			return mValue != null;
		}
	}

	/**
	 * Implement this function to compute the value that will be returned by get()
	 *
	 * @return H
	 */
	abstract protected H initialize();

	public String toString() {
		return "Lazy[" + get().toString() + "]";
	}

}
