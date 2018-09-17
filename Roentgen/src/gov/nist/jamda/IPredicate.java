package gov.nist.jamda;

public interface IPredicate<J> {

	/**
	 * Evaluates J to determine whether the object meets the condition (return true)
	 * or doesn't (return false).
	 * 
	 * @param obj
	 * @return true if meets condition, false otherwise
	 */
	public boolean evaluate(J obj);

}
