package gov.nist.microanalysis.roentgen.physics;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * <p>
 * Description
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author nritchie
 * @version $Rev: $
 */
public class BremsstrahlungXRay extends XRay {

	private final Vector3D mElectronDirection;

	/**
	 * Constructs a BremsstrahlungXRay
	 *
	 * @param energy
	 */
	public BremsstrahlungXRay(final double energy, final Vector3D electronDir) {
		super(energy);
		mElectronDirection = electronDir;
	}

	public Vector3D getElectronDirection() {
		return mElectronDirection;
	}

	@Override
	public double angularDependence(final Vector3D photonDirection) {
		return 1.0;
	}

}
