package gov.nist.microanalysis.roentgen.physics;

/**
 * <p>
 * CODATA 2002 values for some common physical constants in SI units. These are
 * based on the table downloaded from
 * http://physics.nist.gov/cuu/Constants/Table/allascii.txt in 2005. They have
 * not been updated to reflect CODATA 2006 or subsequent.
 * </p>
 * <p>
 * The units with postfix CGS are the constants converted to electrostatic CGS
 * (ESU) units.
 * </p>
 * <p>
 * Please do not misinterpret these values as the latest CODATA values. The
 * latest values can be downloaded from
 * http://physics.nist.gov/cuu/Constants/Table/allascii.txt.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 199 $
 */

public class PhysicalConstants {
   /**
    * AvagadroNumber - Avagadro's number (CODATA 2002)
    */
   static public final double AvagadroNumber = 6.0221415e23; // dimensionless
   /**
    * SpeedOfLight - C, the speed of light in m/s (CODATA 2002)
    */
   static public final double SpeedOfLight = 299792458; // m/s
   static public final double SpeedOfLightCGS = 29979245800.0; // cm/s

   /**
    * ElectronCharge - The charge of a single electron in C (CODATA 2002)
    */
   static public final double ElectronCharge = 1.60217653e-19; // C
   static public final double ElectronChargeCGS = 4.8032068e-10; // C
   /**
    * PlanckConstant - Planck's constant in J s (CODATA 2002)
    */
   static public final double PlanckConstant = 6.6260693e-34; // J s
   static public final double PlanckConstantCGS = 6.6260693e-27; // erg s
   /**
    * PlanckReduced - hBar, Plank's constant over 2pi in J s (CODATA 2002)
    */
   static public final double PlanckReduced = 1.05457266e-34; // J s
   static public final double PlanckReducedCGS = PlanckConstantCGS / (2.0 * Math.PI); // erg
   // s
   /**
    * ElectronMass - The electron mass in kg (CODATA 2002)
    */
   static public final double ElectronMass = 9.1093826e-31; // kg
   static public final double ElectronMassCGS = 9.1093826e-28; // g

   /**
    * ElectronRestMass - The energy equivalent mass of an electron in Joules
    * (CODATA 2002)
    */
   static public final double ElectronRestMass = 8.1871047e-14; // J
   static public final double ElectronRestMassGGS = ElectronMassCGS * SpeedOfLightCGS * SpeedOfLightCGS; // ergs
   static public final double ElectronRestMass_eV = ElectronRestMass / ElectronCharge; // eV
   /**
    * ProtonMass - The proton mass in kg (CODATA 2002)
    */
   static public final double ProtonMass = 1.66053886e-27; // kg
   static public final double ProtonMassCGS = 1.66053886e-24; // g
   /**
    * NeutronMass - The neutron mass in kg (CODATA 2002)
    */
   static public final double NeutronMass = 1.67492728e-27; // kg
   static public final double NeutronMassCGS = 1.67492728e-24; // g
   /**
    * UnifiedAtomicMass - 1 AMU in kg (CODATA 2002)
    */
   static public final double UnifiedAtomicMass = 1.66053886e-27; // kg
   static public final double UnifiedAtomicMassCGS = 1.66053886e-24; // g
   /**
    * PermittivityOfFreeSpace - Permittivity of free space in F/m (CODATA 2002)
    */
   static public final double PermittivityOfFreeSpace = 8.8541878176203898505365630317108e-12; // F/m
   static public final double PermittivityOfFreeSpaceCGS = 1.0 / (4.0 * Math.PI);
   /**
    * PermeabilityOfFreeSpace - Permeability of free space in N/(A^2) (CODATA
    * 2002)
    */
   static public final double PermeabilityOfFreeSpace = 12.566370614359172953850573533118e-7; // N/(A^2)
   static public final double PermeabilityOfFreeSpaceCGS = (4.0 * Math.PI) / Math.pow(SpeedOfLightCGS, 2.0);
   /**
    * BoltzmannConstant - kB, Boltzmann's constant in J/K (CODATA 2002)
    */
   static public final double BoltzmannConstant = 1.3806505e-23; // J/K
   static public final double BoltzmannConstantCGS = 1.3806505e-16; // erg/K

   /**
    * GravitationConstant - G in kg m^2 (CODATA 2002)
    */
   static public final double GravitationConstant = 6.6742e-11; // m^3/(kg s^2)
   static public final double GravitationConstantCGS = 6.6742e-8; // cm^3/(g
   // s^2)
   /**
    * PlanckLength - The Plank length in m (CODATA 2002)
    */
   static public final double PlanckLength = 1.61624e-35; // m
   static public final double PlanckLengthCGS = 1.61624e-33; // cm
   /**
    * PlanckMass - The Plank mass in kg (CODATA 2002)
    */
   static public final double PlanckMass = 2.17645e-8; // kg
   static public final double PlanckMassCGS = 2.17645e-5; // g
   /**
    * PlanckTemperature - The Planck temperature (dimensionless) (CODATA 2002)
    */
   static public final double PlanckTemperature = 1.41679e32; // Dimensionless
   /**
    * PlanckTime - The Planck time in seconds (CODATA 2002)
    */
   static public final double PlanckTime = 5.39121e-44; // s

   /**
    * RydbergEnergy - The Rydberg constant in energy units (Joules) (CODATA
    * 2002)
    */
   static public final double RydbergEnergy = 2.17987209e-18; // Joules
   static public final double RydbergEnergyCGS = 2.17987209e-11; // ergs

   /**
    * BohrRadius - The Bohr radius in meters. (CODATA 2002)
    */
   static public final double BohrRadius = 5.291772083e-11; // Meter
   static public final double BohrRadiusCGS = 5.291772083e-9; // cm

   /**
    * FineStructure - The atomic fine structure constant (dimensionless) (CODATA
    * 2002)
    */
   static public final double FineStructure = 7.297352568e-3; // dimensionless

   /**
    * ClassicalElectronRadius - The classical electron radius in meters (CODATA
    * 2002)
    */
   static public final double ClassicalElectronRadius = 2.817940325e-15; // m
   static public final double ClassicalElectronRadiusCGS = 2.817940325e-13; // cm

   /**
    * IcePoint - The temperature at which water freezes at StandardAtmosphere
    * (CODATA 2002)
    */
   static public final double IcePoint = 273.15; // Kelvin
   /**
    * StandardAtmosphere - The standard atmospheric pressure (CODATA 2002)
    */
   static public final double StandardAtmosphere = 101325.0; // Pascal

}

/*******************************************************************************
 * Fundamental Physical Constants --- Complete Listing From:
 * http://physics.nist.gov/constants Source: Peter J. Mohr and Barry N. Taylor,
 * CODATA Recommended Values of the Fundamental Physical Constants: 2002
 *******************************************************************************/
