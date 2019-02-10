package gov.nist.microanalysis.roentgen.tests;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import gov.nist.microanalysis.roentgen.tests.math.TestIntInterval;
import gov.nist.microanalysis.roentgen.tests.math.TestMultiStepNamedMultivariateJacobianFunction;
import gov.nist.microanalysis.roentgen.tests.math.TestNamedMultivariateJacobianFuntion;
import gov.nist.microanalysis.roentgen.tests.math.TestSafeMultivariateNormalDistribution;
import gov.nist.microanalysis.roentgen.tests.math.TestUncertainValue;
import gov.nist.microanalysis.roentgen.tests.math.TestUncertainValues;
import gov.nist.microanalysis.roentgen.tests.math.TestUtility;
import gov.nist.microanalysis.roentgen.tests.matrixcorrection.KRatioCorrectionModelTest;
import gov.nist.microanalysis.roentgen.tests.matrixcorrection.TwoPointContinuumModelTest;
import gov.nist.microanalysis.roentgen.tests.matrixcorrection.XPPMatrixCorrection2Test;
import gov.nist.microanalysis.roentgen.tests.physics.CompositionTest2;
import gov.nist.microanalysis.roentgen.tests.physics.TestMassAbsorptionCoefficient;
import gov.nist.microanalysis.roentgen.tests.spectrum.TestEDSFittingFilter;
import gov.nist.microanalysis.roentgen.tests.spectrum.TestLinearSpectrumFit;
import gov.nist.microanalysis.roentgen.tests.spectrum.TestLineshapeCalibration;
import gov.nist.microanalysis.roentgen.tests.spectrum.TestReadSpectrum;
import gov.nist.microanalysis.roentgen.tests.spectrum.TestSpectrumFileReader;

@RunWith(Suite.class)
@Suite.SuiteClasses({
  TestIntInterval.class,
  TestMultiStepNamedMultivariateJacobianFunction.class,
  TestNamedMultivariateJacobianFuntion.class,
  TestSafeMultivariateNormalDistribution.class,
  TestUncertainValue.class,
  TestUncertainValues.class,
  TestUtility.class,
  TwoPointContinuumModelTest.class,
  TestMassAbsorptionCoefficient.class,
  TestReadSpectrum.class,
  TestEDSFittingFilter.class,
  TestLinearSpectrumFit.class,
  TestLineshapeCalibration.class,
  TestSpectrumFileReader.class,
  CompositionTest2.class,
  XPPMatrixCorrection2Test.class,
  KRatioCorrectionModelTest.class
})

public class TestSuite {

}

