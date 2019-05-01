/**
 * <p>
 * JUncertainty is a Java library for performing uncertainty calculations on
 * multi-variate measurement models. The strategies implemented are based on the
 * ISO Joint Committee for Guides in Metrology's series of documents JCGM
 * 100:2008, JCGM 101:2008, JCGM 102:2011 which can be found
 * <a href="https://www.bipm.org/en/publications/guides/">here</a>.
 * </p>
 * <p>
 * The general strategy is to propagate uncertainties on input parameters
 * through a measurement model to the output parameters using one of two
 * mechanisms:
 * </p>
 * <ol>
 * 	<li>Computing the Jacobian of the model and performing a Taylor series
 * 		approximation.</li>
 * 	<li>Using a Monte Carlo approach to sample the PDF of the input variables
 * 		propagate these through the model to estimate the PDF of the output
 * 		variables.</li>
 * </ol>
 * <p>
 * Each unique variable (input or output) in a model is identified by a label.
 * The labels allow values to be tracked through the calculation. Collections of
 * variable values and the associated uncertainties are represented by an
 * instance of a class derived from
 * {@link gov.nist.juncertainty.UncertainValuesBase}. Measurement models are
 * represented by instances of classes derived from
 * {@link gov.nist.juncertainty.ExplicitMeasurementModel}. An extension of
 * {@link gov.nist.juncertainty.ExplicitMeasurementModel}
 * ({@link gov.nist.juncertainty.CompositeMeasurementModel}) is provided to
 * facilitate breaking complex measurement models into a series of simpler
 * steps.
 * </p>
 * <p>
 * The library builds on the <a href=
 * "https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html">Apache
 * Commons Math library</a>.
 * </p>
 * 
 * @author Nicholas W. M. Ritchie
 *
 */
package gov.nist.juncertainty;