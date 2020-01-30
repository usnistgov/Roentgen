# Roentgen - Roadmap
Roentgen is divided into three parts
* Reliquary - A scripting (Jython) container application
* ReliquaryUtilities - A handful of helper classes and interface for enhanced scripting
* Roentgen - X-ray microanalysis algorithms

Roentgen is a Java library for X-ray microanalysis
* Roentgen provides many core algorithms necessary to perform common X-ray microanalysis tasks.
   * Matrix correction (including uncertainty propagation)
   *
* Roentgen builds on Apache Math3 for linear algebra, integration, statistical, optimization and other numerical algorithms.
  * Uncertainty calculations are fundamental to Roentgen
* Roentgen provides many basic atomic and X-ray physics classes
* Roentgen provides the XPP algorithm for bulk matrix correction
* Roentgen provides classes for manipulating EDS spectrum data
* Roentgen uses a simple mechanism to allow objects to output their internal state to HTML for reporting
* Roentgen uses Jython to provide a scripting environment for interacting with the library
