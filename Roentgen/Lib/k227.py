# In the 2012 paper, we examine the matrix correction for O K-L3 in K227 relative
# to some common stoichiometric standards and relative to K229, a similar glass.
# Author:  Nicholas Ritchie  3-Sep-2018

import gov.nist.microanalysis.roentgen.matrixcorrection as rmc
import gov.nist.microanalysis.roentgen.math.uncertainty as runc

import com.duckandcover.html.Table as table

std1 = material("Al2O3").asMassFraction()
std2 = material("MgO").asMassFraction()
std3 = material("SiO2").asMassFraction()
std4 = massFraction("K229", { "O":(0.2099,0.002099), "Si":(0.1402,0.001402), "Pb":(0.6498,0.006498) })

unk = massFraction("K227", { "O":(0.1639,0.001639), "Si":(0.0935,0.000935), "Pb":(0.7427,0.007427)})

def uncertainValue(val, sigma=0.0):
    return runc.UncertainValue(val, sigma)

def mcDatum(comp, std, e0, toa=uncertainValue(40.0,0.1), roughness=0.0):
    return rmc.MatrixCorrectionDatum(comp, std, e0, toa, roughness)

cxr = transition("O K-L3")

toa = uncertainValue(40.0, 0.0)

for e0 in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0):
    mcStd1 = mcDatum(std1, True, uncertainValue(e0, 0.0), toa, 0.0)
    mcStd2 = mcDatum(std2, True, uncertainValue(e0, 0.0), toa, 0.0)
    mcStd3 = mcDatum(std3, True, uncertainValue(e0, 0.0), toa, 0.0)
    mcStd4 = mcDatum(std4, True, uncertainValue(e0, 0.0), toa, 0.0)
    
    mcUnk = mcDatum(unk, False, uncertainValue(e0, 0.0), toa, 0.0)
    
    t = table()
    
    for mcStd in (mcStd1, mcStd2, mcStd3, mcStd4):
        xpp = rmc.XPPMatrixCorrection(mcUnk, mcStd, cxr, rmc.XPPMatrixCorrection.defaultVariates())
        results = runc.UncertainValues.propagate(xpp, xpp.buildInput())
        tagZA = rmc.MatrixCorrectionTag(mcUnk, mcStd, cxr)
        val = results.getUncertainValue(tagZA)
        comps = val.getComponentNames()
        for comp in comps:
            t.addRow(table.td(str(e0)), table.td(mcStd.getComposition()), table.td(tagZA), table.td(comp), table.td(val.getComponent(comp)))
    terse(t)