import sys

import java.lang as _jl
import java.util as _ju
import gov.nist.microanalysis.roentgen.physics as _rp
import gov.nist.microanalysis.roentgen.physics.composition as _rpc
import gov.nist.microanalysis.roentgen.math.uncertainty as _runc
import gov.nist.microanalysis.roentgen.matrixcorrection as _rmc

import org.apache.commons.math3.linear as _ml3

print "Loading roentgen.py"

def element(elm):
    if isinstance(elm, _rp.Element):
        return elm
    else:
        return _rp.Element.parse(elm)
    
def family(fam):
    """family(fam)
    Returns the Shell.Principle object associated with the specified family - 'K', 'L', 'M' or 'N' """
    if isinstance(fam, _rp.Shell.Principle):
        return fam
    elif fam=="K":
        return _rp.Shell.Principle.K
    elif fam=="L":
        return _rp.Shell.Principle.L
    elif fam=="M":
        return _rp.Shell.Principle.M
    elif fam=="N":
        return _rp.Shell.Principle.N
    else:
        return None

def massFraction(name, elms):
    """Ex: massFraction("Other",elms = { "Fe": (0.3,0.01), "Mg":(0.2,0.02), "O": 0.5 })
    Constructs a Composition object representing a mass fraction with the specified name and composition.
    elms = { "Fe": (0.3,0.01), "Mg":(0.2,0.02), "O": 0.5 } represents a material with 30+-1% iron, 20+-2% magnesium and 50+-0% oxygen."""
    ell = []
    vals = _ml3.ArrayRealVector(len(elms))
    vars = _ml3.ArrayRealVector(len(elms))
    for i, (elm, v) in enumerate(elms.iteritems()):
        ell.append(_rp.Element.parse(elm))
        if isinstance(v, float):
            vals.setEntry(i, v)
        elif isinstance(v, _runc.UncertainValue):
            vals.setEntry(i,v.doubleValue())
            vars.setEntry(i,v.variance())
        elif isinstance(v, tuple) or isinstance(v, list):
            vals.setEntry(i, v[0])
            if len(v)>1:
                vars.setEntry(i, v[1]*v[1])
        else:
            vals.setEntry(i, v)
    return _rpc.Composition.massFraction(name, ell, vals, vars)
    

def material(chemForm):
    """material(chemForm)
    Create a Composition object representing the chemical formula in chemForm"""    
    if isinstance(chemForm,_rpc.Composition):
        return chemForm
    else:
        return _rpc.Composition.parse(chemForm)


def mcd(mat, std, e0, toa, roughness=0.0):
    """mcd(mat, std, e0, toa, roughness=0.0):
    Create a MatrixCorrectionDatum to represent
    mat: the material
    std=[True,False]: is standard
    e0: beam energy in keV
    toa: take-off angle in degrees
    roughness: in mass thickness g/cm^2"""
    if isinstance(e0, tuple):
        e0= ( uv(e0[0],e0[1]) if len(e0)>1 else uv(e0[0]))
    if isinstance(e0, float):
        e0 = uv(e0)
    if isinstance(toa, tuple):
        toa = (uv(toa[0], toa[1]) if len(toa)>1 else uv(toa[0]))
    if isinstance(toa, float):
        toa = uv(toa)
    return _rmc.MatrixCorrectionDatum(material(mat), True, e0, toa.multiply(_jl.Math.PI/180.), roughness)

def xpp(unk, stds):
    """xpp(unk, stds)
    Returns an initialized LabeledMultivariateJacobianFunction for evaluating the XPP matrix correction model relative
    to the MatrixCorrectionDatum objects in unk and stds.
    unk: a MatrixCorrectionDatum object
    stds: A dictionary of ElementXRaySet to a tuple of ( MatrixCorrectionDatum, UncertainValue(k-ratio) ) object"""
    steps = ju.ArrayList()
    stdMcds, krs = _ju.HashMap(), _ju.HashMap()
    for xrts, ( mcd, kr) in stds.iteritems():
        if isinstance(kr, tuple):
            kr=uv(*kr)
        stdMcds.put(xrts, mcd)
        krs.put( _rmc.KRatioLabel(unk, mcd, xrts, _rmc.KRatioLabel.Method.Measured) , kr)
    xpp = _rmc.XPPMatrixCorrection(unkMcd, stdMcds)
    steps.add(xpp)
    steps.add(_rmc.KRatioCorrectionModel(unkMcd, stdMcds))
    res = _runc.SerialLabeledMultivariateJacobianFunction("Full K-ratio correction", steps)
    res.initializeConstants(unk.getComposition().getValueMap())
    xppInputs = xpp.buildInput()
    fullInputs = runc.UncertainValues.combine(xppInputs, _runc.UncertainValues(krs))
    return (res, fullInputs)


def transition(trStr):
    if isinstance(trStr, _rp.CharacteristicXRay):
        return trStr
    else:
        return _rp.CharacteristicXRay.parse(trStr)

def transitionSet(trs):
    """transitionSet(trs)
    Converts trs into an ElementXRaySet associated with the element elm.
    Example: transitionSet(( "Fe K-L3","Fe K-L2" )) 
    Example: transitionSet(( "Fe K-L3","Ni K-L2" ))"""
    elm=None
    sameElm=True
    if isinstance(trs, ( list, tuple ) ):
        tmp = []
        for tr in trs:
            cxr=transition(tr)
            tmp.append(cxr)
            elm = (cxr.getElement() if not elm else elm)
            if elm<>cxr.getElement():
                sameElm = False
        if sameElm:
            return _rp.XRaySet.ElementXRaySet(tmp)
        else:
            return _rp.XRaySet.CharacteristicXRaySet.build(tmp)
    else:
        return _rp.XRaySet.ElementXRaySet(transition(trs))

def transitions(elm, fam, minWeight=0.001):
    """transitions(elm,fam,[minWeight=0.001])
    Returns a ElementXRaySet containing the transitions associated with the specified element and family.
    Ex: transitions("Fe","K")"""
    return _rp.XRaySet.build(element(elm), family(fam), minWeight)

def uv(val, unc=0.0):
    """uv(x,dx) or uv(x,[0.0])
    Build an uncertain value x plus/minus dx"""
    return _runc.UncertainValue(val, unc)

print "Roentgen scripting initialized..."