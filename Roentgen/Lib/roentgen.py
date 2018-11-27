import sys

import java.lang as _jl
import java.util as _ju
import gov.nist.microanalysis.roentgen.physics as _rp
import gov.nist.microanalysis.roentgen.physics.composition as _rpc
import gov.nist.microanalysis.roentgen.math.uncertainty as _runc

import org.apache.commons.math3.linear as _ml3

print "Loading roentgen.py"

def element(elm):
    if isinstance(elm, Element):
        return elm
    else:
        return _rp.Element.parse(elm)
    
def family(fam):
    """family(fam)
    Returns the Shell.Principle object associated with the specified family - 'K', 'L', 'M' or 'N' """
    from gov.nist.microanalysis.roentgen.physics import Shell
    if isinstance(fam, Shell.Principle):
        return fam
    elif fam=="K":
        return Shell.Principle.K
    elif fam=="L":
        return Shell.Principle.L
    elif fam=="M":
        return Shell.Principle.M
    elif fam=="N":
        return Shell.Principle.N
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
    return _rpc.Composition.parse(chemForm)

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