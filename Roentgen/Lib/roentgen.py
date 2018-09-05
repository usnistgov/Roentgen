import sys

import java.lang as jl
import gov.nist.microanalysis.roentgen.physics as _rp
import gov.nist.microanalysis.roentgen.physics.composition as _rpc
import org.apache.commons.math3.linear as _ml3

print "Loading roentgen.py"

def element(elm):
    from gov.nist.microanalysis.roentgen.physics import Element
    if isinstance(elm, Element):
        return elm
    else:
        return Element.parse(elm)
    
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
        if isinstance(v, tuple) or isinstance(v, list):
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
    from gov.nist.microanalysis.roentgen.physics import CharacteristicXRay
    return CharacteristicXRay.parse(trStr)
    

def transitions(elm, fam, minWeight=0.001):
    """transitions(elm,fam,[minWeight=0.001])
    Returns a ElementXRaySet containing the transitions associated with the specified element and family.
    Ex: transitions("Fe","K")"""
    from gov.nist.microanalysis.roentgen.physics import XRaySet
    elm = element(elm)
    fam = family(fam)
    return XRaySet.build(elm, fam, minWeight)


print "Roentgen scripting initialized..."