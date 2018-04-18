# This script is run each time the application starts to define
# the basic Relinquary functionality.

import jarray
import java.lang as jl
import javax.swing as jxsw

import sys 

import com.duckandcover as dac
import com.duckandcover.html as dach

class PyAction(jxsw.AbstractAction):
    """Defines a Python class the extends javax.swing.AbstractAction using a Python function which takes a single argument, an javax.swing.ActionEvent.
    The constructor takes a name, a function of one argument and an optional javax.swing.Icon"""
    
    def __init__(self, name, func, icon=None):
        """Constructor taking the name of the action, a function implementing the action and an optional javax.swing.Icon"""
        super(PyAction, self).__init__(name, icon)
        self._func = func
        
    def actionPerformed(self, actionEvent):
        self._func(actionEvent)
    
def report(html):
    __worker__.publishHTML(html)
    
def verbose(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(dach.IToHTML.Mode.VERBOSE, rApp.getReport(), dest))            
        else:
            report(obj.toHTML(dach.IToHTML.Mode.VERBOSE))
    else:
        print str(obj)

def normal(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(dach.IToHTML.Mode.NORMAL, rApp.getReport(), dest))            
        else:
            report(obj.toHTML(dach.IToHTML.Mode.NORMAL))
    else:
        print str(obj)

def terse(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(dach.IToHTML.Mode.TERSE, rApp.getReport(), dest))
        else:
            report(obj.toHTML(dach.IToHTML.Mode.TERSE))
    else:
        print str(obj)

# The application 
rApp = dac.Reliquary.instance()
# The JMainFrame containing the main application window
rMainFrame = rApp.getMainFrame()
report("<h4>Reliquary Scriptable Application Container</h4>")
