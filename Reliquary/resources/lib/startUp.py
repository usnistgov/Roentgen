# This script is run each time the application starts to define
# the basic Relinquary functionality.

import jarray
import java.lang as jl
import javax.swing as jxsw

import sys 

import com.duckandcover as dac
import com.duckandcover.html as dach

hVerbose = dach.IToHTML.Mode.VERBOSE
hNormal = dach.IToHTML.Mode.NORMAL
hTerse = dach.IToHTML.Mode.TERSE

class PyAction(jxsw.AbstractAction):
    """Defines a Python class the extends javax.swing.AbstractAction using a Python function which takes a single argument, an javax.swing.ActionEvent.
    The constructor takes a name, a function of one argument and an optional javax.swing.Icon"""
    
    def __init__(self, name, func, icon=None):
        """Constructor taking the name of the action, a function implementing the action and an optional javax.swing.Icon"""
        super(PyAction, self).__init__(name, icon)
        self._func = func
        
    def actionPerformed(self, actionEvent):
        self._func(actionEvent)
    
def report(html, p=True):
    """
    report(html, [p=True])
    Adds text as html to the report document. If p=True wrap the text in <p>...</p>"""
    if not p:
        __worker__.publishHTML(html)
    else:
        __worker__.publishHTML("<p>"+html+"</p>")

def verbose(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(hVerbose, rApp.getReport(), dest))            
        else:
            report(obj.toHTML(hVerbose))
    else:
        print str(obj)

def normal(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(hNormal, rApp.getReport(), dest))            
        else:
            report(obj.toHTML(hNormal))
    else:
        print str(obj)

def terse(obj, dest=None):
    if isinstance(obj, dach.IToHTML):
        if dest and isinstance(obj, dach.IToHTMLExt):
            report(obj.toHTML(hTerse, rApp.getReport(), dest))
        else:
            report(obj.toHTML(hTerse))
    else:
        print str(obj)

# The application 
rApp = dac.Reliquary.instance()
# The JMainFrame containing the main application window
rMainFrame = rApp.getMainFrame()
report("<h4>Reliquary Scriptable Application Container</h4>")
