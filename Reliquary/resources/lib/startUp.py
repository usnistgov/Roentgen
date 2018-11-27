# This script is run each time the application starts to define
# the basic Relinquary functionality.

import jarray
import java.lang as _jl
import javax.swing as _jxsw
import java.util as _ju;

import sys 

import com.duckandcover as _dac
import com.duckandcover.html as _dach

hVerbose = _dach.IToHTML.Mode.VERBOSE
hNormal = _dach.IToHTML.Mode.NORMAL
hTerse = _dach.IToHTML.Mode.TERSE

class PyAction(_jxsw.AbstractAction):
    """Defines a Python class the extends javax.swing.AbstractAction using a Python function which takes a single argument, an javax.swing.ActionEvent.
    The constructor takes a name, a function of one argument and an optional javax.swing.Icon"""
    
    def __init__(self, name, func, icon=None):
        """Constructor taking the name of the action, a function implementing the action and an optional javax.swing.Icon"""
        super(PyAction, self).__init__(name, icon)
        self._func = func
        
    def actionPerformed(self, actionEvent):
        self._func(actionEvent)
    
def report(html, p=True):
    """report(html, [p=True])
    Adds text as html to the report document. If p=True wrap the text in <p>...</p>"""
    if not p:
        __worker__.publishHTML(html)
    else:
        __worker__.publishHTML("<p>"+html+"</p>")

def verboseFunc(obj, dest=None):
    """verboseFunc(obj, dest)
    Converts object to HTML using the IToHTML[Ext] interface with the Mode.VERBOSE"""
    return _dach.HTML.toHTML(obj, hVerbose)
    
def normalFunc(obj, dest=None):
    """normalFunc(obj, dest)
    Converts object to HTML using the IToHTML[Ext] interface with the Mode.NORMAL"""
    return _dach.HTML.toHTML(obj, hNormal)

def terseFunc(obj, dest=None):
    """terseFunc(obj, dest)
    Converts object to HTML using the IToHTML[Ext] interface with the Mode.TERSE"""
    return _dach.HTML.toHTML(obj, hTerse)
    
def terse(obj, dest=None):
    if isinstance(obj, ju.Collection):
        for item in obj:
            terse(item, dest)
    else:
        report(terseFunc(obj), dest)

def normal(obj, dest=None):
    if isinstance(obj, ju.Collection):
        for item in obj:
            normal(item, dest)
    else:
        report(normalFunc(obj), dest)

def verbose(obj, dest=None):
    if isinstance(obj, ju.Collection):
        for item in obj:
            verbose(item, dest)
    else:
        report(verboseFunc(obj, dest), dest)
        
def escape(text):
    if isinstance(text, str) or isinstance(text,unicode):
        return _dach.HTML.escape(text)
    else:
        return _dach.HTML.escape(str(text))

def header(text, level=3):
    report("<h%d>%s</h%d>" % (level, escape(text), level))

def paragraph(text):
    report("<p>"+escape(text)+"</p>")


def tabulate(items, keyFunc=terseFunc, valFunc=terseFunc, p=True):
    rows = ( "<tr><td>" + keyFunc(key)+"</td><td>"+valFunc(val)+"</td></tr>" for key, val in dict(items).iteritems())
    report("<table>"+("".join(rows))+"</table>",p)
    
def ul(items, valFunc=terseFunc, p=True):
    report("<ul>"+("".join(("<li>%s</li>" % (valFunc(s),) for s in items)) )+"</ul>", p)
    
def buildReport(name="Report"):
    return _dach.Report(name)
    

# The application 
rApp = _dac.Reliquary.instance()
# The JMainFrame containing the main application window
rMainFrame = rApp.getMainFrame()
report("<h4>Reliquary Scriptable Application Container</h4>")
