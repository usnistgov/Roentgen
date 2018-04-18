package com.duckandcover;

import java.awt.Component;
import java.awt.Desktop;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;

import javax.swing.AbstractAction;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import com.duckandcover.html.IToHTML;
import com.duckandcover.html.IToHTMLExt;
import com.duckandcover.scripting.JMenuBarBuilder;
import com.duckandcover.scripting.JMenuBuilder;
import com.duckandcover.scripting.JPythonPanel;
import com.duckandcover.scripting.JScriptingPreferencePanel;
import com.duckandcover.swing.JErrorDialog;
import com.duckandcover.swing.JPreferenceDialog;
import com.duckandcover.swing.JPreferencePanel;

/**
 * <p>
 * The primary window for the Relinquary scripting container.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2017
 * </p>
 *
 * @author Nicholas
 * @version $Rev: 312 $
 */
public class JReliquaryFrame extends JFrame {

    private static final long serialVersionUID = 5161432472757718511L;

    private final JSplitPane jSplitPane_Main = new JSplitPane(JSplitPane.VERTICAL_SPLIT, true);
    private final JTabbedPane jTabbedPane_Utility = new JTabbedPane(JTabbedPane.TOP, JTabbedPane.WRAP_TAB_LAYOUT);
    private final JPythonPanel jPythonPanel = new JPythonPanel();

    private final JMenuBarBuilder jMenuBar_MainMenu = buildMainMenu();

    private final List<JPreferencePanel> mPreferencePanels = new ArrayList<>();

    private static final String MAIN_FRAME_POSITION = "MainFrame.Position";
    private static final String MAIN_SPLITTER_POSITION = "MainFrame.Splitter.Divider";

    private void storePreferences() {
	Reliquary.getPreferences().putInt(MAIN_FRAME_POSITION + ".X", getBounds().x);
	Reliquary.getPreferences().putInt(MAIN_FRAME_POSITION + ".Y", getBounds().y);
	Reliquary.getPreferences().putInt(MAIN_FRAME_POSITION + ".W", getBounds().width);
	Reliquary.getPreferences().putInt(MAIN_FRAME_POSITION + ".H", getBounds().height);
	Reliquary.getPreferences().putInt(MAIN_SPLITTER_POSITION, jSplitPane_Main.getDividerLocation());
    }

    private void restorePreferences() {
	final Rectangle rect = new Rectangle(0, 0, 600, 400);
	rect.x = Reliquary.getPreferences().getInt(MAIN_FRAME_POSITION + ".X", rect.x);
	rect.y = Reliquary.getPreferences().getInt(MAIN_FRAME_POSITION + ".Y", rect.y);
	rect.width = Reliquary.getPreferences().getInt(MAIN_FRAME_POSITION + ".W", rect.width);
	rect.height = Reliquary.getPreferences().getInt(MAIN_FRAME_POSITION + ".H", rect.height);
	setBounds(rect);
	final int pos = Reliquary.getPreferences().getInt(MAIN_SPLITTER_POSITION, Integer.valueOf(200));
	jSplitPane_Main.setDividerLocation(pos);
    }

    public JReliquaryFrame() {
	super(Reliquary.APP_NAME + " - " + Reliquary.SLOGAN);
	init();
	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	setResizable(true);
	addWindowListener(new WindowAdapter() {

	    /**
	     * Records the window positions to be restored on restart.
	     *
	     * @param arg0
	     * @see java.awt.event.WindowAdapter#windowClosing(java.awt.event.
	     *      WindowEvent)
	     */
	    @Override
	    public void windowClosing(final WindowEvent arg0) {
		storePreferences();
		Reliquary.getLogger().log(Level.INFO, "Terminating application");
	    }
	});
	restorePreferences();
	addPreferencePanel(new JScriptingPreferencePanel());
	jPythonPanel.loadHTML(Reliquary.getReport());
	final String startUp = JScriptingPreferencePanel.getStartup();
	if (startUp != null) {
	    final File startUpFile = new File(startUp);
	    if (startUpFile.isFile())
		jPythonPanel.setStartupScript(startUpFile);
	}
	jPythonPanel.start();

    }

    private final class ReportOpenInBrowserAction extends AbstractAction {

	public ReportOpenInBrowserAction() {
	    super("Open in Browser");
	}

	private static final long serialVersionUID = 4937508420715721641L;

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    try {
		Desktop.getDesktop().browse(Reliquary.getReport().toURI());
	    } catch (final IOException e) {
		Reliquary.getLogger().log(Level.SEVERE, "Opening report in browser", e);
	    }
	}
    }

    private final class ExitAction extends AbstractAction {

	private static final long serialVersionUID = -9184210420592068241L;

	public ExitAction() {
	    super("Exit");
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    JReliquaryFrame.this.setVisible(false);
	}
    };

    private final class CutAction extends AbstractAction {

	private static final long serialVersionUID = -9184210420592068241L;

	public CutAction() {
	    super("Cut");
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    // MainFrame.this.setVisible(false);
	}
    };

    private final class CopyAction extends AbstractAction {

	private static final long serialVersionUID = -9184210420592068241L;

	public CopyAction() {
	    super("Copy");
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    // MainFrame.this.setVisible(false);
	}
    };

    private final class PasteAction extends AbstractAction {

	private static final long serialVersionUID = -9184210420592068241L;

	public PasteAction() {
	    super("Paste");
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    // MainFrame.this.setVisible(false);
	}
    };

    private final class AboutAction extends AbstractAction {

	private static final long serialVersionUID = -4463170127334246860L;

	private AboutAction() {
	    super("About...");
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    final ImageIcon ii = createImageIcon("Images/Caution_XRays128.png");
	    final String msg = "<html><h2>" + Reliquary.APP_NAME + "</h2>" + Reliquary.SLOGAN + "</html>";
	    JOptionPane.showMessageDialog(JReliquaryFrame.this, msg, "About " + Reliquary.APP_NAME,
		    JOptionPane.INFORMATION_MESSAGE, ii);
	}
    }

    private final class WebSiteAction extends AbstractAction {

	private static final long serialVersionUID = -4463170127334246860L;

	private final String mURL;

	private WebSiteAction(final String menuItem, final String site) {
	    super(menuItem);
	    mURL = site;
	}

	@Override
	public void actionPerformed(final ActionEvent arg0) {
	    try {
		final URL url = new URL(mURL);
		Desktop.getDesktop().browse(url.toURI());
	    } catch (final Exception e) {
		Reliquary.getLogger().log(Level.SEVERE, "Unable to open website: " + mURL.toString(), e);
	    }
	}
    }

    private final class PreferencesAction extends AbstractAction {

	private static final long serialVersionUID = 3386187653218227322L;

	public PreferencesAction() {
	    super("Preferences...");
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
	    final JPreferenceDialog pd = new JPreferenceDialog(JReliquaryFrame.this);
	    for (final JPreferencePanel pp : mPreferencePanels)
		pd.addPanel(pp);
	    pd.setLocationRelativeTo(JReliquaryFrame.this);
	    pd.setVisible(true);
	}

    }

    /** Returns an ImageIcon, or null if the path was invalid. */
    private ImageIcon createImageIcon(final String name) {
	final java.net.URL imgURL = getClass().getResource(name);
	return imgURL != null ? new ImageIcon(imgURL, Reliquary.APP_NAME) : null;
    }

    private JMenuBarBuilder buildMainMenu() {
	final JMenuBarBuilder res = new JMenuBarBuilder();

	final JMenuBuilder file = res.findOrAdd("File", 10);
	file.add(10, new ReportOpenInBrowserAction());
	file.add(20, jPythonPanel.new OpenAction("Open script..."));
	file.add(128, new PreferencesAction());
	file.add(256, new ExitAction());

	final JMenuBuilder edit = res.findOrAdd("Edit", 20);
	edit.add(10, new CutAction());
	edit.add(20, new CopyAction());
	edit.add(30, new PasteAction());

	final JMenuBuilder help = res.findOrAdd("Help", 256);
	help.add(10, new WebSiteAction("Learn about Python 2.7...", "https://docs.python.org/2.7/"));
	help.add(20, new WebSiteAction("Learn about Jython...", "http://www.jython.org/jythonbook/en/1.0/index.html"));
	help.add(256, new AboutAction());
	return res;
    }

    public void addTab(final String name, final Component component) {
	jTabbedPane_Utility.addTab(name, component);
    }

    public void setDisplay(final Component component) {
	jSplitPane_Main.setTopComponent(component);

    }

    private void init() {
	jSplitPane_Main.setBottomComponent(jTabbedPane_Utility);
	jSplitPane_Main.setTopComponent(new JPanel());
	addTab("Report", jPythonPanel);
	add(jSplitPane_Main);
	updateMainMenu();
    }

    public JMenuBarBuilder getMenuBarBuilder() {
	return jMenuBar_MainMenu;
    }

    public void updateMainMenu() {

	final Runnable r = new Runnable() {
	    @Override
	    public void run() {
		setJMenuBar(jMenuBar_MainMenu.build());
		validate();
	    }
	};
	invokeCarefully(r);
    }

    public void invokeCarefully(final Runnable as) {
	if (SwingUtilities.isEventDispatchThread())
	    as.run();
	else
	    try {
		SwingUtilities.invokeAndWait(as);
	    } catch (final Exception e) {
		final Runnable r = new Runnable() {
		    @Override
		    public void run() {
			JErrorDialog.createErrorMessage(JReliquaryFrame.this, Reliquary.APP_NAME, e);
		    }
		};
		invokeCarefully(r);
	    }
    }

    /**
     * Appends an html encoded string to the report.
     *
     * @param text
     */
    public void appendHTMLToReport(final String html) {
	jPythonPanel.append(html);
	jPythonPanel.flush();
    }

    /**
     * Appends a report generated by an object implementing the IToHTML
     * interface to the report.
     *
     * @param html
     *            {@link IToHTML}
     * @param mode
     *            {@link IToHTML.Mode}
     */
    public void appendToReport(final IToHTML html, final IToHTML.Mode mode) {
	appendHTMLToReport(html.toHTML(mode));
    }

    /**
     * Appends a report generated by an object implementing the IToHTMLExt
     * interface to the report.
     *
     * @param html
     *            {@link IToHTMLExt}
     * @param mode
     *            {@link IToHTML.Mode}
     * @param dir
     *            Name of sub-directory into which to copy associated data
     */
    public void appendToReport(final IToHTMLExt html, final IToHTML.Mode mode, final String dir) throws IOException {
	appendHTMLToReport(html.toHTML(mode, Reliquary.getReport(), dir));
    }

    /**
     * Gets the current value assigned to preferencePanels
     *
     * @return Returns the preferencePanels.
     */
    public List<JPreferencePanel> getPreferencePanels() {
	return Collections.unmodifiableList(mPreferencePanels);
    }

    /**
     * Add a JPreferencePanel to the default system preference dialog.
     *
     * @param jpp
     */
    public void addPreferencePanel(final JPreferencePanel jpp) {
	mPreferencePanels.add(jpp);
    }
}
