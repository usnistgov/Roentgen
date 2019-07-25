package com.duckandcover.scripting;

import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.prefs.Preferences;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

import com.duckandcover.Reliquary;
import com.duckandcover.swing.JPreferencePanel;
import com.duckandcover.swing.SwingUtils;
import com.jgoodies.forms.builder.PanelBuilder;
import com.jgoodies.forms.factories.CC;
import com.jgoodies.forms.layout.FormLayout;

/**
 * <p>
 * A panel for configuring scripting related preference values.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class JScriptingPreferencePanel extends JPreferencePanel {

	private static final long serialVersionUID = 1404283203291466195L;

	private final JTextField jTextField_Startup = new JTextField();
	private final JButton jButton_Startup = new JButton(new SelectStartup());
	private final JTextField jTextField_Site = new JTextField();

	private class SelectStartup extends AbstractAction {

		private static final long serialVersionUID = 6759509380749485599L;

		private SelectStartup() {
			super("Select...");
		}

		@Override
		public void actionPerformed(
				ActionEvent arg0
		) {
			final JFileChooser jfc = new JFileChooser();
			jfc.addChoosableFileFilter(new FileNameExtensionFilter("Python script", "py"));
			jfc.setAcceptAllFileFilterUsed(true);
			jfc.setMultiSelectionEnabled(false);
			final String path = Reliquary.getPreferences().get(STARTUP, System.getProperty("user.home"));
			jfc.setCurrentDirectory(new File(path));
			final int res = jfc.showOpenDialog(JScriptingPreferencePanel.this);
			if (res == JFileChooser.APPROVE_OPTION) {
				final String path2 = jfc.getSelectedFile().getPath();
				jTextField_Startup.setText(path2);
				Reliquary.getPreferences().put("Startup Script", path2);
			}
		}
	}

	/**
	 * Constructs a JScriptingPreferencePanel
	 */
	public JScriptingPreferencePanel() {
		super("Scripting", "Configures the startup script and other scripting related preferences");
		init();
	}

	private void init() {
		final FormLayout fl = new FormLayout("5dlu, right:pref, 5dlu, fill:pref:grow, 5dlu, pref ",
				"5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu");
		final PanelBuilder pb = new PanelBuilder(fl, this);
		int yy = 2;
		pb.addSeparator("Startup Script", CC.xyw(2, yy, 5));
		yy += 2;
		pb.add(new JLabel("Python script"), CC.xy(2, yy), jTextField_Startup, CC.xy(4, yy));
		jTextField_Startup.setToolTipText(
				"A script that is run each time the application starts to initialize the scripting environment.");
		pb.add(jButton_Startup, CC.xy(6, yy));
		yy += 2;
		pb.add(new JLabel("Site variable"), CC.xy(2, yy), jTextField_Site, CC.xy(4, yy));
		jTextField_Site.setToolTipText(
				"The variable SITE within the Python scripting environment is set to this text value upon initialization.");
		setBorder(SwingUtils.createDefaultBorder());
	}

	static public final String SITE = "SITE";
	static public final String CLASSPATH = "CLASSPATH";
	static public final String STARTUP = "STARTUP";

	/**
	 * @see com.duckandcover.swing.JPreferencePanel#commit()
	 */
	@Override
	public void commit() {
		final Preferences prefs = Reliquary.getPreferences();
		prefs.put(SITE, jTextField_Site.getText());
		prefs.put(STARTUP, jTextField_Startup.getText());
	}

	public String[] getClassPathItems() throws FileNotFoundException {
		return parseClassPathItems(getClassPath());
	}

	private String[] parseClassPathItems(
			final String cp
	) throws FileNotFoundException {
		final String[] items = cp.split(";");
		final String cwd = System.getProperty("user.dir");
		final ArrayList<String> res = new ArrayList<>();
		final StringBuffer errs = new StringBuffer();
		for (final String item : items) {
			if (item.trim().length() > 0) {
				final String newItem = item.trim().replace("$CWD$", cwd);
				final File test = new File(newItem);
				if (!test.isFile()) {
					if (errs.length() > 0)
						errs.append("\n");
					errs.append("   JAR " + test + " not found.");
				} else
					res.add(newItem);
			}
		}
		if (errs.length() > 0) {
			final String msg = "One or more of the JARs specified in the CLASSPATH are missing.\n" + errs.toString();
			throw new FileNotFoundException(msg);
		}
		return res.toArray(new String[res.size()]);
	}

	/**
	 * @see com.duckandcover.swing.JPreferencePanel#initialize()
	 */
	@Override
	public void initialize() {
		jTextField_Site.setText(getSite());
		jTextField_Startup.setText(getStartup());
	}

	/**
	 * Gets the startup script path
	 *
	 * @return String
	 */
	public static String getStartup() {
		return Reliquary.getPreferences().get(STARTUP, "");
	}

	/**
	 * Gets the value to assign to the variable SITE
	 *
	 * @return String
	 */
	public static String getSite() {
		return Reliquary.getPreferences().get(SITE, "");
	}

	/**
	 * Returns the classpath used to initialize the scripting environment to permit
	 * custom JARs to be loaded.
	 *
	 * @return String
	 */
	public static String getClassPath() {
		return Reliquary.getPreferences().get(CLASSPATH, "");
	}
}
