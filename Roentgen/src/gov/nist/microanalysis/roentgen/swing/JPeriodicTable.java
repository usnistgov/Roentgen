package gov.nist.microanalysis.roentgen.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.SystemColor;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import com.duckandcover.swing.SwingUtils;
import com.google.common.base.Preconditions;
import com.jgoodies.forms.builder.ButtonBarBuilder;

import gov.nist.microanalysis.roentgen.physics.Element;

/**
 * <p>
 * A periodic table control for selecting one or more elements.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2016
 * </p>
 *
 * @author nritchie
 * @version $Rev: 309 $
 */

public class JPeriodicTable extends JComponent {
	private static final long serialVersionUID = 0x1;
	private Dimension mButtonDim;

	private final Set<Element> mSelected = new HashSet<>();
	private final Set<Element> mDisabled = new HashSet<>();
	private Element mDepressedElement;
	private Element mDisplayedElement;
	private Element mSelectedElement; // The
	// last
	// element
	// selected
	private List<Element> mAvailable = Element.range(Element.Hydrogen, Element.Lawrencium);

	public void setLastElement(final Element last) {
		mAvailable = Element.range(Element.Hydrogen, last);
		final Set<Element> sel = trim(mSelected);
		mSelected.clear();
		mSelected.addAll(sel);
	}

	private static final int[] mRow = { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 7, 7, 7, 7, 7, 7, 7,
			7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
			8, 8, 8, 8 };

	private static final int[] mColumn = { 0, 17, 0, 1, 12, 13, 14, 15, 16, 17, 0, 1, 12, 13, 14, 15, 16, 17, 0, 1, 2,
			3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
			16, 17, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
			15, 16, 17, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 };

	private static final int COL_COUNT = 18;
	private static final int ROW_COUNT = 9;

	private transient Vector<ActionListener> actionListeners = new Vector<>();

	protected Rectangle getButtonRect(final Element elm) {
		final int h = mButtonDim.height;
		final int w = mButtonDim.width;
		final int r = (h * mRow[elm.ordinal()]) + ((getHeight() - (ROW_COUNT * h)) / 2);
		final int c = (w * mColumn[elm.ordinal()]) + ((getWidth() - (COL_COUNT * w)) / 2);
		return new Rectangle(c, r, w - 2, h - 2);
	}

	private Element whichElement(final int x, final int y) {
		for (final Element e : mAvailable)
			if (getButtonRect(e).contains(x, y))
				return e;
		return null;
	}

	private Rectangle getElementDisplayRect() {
		return new Rectangle(4 * mButtonDim.width, (mButtonDim.height / 2), 7 * mButtonDim.width,
				(5 * mButtonDim.height) / 2);
	}

	protected void displayPopup(final int x, final int y) {
		final JPopupMenu pm = new JPopupMenu("Periodic table menu");
		{
			final JMenuItem mi = new JMenuItem("Clear all");
			mi.addActionListener(new AbstractAction() {
				static final long serialVersionUID = 0x1;

				@Override
				public void actionPerformed(final ActionEvent e) {
					setAll(false);
				}
			});
			pm.add(mi);
		}
		{
			final JMenuItem mi = new JMenuItem("Select all");
			mi.addActionListener(new AbstractAction() {
				static final long serialVersionUID = 0x1;

				@Override
				public void actionPerformed(final ActionEvent e) {
					setAll(true);
				}
			});
			pm.add(mi);
		}
		this.add(pm);
		{
			final JMenuItem mi = new JMenuItem("Invert selection");
			mi.addActionListener(new AbstractAction() {
				static final long serialVersionUID = 0x1;

				@Override
				public void actionPerformed(final ActionEvent e) {
					final Set<Element> end = new HashSet<>(mAvailable);
					end.removeAll(getSelected());
					setSelected(end, true);
				}
			});
			pm.add(mi);
		}
		this.add(pm);
		pm.show(this, x, y);
	}

	protected void fireSelectionEvent(final Element element) {
		mSelectedElement = element;
		fireActionPerformed(
				new ActionEvent(this, element.getAtomicNumber(), getSelection(element) ? "Selected" : "Deselected"));
	}

	/**
	 * fireActionPerformed - This control fires an event each time an element is
	 * selected or deselected. @param e ActionEvent
	 */
	protected void fireActionPerformed(final ActionEvent e) {
		if (actionListeners.size() > 0)
			for (final Object element : actionListeners)
				((ActionListener) element).actionPerformed(e);
	}

	/**
	 * getLastSelected - Get the last element that was selected or
	 * deselected. @return int
	 */
	public Element getLastSelected() {
		return mSelectedElement;
	}

	/**
	 * JPeriodicTable - Construct a JPeriodicTable object.
	 */
	public JPeriodicTable() {
		super();
		try {
			init();
		} catch (final Exception ex) {
			ex.printStackTrace();
		}
	}

	/**
	 * getSelection - Is the specified element selected in the table? @param element
	 * int - The atomic number of the element. @return boolean - True if selected,
	 * false otherwise.
	 */
	public boolean getSelection(final Element element) {
		return mSelected.contains(element) && mAvailable.contains(element);
	}

	/**
	 * setSelection - Set the selection state of the specified element. @param elm
	 * Element - The element to select @param b boolean - True to select, false to
	 * unselect.
	 */
	public void setSelection(final Element elm, final boolean b) {
		Preconditions.checkArgument(elm != null);
		if (mAvailable.contains(elm))
			if (mSelected.contains(elm) != b) {
				if (b)
					mSelected.add(elm);
				else
					mSelected.remove(elm);
				repaint(mButtonDim.width * mColumn[elm.ordinal()], mButtonDim.height * mRow[elm.ordinal()],
						mButtonDim.width, mButtonDim.height);
			}
	}

	private void init() throws Exception {
		this.setOpaque(true);
		this.setSize(25 * COL_COUNT, 20 * ROW_COUNT);
		this.setMaximumSize(new Dimension(50 * COL_COUNT, 40 * ROW_COUNT));
		this.setMinimumSize(new Dimension(20 * COL_COUNT, 15 * ROW_COUNT));
		this.setPreferredSize(new Dimension(25 * COL_COUNT, 20 * ROW_COUNT));
		this.setToolTipText("Select one or more elements from the periodic table.");
		this.setDoubleBuffered(true);
		this.setOpaque(false);
		this.setLayout(null);
		// this.setBackground(SystemColor.text);
		mButtonDim = new Dimension(getWidth() / COL_COUNT, getHeight() / ROW_COUNT);
		mDepressedElement = null;
		mDisplayedElement = null;

		addMouseListener(new MouseAdapter() {
			@Override
			public void mousePressed(final MouseEvent me) {
				if (me.isPopupTrigger())
					displayPopup(me.getX(), me.getY());
				else if (me.getButton() == MouseEvent.BUTTON1) {
					Element el = whichElement(me.getX(), me.getY());
					if (mDisabled.contains(el))
						el = null;
					mDepressedElement = el;
					if (mDepressedElement != null)
						repaint(getButtonRect(mDepressedElement));
				}
			}

			@Override
			public void mouseReleased(final MouseEvent me) {
				if (me.isPopupTrigger())
					displayPopup(me.getX(), me.getY());
				else if (me.getButton() == MouseEvent.BUTTON1)
					if (mDepressedElement != null) {
						if (whichElement(me.getX(), me.getY()) == mDepressedElement) {
							setSelection(mDepressedElement, !mSelected.contains(mDepressedElement));
							fireSelectionEvent(mDepressedElement);
						}
						mDepressedElement = null;
					}
			}
		});
		addMouseMotionListener(new MouseMotionAdapter() {
			@Override
			public void mouseMoved(final MouseEvent me) {
				final Element el = whichElement(me.getX(), me.getY());
				if (el != mDisplayedElement) {
					mDisplayedElement = el;
					repaint(getElementDisplayRect());
				}
			}
		});
	}

	@Override
	public void paintComponent(final Graphics gr) {
		super.paintComponent(gr);
		final Insets in = getInsets();
		mButtonDim.setSize((getWidth() - (2 * (in.left + in.right))) / COL_COUNT,
				(getHeight() - (2 * (in.top + in.bottom))) / ROW_COUNT);
		FontMetrics fm = gr.getFontMetrics();
		final int h = mButtonDim.height;
		final int w = mButtonDim.width;
		final int dh = (getHeight() - (ROW_COUNT * h)) / 2;
		final int dw = (getWidth() - (COL_COUNT * w)) / 2;
		for (final Element e : mAvailable) {
			int r = (h * mRow[e.ordinal()]) + dh;
			int c = (w * mColumn[e.ordinal()]) + dw;
			if (e == mDepressedElement) {
				++r;
				++c;
			}
			if (mDisabled.contains(e)) {
				gr.setColor(SystemColor.control);
				gr.fillRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
				gr.setColor(SystemColor.controlShadow);
				gr.drawRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
			} else if (mSelected.contains(e)) {
				gr.setColor(SystemColor.controlShadow);
				gr.fillRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
				gr.setColor(SystemColor.controlText);
				gr.drawRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
			} else {
				gr.setColor(SystemColor.control);
				gr.fillRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
				gr.setColor(SystemColor.controlText);
				gr.drawRoundRect(c, r, w - 2, h - 2, w / 10, w / 10);
			}
			final String str = e.getAbbrev();
			if (mDisabled.contains(e))
				gr.setColor(Color.lightGray);
			gr.drawString(str, c + ((w - fm.stringWidth(str)) / 2), r + ((2 * h) / 3));
			if (mDisabled.contains(e))
				gr.setColor(this.getForeground());
		}

		final Rectangle r = getElementDisplayRect();
		gr.clearRect(r.x, r.y, r.width, r.height);
		if (mDisplayedElement != null) {
			if (mDisabled.contains(mDisplayedElement))
				gr.drawString("disabled", (4 * mButtonDim.width) + 2, ((5 * mButtonDim.height) / 2) - 2);
			final Font oldFont = gr.getFont();
			gr.setFont(new Font(oldFont.getName(), oldFont.getStyle(), (3 * oldFont.getSize()) / 2));
			fm = gr.getFontMetrics();
			final int xw = fm.stringWidth("X");
			final String str = mDisplayedElement.toString();
			gr.setColor(SystemColor.textInactiveText);
			gr.drawString(str, r.x + (2 * xw),
					Math.min(r.y + fm.getHeight() + fm.getAscent(), (r.y + r.height) - fm.getDescent()));
			gr.drawString(mDisplayedElement.getAbbrev(), r.x, r.y + fm.getAscent());
			gr.drawString(Integer.toString(mDisplayedElement.getAtomicNumber()), r.x + (4 * xw), r.y + fm.getAscent());
			gr.drawString(Double.toString(mDisplayedElement.getAtomicWeight()), r.x + (8 * xw), r.y + fm.getAscent());
		}
	}

	/**
	 * getSelectedElements - Get a Set object containing a list of selected
	 * element. @return Set&gt;Element&lt; - a set of Element objects
	 */
	public Set<Element> getSelected() {
		return trim(mSelected);
	}

	/**
	 * setSelectedElements - Sets the elements that are selected on the periodic
	 * table control from the list of Element objects in lst. @param lst A
	 * collection of Element objects. @param clear boolean - Determines whether to
	 * clear the table before setting the specified elements.
	 */
	public void setSelected(final Collection<Element> lst, final boolean clear) {
		if (clear)
			mSelected.clear();
		mSelected.addAll(trim(lst));
		repaint();
	}

	private Set<Element> trim(final Collection<Element> elms) {
		final TreeSet<Element> res = new TreeSet<>();
		for (final Element elm : elms)
			if (mAvailable.contains(elm))
				res.add(elm);
		return Collections.unmodifiableSet(res);
	}

	/**
	 * setAll - Select or deselect all elements. @param set boolean - True to set
	 * and false to deselect.
	 */
	public void setAll(final boolean set) {
		if (set)
			mSelected.addAll(mAvailable);
		else
			mSelected.clear();
		repaint();
	}

	/**
	 * setAllExcept - Select or deselect all elements. @param set boolean - True to
	 * set and false to deselect. @param element int - The atomic number of the
	 * element whose state to leave unchanged
	 */
	public void setAllExcept(final boolean set, final Set<Element> element) {
		Preconditions.checkArgument(element != null);
		for (final Element el : mAvailable)
			if (!element.contains(el)) {
				if (set)
					mSelected.add(el);
				else
					mSelected.remove(el);
				repaint(getButtonRect(el));
			}
	}

	/**
	 * removeActionListener - This control fires an event each time an element is
	 * selected or deselected. @param l ActionListener
	 */
	public synchronized void removeActionListener(final ActionListener l) {
		if (actionListeners.contains(l)) {
			final Vector<ActionListener> v = new Vector<>(actionListeners);
			v.removeElement(l);
			actionListeners = v;
		}
	}

	/**
	 * addActionListener - This control fires an event each time an element is
	 * selected or deselected. @param l ActionListener
	 */
	public synchronized void addActionListener(final ActionListener l) {
		if (!actionListeners.contains(l)) {
			final Vector<ActionListener> v = new Vector<>(actionListeners);
			v.addElement(l);
			actionListeners = v;
		}
	}

	/**
	 * setEnabled - Set the enabled state for the button assocaited with the
	 * specified element. @param element int @param enabled boolean
	 */
	public void setEnabled(final Element z, final boolean enabled) {
		Preconditions.checkArgument(z != null);
		if (mDisabled.contains(z) == enabled) {
			if (enabled)
				mDisabled.remove(z);
			else
				mDisabled.add(z);
			repaint(getButtonRect(z));
			if (!enabled)
				mSelected.remove(z);
		}
	}

	public void setEnabled(final Collection<Element> elms, final boolean enabled) {
		for (final Element elm : elms)
			setEnabled(elm, enabled);
	}

	/**
	 * enableAll - Enables (or disables) all buttons. @param enabled boolean
	 */
	public void enableAll(final boolean enabled) {
		for (final Element z : mAvailable)
			setEnabled(z, enabled);
	}

	static public Set<Element> selectElements(final Window window, final String title, final Collection<Element> avail,
			final Collection<Element> selected) {
		final JDialog dia = new JDialog(window);
		dia.setTitle(title);
		dia.setLocationRelativeTo(window);
		final JPeriodicTable jpt = new JPeriodicTable();
		jpt.setBorder(SwingUtils.createEmptyBorder());
		dia.add(jpt, BorderLayout.CENTER);
		final ButtonBarBuilder bbb = new ButtonBarBuilder();
		bbb.addGlue();
		bbb.addButton(new AbstractAction("Ok") {

			private static final long serialVersionUID = -6866060088167807470L;

			@Override
			public void actionPerformed(final ActionEvent arg0) {
				// TODO Auto-generated method stub

			}

		});
		bbb.addRelatedGap();
		bbb.addButton(new AbstractAction("Cancel") {

			private static final long serialVersionUID = -6588697111214385081L;

			@Override
			public void actionPerformed(final ActionEvent arg0) {
				// TODO Auto-generated method stub

			}

		});
		return null;

	}

}
