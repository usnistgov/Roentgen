package com.duckandcover.html;

import java.awt.Desktop;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * <p>
 * Provides a mechanism to combine objects implementing IToHTML and IToHTMLExt
 * into a single enlongated report.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2015
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: $
 */
public class Report implements IToHTMLExt {

	private final String mName;
	private final List<IToHTML> mItems;

	private class ForceMode implements IToHTML, IToHTMLExt {
		
		private final IToHTML mBase;
		private final Mode mMode;
		
		ForceMode(IToHTML base, Mode mode){
			mBase = base;
			mMode = mode;
		}
		

		@Override
		public String toHTML(Mode mode) {
			return mBase.toHTML(mMode);
		}

		@Override
		public String toHTML(Mode mode, File base, String dir) throws IOException {
			if(mBase instanceof IToHTMLExt)
				return ((IToHTMLExt)mBase).toHTML(mMode,base,dir);
			else
				return mBase.toHTML(mMode);
		}
	}
	
	/**
	 * Constructs a HTMLReport
	 */
	public Report(final String name) {
		mName = name;
		mItems = new ArrayList<>();
	}

	public void addHeader(final IToHTML item) {
		mItems.add(Transforms.h3(item));
	}

	public void addHeader(final String html) {
		addHeader(Transforms.createHTML(html));
	}

	public void addSubHeader(final IToHTML item) {
		mItems.add(Transforms.h4(item));
	}

	public void addSubHeader(final String html) {
		addSubHeader(Transforms.createHTML(html));
	}

	public void add(final List<? extends Object> items, final IToHTML.Mode mode) {
		add(Transforms.createList(items), mode);
	}

	public  void add(final List<? extends Object> items) {
		mItems.add(Transforms.createList(items));
	}

	public void addHTML(final String html) {
		mItems.add(Transforms.createHTML(html));
	}
	
	public void addVerbatim(final String text) {
		mItems.add(Transforms.verbatim(text));
	}
	
	public void addThrowable(final Throwable e) {
		final StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		e.printStackTrace(pw);
		mItems.add(Transforms.verbatim(sw.toString(),"Tomato"));
	}

	
	public void addTable(final List<String> colHeaders, final List<String> rowHeaders, final List<List<String>> items) {
		addHTML(HTML.asTable(colHeaders, rowHeaders, items));
	}

	public void addImage(final BufferedImage img, final String caption) {
		add(Transforms.scrollPane(new Image(img, caption)));
	}

	public void add(final Map<? extends Object, ? extends Object> map, final IToHTML.Mode keyMode,
			final IToHTML.Mode valueMode) {
		Table t = new Table();
		for (Map.Entry<? extends Object, ? extends Object> me : map.entrySet())
			t.addRow(Table.td(HTML.toHTML(me.getKey(), keyMode)), Table.td(HTML.toHTML(me.getValue(), valueMode)));
		add(t);
	}

	/**
	 * Adds a block of text. Special characters are escaped to make the text HTML
	 * friendly so no HTML code should be present.
	 *
	 * @param note
	 */
	public void addNote(final String note) {
		addHTML(HTML.p(HTML.escape(note)));
	}

	public void addParagraph(final IToHTML item) {
		add(Transforms.createParagraph(item));
	}

	public void add(final IToHTML item, final IToHTML.Mode mode) {
		mItems.add(new ForceMode(item,mode));
	}

	public void add(final IToHTML item) {
		mItems.add(item);
	}

	public void add(final Table table) {
		mItems.add(Transforms.scrollPane(table));
	}

	@Override
	public String toHTML(final Mode mode) {
		final StringBuffer sb = new StringBuffer();
		for (final IToHTML item : mItems)
			sb.append(HTML.toHTML(item, mode));
		return sb.toString();
	}

	public String highest(final IToHTML item, final Mode mode, final File base, final String dir) throws IOException {
		if (item instanceof IToHTMLExt)
			return ((IToHTMLExt) item).toHTML(mode, base, dir);
		else
			return item.toHTML(mode);
	}

	@Override
	public String toHTML(final Mode mode, final File base, final String dir) throws IOException {
		final File fd = new File(base, dir);
		final StringBuffer sb = new StringBuffer();
		if (fd.isDirectory() || fd.mkdirs()) {
			for (final IToHTML item : mItems)
				sb.append(highest(item, mode, base, dir));
		} else {
			sb.append(HTML.error("Unable to create the report directory. Defaulting to the base report."));
			for (final IToHTML item : mItems)
				sb.append(item.toHTML(mode));

		}
		return sb.toString();
	}

	/**
	 * Writes the report to the specified file.
	 *
	 * @param f    {@link File}
	 * @param mode {@link Mode}
	 * @throws FileNotFoundException
	 */
	public void toFile(final File f, final Mode mode) throws IOException {
		HTML.toFile(this, f, mode, mName);
	}

	/**
	 * Outputs the specified object to an HTML file and opens it in a browser.
	 * 
	 * @param obj
	 * @param mode
	 * @throws IOException
	 */
	public static void dump(Object obj, Mode mode) throws IOException {
		Report r = new Report(obj.toString());
		r.addHTML(HTML.toHTML(obj, mode));
		r.inBrowser(Mode.VERBOSE);
	}
	
	/**
	 * Outputs the specified object to an HTML file and opens it in a browser.
	 * 
	 * @param html
	 * @throws IOException
	 */
	public static void dump(String html) throws IOException {
		Report r = new Report("Dump");
		r.addHTML(html);
		r.inBrowser(Mode.VERBOSE);
	}


	/**
	 * Writes the report to a temporary file and then opens the file in the browser.
	 *
	 * @param mode
	 * @throws IOException
	 */
	public void inBrowser(final Mode mode) throws IOException {
		final File f = File.createTempFile("temp", ".html");
		HTML.toFile(this, f, mode, mName);
		Desktop.getDesktop().browse(f.toURI());
	}

	@Override
	public String toString() {
		return mName;
	}
}
