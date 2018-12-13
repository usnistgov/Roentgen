package com.duckandcover.scripting;

/**
 * <p>
 * Objects representing different types of text output from the ScriptingWorker.
 * </p>
 * <p>
 * Copyright Nicholas W. M. Ritchie 2014-2019
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version $Rev: 293 $
 */
abstract public class ScriptOutput {

    protected final String mValue;

    /**
     * Constructs a ScriptOutput
     */
    protected ScriptOutput(final String text) {
	mValue = text;
    }

    abstract public String asHTML();

    /**
     * <p>
     * A command line
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */
    static public class Command extends ScriptOutput {

	private final int mIndex;

	public Command(final int i, final String val) {
	    super(val);
	    mIndex = i;
	}

	@Override
	public String asHTML() {
	    final StringBuffer res = new StringBuffer();
	    res.append("<p>");
	    final String lineNo = "<span class=\"cmdindex\">" + mIndex + "&gt;&nbsp;</span>";
	    final String[] lines = mValue.split("\n");
	    for (final String line : lines) {
		final String tmp = com.duckandcover.html.HTML.escape(line).replaceAll("\t", "&nbsp;&nbsp;&nbsp;")
			.replaceAll(" ", "&nbsp;");
		res.append(lineNo + "<span class=\"command\">" + tmp + "</span><br/>\n");
	    }
	    res.append("</p>");
	    return res.toString();
	}
    }

    static public class ExecFileCommand extends ScriptOutput {

	private final int mIndex;
	private final String mLink;

	public ExecFileCommand(final int i, final String script, final String link) {
	    super(script);
	    mIndex = i;
	    mLink = link;
	}

	@Override
	public String asHTML() {
	    final StringBuffer res = new StringBuffer();
	    res.append("<p>");
	    final String lineNo = "<span class=\"cmdindex\">" + mIndex + "&gt;&nbsp;</span>";
	    res.append(lineNo);
	    res.append("<span class=\"command\">");
	    res.append("execfile(\"");
	    if (mLink != null)
		res.append("<a href=\"" + mLink + "\" target=\"_blank\">");
	    res.append(mValue);
	    if (mLink != null)
		res.append("</a>");
	    res.append("\")</span></p>");
	    return res.toString();
	}
    }

    /**
     * <p>
     * Output formated as HTML
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */

    static public class HTML extends ScriptOutput {

	public HTML(final String text) {
	    super(text);
	}

	@Override
	public String asHTML() {
	    return mValue;
	}
    }

    /**
     * <p>
     * An error message
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */
    static public class Error extends ScriptOutput {

	public Error(final String text) {
	    super(text);
	}

	@Override
	public String asHTML() {
	    return com.duckandcover.html.HTML.error(com.duckandcover.html.HTML.escape(mValue));
	}
    }

    /**
     * <p>
     * An warning message
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */
    static public class Warning extends ScriptOutput {

	public Warning(final String text) {
	    super(text);
	}

	@Override
	public String asHTML() {
	    return com.duckandcover.html.HTML.warning(com.duckandcover.html.HTML.escape(mValue));
	}
    }

    /**
     * <p>
     * Result of a single command line.
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */
    static public class Result extends ScriptOutput {

	public Result(final String text) {
	    super(text);
	}

	@Override
	public String asHTML() {
	    if (mValue.startsWith("<html>"))
		return mValue.substring(6);
	    else
		return "<p class=\"result\">" + com.duckandcover.html.HTML.escape(mValue) + "</p>";
	}
    }

    static public class EndOfScriptMarker extends ScriptOutput {

	public EndOfScriptMarker() {
	    super("<hr class=\"eos\"/>");
	}

	@Override
	public String asHTML() {
	    return mValue;
	}
    }

    /**
     * <p>
     * Doesn't actually output anything to the report but instead writes the
     * report to disk.
     * </p>
     * <p>
     * Copyright Nicholas W. M. Ritchie 2014-2019
     * </p>
     *
     * @author Nicholas W. M. Ritchie
     * @version $Rev: 293 $
     */
    static public class Flush extends ScriptOutput {

	public Flush() {
	    super("");
	}

	@Override
	public String asHTML() {
	    return mValue;
	}
    }

}
