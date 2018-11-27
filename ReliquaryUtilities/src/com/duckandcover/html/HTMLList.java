package com.duckandcover.html;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

public class HTMLList implements IToHTMLExt {
	
	private final java.util.List<IToHTML> mItems = new ArrayList<>();

	public void add(IToHTML item) {
		mItems.add(item);
	}
	
	public void add(String html) {
		mItems.add(Transforms.createHTML(html));
	}
	
	public void add(Object obj) {
		if(obj instanceof IToHTML)
			mItems.add((IToHTML)obj);
		else 
			mItems.add(Transforms.createHTML(HTML.escape(obj.toString())));
	}
	
	public void addAll(Collection<? extends Object> objs) {
		for(Object obj : objs)
			add(obj);
	}

	@Override
	public String toHTML(Mode mode) {
		StringBuffer sb = new StringBuffer();
		sb.append("<ul>");
		for(IToHTML item : mItems) {
			sb.append("<li>");
			sb.append(item.toHTML(mode));
			sb.append("</li>");
		}
		sb.append("</ul>");
		return sb.toString();
	}

	@Override
	public String toHTML(Mode mode, File base, String dir) throws IOException {
		StringBuffer sb = new StringBuffer();
		sb.append("<ul>");
		for(IToHTML item : mItems) {
			sb.append("<li>");
			sb.append(Transforms.promote(item).toHTML(mode, base, dir));
			sb.append("</li>");
		}
		sb.append("</ul>");
		return sb.toString();
	}
}
