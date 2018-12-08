package com.duckandcover.scripting;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;

/**
 * 
 * 
 * 
 * @author Nicholas
 *
 */
public class ScriptClassLoader extends URLClassLoader {

	private static URL[] extURL() throws MalformedURLException {
		final String jcp = "D:\\Users\\Nicholas\\git\\Roentgen\\Reliquary";// System.getProperty("user.dir");
		if (jcp != null) {
			final File fjcp = new File(jcp,"Lib");
			File[] files = fjcp.listFiles();
			ArrayList<URL> urls = new ArrayList<URL>();
			for(File f : files) 
				if(f.isFile() && f.getName().toLowerCase().endsWith(".jar")) {
					URL url = f.toURI().toURL();
					urls.add(url);
					System.out.println(url);
				}
			if(urls.size()>0)
				return urls.toArray(new URL[urls.size()]);
		}
		return new URL[0];
	}

	public ScriptClassLoader(ClassLoader parent) throws MalformedURLException {
		super(extURL(), parent);
	}
	
	@Override
	public Class<?> findClass(String cls) throws ClassNotFoundException {
		System.out.println("Finding "+cls);
		Class<?> res= super.findClass(cls);
		System.out.println("Result is "+res);
		return res;
	}

	@Override
	public String toString() {
		return "ScriptClassLoader[...]";
	}
}
