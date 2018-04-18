package com.duckandcover.scripting;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;

import javax.swing.SwingWorker;

import org.python.util.InteractiveConsole;

import com.duckandcover.Reliquary;
import com.duckandcover.html.HTML;

/**
 * <p> A SwingWorker thread which runs a scripting engine and provides
 * mechanisms for queuing tasks on the scripting engine. </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version 1.0
 */
public class ScriptingWorker
   extends
   SwingWorker<Object, ScriptOutput> {


   private InteractiveConsole mScripter;
   private File mStartup;
   private int mCmdIndex = 0;

   private static final String QUIT = "!QUIT!";
   private static final String TERMINATED = "terminated";

   private final LinkedBlockingQueue<Object> mCommands = new LinkedBlockingQueue<>();

   private final JPythonPanel mPanel;

   private class Assignment {
      private final String mName;
      private final Object mValue;

      private Assignment(final String name, final Object val) {
         mName = name;
         mValue = val;
      }
   }

   abstract private class PyOutStream
      extends
      Writer {

      private boolean mClosed;
      private final StringBuffer mBuffer = new StringBuffer();

      private PyOutStream() {
         mClosed = false;
      }

      abstract void send(String msg);

      @Override
      public void close()
            throws IOException {
         flush();
         mClosed = true;
      }

      @Override
      public void flush()
            throws IOException {
         synchronized(mBuffer) {
            if(mBuffer.length() > 0) {
               send(mBuffer.toString());
               mBuffer.setLength(0);
            }
         }
      }

      @Override
      public void write(final char[] cbuf, final int off, final int len)
            throws IOException {
         if(!mClosed)
            synchronized(mBuffer) {
               for(int i = 0; i < len; ++i)
                  if(cbuf[i + off] != '\n')
                     mBuffer.append(cbuf[i + off]);
                  else {
                     send(mBuffer.toString());
                     mBuffer.setLength(0);
                  }
            }
      }
   }

   /**
    * Constructs a ScriptingWorker
    */
   public ScriptingWorker(final JPythonPanel panel) {
      mPanel = panel;
   }

   public void terminate() {
      if(mScripter != null)
         mScripter.set(TERMINATED, Boolean.TRUE);
   }

   @Override
   protected void process(final List<ScriptOutput> chunks) {
      for(final ScriptOutput tmp : chunks)
         process(tmp);
      mPanel.flush();
   }

   private void process(final ScriptOutput chunk) {
      if(chunk instanceof ScriptOutput.Flush)
         mPanel.flush();
      else
         mPanel.append(chunk.asHTML());
   }

   /**
    * Assign the @param name @param value @throws InterruptedException
    */
   public void put(final String name, final Object value) {
      mCommands.offer(new Assignment(name, value));
   }

   /**
    * Clears the assignment of the variable called 'name'. @param name @throws
    * InterruptedException
    */
   public void clear(final String name) {
      put(name, null);
   }

   /**
    * Excute a single command string using the scripting engine. @param
    * cmd @throws InterruptedException
    */
   public void execute(final String cmd) {
      mCommands.offer(cmd);
   }

   /**
    * Excute the commands in a file using the scripting engine. @param
    * file @throws InterruptedException
    */
   public void execute(final File file) {
      mCommands.offer(file);
   }

   private void execNow(final InputStream inputStream)
         throws IOException {
      final ByteArrayOutputStream result = new ByteArrayOutputStream();
      final byte[] buffer = new byte[1024];
      int length;
      while((length = inputStream.read(buffer)) != -1) {
         result.write(buffer, 0, length);
      }
      final String script = result.toString(StandardCharsets.US_ASCII.name()) + "\n";
      mScripter.exec(script);
   }

   public void setStartupScript(final File startup) {
      mStartup = startup;
   }

   private String elapseToString(final long ms) {
      final int hr = (int) (ms / 3600000), min = (int) ((ms / 60000) % 60), sec = (int) ((ms / 1000) % 60),
            milli = (int) (ms % 1000);
      final StringBuffer sb = new StringBuffer();
      if(hr > 0)
         sb.append(Integer.toString(hr) + " hr ");
      if((hr > 0) || (min > 0))
         sb.append(Integer.toString(min) + " min ");
      sb.append(Integer.toString(sec) + "." + Integer.toString(milli) + " s");
      return sb.toString();
   }

   private String removeComment(final String cmd) {
      final int p = cmd.indexOf('#');
      return p != -1 ? cmd.substring(0, p) : cmd;
   }

   /**
    * Output from the scripting thread happens two different ways. Most of the
    * output go through the stdout and stderr mechanism. Error output is
    * generated directly in the doInBackground method.
    *
    * @return null
    * @throws Exception
    * @see javax.swing.SwingWorker#doInBackground()
    */
   @Override
   protected Object doInBackground()
         throws Exception {
      final String oldName = Thread.currentThread().getName();
      try {
         Thread.currentThread().setName("Jython Interactive Console");
         Properties props = new Properties();
         final String path = System.getProperty("user.dir");
         props.put("python.path", path);
         props.put("python.cachedir.skip", "false");
         props.put("python.cachedir", path + "\\jython_cache");

         mScripter = new InteractiveConsole();
         mScripter.set("__worker__", this);
         // mScripter.exec("import sys");
         // mScripter.exec("sys.path.append(\"" + path + "\")");
         // mScripter.exec("sys.packageManager.addJarDir(\"" + path + "\",True)");
         mScripter.setOut(new PyOutStream() {
            @Override
            void send(final String msg) {
               publish(new ScriptOutput.Result(msg));
            }
         });
         mScripter.setErr(new PyOutStream() {
            @Override
            void send(final String msg) {
               publish(new ScriptOutput.Error(msg));
            }
         });
         try {
            ClassLoader classLoader = ScriptingWorker.class.getClassLoader();
            final URL url = classLoader.getResource("lib/startUp.py");
            if(url != null)
               try (InputStream inputStream = url.openStream()) {
                  execNow(inputStream);
               }
            else
               publish(new ScriptOutput.Error("Can't locate startUp.py"));
         }
         catch(final Throwable e1) {
            publish(new ScriptOutput.Error(e1.toString()));
            Reliquary.getLogger().log(Level.WARNING, "Error executing script: Internal startUp.py script from resource.", e1);
         }
         if(mStartup != null) {
            final String pyPath = mStartup.getCanonicalPath().replace("\\", "/");
            try {
               mScripter.execfile(pyPath);
            }
            catch(final Throwable e) {
               publish(new ScriptOutput.Error(e.toString()));
               Reliquary.getLogger().log(Level.WARNING, "Error executing script: " + pyPath, e);
            }
         }
         // publish(new ScriptOutput.HTML("<p>Ready...</p>"));
         while(true) {
            if(mCommands.contains(QUIT))
               break;
            final Object obj = mCommands.take();
            if(obj instanceof File) {
               final File file = (File) obj;
               String dupName = null;
               try {
                  // Create a replica of the script and then link to it.
                  final File dup = File.createTempFile("script", ".py", Reliquary.getReport().getParentFile());
                  Files.copy(file.toPath(), dup.toPath(), StandardCopyOption.REPLACE_EXISTING);
                  dup.setReadOnly();
                  dupName = dup.getName();
               }
               catch(final Throwable e) {
                  Reliquary.getLogger().log(Level.WARNING, "Can't execute script: " + file.getPath(), e);
               }
               final String pyPath = file.getCanonicalPath().replace("\\", "/");
               publish(new ScriptOutput.ExecFileCommand(++mCmdIndex, pyPath, dupName));
               final long ctm = System.currentTimeMillis();
               try {
                  mScripter.set(TERMINATED, Boolean.FALSE);
                  mScripter.push("execfile('" + pyPath + "')");
                  if(mScripter.get(TERMINATED).__nonzero__())
                     publish(new ScriptOutput.Error("Warning: Script execution may have been terminated prematurely."));
               }
               catch(final Throwable e) {
                  publish(new ScriptOutput.Error(e.toString()));
               }
               final long elapse = System.currentTimeMillis() - ctm;
               if(elapse > (10 * 1000)) // 10 s
                  publish(new ScriptOutput.HTML(HTML.br()
                        + HTML.p(HTML.i("Elapse:") + "&nbsp;" + elapseToString(elapse) + HTML.br())));
               publish(new ScriptOutput.EndOfScriptMarker());
            } else if(obj instanceof String) {
               final String cmds = ((String) obj);
               if(cmds.length() > 0) {
                  publish(new ScriptOutput.Command(++mCmdIndex, cmds));
                  try {
                     mScripter.set(TERMINATED, Boolean.FALSE);
                     final String[] cmdLines = cmds.split(System.lineSeparator());
                     boolean inMultiLine = false;
                     for(String cmd : cmdLines) {
                        cmd = removeComment(cmd);
                        if(inMultiLine && (cmd.length() > 0) && (!Character.isWhitespace(cmd.charAt(0))))
                           mScripter.push("");
                        if(cmd.trim().length() > 0)
                           inMultiLine = mScripter.push(cmd);
                     }
                     if(inMultiLine)
                        mScripter.push("");
                     if(mScripter.get(TERMINATED).__nonzero__())
                        publish(new ScriptOutput.Error("Warning: Command execution may have been terminated prematurely."));
                  }
                  catch(final Throwable e) {
                     publish(new ScriptOutput.Error(e.toString()));
                  }
               }
            } else if(obj instanceof Assignment) {
               final Assignment ass = (Assignment) obj;
               try {
                  mScripter.set(ass.mName, ass.mValue);
               }
               catch(final Throwable e) {
                  publish(new ScriptOutput.Error(e.getMessage()));
               }
            }
            publish(new ScriptOutput.Flush());
         }
      }
      catch(final Throwable t) {
         publish(new ScriptOutput.Error(t.getMessage()));
      }
      Thread.currentThread().setName(oldName);
      return null;
   }

   public void publishHTML(final String html) {
      publish(new ScriptOutput.HTML(html));
   }
}
