package com.duckandcover.swing;

import java.io.Serializable;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

import com.jgoodies.forms.builder.PanelBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

/**
 * <p>
 * A progress panel with a progress bar that scrolls while a computation
 * returning an object of type T is performed in a background thread. When the
 * computation is done, the result of type T can be accessed using
 * </p>
 *
 * @author Nicholas W. M. Ritchie
 * @version 1.0
 */
public class JWizardProgressPanel<W extends Serializable, T>
   extends
   JWizardPanel<W> {

   private static final long serialVersionUID = 4976247586981195458L;

   public abstract class ProgressWorker<WW, TT>
      extends
      SwingWorker<TT, Integer> {

      abstract public boolean initialize(WW config);

      @Override
      protected void process(final List<Integer> chunks) {
         final int val = Integer.MIN_VALUE;
         for(int i : chunks)
            if(i > val)
               i = val;
         jProgressBar_Progress.setValue(val);
         jLabel_Progress.setText(format(val) + mOutOf);
         if(mProgress >= jProgressBar_Progress.getMaximum())
            getWizardDialog().enableFinish(true);
      }

   }

   private final JProgressBar jProgressBar_Progress = new JProgressBar();
   private final JLabel jLabel_Progress = new JLabel("");
   private final ProgressWorker<W, T> mThread;
   private T mResults;
   private String mMessage;
   private String mOutOf;
   private final int mProgress = 0;

   private String format(final int val) {
      final StringBuffer res = new StringBuffer();
      final String tmp = Integer.toString(val);
      final int len = tmp.length();
      for(int i = 0; i < len; ++i) {
         if((len - i) % 3 == 0)
            res.append('\u2009');
         res.append(tmp.charAt(i));
      }
      return res.toString();
   }

   /**
    * <p>
    * Constructs a JWizardProgressPanel to perform and track the progress of a
    * task with a long duration.
    * </p>
    * .
    * <p>
    * The thread should implement {@link SwingWorker} methods
    * 'process(java.util.List&lt;Integer&gt; chunks)' to update this panel using
    * the setProgress method and 'done()' to set 'enableFinish(true)'.
    * </p>
    *
    * @param wiz
    * @param message
    * @param thread
    */
   public JWizardProgressPanel(final JWizardDialog<W> wiz, final String message, final ProgressWorker<W, T> sw) {
      super(wiz);
      mMessage = message;
      mThread = sw;
      jLabel_Progress.setHorizontalAlignment(SwingConstants.CENTER);
      initialize();
   }

   private void initialize() {
      final FormLayout layout = new FormLayout("200dlu", "pref, 2dlu, 20dlu, 5dlu, pref");
      setLayout(layout);
      final PanelBuilder pb = new PanelBuilder(layout, this);
      final CellConstraints cc = new CellConstraints();
      pb.addSeparator("Progress", cc.xy(1, 1));
      setRange(0, 100);
      pb.add(jProgressBar_Progress, cc.xy(1, 3));
      pb.add(jLabel_Progress, cc.xy(1, 5));
   }

   @Override
   public void onShow() {
      getWizardDialog().setMessageText(mMessage);
      getWizardDialog().enableFinish(false);
      if(mThread != null) {
         mThread.initialize(getWizardDialog().getResult());
         mThread.execute();
      } else {
         getWizardDialog().setMessageText("Nothing to do!");
      }
   }

   public void setRange(final int min, final int max) {
      jProgressBar_Progress.setMinimum(min);
      jProgressBar_Progress.setMaximum(max);
      final int val = Math.max(min, Math.min(max, jProgressBar_Progress.getValue()));
      if(jProgressBar_Progress.getValue() != val)
         jProgressBar_Progress.setValue(val);
      mOutOf = " out of " + format(jProgressBar_Progress.getMaximum() - jProgressBar_Progress.getMinimum());
   }

   /**
    * A thread safe mechanism to update the message displayed in the
    * {@link JWizardDialog}.
    *
    * @param msg
    */
   public void setMessage(final String msg) {
      mMessage = msg;
      SwingUtilities.invokeLater(new Runnable() {
         @Override
         public void run() {
            getWizardDialog().setMessageText(mMessage);
         }
      });
   }

   @Override
   public void onHide() {
      mResults = null;
      if(mThread != null)
         if(mThread.isDone())
            try {
               mResults = mThread.get();
            }
            catch(final Exception e) {
               mResults = null;
            }
         else
            mThread.cancel(true);
   }

   /**
    * Returns the calculation results if the calculation ended without error or
    * null otherwise.
    *
    * @return T
    */
   public T getResults() {
      return mResults;
   }
}
