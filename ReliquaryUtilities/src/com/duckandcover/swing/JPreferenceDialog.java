package com.duckandcover.swing;

import java.awt.Font;
import java.awt.Frame;
import java.awt.HeadlessException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.border.Border;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.DefaultTreeSelectionModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import com.jgoodies.forms.builder.ButtonBarBuilder;
import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

/**
 * <p>
 * A mechanism for displaying user preferences in a unified fashion.
 * </p>
 * <p>
 * Libraries can add preference panels using the {@link JPreferencePanel} to
 * allow the user to specify library specific options.
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class JPreferenceDialog
   extends
   JDialog {

   private static final long serialVersionUID = 0x95962961118L;

   private final JTree mPreferenceTree = new JTree();
   private DefaultMutableTreeNode mRoot;
   private final JLabel mMessageLabel = new JLabel();
   private final JLabel mBannerLabel = new JLabel();
   private JScrollPane mContentPane = null;

   /**
    * Constructs a PreferenceDialog
    *
    * @throws HeadlessException
    */
   public JPreferenceDialog()
         throws HeadlessException {
      super();
      try {
         setDefaultCloseOperation(DISPOSE_ON_CLOSE);
         initialize();
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   private void initialize() {
      setLayout(new FormLayout("5dlu, 150dlu, 3dlu, 300dlu, 5dlu", "5dlu, pref, 5dlu, 200dlu, 5dlu, pref, 5dlu, pref, 5dlu"));
      final CellConstraints cc = new CellConstraints();
      final Border bdr = BorderFactory.createCompoundBorder(SwingUtils.createDefaultBorder(), SwingUtils.createEmptyBorder());

      mRoot = new DefaultMutableTreeNode("Root");
      mPreferenceTree.setModel(new DefaultTreeModel(mRoot));
      mPreferenceTree.setRootVisible(false);
      mPreferenceTree.setExpandsSelectedPaths(true);
      final TreeSelectionModel tsm = new DefaultTreeSelectionModel();
      tsm.setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
      mPreferenceTree.setSelectionModel(tsm);
      tsm.addTreeSelectionListener(new TreeSelectionListener() {
         @Override
         public void valueChanged(final TreeSelectionEvent e) {
            final DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) mPreferenceTree.getLastSelectedPathComponent();
            if(dmtn != null) {
               final Object obj = dmtn.getUserObject();
               mBannerLabel.setText(obj.toString());
               if(obj instanceof JPreferencePanel) {
                  final JPreferencePanel pp = (JPreferencePanel) obj;
                  pp.initialize();
                  if(mContentPane != null) {
                     remove(mContentPane);
                     mContentPane = null;
                  }
                  mMessageLabel.setText(pp.getDescription());
                  mContentPane = new JScrollPane(pp);
                  final CellConstraints cc = new CellConstraints();
                  pp.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
                  add(mContentPane, cc.xy(4, 4, CellConstraints.FILL, CellConstraints.FILL));
               } else
                  mMessageLabel.setText("Unexpected panel");
            }
         }
      });

      final JScrollPane sp = new JScrollPane(mPreferenceTree);
      add(sp, cc.xywh(2, 2, 1, 5));
      {
         mBannerLabel.setText("User preferences");
         final Font def = mBannerLabel.getFont();
         mBannerLabel.setFont(new Font(def.getName(), def.getStyle(), 2 * def.getSize()));
         {
            URL url = null;
            if(System.currentTimeMillis() % 2 == 0)
               url = JPreferenceDialog.class.getResource("femaleScientist64.png");
            else
               url = JPreferenceDialog.class.getResource("maleScientist64.png");
            if(url != null) {
               final ImageIcon ii = new ImageIcon(url);
               if(ii != null)
                  mBannerLabel.setIcon(ii);
            }
         }
         mBannerLabel.setBorder(bdr);
         add(mBannerLabel, cc.xy(4, 2));
      }
      {
         final JPanel msgPanel = new JPanel(new FormLayout("5dlu, left:pref, 5dlu", "1dlu, center:15dlu, 1dlu"));
         msgPanel.setBorder(bdr);
         msgPanel.add(mMessageLabel, cc.xy(2, 2));
         add(msgPanel, cc.xyw(4, 6, 1));
      }
      {
         final JButton ok = new JButton("OK");
         ok.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(final ActionEvent e) {
               if(commit(mRoot))
                  setVisible(false);
            }
         });
         final JButton cancel = new JButton("Cancel");
         cancel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(final ActionEvent e) {
               setVisible(false);
            }
         });
         final JButton apply = new JButton("Apply");
         apply.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(final ActionEvent e) {
               commit(mRoot);
            }
         });
         final ButtonBarBuilder bbb = new ButtonBarBuilder();
         bbb.addGlue();
         bbb.addButton(ok, cancel, apply);
         add(bbb.build(), cc.xyw(2, 8, 3));
         getRootPane().setDefaultButton(ok);
         setResizable(false);
      }
   }

   private boolean commit(final DefaultMutableTreeNode node) {
      boolean ok = true;
      final Object obj = node.getUserObject();
      if(obj instanceof JPreferencePanel) {
         final JPreferencePanel pp = (JPreferencePanel) obj;
         if(pp.permitExit())
            pp.commit();
         else
            ok = false;
      }
      for(int cx = node.getChildCount() - 1; cx >= 0; --cx)
         ok &= commit((DefaultMutableTreeNode) node.getChildAt(cx));
      return ok;
   }

   /**
    * Adds a new PreferencePanel to the parent TreeNode. If the PreferencePanel
    * has children then these children are added recursively.
    *
    * @param parent
    * @param pp
    * @return The DefaultMutableTreeNode representing the child PreferencePanel
    */
   private DefaultMutableTreeNode addNode(final DefaultMutableTreeNode parent, final JPreferencePanel pp) {
      final DefaultMutableTreeNode newNode = new DefaultMutableTreeNode(pp);
      for(final JPreferencePanel ch : pp.getChildren())
         addNode(newNode, ch);
      parent.add(newNode);
      return newNode;
   }

   /**
    * Recursively checks child tree nodes to find the parent. When the parent is
    * found the child is added to the parent.
    */
   private void addNode(final DefaultMutableTreeNode base, final JPreferencePanel parent, final JPreferencePanel child, final boolean selected) {
      final DefaultMutableTreeNode cn = findPanel(base, parent);
      if(cn != null) {
         final DefaultMutableTreeNode newNode = addNode(cn, child);
         parent.addChildPanel(child);
         mPreferenceTree.setModel(new DefaultTreeModel(mRoot));
         if(selected)
            mPreferenceTree.setSelectionPath(new TreePath(newNode.getPath()));
      }
   }

   /**
    * Adds a child PreferencePanel to the specified parent PreferencePanel.
    *
    * @param parent
    * @param child
    */
   public void addPanel(final JPreferencePanel parent, final JPreferencePanel child, final boolean selected) {
      addNode(mRoot, parent, child, selected);
   }

   /**
    * Adds a preferences panel at the base level.
    *
    * @param pp
    */
   public void addPanel(final JPreferencePanel pp) {
      addNode(mRoot, pp);
      mPreferenceTree.setModel(new DefaultTreeModel(mRoot));
      for(int row = 0; row < mPreferenceTree.getRowCount(); ++row)
         mPreferenceTree.expandRow(row);

   }

   /**
    * Constructs a PreferenceDialog
    *
    * @param owner
    * @throws HeadlessException
    */
   public JPreferenceDialog(final Frame owner)
         throws HeadlessException {
      super(owner, "Preferences", true);
      try {
         initialize();
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   /**
    * Constructs a PreferenceDialog
    *
    * @param owner
    * @param modal
    * @throws HeadlessException
    */
   public JPreferenceDialog(final Frame owner, final boolean modal)
         throws HeadlessException {
      super(owner, "Preferences", modal);
      try {
         initialize();
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   /**
    * Constructs a PreferenceDialog
    *
    * @param owner
    * @param title
    * @param modal
    * @throws HeadlessException
    */
   public JPreferenceDialog(final Frame owner, final String title, final boolean modal)
         throws HeadlessException {
      super(owner, title, modal);
      try {
         initialize();
         pack();
      }
      catch(final Exception ex) {
         ex.printStackTrace();
      }
   }

   public void setMessage(final String msg) {
      mMessageLabel.setText(msg);
   }

   private DefaultMutableTreeNode findPanel(final DefaultMutableTreeNode root, final JPreferencePanel pp) {
      if(root.getUserObject() == pp)
         return root;
      for(int cx = root.getChildCount() - 1; cx >= 0; --cx) {
         final DefaultMutableTreeNode node = findPanel((DefaultMutableTreeNode) root.getChildAt(cx), pp);
         if(node != null)
            return node;
      }
      return null;
   }

   /**
    * Provide a new label for this page. This label is used in the header and in
    * the tree view.
    *
    * @param label
    */
   void updateHeader(final String label) {
      mBannerLabel.setText(label);
      final JTree prefTree = mPreferenceTree;
      final TreePath sel = prefTree.getSelectionPath();
      prefTree.setModel(new DefaultTreeModel(mRoot));
      prefTree.setSelectionPath(sel);
   }
}
