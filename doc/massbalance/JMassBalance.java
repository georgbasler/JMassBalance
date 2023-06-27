/**
 * 
 */
package massbalance;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

/**
 * Graphical User Interface for mass-balanced randomization.
 * 
 * @author Georg Basler
 *
 */
public class JMassBalance extends JPanel implements ActionListener {
	private static final long serialVersionUID = -3531306005220407446L;
	JFrame frame;
	Thread randomize, properties;
	JLabel outputDirLabel, inputDirLabel;
	JTextField networkFileField, outputDirField, nameField, inputDirField;
	JButton networkFileButton, outputDirButton, nameButton, inputDirButton, startButton, stopButton, startButton2, stopButton2;
	JRadioButton massbalancedButton, switchButton, massbalancedButton2, switchButton2;
	
	// randomization components
	JCheckBox fix, compartments, reversible, iterative, strict;
	JSpinner fromSpinner, toSpinner, fromSpinner2, toSpinner2, depthSpinner, probSpinner;
	
	// properties components
	JCheckBox pathLength, clustering, assortativity, degrees, weights, scopes, nCycles;
	JSpinner nCyclesSpinner, scopesSpinner;
	
	JTextArea console;
	boolean randomizeStarted = false, propertiesStarted = false;
	
	/**
	 * GUI main application.
	 * 
	 * @param args
	 */
    public static void main(String[] args) {
        // Schedule a job for the event dispatch thread:
        // creating and showing this application's GUI.
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                // Turn off metal's use of bold fonts
                UIManager.put("swing.boldMetal", Boolean.FALSE); 
                createAndShowGUI();
            }
        });
    }
    
    /**
     * Create the GUI and show it. For thread safety, this method should be invoked
     * from the event dispatch thread.
     */
    private static void createAndShowGUI() {
        // Create and set up the window.
        JFrame frame = new JFrame("JMassBalance");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setPreferredSize(new Dimension(550,650));
        frame.setResizable(false);
        frame.setIconImage(new ImageIcon("icon.jpg").getImage());
        
        // Create and set up the content pane.
        JMassBalance mainPanel = new JMassBalance(frame);
        mainPanel.setOpaque(true); // content panes must be opaque
        frame.setContentPane(mainPanel);
        
//        Dimension preferredSize = mainPanel.getPreferredSize();
//        frame.setPreferredSize(new Dimension(preferredSize.width, preferredSize.height));
        
        // Display the window
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
	
	/**
	 * Constructor creating the components of the GUI.
	 * 
	 * Component hierarchy:
	 * 
	 * JMassBalance (JPanel)
	 * `-- splitPane (JSplitPane)
	 *     |-- tabbedPane (JTabbedPane)
	 *     |   |-- randomizationPanel (JPanel)
	 *     |   |   |-- topPanel (JPanel)
	 *     |   |   |-- middlePanel (JPanel)
	 *     |   |   `-- lowerPanel (JPanel)     
	 *     networkFileField (JTextField)
	 *     outputDirField (JTextField)
	 *     |   |   `-- rTypePanel (JPanel)
	 *     |   |       |-- massbalancedButton (JButton)
	 *     |   |       `-- switchButton (JButton)
	 *     |   `-- propertiesPanel (JPanel)
	 *     |       |-- nameFileField (JTextField)
	 *     |       `-- inputDirField (JTextField)
	 *     `-- consolePane (JScrollPane)
	 *         `-- console (JTextArea)
	 */
	public JMassBalance(JFrame frame) {
		// layout for the single main component (splitPane)
		super(new GridLayout(1,1));
		this.frame = frame;
		
        // Create the console first, because the action listeners may need to refer to it
        console = new JTextArea();
        console.setMargin(new Insets(5,5,5,5));
        console.setEditable(false);
        console.setLineWrap(true);
        console.setWrapStyleWord(true);
        
        // override an output stream for updating the console
		OutputStream redirect = new OutputStream() {
			@Override
			public void write(int b) throws IOException {
				updateConsole(String.valueOf((char) b));
			}
			@Override
			public void write(byte[] b, int off, int len) throws IOException {
				updateConsole(new String(b, off, len));
			}
			@Override
			public void write(byte[] b) throws IOException {
				write(b, 0, b.length);
			}
		};
		// set the standard out and error to the redirecting output stream
		System.setOut(new PrintStream(redirect, true));
		System.setErr(new PrintStream(redirect, true));
        JScrollPane consolePane = new JScrollPane(console);
        
        // create the randomization and properties tab panels
        JPanel randomizationPanel = createRandomizationPanel();
        JPanel propertiesPanel = createPropertiesPanel();
        // set both panel's preferred sizes to the one of the larger
//        int height = Math.max(randomizationPanel.getPreferredSize().height, propertiesPanel.getPreferredSize().height);
//        int width = Math.max(randomizationPanel.getPreferredSize().width, propertiesPanel.getPreferredSize().width);
//        randomizationPanel.setPreferredSize(new Dimension(width, height));
//        propertiesPanel.setPreferredSize(new Dimension(width, height));
        
        // create the tabbed pane from the randomization panel and the properties panel
        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.addTab("Randomization", randomizationPanel);
        tabbedPane.addTab("Properties", propertiesPanel);
        
        // add the components to this panel on a split pane
        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, tabbedPane, consolePane);
        // the split pane should be expandable, but not draggable
        splitPane.setOneTouchExpandable(true);
        randomizationPanel.setMinimumSize(randomizationPanel.getMaximumSize());
        propertiesPanel.setMinimumSize(propertiesPanel.getMaximumSize());
        
        // this is not too clean: the correct height should be used here with a small offset
        splitPane.setDividerLocation((int)(tabbedPane.getLayout().preferredLayoutSize(tabbedPane).height*1.5));
        add(splitPane);

        updateParameterButtons();
        updatePropertyButtons();
	}
	
	/**
	 * Creates the GUI components of the randomization tab.
	 * 
	 * @return The JPanel containing all GUI components.
	 */
	private JPanel createRandomizationPanel() {
		JPanel randomizationPanel = new JPanel();
//		randomizationPanel.setLayout(new BoxLayout(randomizationPanel, BoxLayout.Y_AXIS));
        
        JPanel topPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        
        // create the input file and output directory choosers with action listeners
        JLabel networkFileLabel = new JLabel("Network file ");
        outputDirLabel = new JLabel("Output directory ");
        networkFileField = new JTextField(25);
        outputDirField = new JTextField(25);
        networkFileButton = new JButton("Browse...");
        networkFileButton.addActionListener(this);
        outputDirButton = new JButton("Browse...");
        outputDirButton.addActionListener(this);
        // add the file file choosers to the randomization panel
        gbc.gridx = 0; gbc.gridy = 0;
        topPanel.add(networkFileLabel, gbc);
        gbc.gridx = 1; gbc.gridy = 0;
        topPanel.add(networkFileField, gbc);
        gbc.gridx = 2; gbc.gridy = 0;
        topPanel.add(networkFileButton, gbc);
        gbc.gridx = 0; gbc.gridy = 1;
        topPanel.add(outputDirLabel, gbc);
        gbc.gridx = 1; gbc.gridy = 1;
        topPanel.add(outputDirField, gbc);
        gbc.gridx = 2; gbc.gridy = 1;
        topPanel.add(outputDirButton, gbc);
        randomizationPanel.add(topPanel);
        
        JPanel middlePanel = new JPanel(new GridLayout(1,2));
        massbalancedButton = new JRadioButton("Mass-balanced");
        switchButton = new JRadioButton("Switch");
        fromSpinner = new JSpinner(new SpinnerNumberModel(0, 0, 9999, 1));
        toSpinner = new JSpinner(new SpinnerNumberModel(99, 0, 9999, 1));
        middlePanel.add(createLeftPanel("Randomization", "Networks to generate", massbalancedButton, switchButton, fromSpinner, toSpinner));
        
        // create the check boxes and text fields for the parameters
        JPanel paramPanel = new JPanel(new GridLayout(7,1));
        paramPanel.setBorder(BorderFactory.createTitledBorder("Parameters"));
        fix = new JCheckBox("Fix unbalanced reactions");
        compartments = new JCheckBox("Use compartments");
        reversible = new JCheckBox("All reactions reversible");
        iterative = new JCheckBox("Iterative randomization");
        strict = new JCheckBox("Strict switch randomization");
        // add the listeners
        fix.setSelected(true);
        iterative.setSelected(true);
        strict.setSelected(true);
        fix.addActionListener(this);
        compartments.addActionListener(this);
        reversible.addActionListener(this);
        iterative.addActionListener(this);
        strict.addActionListener(this);
        
        // randomization depth
        JPanel depthPanel = new JPanel(new BorderLayout());
        JLabel depthLabel = new JLabel("Depth ");
        depthSpinner = new JSpinner(new SpinnerNumberModel(1, 0, 100, 0.1));
        depthPanel.add(depthLabel, BorderLayout.WEST);
        depthPanel.add(depthSpinner, BorderLayout.EAST);
        // randomization probability
        JPanel probPanel = new JPanel(new BorderLayout());
        JLabel probLabel = new JLabel("Probability ");
        probSpinner = new JSpinner(new SpinnerNumberModel(1, 0, 1, 0.1));
        // set the preferred size to the preferred size of the larger spinner
        probSpinner.setPreferredSize(depthSpinner.getPreferredSize());
        probPanel.add(probLabel, BorderLayout.WEST);
        probPanel.add(probSpinner, BorderLayout.EAST);
        // add the check boxes and text fields to the parameters panel 
        paramPanel.add(fix);
        paramPanel.add(compartments);
        paramPanel.add(reversible);
        paramPanel.add(iterative);
        paramPanel.add(strict);
        paramPanel.add(depthPanel);
        paramPanel.add(probPanel);
        middlePanel.add(paramPanel);
        randomizationPanel.add(middlePanel);
        
        // create the Start and Stop boxes
        JPanel lowerPanel = new JPanel();
        startButton = new JButton("Start");
        startButton.addActionListener(this);
        stopButton = new JButton("Stop");
        stopButton.addActionListener(this);
        lowerPanel.add(startButton);
        lowerPanel.add(stopButton);
        
        randomizationPanel.add(lowerPanel);
        
        return randomizationPanel;
	}
	
	/**
	 * Creates the GUI components of the properties tab.
	 * 
	 * @return The JPanel containing the components of the properties tab.
	 */
	private JPanel createPropertiesPanel() {
        JPanel propertiesPanel = new JPanel();
//        propertiesPanel.setLayout(new BoxLayout(propertiesPanel, BoxLayout.Y_AXIS));
        
        // create the file panel
        JPanel topPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        
        // create the network name and input directoy choosers with action listeners
        JLabel nameLabel = new JLabel("Network name ");
        inputDirLabel = new JLabel("Input directory ");
        inputDirLabel.setPreferredSize(outputDirLabel.getPreferredSize());
        nameField = new JTextField(25);
        inputDirField = new JTextField(25);
        nameButton = new JButton("Browse...");
        nameButton.addActionListener(this);
        inputDirButton = new JButton("Browse...");
        inputDirButton.addActionListener(this);
        // add the file file choosers to the randomization panel
        gbc.gridx = 0; gbc.gridy = 0;
        topPanel.add(nameLabel, gbc);
        gbc.gridx = 1; gbc.gridy = 0;
        topPanel.add(nameField, gbc);
        gbc.gridx = 2; gbc.gridy = 0;
        topPanel.add(nameButton, gbc);
        gbc.gridx = 0; gbc.gridy = 1;
        topPanel.add(inputDirLabel, gbc);
        gbc.gridx = 1; gbc.gridy = 1;
        topPanel.add(inputDirField, gbc);
        gbc.gridx = 2; gbc.gridy = 1;
        topPanel.add(inputDirButton, gbc);
        propertiesPanel.add(topPanel);
        
        JPanel middlePanel = new JPanel(new GridLayout(1,2));
        massbalancedButton2 = new JRadioButton("Mass-balanced");
        switchButton2 = new JRadioButton("Switch");
        fromSpinner2 = new JSpinner(new SpinnerNumberModel(0, 0, 9999, 1));
        toSpinner2 = new JSpinner(new SpinnerNumberModel(99, 0, 9999, 1));
        middlePanel.add(createLeftPanel("Randomized networks", "Networks to analyze", massbalancedButton2, switchButton2, fromSpinner2, toSpinner2));
        
        // create the check boxes and text fields with action listeners for the parameters
        JPanel paramPanel = new JPanel(new GridLayout(7,1));
        paramPanel.setBorder(BorderFactory.createTitledBorder("Properties"));
        
        // create the check boxes and spinner for the properties to calculate
        pathLength = new JCheckBox("Average path length");
        clustering  = new JCheckBox("Clustering coefficient");
        assortativity = new JCheckBox("Assortativities");
        degrees = new JCheckBox("Degree distribution");
        weights = new JCheckBox("Weights distribution");
        // scopes panel
        JPanel scopesPanel = new JPanel(new BorderLayout());
        scopes = new JCheckBox("Scope distribution, seed");
        scopesSpinner = new JSpinner(new SpinnerNumberModel(4, 1, 256, 1));
        scopesSpinner.setEnabled(false);
        scopesPanel.add(scopes, BorderLayout.WEST);
        scopesPanel.add(scopesSpinner, BorderLayout.EAST);
        // n-cycles panel
        JPanel nCyclesPanel = new JPanel(new BorderLayout());
        nCycles = new JCheckBox("Cycles of length");
        nCyclesSpinner = new JSpinner(new SpinnerNumberModel(2, 2, 20, 1));
        // set the preferred size to the preferred size of the larger spinner
        nCyclesSpinner.setPreferredSize(scopesSpinner.getPreferredSize());
        nCyclesSpinner.setEnabled(false);
        nCyclesPanel.add(nCycles, BorderLayout.WEST);
        nCyclesPanel.add(nCyclesSpinner, BorderLayout.EAST);
        // add the listeners
        pathLength.addActionListener(this);
        clustering.addActionListener(this);
        assortativity.addActionListener(this);
        degrees.addActionListener(this);
        weights.addActionListener(this);
        scopes.addActionListener(this);
        nCycles.addActionListener(this);
        
        // add the check boxes and text fields to the parameters panel 
        paramPanel.add(pathLength);
        paramPanel.add(clustering);
        paramPanel.add(assortativity);
        paramPanel.add(degrees);
        paramPanel.add(weights);
        paramPanel.add(scopesPanel);
        paramPanel.add(nCyclesPanel);
        
        middlePanel.add(paramPanel);
        propertiesPanel.add(middlePanel);
        
        // create the Start and Stop boxes
        JPanel lowerPanel = new JPanel();
        startButton2 = new JButton("Start");
        startButton2.addActionListener(this);
        stopButton2 = new JButton("Stop");
        stopButton2.addActionListener(this);
        lowerPanel.add(startButton2);
        lowerPanel.add(stopButton2);
        
        propertiesPanel.add(lowerPanel);
        
        return propertiesPanel;
	}
	
	/**
	 * Creates the left-side panel, containing the randomization type radio buttons
	 * and the network range spinners.
	 * 
	 * @param randomizationLabel
	 * @param fromToLabel
	 * @param massbalancedButton
	 * @param switchButton
	 * @param fromSpinner
	 * @param toSpinner
	 * @return The JPanel containing the the randomization type radio buttons
	 * and the network range spinners.
	 */
	private JPanel createLeftPanel(String randomizationLabel, String fromToLabel, JRadioButton massbalancedButton, JRadioButton switchButton, JSpinner fromSpinner, JSpinner toSpinner) {
		
        JPanel leftPanel = new JPanel(new GridLayout(2,1));
        
        // create the radio buttons with for the randomization type
        JPanel rTypePanel = new JPanel(new GridLayout(2,1));
        rTypePanel.setBorder(BorderFactory.createTitledBorder(randomizationLabel));
        massbalancedButton.setSelected(true);
        massbalancedButton.addActionListener(this);
        switchButton.addActionListener(this);
        ButtonGroup rGroup = new ButtonGroup();
        rGroup.add(massbalancedButton);
        rGroup.add(switchButton);
        // add the radio buttons to the randomization panel
        rTypePanel.add(massbalancedButton);
        rTypePanel.add(switchButton);
        leftPanel.add(rTypePanel);
        
        // create the text fields for the index range of networks to generate
        JPanel fromToPanel = new JPanel();
//        fromToPanel.setLayout(new BoxLayout(fromToPanel, BoxLayout.X_AXIS));
        fromToPanel.setBorder(BorderFactory.createTitledBorder(fromToLabel));
        JLabel fromLabel = new JLabel("From ");
        JLabel toLabel = new JLabel(" to ");
        fromToPanel.add(fromLabel);
        fromToPanel.add(fromSpinner);
        fromToPanel.add(toLabel);
        fromToPanel.add(toSpinner);
        leftPanel.add(fromToPanel);
        
        return leftPanel;
	}
	
	/**
	 * Event handling for this ActionListener. ActionEvents are generated from the
	 * GUI components and trigger this method to be executed. Depending on the source
	 * GUI component, the appropriate action is performed. 
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		
		// perform the appropriate action depending on the source of the event
		
		// Randomization tab
		
		if (e.getSource() == networkFileButton) {
			// open the network file chooser
			JFileChooser fileChooser = new JFileChooser(networkFileField.getText());
			int returnVal = fileChooser.showOpenDialog(JMassBalance.this);
            if (returnVal == JFileChooser.APPROVE_OPTION)
                networkFileField.setText(fileChooser.getSelectedFile().getAbsolutePath());
            
		} else if (e.getSource() == outputDirButton) {
			// open the output directory chooser
			JFileChooser fileChooser = new JFileChooser(outputDirField.getText());
			fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = fileChooser.showOpenDialog(JMassBalance.this);
            if (returnVal == JFileChooser.APPROVE_OPTION)
                outputDirField.setText(fileChooser.getSelectedFile().getAbsolutePath());
            
		} else if (e.getSource() == massbalancedButton || e.getSource() == switchButton) {
			// enable/disable parameters according to selected randomization type
			updateParameterButtons();
            
		} else if (e.getSource() == startButton) {
			// extract the parameters from the GUI components
			String networkFile = networkFileField.getText();
			String outputDir = outputDirField.getText();
			String name = "";
			int from = (Integer)fromSpinner.getValue();
			int to = (Integer)toSpinner.getValue();
			double depth = (Double)depthSpinner.getValue();
			double prob = (Double)probSpinner.getValue();
			
			// check parameters for validity
			// JSpinner values depth and prob are enforced valid by the GUI
			if (!new File(networkFile).isFile()) {
				JOptionPane.showMessageDialog(frame, "Error opening "+networkFile+".", "Invalid network file", JOptionPane.ERROR_MESSAGE);
				return;
			}
			if (new File(outputDir).isFile()) {
				JOptionPane.showMessageDialog(frame, "Not a directory: "+outputDir+".", "Invalid output directory", JOptionPane.ERROR_MESSAGE);
				return;
			}
			if (from > to) {
				JOptionPane.showMessageDialog(frame, "First index (from) must not be larger than second index (to).", "Invalid network range", JOptionPane.ERROR_MESSAGE);
				return;
			}
			
			// set the network name and input directory of the properties tab accordingly
			File file = new File(networkFile);
			try {
				name = Utilities.getVersion(file.getAbsolutePath());
			} catch (IOException e1) {
				e1.printStackTrace();
				return;
			}
			nameField.setText(name);
			inputDirField.setText(outputDir);
			
			// set the instance parameters
			Randomize r = new Randomize();
			r.networkFile = networkFile;
			r.outputDir = outputDir;
			r.massbalance = massbalancedButton.isSelected();
			r.switchRand = switchButton.isSelected();
			r.from = from;
			r.to = to;
			r.nofix = !fix.isSelected();
			r.compartments = compartments.isSelected();
			r.reversible = reversible.isSelected();
			r.noniterative = !iterative.isSelected();
			r.strict = strict.isSelected();
			r.depth = (float)depth;
			r.p = (float)prob;
			
			// clear the console
			console.setText("");
			
			// start randomization in a concurrent thread
			randomizeStarted = true;
			randomize = new Thread(r);
			randomize.setPriority(Thread.currentThread().getPriority()-1);
			randomize.start();
			
		} else if (e.getSource() == stopButton) {
			if (randomizeStarted && randomize != null) {
				randomizeStarted = false;
				console.append("\nStopping...");
				// interrrupt the randomization thread
				randomize.interrupt();
			}

		// Properties tab
			
		} else if (e.getSource() == nameButton) {
			// open the network file chooser
			JFileChooser fileChooser = new JFileChooser(nameField.getText());
			int returnVal = fileChooser.showOpenDialog(JMassBalance.this);
            if (returnVal == JFileChooser.APPROVE_OPTION)
            	nameField.setText(fileChooser.getSelectedFile().getAbsolutePath());
            
		} else if (e.getSource() == inputDirButton) {
			// open the output directory chooser
			JFileChooser fileChooser = new JFileChooser(inputDirButton.getText());
			fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int returnVal = fileChooser.showOpenDialog(JMassBalance.this);
            if (returnVal == JFileChooser.APPROVE_OPTION)
            	inputDirField.setText(fileChooser.getSelectedFile().getAbsolutePath());
            
		} else if (e.getSource() == massbalancedButton2 || e.getSource() == switchButton2) {
			// enable/disable parameters according to selected randomization type
			updatePropertyButtons();
			
		} else if (e.getSource() == scopes) {
			scopesSpinner.setEnabled(scopes.isSelected());

		} else if (e.getSource() == nCycles) {
			nCyclesSpinner.setEnabled(nCycles.isSelected());
			
		} else if (e.getSource() == startButton2) {
			
			// get version name
			
			// extract the parameters from the GUI components
			String name = nameField.getText();
			String inputDir = inputDirField.getText();
			int from = (Integer)fromSpinner2.getValue();
			int to = (Integer)toSpinner2.getValue();
			int seedSize = (Integer)scopesSpinner.getValue();
			int cycleLength = (Integer)nCyclesSpinner.getValue();
			
			// check parameters for validity
			// JSpinner values seedSize and cycleLength are enforced valid by the GUI
			File file = new File(nameField.getText());
			if (file.exists()) {
				try {
					name = Utilities.getVersion(file.getAbsolutePath());
				} catch (IOException e1) {
					// try to use the text in the name field
				}
			}
			if (!new File(inputDir).isDirectory()) {
				JOptionPane.showMessageDialog(frame, "Error accessing "+inputDir+".", "Invalid input directory", JOptionPane.ERROR_MESSAGE);
				return;
			}
			if (from > to) {
				JOptionPane.showMessageDialog(frame, "First index (from) must not be larger than second index (to).", "Invalid network range", JOptionPane.ERROR_MESSAGE);
				return;
			}
			
			// set the instance parameters
			Properties p = new Properties();
			p.version = name;
			p.inputDir = inputDir;
			p.massbalance = massbalancedButton2.isSelected();
			p.switchRand = switchButton2.isSelected();
			p.from = from;
			p.to = to;
			p.pathLength = pathLength.isSelected();
			p.cluster = clustering.isSelected();
			p.assortativity = assortativity.isSelected();
			p.degrees = degrees.isSelected();
			p.weights = weights.isSelected();
			p.scopes = scopes.isSelected();
			p.cycles = nCycles.isSelected();
			p.seedSizes = new int[]{seedSize};
			p.nCycles = new ArrayList<Integer>(1);
			p.nCycles.add(cycleLength);
			
			// clear the console
			console.setText("");
			
			// start calculating the properties in a concurrent thread
			propertiesStarted = true;
			properties = new Thread(p);
			properties.setPriority(Thread.currentThread().getPriority()-1);
			properties.start();
			
		} else if (e.getSource() == stopButton2) {
			if (propertiesStarted && properties != null) {
				propertiesStarted = false;
				console.append("\nStopping...");
				// interrrupt the randomization thread
				properties.interrupt();
			}
			
		}
	}
	
	/**
	 * Enables/disables the randomization parameter buttons according
	 * to the currently selected randomization type.
	 */
	private void updateParameterButtons() {
		boolean massbalanced = massbalancedButton.isSelected();
		iterative.setEnabled(massbalanced);
		strict.setEnabled(!massbalanced);
		depthSpinner.setEnabled(massbalanced);
		probSpinner.setEnabled(massbalanced);
	}
	
	/**
	 * Enables/disables the properties buttons according
	 * to the currently selected randomization type.
	 */
	private void updatePropertyButtons() {
		boolean massbalanced = massbalancedButton2.isSelected();
		degrees.setEnabled(massbalanced);
		weights.setEnabled(massbalanced);
	}
	
	/**
	 * Updates the console from the redirected standard and error
	 * output streams.
	 * 
	 * @param text
	 */
	private void updateConsole(final String text) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				console.append(text);
			}
		});
	}
}
