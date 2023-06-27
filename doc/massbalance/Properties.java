/**
 * 
 */
package massbalance;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.sbml.libsbml.ListOfReactions;
import org.sbml.libsbml.ListOfSpecies;
import org.sbml.libsbml.ListOfSpeciesTypes;
import org.sbml.libsbml.Model;
import org.sbml.libsbml.Reaction;
import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.Species;
import org.sbml.libsbml.SpeciesReference;

/**
 * Class for parsing, graph creation and calculating properties of randomized networks.
 * 
 * @author Georg Basler
 *
 */
public class Properties extends Thread {

	/**
	 * Static initializer for loading the configuration parameters.
	 */
	static {
		Utilities.loadConfiguration();
	}
	
	static final String OPTIONS =
		"massbalance\t\tCalculate properties of mass-balanced randomized networks.\n"
		+"switch\t\t\tCalculate properties of switch randomized networks.\n"
		+"Properties:\n"
		+"matrix[-labelled][-sorted][-reversible]\tPrint the stoichiometric matrices.\n"
		+"massMatrix[-labelled]\tPrint a matrix of mass vectors.\n"
		+"degrees\t\t\tPrint the degree distributions.\n"
		+"weights\t\t\tPrint the weight distributions (stoichiometric coefficients).\n"
		+"scopes[=seedsize]\tCalculate the scope size distributions with the given seed size (default: 1, 2, 4, 8, 16, and 32).\n"
		+"pathLength\t\tCalculate the characteristic path lengths.\n"
		+"clustering\t\tCalculate the clustering coefficients.\n"
		+"n-cycles*\t\tCalculate the number of cycles of length n.\n"
		+"path[-directed]=compounds...\tTest whether an (un-)directed path consisting of the given compounds exists.\n"
		+"connected=compounds...\tTest whether the compounds are connected by (un-)directed paths.\n"
		+"deltaGf\t\t\tPrint the standard Gibbs free energy of all compounds (requires a deltaG annotation file).\n"
		+"deltaGr\t\t\tPrint the standard Gibbs free energy distributions of reactions (requires a deltaG annotation file).\n"
		+"deltaGn\t\t\tPrint the sum of standard Gibbs free energy distributions of reactions (requires a deltaG annotation file).\n"
		+"assortativity\t\tCalculate the in-degree-in-degree and out-degree-out-degree correlation coefficients of compounds.\n"
		+"transitionDegree\tCalculate the degree of the network node in the transition graph, i.e., the number of possible unique substitutions.\n"
		+"matrixSubstitutions[-reversible]\tPrint all possible substitutions of the stoichiometric matrix to a tab-separated file.\n"
		+"equationSubstitutions\tPrint all possible substitutions as reaction equations to a tab-separated file.\n"
		+"centrality[-reversible][-double],d=[0,1]\tCalculate the global reaction knockout centralities. Reversible: consider reversible reactions as one; double: double precision; d: damping factor.\n"
		+"localCentrality[-reversible]\tCalculate the local reaction knockout centralities. Reversible: consider reversible reactions as one.\n"
		+"knockoutSet[,e=]0,1]],reaction\tDetermine the set of reactions affected by the given reaction's knockout up to the specified strength, e.\n"
		+"D=vertices...\tPreserve the given reactions/compounds.\n"
		+"write\t\t\tExport the networks to a tab-delimited text file.\n"
		+"write2\t\t\tExport the networks to a tab-delimited text file (alternative format).\n"
		+"from-to\t\t\tIndices of randomized networks to read and parallelization. (default: 0-999).\n"
		+"See Documentation for details.";
	
	static final String USAGE = "Usage: <name> <inputDir> [options]\nOptions:\n"
		+"<name>\t\t\tThe network name generated during randomization.\n"
		+"<inputDir>\t\tDirectory specified during randomization.\n"
		+"Options:\n"+OPTIONS;
	
	private static final int NUMREACTIONS = 0, TRANSPORTER = 1, DUPLICATES = 2, CYCLIC = 3, ZERODEGREE = 4, UNANNOTATED = 5, UNBALANCED = 6;
	
	protected String inputDir = "", version = "";
	protected boolean append = true, massbalance = false, switchRand = false, matrix = false, matrixReversible = false, matrixLabelled = false, matrixSorted = false, massMatrix = false, massMatrixLabelled = false, degrees = false, reactionDegrees = false, weights = false, scopes = false, pathLength = false, cycles = false, isPath = false, directed = false, connected = false, cluster = false, deltaGf = false, deltaGn = false, deltaGr = false, assortativity = false, transitionDegree = false, matrixSubstitutions = false, equationSubstitutions = false, matrixSubstitutionsReversible = false, centrality = false, centrality_reversible = false, centrality_double = false, localCentrality = false, localCentrality_reversible = false, knockoutSet = false, preserve = false, write = false, write2 = false;
	protected int from = 0, to = 999;
	protected int[] seedSizes = {1, 2, 4, 8, 16, 32};
	protected String knockout = null;
	protected double knockoutThreshold = 1d, dampingFactor = 0.85d;
	protected ArrayList<Integer> nCycles = new ArrayList<Integer>();
	protected HashMap<ArrayList<Integer>, EquivalenceClass> classes = null;
	protected ArrayList<String> path = null, intermediaries = null, compartments = null;
	protected ArrayList<Boolean> reversible = null;
	protected HashSet<String> preserved;
	
	/**
	 * Method for calculating graph properties. Called from the command line using 'java -cp massbalance.jar massbalance/Properties'.<br>
	 * <br>
     * <strong>Usage:</strong><br>
     * <code>name inputDir [options]</code><br>
	 * <br>
	 * <table>
	 * <tr><td><code>inputDir</code></td>		<td>Input directory specified during randomization.</td></tr>
	 * <tr><td><code>name</code></td>		<td>The network name generated during randomization.</td></tr>
	 * </table>
	 * <br>
	 * <strong>Options:</strong>
	 * <table>
	 * <tr><td><code>massbalance</code></td>		<td>Calculate properties of mass-balanced randomized networks.</td></tr>
	 * <tr><td><code>switch</code></td>			<td>Calculate properties of switch randomized networks.</td></tr>
	 * </table>
	 * <br>
	 * <strong>Properties:</strong>
	 * <table>
	 * <tr><td><code>matrix[-labelled][-sorted][-reversible]</code></td>			<td>Print the stoichiometric matrices.</td></tr>
	 * <tr><td><code>massMatrix[-labelled]</code></td>			<td>Print a matrix of mass vectors.</td></tr>
	 * <tr><td><code>degrees</code></td>			<td>Print the degree distributions.</td></tr>
	 * <tr><td><code>weights</code></td>			<td>Print the weight distributions (stoichiometric coefficients).</td></tr>
	 * <tr><td><code>scopes[=seedsize]</code></td>	<td>Calculate the scope size distributions with the given seed size (default: 1, 2, 4, 8, 16, and 32).</td></tr>
	 * <tr><td><code>pathLength</code></td>		<td>Calculate the characteristic path lengths.</td></tr>
	 * <tr><td><code>clustering</code></td>		<td>Calculate the clustering coefficients.</td></tr>
	 * <tr><td><code>n-cycles*</code></td>		<td>Calculate the number of cycles of length n.</td></tr>
	 * <tr><td><code>path[-directed]=compounds...</code></td>	<td>Test whether an (un-)directed path consisting of the given compounds exists.</td></tr>
	 * <tr><td><code>connected=compounds...</code></td>		<td>Test whether the compounds are connected by (un-)directed paths.</td></tr>
	 * <tr><td><code>deltaGf</code></td>			<td>Print the standard Gibbs free energy of all compounds (requires a deltaG annotation file).</td></tr>
	 * <tr><td><code>deltaGr</code></td>			<td>Print the standard Gibbs free energy distributions of reactions (requires a deltaG annotation file).</td></tr>
	 * <tr><td><code>deltaGn</code></td>			<td>Print the sum of standard Gibbs free energy distributions of reactions (requires a deltaG annotation file).</td></tr>
	 * <tr><td><code>assortativity</code></td>		<td>Calculate the in-degree-in-degree and out-degree-out-degree correlation coefficients of compounds.</td></tr>
	 * <tr><td><code>transitionDegree</code></td>	<td>Calculate the degree of the network node in the transition graph, i.e., the number of possible unique substitutions.</td></tr>
	 * <tr><td><code>matrixSubstitutions[-reversible]</code></td>	<td>Print all possible substitutions of the stoichiometric matrix to a tab-separated file.</td></tr>
	 * <tr><td><code>equationSubstitutions</code></td>	<td>Print all possible substitutions as reaction equations to a tab-separated file.</td></tr>
	 * <tr><td><code>centrality[-reversible][-double],d=[0,1]</code></td>	<td>Calculate the global reaction knockout centralities. Reversible: consider reversible reactions as one; double: double precision; d: damping factor.</td></tr>
	 * <tr><td><code>localCentrality[-reversible]</code></td>	<td>Calculate the local reaction knockout centralities. Reversible: consider reversible reactions as one.</td></tr>
	 * <tr><td><code>knockoutSet[,e=]0,1]],reaction</code></td>	<td>Determine the set of reactions affected by the given reaction's knockout up to the specified strength, e.</td></tr>
	 * <tr><td><code>D=vertices...</code></td><td>	Preserve the given reactions/compounds.</td></tr> 
	 * <tr><td><code>write</code></td>			<td>Export the networks to a tab-delimited text file.</td></tr>
	 * <tr><td><code>write2</code></td>			<td>Export the networks to a tab-delimited text file (alternative format).</td></tr>
	 * <tr><td><code>from-to</code></td>			<td>Indices of randomized networks to read and parallelization. (default: 0-999).</td></tr>
	 * </table>
	 */
	public static void main(String[] args) throws InterruptedException {
		
		Properties properties = new Properties();
		properties.parseArgs(args);
		
		// run randomization
		Thread t = new Thread(properties);
		t.start();
		
	}
	
	public void parseArgs(String[] args) {		
		// argument check
		if (args.length < 2) {
			System.out.println(USAGE);
			System.exit(-1);
		}
		if (args[1].endsWith(File.separator))
			args[1] = args[0].substring(0, args[1].length()-File.separator.length());
		if (!(new File(args[1]).isDirectory())) {
			System.err.println("Not a directory: "+args[1]);
			System.out.println(USAGE);
			System.exit(-1);
		}
		
		int i=0;
		version = args[i++];
		inputDir = args[i++];
		
		try {
			while (args.length > i) {
				if (args[i].equals("massbalance"))
					massbalance = true;
				else if (args[i].equals("switch"))
					switchRand = true;
				else if (args[i].startsWith("matrix") && !args[i].startsWith("matrixSubstitutions")) {
					matrix = true;
					if (args[i].contains("-reversible"))
						matrixReversible = true;
					if (args[i].contains("-labelled"))
						matrixLabelled = true;
					if (args[i].contains("-sorted"))
						matrixSorted = true;
				} else if (args[i].startsWith("massMatrix")) {
					massMatrix = true;
					if (args[i].contains("-labelled"))
						massMatrixLabelled = true;
				} else if (args[i].equals("degrees"))
					degrees = true;
				else if (args[i].equals("reactiondegrees"))
					reactionDegrees = true;
				else if (args[i].equals("weights"))
					weights = true;
				else if (args[i].startsWith("scopes")) {
					scopes = true;
					if (args[i].contains("="))
						seedSizes = new int[] {Integer.parseInt(args[i].substring(args[i].indexOf("=")+1, args[i].length()))};
				} else if (args[i].equals("pathLength"))
					pathLength = true;
				else if (args[i].equals("clustering"))
					cluster = true;
				else if (args[i].equals("deltaGf"))
					deltaGf= true;								
				else if (args[i].equals("deltaGr"))
					deltaGr = true;				
				else if (args[i].equals("deltaGn"))
					deltaGn = true;
				else if (args[i].equals("assortativity"))
					assortativity = true;
				else if (args[i].equals("transitionDegree"))
					transitionDegree = true;
				else if (args[i].startsWith("matrixSubstitutions")) {
					matrixSubstitutions = true;
					if (args[i].contains("-reversible"))
						matrixSubstitutionsReversible = true;
				} else if (args[i].startsWith("equationSubstitutions")) {
						equationSubstitutions = true;
				} else if (args[i].startsWith("centrality")) {
					centrality = true;
					if (args[i].contains("-double"))
						centrality_double = true;
					if (args[i].contains("-reversible"))
						centrality_reversible = true;
					if (args[i].contains("."))
						dampingFactor = Double.parseDouble(args[i].substring(args[i].indexOf(".")-1, args[i].length()));
				} else if (args[i].startsWith("localCentrality")) {
					localCentrality = true;
					if (args[i].contains("-double"))
						localCentrality_reversible = true;
				} else if (args[i].startsWith("knockoutSet,")) {
					knockoutSet = true;
					// parse the centrality threshold
					int init=13;
					if (args[i].substring(13, 15).equals("e=")) {
						init = args[i].indexOf(",", 16)+1;
						knockoutThreshold = Double.parseDouble(args[i].substring(15, init-1));
					}
					// parse the reaction name
					knockout = args[i].substring(init);
					// parse the reaction strings into the ArrayList
//					knockouts = new ArrayList<String>();
//					int end = args[i].indexOf(",", init);
//					for (; end!=-1; end=args[i].indexOf(",", init)) {
//						knockouts.add(args[i].substring(init, end));
//						init = end+1;
//					}
//					knockouts.add(args[i].substring(init));
				} else if (args[i].equals("write"))
					write = true;
				else if (args[i].equals("write2"))
					write2 = true;				
				else if (args[i].startsWith("path")) {
					isPath = true;
					if (args[i].substring(10, 19).equals("-directed"))
						directed = true;
					path = new ArrayList<String>();
					for (int j=args[i].indexOf("=")+1; j<args[i].length();) {
						int end = (args[i].indexOf(",", j) > -1 ? args[i].indexOf(",", j) : args[i].length());
						path.add(args[i].substring(j, end));
						j = end+1;
					}
				} else if (args[i].startsWith("connected")) {
					connected = true;
					intermediaries = new ArrayList<String>();
					reversible = new ArrayList<Boolean>();
					compartments = new ArrayList<String>();
					for (int j=args[i].indexOf("=")+1; j<args[i].length();) {
						// get the index of the next irreversible character '>' or reversible character '='
						int irrev = (args[i].indexOf(">", j) > -1 ? args[i].indexOf(">", j) : args[i].length());
						int rev = (args[i].indexOf("=", j) > -1 ? args[i].indexOf("=", j) : args[i].length());
						int end = (irrev < rev ? irrev : rev);
						// set this interaction to reversible or irreversible
						if (irrev < args[i].length() || rev < args[i].length())
							reversible.add(rev < irrev);
						// extract the compartment and add the compound
						if (args[i].substring(j, end).contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
							intermediaries.add(args[i].substring(j, args[i].indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT, j)));
							compartments.add(args[i].substring(args[i].indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT, j)+1, end));
						} else {
							intermediaries.add(args[i].substring(j, end));
							compartments.add(null);
						}
						j = end+1;
					}
//					String[] tcaIntermediaries = new String[]{"accoa", "cit", "icit", "akg", "succoa", "succ", "fum", "mal-L", "oaa", "cit"};
//					boolean[] reversible = new boolean[]{false, true, true, false, true, false, true, true, false};
				} else if (args[i].endsWith("cycles")) {
					// parse the cycle-lengths
					nCycles.add(Integer.parseInt(args[i].substring(0, 1)));
					cycles = true;
				} else if (args[i].startsWith("D=")) {
					preserve = true;
					preserved = new HashSet<String>();
					String vertexNames = args[i].substring(args[i].indexOf("=")+1);
					for (int j=0; j<vertexNames.length();) {
						// get the index of the next separator ','
						int end = (vertexNames.indexOf(",",j) > -1 ? vertexNames.indexOf(",",j) : vertexNames.length());
						// extract vertex names
						preserved.add(vertexNames.substring(j, end));
						j = end+1;
					}				
				} else if (args[i].contains("-") && Character.isDigit(args[i].charAt(0)) && Character.isDigit(args[i].charAt(args[i].length()-1))) {
					// parse the start and end indices to use
					from = Integer.parseInt(args[i].substring(0, args[i].indexOf("-")));
					to = Integer.parseInt(args[i].substring(args[i].indexOf("-")+1, args[i].length()));
				} else {
					System.err.println("Unrecognized argument: "+args[i]);
					System.out.println(USAGE);
					System.exit(-1);
				}
				i++;
			}
			if (switchRand && !massbalance && degrees) {
				System.out.println("Degrees not applicable to switch randomization.");
				System.exit(-1);
			}
			if (switchRand && !massbalance && weights) {
				System.out.println("Weights not implemented for switch randomization.");
				System.exit(-1);
			}
		} catch (NumberFormatException e) {
			System.err.println("Unrecognized number format: "+args[i]);
			System.out.println(USAGE);
			System.exit(-1);
		}
		
	}
	
	@Override
	public synchronized void run() {
		long start = System.currentTimeMillis();
		long time = 0;
		
		MetabolicGraph originalGraph = null;
		String graphPrefix = inputDir+File.separator+version+"-graphs"+File.separator+version;
		String matrixPrefix = inputDir+File.separator+version+"-matrices"+File.separator+version;
		String outputPrefix = inputDir+File.separator+version;
		HashSet<Vertex> preservedV = new HashSet<Vertex>();
		
		// add the index range to the output file name
		outputPrefix += "." + from + "-" + to;
		
		try {
			// allow to interrupt this thread
			wait(1);
			
			System.out.println("Writing properties..."); 
			
			if (transitionDegree || matrixSubstitutions || equationSubstitutions) {
				// read equivalence classes from file
				String classesFile = inputDir+File.separator+version+".classes";
				System.out.print("Reading classes from "+classesFile+". ");
				time = System.currentTimeMillis();
				classes = Equivalence.readClasses(classesFile);
				System.out.println(" ("+(System.currentTimeMillis()-time)/1000 + " seconds).");
			}
			
			if (from == 0)
				System.out.println("Original graph:");
			
			if (from == 0 || centrality) {
				// read the original graph
				ObjectInputStream graphReader = new ObjectInputStream(new FileInputStream(graphPrefix));
				originalGraph = (MetabolicGraph)graphReader.readObject();
				graphReader.close();
			}
			
			if (from == 0) {
				// print the connected components and reaction types
				originalGraph.connectedComponents(true);
				Utilities.countReactions(originalGraph, true);
				
				// first get the set of preserved vertices
				if (preserve) {
					for (String name : preserved) {
						String vertex = name, compartment = null;
						if (name.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
							vertex = name.substring(0, name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
							compartment = name.substring(name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1);
						}
						Vertex v = originalGraph.getReaction(vertex);
						if (v == null)
							v = originalGraph.getCompound(vertex, compartment);
						preservedV.add(v);
					}
				}
				
				if (matrix) {
					if (!new File(inputDir+File.separator+version+"-matrices").exists())
						new File(inputDir+File.separator+version+"-matrices").mkdir();
//					graph.adjacencyMatrix(matrixPrefix+".matrix");
					originalGraph.stoichiometricMatrix(matrixPrefix+".stmatrix", matrixReversible, matrixLabelled, matrixSorted, true);
				}
				if (massMatrix) {
					if (!new File(inputDir+File.separator+version+"-matrices").exists())
						new File(inputDir+File.separator+version+"-matrices").mkdir();
					originalGraph.massMatrix(matrixPrefix+".massMatrix", massMatrixLabelled);
				}
				if (massbalance && degrees) {
					originalGraph.compoundDegrees(outputPrefix+".massbalance.degrees", false);
					if (!originalGraph.isReversible()) {
						originalGraph.compoundInOutDegrees(outputPrefix+".massbalance.degrees.in", true, false);
						originalGraph.compoundInOutDegrees(outputPrefix+".massbalance.degrees.out", false, false);
					}
				}
				if (massbalance && reactionDegrees)
					originalGraph.reactionDegrees(outputPrefix+".massbalance.reactiondegrees", false);
				if (massbalance && weights)
					originalGraph.weights(outputPrefix+".massbalance.weights", false);
				if ((massbalance || switchRand) && scopes) {
					if (massbalance)
						originalGraph.randomScopeSizes(seedSizes, outputPrefix+".massbalance.scopes", (switchRand ? outputPrefix+".switch.scopes" : null), false);
					else
						originalGraph.randomScopeSizes(seedSizes, outputPrefix+".switch.scopes", null, false);
				}
				if (pathLength) {
					float diam = originalGraph.characteristicPathLength();
					System.out.println("Characteristic path length: "+diam);
					if (massbalance) {
						BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.pathLength", false));
						writer.write(String.valueOf(diam)+"\n");
						writer.close();
					}
					if (switchRand) {
						BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".switch.pathLength", false));
						writer.write(String.valueOf(diam)+"\n");
						writer.close();
					}
				}
				if (cycles) {
					// get the strongly connected component
					StrongConnectivityInspector<Vertex, DefaultWeightedEdge> connectivity = new StrongConnectivityInspector<Vertex, DefaultWeightedEdge>(originalGraph);
					List<Set<Vertex>> components = connectivity.stronglyConnectedSets();
					// count the n-cycles
					for (Iterator<Integer> nCyclesIt = nCycles.iterator(); nCyclesIt.hasNext();) {
						int n = nCyclesIt.next();
						time = System.currentTimeMillis();
						int cycleCount = originalGraph.cycleCount(components, n, true, n==nCycles.get(0));
						System.out.println("Counting of "+cycleCount+" "+n+"-cycles in the original graph took "+(System.currentTimeMillis()-time)+" ms.");
						time = System.currentTimeMillis();
						if (massbalance) {
							BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".massbalance."+n+"cycles", false));
							writer.write(String.valueOf(cycleCount)+"\n");
							writer.close();
						}
						if (switchRand) {
							BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".switch."+n+"cycles", false));
							writer.write(String.valueOf(cycleCount)+"\n");
							writer.close();						
						}
					}
				}
				if (isPath) {
					boolean hasPath = originalGraph.isPath(path, !directed);
					System.out.println((directed?"Directed":"Undirected")+" path "+(hasPath?"":"not ")+"found in "+version+".");
					// write 1 if the path exist, 0 otherwise to the original, massbalance and/or switch files
					for (int j=0; j<1+(massbalance && switchRand ? 1:0); j++) {
						BufferedWriter writer = new BufferedWriter(new FileWriter((j==0 && massbalance ? outputPrefix+".massbalance" : (j==1 || switchRand ? outputPrefix+".switch" : inputDir+File.separator+version))+".path", false));
						writer.write((hasPath ? "1":"0")+"\n");
						writer.close();
					}
				}
				if (connected) {
					// check if the network contains the specified path
					boolean[] isConnected = originalGraph.isConnected(intermediaries, compartments, reversible);
					// write the connections to the original, massbalance and/or switch files
					for (int j=0; j<1+(massbalance && switchRand ? 1:0); j++) {
						BufferedWriter writer = new BufferedWriter(new FileWriter((j==0 && massbalance ? outputPrefix+".massbalance" : (j==1 || switchRand ? outputPrefix+".switch" : inputDir+File.separator+version))+".connected", false));
						for (int k=0; k<isConnected.length; k++) {
							if (k<isConnected.length-1)
								writer.write((isConnected[k]?"1":"0")+"\t");
							else
								writer.write((isConnected[k]?"1":"0")+"\n");
						}
						writer.close();
					}
				}
				if (cluster) {
					double clustering = originalGraph.clusteringCoefficientUndirected();
					System.out.println("Clustering coefficient: "+clustering);
					if (massbalance) {
						BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.clustering", false));
						writer.write(clustering+"\n");
						writer.close();
					}
					if (switchRand) {
						BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".switch.clustering", false));
						writer.write(clustering+"\n");
						writer.close();
					}
				}
				if (deltaGf) {
					// write the deltaGf of all compounds
					BufferedWriter writer = new BufferedWriter(new FileWriter(inputDir+File.separator+version+".deltaGf", false));
					Collection<Vertex> compoundsSet = originalGraph.getCompounds();
					// sort the compounds by name and compartment
					ArrayList<Vertex> compoundsList = new ArrayList<Vertex>(compoundsSet);
			    	Collections.sort(compoundsList, new VertexComparator());
					for (Vertex compound : compoundsList)
						writer.write(compound.getName()+(originalGraph.hasCompartments()?Utilities.DELIMITER_COMPOUND_COMPARTMENT+compound.getCompartment():"")+"\tD="+String.valueOf(compound.getDeltaGf())+"\n");
					writer.close();
				}
				if (deltaGn) {
					// write the deltaGn distributions to the original, massbalance and/or switch files
					for (int j=0; j<1+(massbalance && switchRand ? 1:0); j++) {
						BufferedWriter writer = new BufferedWriter(new FileWriter((j==0 && massbalance ? outputPrefix+".massbalance" : (j==1 || switchRand ? outputPrefix+".switch" : inputDir+File.separator+version))+".deltaGn", false));
						writer.write(originalGraph.getDeltaGn().toString()+"\n");
						writer.close();
					}
				}
				if (deltaGr) {
					// write the deltaGr distributions to the original, massbalance and/or switch files
					for (int j=0; j<1+(massbalance && switchRand ? 1:0); j++) {
						BufferedWriter writer = new BufferedWriter(new FileWriter((j==0 && massbalance ? outputPrefix+".massbalance" : (j==1 || switchRand ? outputPrefix+".switch" : inputDir+File.separator+version))+".deltaGr", false));
						for (Iterator<Vertex> reactionIt = originalGraph.getReactions().iterator(); reactionIt.hasNext(); ) {
							Double deltaG = originalGraph.getDeltaGr(reactionIt.next(), true);
							if (deltaG == null)
								writer.write("NaN");
							else
								writer.write(String.valueOf(deltaG));
							if (reactionIt.hasNext())
								writer.write("\t");
						}
						writer.write("\n");
						writer.close();
					}
				}
				if (assortativity) {
					double[] assort = originalGraph.assortativities();
					System.out.println("Assortativity coefficients: "+assort[0]+", "+assort[1]);
					// write the assortativities to the original, massbalance and/or switch files
					for (int j=0; j<1+(massbalance && switchRand ? 1:0); j++) {
						BufferedWriter writer = new BufferedWriter(new FileWriter((j==0 && massbalance ? outputPrefix+".massbalance" : (j==1 || switchRand ? outputPrefix+".switch" : inputDir+File.separator+version))+".assortativity", false));
						writer.write(assort[0]+"\t"+assort[1]+"\n");
						writer.close();
					}
				}
				if (transitionDegree) {
					// calculate transition degree
					int transDegree = originalGraph.getTransitionDegree(classes, preservedV);
					System.out.println("Transition degree: "+transDegree);
					BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.transitiondegree", false));
					writer.write(String.valueOf(transDegree)+"\n");
					writer.close();
				}
				if (matrixSubstitutions) {
					// write the stoichiometric matrix of the original network, exclude reversible reactions
					originalGraph.stoichiometricMatrix(inputDir+File.separator+version+".stmatrix", matrixSubstitutionsReversible, true, true, true);
					// write all possible matrix substitutions
					originalGraph.matrixSubstitutions(classes, preservedV, inputDir+File.separator+version+".matrixsubstitutions", true, matrixSubstitutionsReversible);
				}
				if (equationSubstitutions) {
					// write all possible equation substitutions
					originalGraph.equationSubstitutions(originalGraph, classes, preservedV, inputDir+File.separator+version+".equationsubstitutions", true);
				}
				if (localCentrality) {
					// calculate centralities for all forward and reversed reactions
					ArrayList<Double> centralities = new ArrayList<Double>(originalGraph.getReactions().size()*2);
					BufferedWriter writer = new BufferedWriter(new FileWriter(inputDir+File.separator+version+".centrality.local"+(localCentrality_reversible?".reversible":""), false));
					// sort by reactions
			    	ArrayList<Vertex> reactionsList = new ArrayList<Vertex>(originalGraph.getReactions());
			    	Collections.sort(reactionsList, new VertexComparator());
					
					for (int j=0; j<reactionsList.size(); j++) {
						double value = originalGraph.localEssentiality(reactionsList.get(j), localCentrality_reversible);
						centralities.add(value);
						writer.write(reactionsList.get(j).getName()+"\t"+value+"\n");
						if (!localCentrality_reversible && reactionsList.get(j).reversedReaction() != null) {
							value = originalGraph.localEssentiality(reactionsList.get(j).reversedReaction(), false);
							centralities.add(value);
							writer.write(reactionsList.get(j).getName()+"$rev\t"+value+"\n");
						}
					}
					writer.close();
					// write to output files
					if (massbalance) {
						BufferedWriter writer2 = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.centrality.local"+(localCentrality_reversible?".reversible":""), false));
						for (int j=0; j<centralities.size(); j++) {
							writer2.write(centralities.get(j).toString());
							if (j < centralities.size()-1)
								writer2.write("\t");
						}
						writer2.write("\n");
						writer2.close();
					}
					if (switchRand) {
						BufferedWriter writer2 = new BufferedWriter(new FileWriter(outputPrefix+".switch.centrality.local"+(localCentrality_reversible?".reversible":""), false));
						for (int j=0; j<centralities.size(); j++) {
							writer2.write(centralities.get(j).toString());
							if (j < centralities.size()-1)
								writer2.write("\t");
						}
						writer2.write("\n");
						writer2.close();
					}
				}
				if (centrality) {
					// write reaction names and centralities
					TreeMap<Vertex, Number> centralities = originalGraph.reactionCentralities(null, dampingFactor, centrality_double, centrality_reversible);
					BufferedWriter writer = new BufferedWriter(new FileWriter(inputDir+File.separator+version+".centrality"+(centrality_reversible?".reversible":"")+".d"+dampingFactor, false));
					for (Vertex reaction : centralities.keySet())
						writer.write(reaction.getName()+(reaction.isReversed()?"$rev":"")+"\t"+(centralities.get(reaction) == null ? "NA" : String.valueOf(centralities.get(reaction)))+"\n");
					writer.close();
					// write centralities of the forward reactions to the files for randomized distributions 
					if (massbalance) {
						BufferedWriter mbWriter = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.centrality"+(centrality_reversible?".reversible":"")+".d"+dampingFactor, false));
						for (Vertex reaction : centralities.keySet())
							mbWriter.write((centralities.get(reaction) == null ? "NA" : centralities.get(reaction).toString())+"\t");
						mbWriter.write("\n");
						mbWriter.close();
					}
					if (switchRand) {
						BufferedWriter swWriter = new BufferedWriter(new FileWriter(outputPrefix+".switch.centrality"+(centrality_reversible?".reversible":"")+".d"+dampingFactor, false));
						for (Vertex reaction : centralities.keySet())
							swWriter.write((centralities.get(reaction) == null ? "NA" : centralities.get(reaction).toString())+"\t");
						swWriter.write("\n");
						swWriter.close();
					}
				}
				if (knockoutSet) {
					Vertex knockoutReaction = originalGraph.getReaction(knockout);
					if (knockoutReaction == null)
						System.err.println("Reaction "+knockout+" not found in network.");
					else {
						HashMap<Vertex, Double[]> reactions = originalGraph.knockoutSet(knockoutReaction, knockoutThreshold);
						BufferedWriter writer = new BufferedWriter(new FileWriter(inputDir+File.separator+version+"."+knockout+".knockout", false));
						for (Vertex reaction : reactions.keySet()) {
							writer.write(reaction.getName()+(reaction.isReversed()?"$rev":"")+(originalGraph.hasCompartments()?Utilities.DELIMITER_COMPOUND_COMPARTMENT+reaction.getCompartment():"")+"\t"+reactions.get(reaction)[0]+"\t"+reactions.get(reaction)[1]+"\n");
						}
						writer.close();
					}
				}
				
				if (write)
					originalGraph.write(graphPrefix+".network", false);
				if (write2)
					originalGraph.write2(graphPrefix+".network2.", false);
				
			} else {
				// if only a subset is calculated any previous output files with the same
				// name have to be deleted. If from==0, this is done in the original run above.
				append = false;
			}
			
			int decimal = 1, fifth = 1;
			int run = 1, runs = (to-from+1) * ((massbalance && switchRand) ? 2 : 1);
			for (int index=from; index<=to && (massbalance || switchRand); index++, run++) {
				boolean first = (index == from);
				
				MetabolicGraph graph = null;
				String prefix = null, graphPrefixR = null;
				
				// allow to interrupt this thread
				wait(1);
				
				if (massbalance) {
					// read the massbalance randomized graph
					if (first)
						System.out.println("Mass-balanced randomized graph ["+index+"]:");
					ObjectInputStream reader = new ObjectInputStream(new FileInputStream(graphPrefix+".massbalance."+index));
					graph = (MetabolicGraph)reader.readObject();
					reader.close();
					matrixPrefix = inputDir+File.separator+version+"-matrices"+File.separator+version+".massbalance.stmatrix."+index;
					graphPrefixR = graphPrefix+".massbalance";
					prefix = outputPrefix+".massbalance";
				} else if (switchRand) {
					// read the switch randomized graph
					if (first)
						System.out.println("Switch randomized graph ["+index+"]:");
					ObjectInputStream reader = new ObjectInputStream(new FileInputStream(graphPrefix+".switch."+index));
					graph = (MetabolicGraph)reader.readObject();
					reader.close();
					matrixPrefix = inputDir+File.separator+version+"-matrices"+File.separator+version+".switch.stmatrix."+index;
					graphPrefixR = graphPrefix+".switch";
					prefix = outputPrefix+".switch";
				}

				if (first) {
					// print the connected components and reaction types
					graph.connectedComponents(true);
					Utilities.countReactions(graph, false);
				}
				
				// first get the set of preserved vertices
				if (preserve) {
					for (String name : preserved) {
						String compartment = null;
						if (name.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
							name = name.substring(0, name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
							compartment = name.substring(name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
						}
						Vertex v = graph.getCompound(name, compartment);
						if (v == null)
							v = graph.getReaction(name);
						preservedV.add(v);
					}
				}
				
				if (matrix)
					graph.stoichiometricMatrix(matrixPrefix, matrixReversible, matrixLabelled, matrixSorted, false);
				if (massbalance && degrees) {
					graph.compoundDegrees(prefix+".degrees", append);
					if (!graph.isReversible()) {
						graph.compoundInOutDegrees(prefix+".degrees.in", true, append);
						graph.compoundInOutDegrees(prefix+".degrees.out", false, append);
					}
				}
				if (massbalance && reactionDegrees)
					graph.reactionDegrees(prefix+".reactiondegrees", append);
				if (massbalance && weights)
					graph.weights(prefix+".weights", append);
				if (scopes)
					graph.randomScopeSizes(seedSizes, prefix+".scopes", null, append);
				
				// allow to interrupt this thread
				wait(1);
				
				if (pathLength) {
					float diam = graph.characteristicPathLength();
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".pathLength", append));
					writer.write(String.valueOf(diam)+"\n");
					writer.close();
				}
				
				// allow to interrupt this thread
				wait(1);
				
				if (cycles) {
					StrongConnectivityInspector<Vertex, DefaultWeightedEdge> connectivity = new StrongConnectivityInspector<Vertex, DefaultWeightedEdge>(graph);
					List<Set<Vertex>> components = connectivity.stronglyConnectedSets();
					for (Iterator<Integer> nCyclesIt = nCycles.iterator(); nCyclesIt.hasNext();) {
						int n = nCyclesIt.next();
						int cycleCount = graph.cycleCount(components, n, first, first && n==nCycles.get(0));
						BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+"."+n+"cycles", append));
						writer.write(String.valueOf(cycleCount)+"\n");
						writer.close();
					}
				}
				if (isPath) {
					boolean hasPath = graph.isPath(path, !directed);
					// write 1 if the path exist, 0 otherwise to the massbalance or switch file
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".path", append));
					writer.write((hasPath ? "1":"0")+"\n");
					writer.close();
				}
				if (connected) {
					// check if the network contains the specified path
					boolean[] isConnected = graph.isConnected(intermediaries, compartments, reversible);
					// write the connections to the massbalance or switch file
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".connected", append));
					for (int k=0; k<isConnected.length; k++) {
						if (k<isConnected.length-1)
							writer.write((isConnected[k]?"1":"0")+"\t");
						else
							writer.write((isConnected[k]?"1":"0")+"\n");
					}
					writer.close();
				}
				if (cluster) {
					double clustering = graph.clusteringCoefficientUndirected();
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".clustering", append));
					writer.write(clustering+"\n");
					writer.close();
				}
				if (deltaGn) {
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".deltaGn", append));
					writer.write(graph.getDeltaGn().toString()+"\n");
					writer.close();
				}
				if (deltaGr) {
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".deltaGr", append));
					for (Iterator<Vertex> reactionIt = graph.getReactions().iterator(); reactionIt.hasNext(); ) {
						Double deltaG = graph.getDeltaGr(reactionIt.next(), true);
						if (deltaG == null)
							writer.write("NaN");
						else
							writer.write(String.valueOf(deltaG));
						if (reactionIt.hasNext())
							writer.write("\t");
					}
					writer.write("\n");
					writer.close();
				}
				if (assortativity) {
					double[] assort = graph.assortativities();
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".assortativity", append));
					writer.write(assort[0]+"\t"+assort[1]+"\n");
					writer.close();
				}
				if (localCentrality) {
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".centrality.local"+(localCentrality_reversible?".reversible":""), append));
					ArrayList<Double> centralities = new ArrayList<Double>(graph.getReactions().size()*2);
					// sort by reactions
			    	ArrayList<Vertex> reactionsList = new ArrayList<Vertex>(graph.getReactions());
			    	Collections.sort(reactionsList, new VertexComparator());
					
					for (int j=0; j<reactionsList.size(); j++) {
						centralities.add(graph.localEssentiality(reactionsList.get(j), localCentrality_reversible));
						if (!localCentrality_reversible && reactionsList.get(j).reversedReaction() != null)
							centralities.add(graph.localEssentiality(reactionsList.get(j).reversedReaction(), false));
					}
					for (int j=0; j<centralities.size(); j++) {
						writer.write(centralities.get(j).toString());
						if (j < centralities.size()-1)
							writer.write("\t");
					}
					writer.write("\n");
					writer.close();
				}
				
				// allow to interrupt this thread
				wait(1);
				
				if (centrality) {
					TreeMap<Vertex, Number> centralities = graph.reactionCentralities(originalGraph, dampingFactor, centrality_double, centrality_reversible);
					BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+".centrality"+(centrality_reversible?".reversible":"")+".d"+dampingFactor, append));
					for (Iterator<Vertex> it = centralities.keySet().iterator(); it.hasNext(); ) {
						Vertex reaction = it.next();
						writer.write(centralities.get(reaction) == null ? "NA" : String.valueOf(centralities.get(reaction)));
						if (it.hasNext())
							writer.write("\t");
					}
					writer.write("\n");
					writer.close();
				}
				
				// allow to interrupt this thread
				wait(1);
				
				if (transitionDegree) {
					int transDegree = graph.getTransitionDegree(classes, preservedV);
					BufferedWriter writer = new BufferedWriter(new FileWriter(outputPrefix+".massbalance.transitiondegree", append));
					writer.write(String.valueOf(transDegree)+"\n");
					writer.close();
				}
				
				if (write)
					graph.write(graphPrefixR+".network."+index, false);
				if (write2)
					graph.write2(graphPrefixR+".network2."+index, false);

				// print ecoli reactions of the tca cycle
//				System.out.println("*********************");
//				Utilities.printTCACycleReactions(graph);
//				for (Vertex neighbour : Graphs.neighborListOf(graph, graph.getCompound("succ", null))) {
//					Utilities.printReaction(graph, neighbour);
//				}
				
				// progress output
				if (runs>= 10 && (run*10/decimal) >= runs) {
					System.out.print((10*decimal) + "%");
					decimal++;
					fifth++;
				} else if (runs>=50 && (run*50/fifth) >= runs) {
					System.out.print(".");
					fifth++;
				}
				
				// always append from the second run
				append = true;
				
				// after massbalance continue with switch randomization
				if (index == to) {
					if (massbalance && switchRand) {
						massbalance = false;
						index = from-1;
						append = (from == 0);
					}
				}
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			System.out.println("execution stopped.");
		}
		
		System.out.println("Total elapsed: " + ((float)(System.currentTimeMillis()-start)/60000f + " minutes."));
	}
	
	/**
	 * Creates a MetabolicGraph object from a tab-delimited network file. The format is freely configurable
	 * by defining the delimiters in the included configuration file.
	 * 
	 * @param reactionsFile Tab-delimited network file.
	 * @param logFile Output log-file.
	 * @param masses HashMap of compound names and their mass vectors.
	 * @param reversible If true, every reaction is considered reversible, irrespective of the annotation.
	 * @param compartments If true, compartment delimiters are used to create a compartment-specific network.
	 * @param fixBalance If true, hydrogen-unbalanced reactions are fixed either by phosphate of hydrogen atoms.
	 * @return The created MetabolicGraph.
	 * @throws IOException
	 */
	public static MetabolicGraph createGraph(String reactionsFile, String compoundsFile, String logFile, HashMap<String, int[]> masses, boolean reversible, boolean compartments, boolean fixBalance) throws IOException {
		boolean debug = false;
		BufferedReader reactionsReader = new BufferedReader(new FileReader(reactionsFile));
		BufferedReader deltaGReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(logFile, false));
		int[] counters = new int[7];
		int compounds = 0, validCompounds = 0, fixedReactions = 0, deltaGValues = 0;
		HashMap<Vertex, Integer> hydrogenImbalanced = new HashMap<Vertex, Integer>();
		HashMap<String, Double> deltaG = new HashMap<String, Double>();
		String line;
		String compartmentR = null;
		
		String version = Utilities.getVersion(reactionsFile);
		MetabolicGraph graph = new MetabolicGraph(version, masses.size(), reversible, compartments);
		infoWriter.write("# "+version+"\n");
		
		// extract deltaG values from the compounds file
		while ((line = deltaGReader.readLine()) != null) {

			if (line.startsWith("#"))
				continue;
			
			String compoundName = line.substring(0, line.indexOf("\t")).trim();
			ArrayList<String> tokens = Utilities.parseTokens(line.substring(line.indexOf("\t")+1), "\t", false);
			Double deltaGf;
			for (String token : tokens) {
				if (token.startsWith("D=")) {
					try {
						deltaGf = new Double(Double.parseDouble(token.substring(2)));
						deltaG.put(compoundName, deltaGf);
						break;
					} catch (NumberFormatException e) {
						// no double
					}
				}
			}
		}
		deltaGReader.close();
		
		while ((line = reactionsReader.readLine()) != null) {

			if (line.startsWith("#"))
				continue;
			
			// parse the reaction equation
			String equation = line.substring(line.indexOf("\t")+1, line.length()).trim();
			// parse the reaction compartment
			if (equation.charAt(0) == Utilities.DELIMITER_REACTION_COMPARTMENT_START && equation.indexOf(Utilities.DELIMITER_REACTION_COMPARTMENT_END) != -1) {
				compartmentR = (compartments ? equation.substring(1, equation.indexOf(Utilities.DELIMITER_REACTION_COMPARTMENT_END)) : null);
				equation = equation.substring(equation.indexOf(Utilities.DELIMITER_REACTION_COMPARTMENT_END)+1, equation.length());
			}
			
			boolean newToken = false, left = true;
			double coefficient = 1;
			
			// set reversibility if the graph allows irreversible reactions
			boolean reversibleReaction = reversible || equation.contains(Utilities.DELIMITER_EQUALS);
			
			// add reaction vertex
			counters[NUMREACTIONS]++;
			String reactionName = line.substring(0, line.indexOf("\t")).trim();
			Vertex reaction = new Vertex(reactionName, reversibleReaction);
			graph.addVertex(reaction);
			if (!reactionName.equals(reaction.getName()))
				infoWriter.write("Replaced duplicate reaction name "+reactionName+" by "+reaction.getName()+"\n");
			
			if (debug)
				System.out.println(reaction.getName()+"\t"+equation);
			
			for (int i=0, j=0; i<equation.length(); i++) {
				boolean lastToken = false;
				int increment = 0;
				String compartmentC = null;
				
				if (i == equation.length()-1)
					lastToken = true;
				
				if (equation.charAt(i) == Utilities.DELIMITER_COEFFICIENT_START) {
					// parse the stoichiometric coefficient
					coefficient = Utilities.parseCoefficient(equation.substring(i, equation.indexOf(Utilities.DELIMITER_COEFFICIENT_END, i+1)+1));
					if (coefficient <= 0)
						coefficient = 1;
					i = equation.indexOf(Utilities.DELIMITER_COEFFICIENT_END, i);
					j = i+1;
					
				} else if (lastToken && newToken || (increment = Utilities.isDelimiter(equation, i)) != -1) {
					// add the previous token
					String name = equation.substring(j, lastToken ? i+1 : i).trim();
					
					if (name.length() > 0) {
					
						if (name.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
							compartmentC = compartments ? name.substring(name.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1, name.length()) : null;
							name = name.substring(0, name.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
						}
						// use the reaction's compartment if the compound has no individual compartment 
						String compartment = compartmentC != null ? compartmentC : compartmentR;
						if (compartment == null)
							compartment = Utilities.DEFAULT_COMPARTMENT;
						Vertex compound = graph.getCompound(name, compartment);
						
						if (compound == null) {
							// add a new compound vertex
							compounds++;
							int[] mass = masses.get(name);
							compound = new Vertex(name, mass, compartment);
							if (deltaG.containsKey(name)) {
								compound.setDeltaGf(deltaG.get(name));
								deltaGValues++;
							}
							graph.addVertex(compound);
							if (mass != null)
								validCompounds++;
							else
								infoWriter.write("Compound from "+reaction.getName()+" not annotated in compounds file: "+name+". Adding null mass vector.\n");
						}
						// add the edge
						if (left) {
							DefaultWeightedEdge edge = graph.addEdge(compound, reaction);
							if (edge == null) {
								edge = graph.getEdge(compound, reaction);
								coefficient += graph.getEdgeWeight(edge);
							}
							graph.setEdgeWeight(edge, coefficient);
						} else {
							DefaultWeightedEdge edge = graph.addEdge(reaction, compound);
							if (edge == null) {
								edge = graph.getEdge(reaction, compound);
								coefficient += graph.getEdgeWeight(edge);
							}
							graph.setEdgeWeight(edge, coefficient);
						}
						if (debug)
							System.out.println(compound.getName()+" "+Utilities.DELIMITER_COEFFICIENT_START+(left ? graph.getEdgeWeight(graph.getEdge(compound, reaction)) : graph.getEdgeWeight(graph.getEdge(reaction, compound)))+Utilities.DELIMITER_COEFFICIENT_END);
					}
					
					// parse the equation side delimiter
					if (Utilities.isSideDelimiter(equation, i) != -1)
						left = false; 
					
					coefficient = 1;
					newToken = false;
					
					if (increment > 0)
						i += (increment-1);
					
				} else if (newToken == false) {
					// beginning of a new token
					newToken = true;
					j = i;
				}
			}
			// post-process the created reaction
			postProcessReaction(graph, reaction, hydrogenImbalanced, counters, fixBalance, infoWriter);
			
			if (debug) {
				if (reversibleReaction)
					System.out.println("Reversed reaction created.");
				System.out.println();
			}
		}
		
		// remove isolated compounds (left over from reaction removal)
		Vertex[] compound = new Vertex[graph.getCompounds().size()];
		graph.getCompounds().toArray(compound);
		for (Vertex c : compound) {
			if (graph.inDegreeOf(c) == 0 && graph.outDegreeOf(c) == 0) {
				graph.removeVertex(c);
				if (c.getMass() != null)
					validCompounds--;
				compounds--;
			}
		}
		
		// fix hydrogen unbalanced reactions
		for (Iterator<Vertex> it = hydrogenImbalanced.keySet().iterator(); it.hasNext();) {
			Vertex imbalanced = it.next();
			int added = graph.fixBalance(imbalanced, hydrogenImbalanced.get(imbalanced), true, true, infoWriter);
			if (added != -1)
				fixedReactions++;
			compounds += added;
			validCompounds += added;
		}
		
		// validation
		if (counters[NUMREACTIONS] != graph.getReactions().size())
			throw new RuntimeException("reactions: "+counters[NUMREACTIONS]+", getReactions().size(): "+graph.getReactions().size()+", "+counters[DUPLICATES]+" duplicates removed.");
		if (compounds != graph.getCompounds().size())
			throw new RuntimeException("compunds: "+compounds+", getCompounds().size(): "+graph.getCompounds().size());
		
		float unbalancedRatio = (counters[NUMREACTIONS]>0 ? (float)counters[UNBALANCED]*100/counters[NUMREACTIONS] : 100);
		float fixedRatio = (counters[UNBALANCED]>0 ? (float)fixedReactions*100/counters[UNBALANCED] : 100);
		float remainingUnbalancedRatio = (counters[NUMREACTIONS]>0?(float)(counters[UNBALANCED]-fixedReactions)*100/counters[NUMREACTIONS]: 100);
		float annotatedRatio = (compounds>0?(float)(validCompounds)*100/compounds: 100);
		
		String compartmentString = "";
		if (compartments) {
			Set<String> compartmentStrings = graph.getCompartments().keySet();
			for (Iterator<String> it = compartmentStrings.iterator(); it.hasNext();) {
				compartmentString += it.next();
				if (it.hasNext())
					compartmentString += " ";
			}
		}
		
		String message = compounds+" compound vertices"+(compartments?" in compartments":"")+", "+validCompounds+" annotated with mass ("+annotatedRatio+"%), "+deltaGValues+" annotated with deltaG, "+counters[NUMREACTIONS]+" reactions after removing "+counters[DUPLICATES]+" duplicate, "+counters[CYCLIC]+" cyclic and "+counters[ZERODEGREE]+" zero-degree reactions.\n"
			+(compartments ? "Compartments: "+compartmentString+". "+counters[TRANSPORTER]+" transport reactions.\n":"")
			+counters[UNBALANCED]+" unbalanced reactions ("+unbalancedRatio+"%), "+fixedReactions+" fixed ("+fixedRatio+"%), "+remainingUnbalancedRatio+"% remain unbalanced, "
			+counters[UNANNOTATED]+" unannotated ("+(counters[NUMREACTIONS]>0?counters[UNANNOTATED]*100/counters[NUMREACTIONS]:100)+"%).\n"
			+"Graph is "+(reversible? "" : "not ")+"reversible, has "+graph.vertexSet().size()+" vertices and "+graph.edgeSet().size()+" edges.\n";
		
		System.out.print(message);
		infoWriter.write(message);
		
		infoWriter.close();
		reactionsReader.close();
		
		return graph;
	}
	
	/**
	 * Parses a BioCyc flat file and creates a MetabolicGraph with weights according to the stoichiometric
	 * coefficients. If reversible==true, for every reaction a reversed reaction is added which converts
	 * the products into the substrates.
	 * 
	 * When a stoichiometric coefficient cannot be parsed as an integer, then it is substituted by 1.
	 * A warning is printed and invalidCoefficients is incremented.
	 * 
	 * @param reactionsFile
	 * @param reversible
	 * @throws IOException
	 */
	public static MetabolicGraph createBiocycGraph(String reactionsFile, String pathwaysFile, String logFile, HashMap<String, int[]> masses, boolean reversible, boolean compartments, boolean fixBalance) throws IOException {
		BufferedReader reactionsReader = new BufferedReader(new FileReader(reactionsFile));
		BufferedReader preParser = new BufferedReader(new FileReader(reactionsFile));
		BufferedReader pathwaysReader = null;
		if (pathwaysFile != null && new File(pathwaysFile).exists())
			pathwaysReader = new BufferedReader(new FileReader(pathwaysFile));
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(logFile, false));		
		String line;
		Vertex reaction = null;
		boolean left = false, right = false, smallMoleculeReaction = false;
		int[] counters = new int[7];
		int direction = 0, compounds = 0, validCompounds = 0, invalidCoefficients = 0, fixedReactions = 0;
		DefaultWeightedEdge newEdge = null;
		HashMap<Vertex, Integer> hydrogenImbalanced = new HashMap<Vertex, Integer>();
		HashMap<String, Integer> reversibilities = new HashMap<String, Integer>();
		HashMap<String, String> compartmentsMap = new HashMap<String, String>();
		
		String version = Utilities.getVersion(reactionsFile);
		MetabolicGraph graph = new MetabolicGraph(version, masses.size(), false, compartments);
		infoWriter.write("# "+version+"\n");
		
		// Read the global reaction directions from the pathways file. Reactions are considered reversible,
		// if both directions occur in any pathway. This information is only used, if there is no REACTION-DIRECTION
		// attribute in the reactions file, which is the case e.g. for AraCyc/TAIR up to version 6.0.
		if (!reversible && pathwaysReader != null) {
			while ((line = pathwaysReader.readLine()) != null) {
				
				if (line.startsWith("#"))
					continue;
				
				// parse the reaction direction: -1, 0, 1 for left-to-right, reversible, and right-to-left
				if (line.length() >= 19 && line.substring(0, 19).equals("REACTION-LAYOUT - (")) {
					
					int end = line.indexOf('(', line.indexOf('(')+1);
					String reactionName = line.substring(19, end).trim();
					end = line.indexOf(":DIRECTION")+10;
					
					if (end > 0) {
						// skip if the reaction is already reversible
						if (!reversibilities.containsKey(reactionName) || reversibilities.get(reactionName).intValue() != 0) {
							// store the direction: set to reversible if the opposite direction occurred before
							if (line.substring(end, end+5).trim().equals(":L2R")) {
								if (reversibilities.containsKey(reactionName) && reversibilities.get(reactionName).intValue() >= 0)
									reversibilities.put(reactionName, new Integer(0));
								else
									reversibilities.put(reactionName, new Integer(-1));
							} else if (line.substring(end, end+5).trim().equals(":R2L")) {
								if (reversibilities.containsKey(reactionName) && reversibilities.get(reactionName).intValue() <= 0)
									reversibilities.put(reactionName, new Integer(0));
								else
									reversibilities.put(reactionName, new Integer(1)); 
							}
						}
					}
				}
			}
		}
		
		// Pre-parse the network file to extract reversibilities and compartments.
		// Allows for sequential creation of the graph in the next loop.
		String rName = null, cName = null;
		while ((line = preParser.readLine()) != null) {
			
			if (line.startsWith("#"))
				continue;
			
			if (line.length() >= 12 && line.substring(0, 12).equals("UNIQUE-ID - ")) {
				rName = line.substring(12);
			} else if (!reversible && rName != null && line.length() >= 21 && line.substring(0, 21).equals("REACTION-DIRECTION - ")) {
				// overwrite the reversibilities from the pathways file
				if (line.substring(21, line.length()).trim().equals("REVERSIBLE"))
					reversibilities.put(rName, new Integer(0));
				else if (line.substring(21, line.length()).trim().equals("LEFT-TO-RIGHT"))
					reversibilities.put(rName, new Integer(-1));
				else if (line.substring(21, line.length()).trim().equals("RIGHT-TO-LEFT"))
					reversibilities.put(rName, new Integer(1));
			} else if (line.length() >= 7 && line.substring(0, 7).equals("LEFT - ")) {
				left = true;
				cName = line.substring(7);
				// remove the pipes				
				if (cName.charAt(0) == '|' && cName.charAt(cName.length()-1) == '|')
					cName = cName.substring(1, cName.length()-1);
				
			} else if (line.length() >= 8 && line.substring(0, 8).equals("RIGHT - ")) {
				left = false;
				cName = line.substring(8);
				// remove the pipes				
				if (cName.charAt(0) == '|' && cName.charAt(cName.length()-1) == '|')
					cName = cName.substring(1, cName.length()-1);
				
			} else if (compartments && line.length() >= 15 && line.substring(0, 15).equals("^COMPARTMENT - ")) {
				// add the compound's compartment
				compartmentsMap.put(left ? cName+"$"+rName : rName+"$"+cName, line.substring(15));
			} else if (line.length() >= 2 && line.substring(0, 2).equals("//")) {
				rName = null;
				cName = null;
			}
		}
		preParser.close();
		
		left = false;
		right = false;
		// parse the reactions and create the graph
		while ((line = reactionsReader.readLine()) != null) {
			
			if (line.startsWith("#"))
				continue;
			
			// add reaction as vertex
			if (line.length() >= 12 && line.substring(0, 12).equals("UNIQUE-ID - ")) {
				counters[NUMREACTIONS]++;
				String reactionName = line.substring(12);
				direction = reversibilities.containsKey(reactionName) ? reversibilities.get(reactionName).intValue() : 0;
				reaction = new Vertex(reactionName, reversible || direction == 0);
				graph.addVertex(reaction);
				if (!reactionName.equals(reaction.getName()))
					infoWriter.write("Replaced duplicate reaction name "+reactionName+" by "+reaction.getName()+"\n");
				smallMoleculeReaction = false;
			}
			
			if (line.length() == 32 && line.substring(0, 7).equals("TYPES - Small-Molecule-Reactions"))
				smallMoleculeReaction = true;
			else if (line.length() >= 7 && line.substring(0, 7).equals("LEFT - "))
				left = true;
			else if (line.length() >= 8 && line.substring(0, 8).equals("RIGHT - "))
				right = true;
			else if (line.length() >= 15 && line.substring(0, 15).equals("^COEFFICIENT - ")) {
				// set weight to stoichiometric coefficient
				int coefficient = 0;
				try {
					coefficient = Integer.parseInt(line.substring(15));
				} catch (NumberFormatException e) {
					// exception is handled in finally-clause
				} finally {
					if (coefficient < 1) {
						infoWriter.write("Invalid coefficient in reaction " + reaction.getName() + ", substituted by 1.\n");
						coefficient = 1;
						invalidCoefficients++;
					}
				}
				// increase the weight by the parsed coefficient-1
				graph.setEdgeWeight(newEdge, graph.getEdgeWeight(newEdge)+coefficient-1);
			} else if (line.length() >= 2 && line.substring(0, 2).equals("//")) {
				// post-process the created reaction
				postProcessReaction(graph, reaction, hydrogenImbalanced, counters, fixBalance, infoWriter);
			}
			
			if (left || right) {
				String name = (left? line.substring(7) : line.substring(8));
				String compartment = null;
				
				// remove the pipes				
				if (name.charAt(0) == '|' && name.charAt(name.length()-1) == '|')
					name = name.substring(1, name.length()-1);
				
				int[] mass = masses.get(name);
				if (compartments) {
					// here we rely on unique reaction names provided in reactions.dat, otherwise the key might be different from reaction.getName()
					compartment = compartmentsMap.get(left ? name+"$"+reaction.getName() : reaction.getName()+"$"+name);
					if (compartment == null)
						compartment = Utilities.DEFAULT_COMPARTMENT;
				}
				
				// add reactant as vertex if not already present
				Vertex compound = graph.getCompound(name, compartment);
				if (compound == null) {
					compounds++;
					compound = new Vertex(name, mass, compartment);
					graph.addVertex(compound);
					if (mass != null)
						validCompounds++;
					else {
						if (smallMoleculeReaction)
							infoWriter.write("Compound from Small-Molecule-Reaction "+reaction.getName()+" not annotated in compounds.dat: "+name+". Adding null mass vector.\n");							
						else
							infoWriter.write("Compound from "+reaction.getName()+" not annotated in compounds.dat: "+name+". Adding null mass vector.\n");
					}
				}
				
				if (direction == 1) {
					// switch the reaction direction
					left = !left;
					right = !right;
				}
				
				// create an edge in the respective direction (reversible edges are added implicitly)
				newEdge = (left ? graph.addEdge(compound, reaction) : graph.addEdge(reaction, compound));
				
				if (newEdge == null) {
					// if the edge already exists increment its coefficient					
					newEdge = left ? graph.getEdge(compound, reaction) : graph.getEdge(reaction, compound); 
					graph.setEdgeWeight(newEdge, graph.getEdgeWeight(newEdge)+1);
				}
				
				left = false;
				right = false;
			}
		}
		
		// remove isolated compounds (left over from reaction removal)
		Vertex[] compound = new Vertex[graph.getCompounds().size()];
		graph.getCompounds().toArray(compound);
		for (Vertex c : compound) {
			if (graph.inDegreeOf(c) == 0 && graph.outDegreeOf(c) == 0) {
				graph.removeVertex(c);
				if (c.getMass() != null)
					validCompounds--;
				compounds--;
			}
		}
		
		// fix hydrogen unbalanced reactions
		for (Iterator<Vertex> it = hydrogenImbalanced.keySet().iterator(); it.hasNext();) {
			Vertex imbalanced = it.next();
			int added = graph.fixBalance(imbalanced, hydrogenImbalanced.get(imbalanced), true, true, infoWriter);
			if (added != -1)
				fixedReactions++;
			compounds += added;
			validCompounds += added;
		}
		
		// validation
		if (counters[NUMREACTIONS] != graph.getReactions().size())
			throw new RuntimeException("reactions: "+counters[NUMREACTIONS]+", getReactions().size(): "+graph.getReactions().size()+", "+counters[DUPLICATES]+" duplicates removed.");
		if (compounds != graph.getCompounds().size())
			throw new RuntimeException("compunds: "+compounds+", getCompounds().size(): "+graph.getCompounds().size());
		
		// use reactions and compounds to test the code
//		if (compounds != graph.getCompounds().size())
//			throw new RuntimeException("Wrong number of compounds: "+compounds+" != "+graph.getCompounds().size());
//		if (reactions != graph.getReactions().size())
//			throw new RuntimeException("Wrong number of reactions: "+reactions+" != "+graph.getReactions().size());
		
		float unbalancedRatio = (counters[NUMREACTIONS]>0 ? (float)counters[UNBALANCED]*100/counters[NUMREACTIONS] : 100);
		float fixedRatio = (counters[UNBALANCED]>0 ? (float)fixedReactions*100/counters[UNBALANCED] : 100);
		float remainingUnbalancedRatio = (counters[NUMREACTIONS]>0?(float)(counters[UNBALANCED]-fixedReactions)*100/counters[NUMREACTIONS]: 100);
		float annotatedRatio = (compounds>0?(float)(validCompounds)*100/compounds: 100);
		
		counters[NUMREACTIONS] = graph.getReactions().size();
		compounds = graph.getCompounds().size();
		
		String compartmentString = "";
		if (compartments) {
			Set<String> compartmentStrings = graph.getCompartments().keySet();
			for (Iterator<String> it = compartmentStrings.iterator(); it.hasNext();) {
				compartmentString += it.next();
				if (it.hasNext())
					compartmentString += " ";
			}
		}
		
		String message = compounds+" compound vertices"+(compartments?" in compartments":"")+", "+validCompounds+" annotated with mass ("+annotatedRatio+"%), "+counters[NUMREACTIONS]+" reactions after removing "+counters[DUPLICATES]+" duplicate, "+counters[CYCLIC]+" cyclic and "+counters[ZERODEGREE]+" zero-degree reactions.\n"
			+(compartments ? "Compartments: "+compartmentString+". "+counters[TRANSPORTER]+" transport reactions.\n":"")
			+counters[UNBALANCED]+" unbalanced reactions ("+unbalancedRatio+"%), "+fixedReactions+" fixed ("+fixedRatio+"%), "+remainingUnbalancedRatio+"% remain unbalanced, "
			+counters[UNANNOTATED]+" unannotated ("+(counters[NUMREACTIONS]>0?counters[UNANNOTATED]*100/counters[NUMREACTIONS]:100)+"%).\n"
			+"Graph is "+(reversible? "" : "not ")+"reversible, has "+graph.vertexSet().size()+" vertices and "+graph.edgeSet().size()+" edges.\n";
		
		System.out.print(message);
		infoWriter.write(message);
		
		infoWriter.close();
		reactionsReader.close();
		
		return graph;
	}
	
	/**
	 * Parses an SBML file and creates a MetabolicGraph with weights according to the stoichiometric
	 * coefficients. if reversible==true, for every reaction a reversed reaction is added which converts
	 * the products into the substrates.
	 * 
	 * When a stoichiometric coefficient cannot be parsed as an integer, then it is substituted by 1.
	 * A warning is printed and invalidCoefficients is incremented.
	 * 
	 * @param sbmlFile The SBML file name.
	 * @param logFile Log file.
	 * @param masses HashMap of compound names and their mass vectors.
	 * @param reversible If true, every reaction is considered reversible, irrespective of the annotation.
	 * @param compartments If true, a compartment-specific network is created.
	 * @param fixBalance If true, hydrogen-unbalanced reactions are fixed either by phosphate of hydrogen atoms.
	 * @return The created MetabolicGraph.
	 * @throws IOException
	 */
	public static MetabolicGraph createSBMLGraph(String sbmlFile, String logFile, HashMap<String, int[]> masses, boolean reversible, boolean compartments, boolean fixBalance) throws IOException {
		boolean useTypes = false;
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(logFile, false));
		int[] counters = new int[7];
		int compounds = 0, validCompounds = 0, fixedReactions = 0;
		HashMap<Vertex, Integer> hydrogenImbalanced = new HashMap<Vertex, Integer>();
		String compartment = null;
		SBMLDocument sbmlDoc = Utilities.initSBML(sbmlFile);
		Model sbmlModel = sbmlDoc.getModel();
		String version = Utilities.getVersion(sbmlFile);
		
		// initialize the graph object
		MetabolicGraph graph = new MetabolicGraph(version, masses.size(), reversible, compartments);
		infoWriter.write("# "+version+"\n");
		
		// The SpeciesType object class is only available in SBML Level 2 Versions 2-4. It is not available in Level 1 nor Level 3.
		ListOfSpeciesTypes speciesTypes = null;
		if (sbmlModel.getLevel() == 2 && sbmlModel.getVersion() >= 2 && sbmlModel.getVersion() <= 4) {
			speciesTypes = sbmlModel.getListOfSpeciesTypes();
			// check if we have SpeciesTypes as in YeastNet (Herrgard et al., Nat Biotechnol 26(10), 1155-1160 Oct (2008))
			if (speciesTypes.size() > 0)
				useTypes = true;
		}
		
		// parse the network
		ListOfReactions reactions = sbmlModel.getListOfReactions();
		for (int i=0; i<reactions.size(); i++) {
			Reaction reactionElement = reactions.get(i);
			Vertex reaction = new Vertex(reactionElement.getId(), reversible || reactionElement.getReversible());
			graph.addVertex(reaction);
			counters[NUMREACTIONS]++;
			if (!reactionElement.getId().equals(reaction.getName()))
				infoWriter.write("Replaced duplicate reaction name "+reactionElement.getId()+" by "+reaction.getName()+"\n");
			
			ListOfSpecies speciesList = sbmlModel.getListOfSpecies();
			// the following call requires SBML Level 3
//			compartment = (compartments ? reactionElement.getCompartment() : null);
			
			// create the substrates and products
			long degree = reactionElement.getNumReactants()+reactionElement.getNumProducts();
			for (int j=0; j<degree; j++) {
				boolean left = j<reactionElement.getNumReactants();
				SpeciesReference speciesRef = (left ? reactionElement.getReactant(j) : reactionElement.getProduct(j-reactionElement.getNumReactants()));
				String id = speciesRef.getSpecies();
				Species species = speciesList.get(id);
				
				// skip substrates and products with boundaryCondition="true"
				if (species.getBoundaryCondition())
					continue;
				
				double coefficient = speciesRef.getStoichiometry();
				if (coefficient < 0)
					throw new RuntimeException("Invalid stoichiometric coefficient: "+coefficient+" for "+(left?"substrate":"product")+id+" in reaction "+reactionElement.getId()+" (must be larger than 0).");
				
				// YeastNet format: get the id from species type
				if (useTypes)
					id = speciesTypes.get(species.getSpeciesType()).getId();
				
				compartment = (compartments ? species.getCompartment() : null);
				
				Vertex compound = graph.getCompound(id, compartment);
				if (compound == null) {
					// add a new compound vertex
					compounds++;
					int[] mass = masses.get(id);
					compound = new Vertex(id, mass, compartment);
					graph.addVertex(compound);
					if (mass != null)
						validCompounds++;
					else
						infoWriter.write("Compound "+id+" from "+reaction.getName()+" has no mass. Adding null mass vector.\n");
				}
				
				// add the edge
				if (left) {
					DefaultWeightedEdge edge = graph.addEdge(compound, reaction);
					if (edge == null) {
						edge = graph.getEdge(compound, reaction);
						coefficient += graph.getEdgeWeight(edge);
					}
					graph.setEdgeWeight(edge, coefficient);
				} else {
					DefaultWeightedEdge edge = graph.addEdge(reaction, compound);
					if (edge == null) {
						edge = graph.getEdge(reaction, compound);
						coefficient += graph.getEdgeWeight(edge);
					}
					graph.setEdgeWeight(edge, coefficient);
				}
			}
			
			// post-process the created reaction
			postProcessReaction(graph, reaction, hydrogenImbalanced, counters, fixBalance, infoWriter);
		}
		
		// remove isolated compounds (left over from reaction removal)
		Vertex[] compound = new Vertex[graph.getCompounds().size()];
		graph.getCompounds().toArray(compound);
		for (Vertex c : compound) {
			if (graph.inDegreeOf(c) == 0 && graph.outDegreeOf(c) == 0) {
				graph.removeVertex(c);
				if (c.getMass() != null)
					validCompounds--;
				compounds--;
			}
		}
		
		// fix hydrogen unbalanced reactions
		for (Iterator<Vertex> it = hydrogenImbalanced.keySet().iterator(); it.hasNext();) {
			Vertex imbalanced = it.next();
			int added = graph.fixBalance(imbalanced, hydrogenImbalanced.get(imbalanced), true, true, infoWriter);
			if (added != -1)
				fixedReactions++;
			compounds += added;
			validCompounds += added;
		}
		
		// validation
		if (counters[NUMREACTIONS] != graph.getReactions().size())
			throw new RuntimeException("reactions: "+counters[NUMREACTIONS]+", getReactions().size(): "+graph.getReactions().size()+", "+counters[DUPLICATES]+" duplicates removed.");
		if (compounds != graph.getCompounds().size())
			throw new RuntimeException("compunds: "+compounds+", getCompounds().size(): "+graph.getCompounds().size());
		
		float unbalancedRatio = (counters[NUMREACTIONS]>0 ? (float)counters[UNBALANCED]*100/counters[NUMREACTIONS] : 100);
		float fixedRatio = (counters[UNBALANCED]>0 ? (float)fixedReactions*100/counters[UNBALANCED] : 100);
		float remainingUnbalancedRatio = (counters[NUMREACTIONS]>0?(float)(counters[UNBALANCED]-fixedReactions)*100/counters[NUMREACTIONS]: 100);
		float annotatedRatio = (compounds>0?(float)(validCompounds)*100/compounds: 100);
		
		String compartmentString = "";
		if (compartments) {
			Set<String> compartmentStrings = graph.getCompartments().keySet();
			for (Iterator<String> it = compartmentStrings.iterator(); it.hasNext();) {
				compartmentString += it.next();
				if (it.hasNext())
					compartmentString += " ";
			}
		}
		
		String message = compounds+" compound vertices"+(compartments?" in compartments":"")+", "+validCompounds+" annotated with mass ("+annotatedRatio+"%), "+counters[NUMREACTIONS]+" reactions after removing "+counters[DUPLICATES]+" duplicate, "+counters[CYCLIC]+" cyclic and "+counters[ZERODEGREE]+" zero-degree reactions.\n"
			+(compartments ? "Compartments: "+compartmentString+". "+counters[TRANSPORTER]+" transport reactions.\n":"")
			+counters[UNBALANCED]+" unbalanced reactions ("+unbalancedRatio+"%), "+fixedReactions+" fixed ("+fixedRatio+"%), "+remainingUnbalancedRatio+"% remain unbalanced, "
			+counters[UNANNOTATED]+" unannotated ("+(counters[NUMREACTIONS]>0?counters[UNANNOTATED]*100/counters[NUMREACTIONS]:100)+"%).\n"
			+"Graph is "+(reversible? "" : "not ")+"reversible, has "+graph.vertexSet().size()+" vertices and "+graph.edgeSet().size()+" edges.\n";
		
		System.out.print(message);
		infoWriter.write(message);
		
		infoWriter.close();
		
		return graph;
	}
	
	/**
	 * Private method for post-processing created reactions:
	 * 
	 * 1. remove (irreversible) duplicates of reactions
	 * 2. remove zero-degree reactions
	 * 3. remove cyclic reactions
	 * 4. add unbalanced reactions to HashMap<Vertex, Integer> hydrogenImbalanced
	 * 
	 * Returns the updated counters of reactions, transporters, duplicates, cyclic
	 * reactions, zero-degree, unannotated, and unbalanced reactions in the network.
	 * 
	 * @param graph
	 * @param reaction
	 * @param numReactions
	 * @param fixBalance
	 * @param infoWriter
	 * @return
	 * @throws IOException
	 */
	private static int[] postProcessReaction(MetabolicGraph graph, Vertex reaction, HashMap<Vertex, Integer> hydrogenImbalanced, int[] counters, boolean fixBalance, BufferedWriter infoWriter) throws IOException {
		
		// fix non-integer stoichiometry
		BigDecimal fix = graph.fixStoichiometry(reaction);
		if (fix.compareTo(new BigDecimal("1")) != 0)
			infoWriter.write("Fixed stoichiometry of reaction "+reaction.getName()+" by factor "+fix.floatValue()+".\n");
		
		// check if there is a duplicate
		Vertex duplicate = graph.isDuplicate(reaction);
			
		if (duplicate != null) {
			counters[DUPLICATES]++;
			counters[NUMREACTIONS]--;
			// keep the reversible reaction if one of the duplicates is irreversible
			if (reaction.reversedReaction() != null && duplicate.reversedReaction() == null) {
				graph.removeVertex(duplicate);
				hydrogenImbalanced.remove(duplicate);
				infoWriter.write("Duplicate reactions: "+reaction.getName()+" and "+duplicate.getName()+", removing irreversible "+duplicate.getName()+"\n");
			} else {
				graph.removeVertex(reaction);
				infoWriter.write("Duplicate reactions: "+reaction.getName()+" and "+duplicate.getName()+", removing "+((reaction.reversedReaction()==null && duplicate.reversedReaction()!=null)?"irreversible ":"")+reaction.getName()+"\n");
				return counters;
			}
		}
		
		// check zero-degree, cyclic and balance
		if (graph.inDegreeOf(reaction) == 0 && graph.outDegreeOf(reaction) == 0) {
			// remove reactions with no substrates and no products
			// note that biomass/export reactions may have no products
			counters[ZERODEGREE]++;
			counters[NUMREACTIONS]--;
			infoWriter.write("Removing zero-degree reaction: "+reaction.getName()+"\n");
			graph.removeVertex(reaction);
			
		} else if (graph.isCyclic(reaction, false)) {
			// if the created reaction is a cyclic reaction, remove it
			counters[CYCLIC]++;
			counters[NUMREACTIONS]--;
			infoWriter.write("Removing cyclic reaction: "+reaction.getName()+"\n");
			graph.removeVertex(reaction);
			
		// don't modify one-sided reactions (e.g. import/export reactions)
		} else if (graph.inDegreeOf(reaction) > 0 && graph.outDegreeOf(reaction) > 0) {
			// check for transporter
			if (graph.getCompartments(reaction).keySet().size() > 1) {
				counters[TRANSPORTER]++;
				String message = "Reaction "+reaction.getName()+" is a transporter through compartments ";
				for (String c : graph.getCompartments(reaction).keySet())
					message += c+" ";
				infoWriter.write(message+"\n");
			}
			
			// check the atom balance of the reaction, in particular hydrogen
			int[] imbalance = graph.balance(reaction);
			if (imbalance == null) {
				counters[UNANNOTATED]++;
				infoWriter.write("Reaction "+reaction.getName()+" contains unannotated compound(s).\n");
			} else {
				boolean balanced = true;
				int hydrogenImbalance = 0;
				for (int k=0; k<imbalance.length; k++) {
					if (Utilities.ELEMENTS[k].equals("H"))
						hydrogenImbalance = imbalance[k];
					else if (imbalance[k] != 0)
						balanced = false;
				}
				// if the reaction is hydrogen unbalanced, tag it for fixing
				if (!balanced || hydrogenImbalance != 0) {
					counters[UNBALANCED]++;
					String message = "";
					for (int k=0; k<imbalance.length; k++)
						if (imbalance[k] != 0)
							message += (Utilities.ELEMENTS[k]+"["+imbalance[k]+"] ");
					infoWriter.write("Reaction "+reaction.getName()+" is unbalanced: "+message+"\n");
				}
				if (fixBalance && balanced && hydrogenImbalance != 0)
					hydrogenImbalanced.put(reaction, hydrogenImbalance);
			}
		}
		
		return counters;
	}
	
	/**
	 * Returns the type of phosphate represented by the mass or -1, if the mass does not represent
	 * any phosphate form. The returned type corresponds to the number of hydrogen atoms:
	 * 0: phosphate, 1: hydrogen phosphate, 2: dihydrogen phosphate, 3: phosphoric acid.
	 * 
	 * @param mass
	 * @return The type of phosphate or -1, if the mass does not represent any phosphate form.
	 */
	public static int isPhosphate(int[] mass) {
		int type = -1;
		
		if (mass == null)
			return -1;
		
		if (mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: "+mass.length+".");
		
		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			if (Utilities.ELEMENTS[i].equals("H")) {
				type = mass[i];
			} else if (Utilities.ELEMENTS[i].equals("O")) {
				if (mass[i] != 4)
					return -1;
			} else if (Utilities.ELEMENTS[i].equals("P")) {
				if (mass[i] != 1)
					return -1;
			} else if (mass[i] != 0) {
				return -1;
			}
		}
		
		if (type < 0 || type > 3)
			return -1;
		else
			return type;
	}
	
	/**
	 * Determines whether the given mass vector represents
	 * hydrogen (H).
	 * 
	 * @param mass
	 * @return true, if the mass vector represents hydrogen, false otherwise.
	 */
	public static boolean isHydrogen(int[] mass) {
		boolean hydrogen = false;
		
		if (mass == null)
			return false;
		
		if (mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: "+mass.length+".");
		
		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			if (Utilities.ELEMENTS[i].equals("H")) {
				if (mass[i] == 1)
					hydrogen = true;
			} else if (mass[i] != 0) {
				return false;
			}
		}
		
		return hydrogen;
	}
	
	/**
	 * Compares two graph objects for equality. Equal means both graphs have reactions with 
	 * the same names, compounds with the same names and masses, and edges between equally
	 * named vertices.
	 * 
	 * @param graph1
	 * @param graph2
	 * @return True if the graphs are equal, false otherwise.
	 */
	public static boolean compare(MetabolicGraph graph1, MetabolicGraph graph2) {
		
		if (graph1.isReversible() != graph2.isReversible()) {
			System.out.println("Reversiblitites don't match.");
			return false;
		}
		
		if (!graph1.version.equals(graph2.version)) {
			System.out.println("Versions don't match.");
			return false;
		}
		
		if (graph1.vertexSet().size() != graph2.vertexSet().size()) {
			System.out.println("Different number of vertices.");
			return false;
		}

		if (graph1.edgeSet().size() != graph2.edgeSet().size()) {
			System.out.println("Different number of edges.");
			return false;
		}
		
		if (graph1.getCompounds().size() != graph2.getCompounds().size()) {
			System.out.println("Different number of compounds.");
			return false;
		}
		
		if (graph1.getReactions().size() != graph2.getReactions().size()) {
			System.out.println("Different number of reactions.");
			return false;
		}
		
		// compare vertices
		for (Iterator<Vertex> vertices = graph1.vertexSet().iterator(); vertices.hasNext();) {
			Vertex vertex = vertices.next();
			if (vertex.getType() == Vertex.REACTION) {
				if (vertex.isReversed()) {
					if (graph2.getReaction(vertex.getName()).reversedReaction() == null) {
						System.out.println("Reversed vertex not found in graph 2: "+vertex.getName());
						return false;
					}
				} else if (graph2.getReaction(vertex.getName()) == null) {
					System.out.println("Vertex not found in graph 2: "+vertex.getName());
					return false;
				}
			} else if (vertex.getType() == Vertex.COMPOUND) {
				if (graph2.getCompound(vertex.getName(), vertex.getCompartment()) == null || !Arrays.equals(vertex.getMass(), graph2.getCompound(vertex.getName(), vertex.getCompartment()).getMass())) {
					System.out.println("Vertex not found/different in graph 2: "+vertex.getName());
					return false;
				}
			}
		}
		
		// compare edges
		for (Iterator<DefaultWeightedEdge> edges = graph1.edgeSet().iterator(); edges.hasNext();) {
			DefaultWeightedEdge edge = edges.next();
			String sourceName = graph1.getEdgeSource(edge).getName();
			String sourceCompartment = graph1.getEdgeSource(edge).getCompartment();
			String targetName = graph1.getEdgeTarget(edge).getName();
			String targetCompartment = graph1.getEdgeTarget(edge).getCompartment();
			boolean substrate = graph1.getEdgeSource(edge).getType() == Vertex.COMPOUND;
			Vertex reaction1 = substrate ? graph1.getEdgeTarget(edge) : graph1.getEdgeSource(edge);
			Vertex reaction2;
			if (reaction1.isReversed()) {
				reaction2 = substrate ? graph2.getReaction(targetName).reversedReaction() : graph2.getReaction(sourceName).reversedReaction();
			} else {
				reaction2 = substrate ? graph2.getReaction(targetName) : graph2.getReaction(sourceName);
			}
			if (substrate && graph2.getEdge(graph2.getCompound(sourceName, sourceCompartment), reaction2) == null) {
				System.out.println("Substrate edge not found in graph 2: "+sourceName+" ---> "+reaction2.getName());
				return false;
			} else if (!substrate && graph2.getEdge(reaction2, graph2.getCompound(targetName, targetCompartment)) == null) {
				System.out.println("Product edge not found in graph 2: "+reaction2.getName()+" ---> "+targetName);
				return false;
			}
		}
		
		return true;
	}
	
	/**
	 * The following static block is needed in order to load the
	 * libSBML Java interface library when the application starts.
	 * If it fails, Utilities.SMBL is set to false. The error messages
	 * are only printed if a method requiring SBML is called.
	 */
	static {
		String varname;
		String shlibname;

		if (System.getProperty("os.name").startsWith("Mac OS")) {
			varname = "DYLD_LIBRARY_PATH"; // We're on a Mac.
			shlibname = "libsbmlj.jnilib and/or libsbml.dylib";
		} else {
			varname = "LD_LIBRARY_PATH"; // We're not on a Mac.
			shlibname = "libsbmlj.so and/or libsbml.so";
		}

		try {
			System.loadLibrary("sbmlj");
			// For extra safety, check that the jar file is in the classpath.
			Class.forName("org.sbml.libsbml.libsbml");
			Utilities.sbml = true;
			
		} catch (UnsatisfiedLinkError e) {
			Utilities.sbmlThrowable = e;
			Utilities.sbmlMessage = "Error encountered while attempting to load libSBML:\n"
				+ "Please check the value of your " + varname
				+ " environment variable and/or"
				+ " your 'java.library.path' system property"
				+ " (depending on which one you are using) to"
				+ " make sure it lists all the directories needed to"
				+ " find the " + shlibname + " library file and the"
				+ " libraries it depends upon (e.g., the XML parser).";
		} catch (ClassNotFoundException e) {
			Utilities.sbmlThrowable = e;
			Utilities.sbmlMessage = "Error: unable to load the file 'libsbmlj.jar'."
					+ " It is likely that your -classpath command line "
					+ " setting or your CLASSPATH environment variable "
					+ " do not include the file 'libsbmlj.jar'.";
		} catch (SecurityException e) {
			Utilities.sbmlThrowable = e;
			Utilities.sbmlMessage = "Error encountered while attempting to load libSBML:\n"
				+ "Could not load the libSBML library files due to a"
				+ " security exception.";
		}
	}
}
