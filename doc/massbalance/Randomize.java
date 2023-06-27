/**
 * 
 */
package massbalance;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Main class for randomizing a MetabolicGraph using the pre-calculated mass equivalence classes.
 * 
 * @author Georg Basler
 *
 */
public class Randomize extends Thread {
	
	/**
	 * Static initializer for loading the configuration parameters.
	 */
	static {
		Utilities.loadConfiguration();
	}
	
	static final String OPTIONS = "nofix\t\tDo not fix unbalanced reactions.\n"
		+"compartments\tGenerate a network with compartments.\n"
		+"reversible\tConsider every reaction reversible.\n"
		+"massbalance\tRun mass-balanced randomization.\n"
		+"noniterative\tReactions are chosen noniteratively for mass-balanced randomization.\n"
		+"depth=[0,...]\tMass-balanced randomization depth (0: no randomization, 1: default)\n"
		+"p=[0,1]\t\tProbability of mass-balanced randomizing a chosen reaction (1: default).\n"
		+"substitutability Print sample space and substitutability class size distributions.\n"
		+"switch\t\tRun switch randomization.\n"
		+"strict\t\tRun strict switch randomization.\n"
		+"D=vertices...\tPreserve the given reactions/compounds.\n"		
		+"from-to\t\tIndices of randomized networks to generate and parallelization. (default: 0-999).\n"
		+"continue\tContinue randomization on existing randomized object files.\n"		
		+"See Documentation for details.";
	static final String USAGE = "Usage: <networkFile> <outputDir> [options]\n"
		+"<networkFile>\tNetwork file as SBML, BioCyc flat file, or tab-separated list of reactions.\n"
		+"<outputDir>\tOutput directory for the generated networks and log-files.\n"
		+"Options:\n"+OPTIONS;

	// global HashMap for storing the compounds in equivalence classes
	static HashMap<ArrayList<Integer>, EquivalenceClass> classes;
	
	// randomizatin parameters, either parsed from the command line or passed by the GUI
	protected String networkFile = "", pathwayFile = "", outputDir = ""; 
	protected boolean reversible = false, nofix = false, compartments = false, continueRand = false, switchRand = false, massbalance = false, noniterative = false, strict = false, deltaG = false, substitutability = false, preserve = false;
	protected int from = 0, to = 999;
	protected float p = 1f, depth = 1f;
	protected HashSet<String> preserved;
	
	/**
	 * 
	 * Main method called from the command line.<br>
	 * <br>
	 * <strong>Usage:</strong><br>
	 * <code>networkFile outputDir [options]</code> (see below)<br>
	 * <strong>or</strong><br>
	 * <code>properties [args]</code> (calls {@link Properties#main(String[]) MetabolicGraphs.main(args)})<br>
	 * <strong>or</strong><br>
	 * <code>utilities [args])</code> (calls {@link Utilities#main(String[]) Utilities.main(args)})<br>
	 * <br>
	 * <table>
	 * <tr><td><code>networkFile</code</td><td>	Network file as SBML, BioCyc flat file, or tab-separated list of reactions.</td></tr>
	 * <tr><td><code>outputDir</code</td><td>	Output directory for the generated networks and log-files.</td></tr>
	 * </table>
	 * <br>
	 * <strong>Options:</strong>
	 * <table>
	 * <tr><td><code>nofix</code</td><td>		Do not fix unbalanced reactions.</td></tr>
	 * <tr><td>c<code>ompartments</code</td><td>	Generate a network with compartments.</td></tr>
	 * <tr><td><code>reversible</code</td><td>	Consider every reaction reversible.</td></tr>
	 * <tr><td><code>massbalance</code</td><td>	Run mass-balanced randomization.</td></tr>
	 * <tr><td><code>noniterative</code</td><td>	Reactions are chosen noniteratively for mass-balanced randomization.</td></tr>
	 * <tr><td><code>depth=[0,...]</code</td><td>	mass-balanced randomization depth (0: no randomization, 1: default).</td></tr>
	 * <tr><td><code>p=[0,1]</code</td><td>		Probability of mass-balanced randomizing a chosen reaction (1: default).</td></tr>
	 * <tr><td><code>substitutability</code</td><td> Print sample space and substitutability class size distributions.</td></tr>
	 * <tr><td><code>switch</code</td><td>		Run switch randomization.</td></tr>
	 * <tr><td><code>strict</code</td><td>		Run strict switch randomization.</td></tr>
	 * <tr><td><code>D=vertices...</code></td><td>	Preserve the given reactions/compounds.</td></tr> 
	 * <tr><td>f<code>rom-to</code</td><td>		Indices of randomized networks to generate and parallelization. (default: 0-999).</td></tr>
	 * <tr><td><code>continue</code</td><td>	Continue randomization on existing randomized object files.</td></tr>
	 * </table>
	 * <br>
	 */
	public static void main(String[] args) throws InterruptedException {
		// user-friendly way to call MetabolicGraphs.main by providing the 'properties' argument
		if (args.length > 0) {
			if (args[0].equals("properties")) {
				Properties.main(Arrays.copyOfRange(args, 1, args.length));
				System.exit(0);
			}
			if (args[0].equals("utilities")) {
				Utilities.main(Arrays.copyOfRange(args, 1, args.length));
				System.exit(0);
			}
		}
		
		Randomize randomize = new Randomize();
		randomize.parseArgs(args);
		
		// run randomization
		Thread t = new Thread(randomize);
		t.start();
		
	}
	
	public void parseArgs(String[] args) {
		boolean indexSpecified = false;
		
		// argument check
		if (args.length < 2) {
			System.out.println(USAGE);
			System.exit(-1);
		}
		
		// argument assignment
		int i=0;
		networkFile = args[i++];
		outputDir = args[i++];
		if (outputDir.endsWith(File.separator))
			outputDir = outputDir.substring(0, outputDir.length()-File.separator.length());
		
		if (!new File(networkFile).exists() || new File(networkFile).isDirectory()) {
			System.err.println("File not found: "+networkFile);
			System.out.println(USAGE);
			System.exit(-1);
		}
		
		// parse optional arguments
		while (args.length > i) {
			if (args[i].equals("continue"))
				continueRand = true;
			else if (args[i].equals("switch"))
				switchRand = true;
			else if (args[i].equals("massbalance"))
				massbalance = true;
			else if (args[i].equals("noniterative"))
				noniterative = true;
			else if (args[i].equals("strict"))
				strict = true;
			else if (args[i].equals("reversible"))
				reversible = true;
			else if (args[i].equals("nofix"))
				nofix = true;
			else if (args[i].equals("compartments"))
				compartments = true;
			else if (args[i].equals("substitutability"))
				substitutability = true;
			else if (!indexSpecified && args[i].contains("-") && Character.isDigit(args[i].charAt(0)) && Character.isDigit(args[i].charAt(args[i].length()-1))) {
				// parse the start and end indices to use
				from = Integer.parseInt(args[i].substring(0, args[i].indexOf("-")));
				to = Integer.parseInt(args[i].substring(args[i].indexOf("-")+1, args[i].length()));
				indexSpecified = true;
			} else if (args[i].startsWith("p=")) {
				// parse the perturbation probability p
				p = Float.parseFloat(args[i].substring(2, args[i].length()));
			} else if (args[i].startsWith("depth=")) {
				// parse the randomization depth
				depth = Float.parseFloat(args[i].substring(6, args[i].length()));
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
			} else {
				System.err.println("Unrecognized argument: "+args[i]+(indexSpecified && args[i].contains("-")?" (only one index argument allowed)":""));
				System.out.println(USAGE);
				System.exit(-1);
			}
			i++;
		}
	}
	
	/**
	 * Runs the randomization.
	 * 
	 * The parameters are either passed by the GUI ({@link JMassBalance}),
	 * or parsed from the command line ({@link Randomize#main(String[])}).
	 */
	@Override
	public synchronized void run() {
		long start = System.currentTimeMillis(), time = System.currentTimeMillis(), switchTime = 0, massBalanceTime = 0;
		String graphLog = "", randomizeLog = "", version = "", substitutabilityFile;
		HashSet<Vertex> preservedV = new HashSet<Vertex>();
		
		try {
			
			// prepare the file names from the input directory
			version = Utilities.getVersion(networkFile);
			graphLog = outputDir+File.separator+version+".graph.log";
			new File(graphLog).delete();
			randomizeLog = outputDir+File.separator+version+".randomize.log";
			new File(randomizeLog).delete();
			
			HashMap<String, int[]> massesMap;
			substitutabilityFile = outputDir+File.separator+version+".substitutability";
			
			// allow to interrupt this thread
			wait(1);
			
			String classesFile = outputDir+File.separator+version+".classes";
			if (new File(classesFile).exists()) {
				// read equivalence classes from file and create the mass vectors
				System.out.print("Reading classes from "+classesFile+". ");
				time = System.currentTimeMillis();
				classes = Equivalence.readClasses(classesFile);
				System.out.println(" ("+(System.currentTimeMillis()-time)/1000 + " seconds).");
				
				massesMap = Equivalence.getMasses(classes);
			
			} else {
				if (continueRand) {
					System.err.println("Missing classes file: "+classesFile);
					System.exit(-1);
				}
				// create new mass equivalence classes
				String[] arguments = {networkFile, outputDir};
				Equivalence.main(arguments);
				classes = Equivalence.classes;
				massesMap = Equivalence.massesMap;
			}
			// check classes for consistency
//			Equivalence.consistency(classes, massesMap.size(), true);
			
			String graphDir = outputDir+File.separator+version+"-graphs";
			new File(graphDir).mkdir();
			String graphPrefix = graphDir+File.separator+version;
			String newGraphPrefix = null;
			int newGraphIndex = 1;
			
			// create the graph from the original network
			MetabolicGraph graph = null;
			if (continueRand) {
				// read the original graph
				ObjectInputStream graphReader = new ObjectInputStream(new FileInputStream(graphPrefix));
				graph = (MetabolicGraph)graphReader.readObject();
				graphReader.close();
				// create a new graph directory
				while (new File(graphDir+"-"+newGraphIndex).exists())
					newGraphIndex++;
				newGraphPrefix = graphDir+"-"+newGraphIndex+File.separator+version;
				new File(graphDir+"-"+newGraphIndex).mkdir();
				// check if we have deltaG values
				for (Vertex compound : graph.getCompounds()) {
					if (compound.getDeltaGf() != null) {
						deltaG = true;
						break;
					}
				}
			} else {
				// create the original graph
				String deltaGFile;
				
				if (networkFile.endsWith(".xml")) {
					graph = Properties.createSBMLGraph(networkFile, graphLog, massesMap, reversible, compartments, !nofix);
					deltaGFile = networkFile.replace(".xml", ".deltaG");
				} else if (networkFile.endsWith(".dat")) {
					pathwayFile = Utilities.getPathwaysFile(networkFile);
					graph = Properties.createBiocycGraph(networkFile, pathwayFile, graphLog, massesMap, reversible, compartments, !nofix);
					deltaGFile = networkFile.replace(".dat", ".deltaG");
				} else {
					graph = Properties.createGraph(networkFile, Utilities.getCompoundsFile(networkFile), graphLog, massesMap, reversible, compartments, !nofix);
					deltaGFile = networkFile+".deltaG";
				}
				// read deltaG if there is a file
				try {
					graph.readDeltaG(deltaGFile, graphLog);
					deltaG = true;
				} catch (FileNotFoundException e) {
					// no deltaG file found
				}
			}
			
			// get the set of preserved vertices
			if (preserve) {
				for (String name : preserved) {
					String vertex = name, compartment = null;
					if (name.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
						vertex = name.substring(0, name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
						compartment = name.substring(name.lastIndexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1);
					}
					Vertex v = graph.getCompound(vertex, compartment);
					if (v == null)
						v = graph.getReaction(vertex);
					preservedV.add(v);
				}
			}
			
//			graph.write2(graphPrefix+".network", false);

			// remove to increase performance
			System.out.println("Characteristic path length: "+graph.characteristicPathLength()+", clustering coefficient: "+graph.clusteringCoefficientUndirected());
			
			if (!continueRand && from == 0) {
				ObjectOutputStream writer = new ObjectOutputStream(new FileOutputStream(graphPrefix));
				writer.writeObject(graph);
				writer.close();
				if (substitutability) {
					// print sample space and substitutability classes
					graph.substitutability(classes, substitutabilityFile);
				}
//				int deadEnds = graph.compoundDeadEnds(outputDir+File.separator+version+".deadends");
//				System.out.println("Found "+deadEnds+" dead-end compounds.");
			}
			
			// randomized graphs
//			BigDecimal modifiedRatio = new BigDecimal("0");
			int decimal = 1, fifth = 1;
			int run = 1, runs = to-from+1;
			boolean progress = false;
			for (int index=from; (massbalance || switchRand) && index<=to; index++, run++) {
				
				// allow to interrupt this thread
				wait(1);
				
				if (switchRand) {
					MetabolicGraph switchGraph = null;
					// switch randomization
					if (continueRand) {
						ObjectInputStream reader = new ObjectInputStream(new FileInputStream(graphPrefix+".switch."+index));
						switchGraph = (MetabolicGraph)reader.readObject();
						reader.close();
					} else {
						switchGraph = graph.clone();
						if (index == from) {
							if (!switchGraph.compare(graph))
								throw new RuntimeException("Unidentical graph clone.");
						}
					}
					time = System.currentTimeMillis();
					switchGraph.switchRandomize(strict);
					switchTime += (System.currentTimeMillis()-time);
					// serialize graph
					ObjectOutputStream graphWriter = new ObjectOutputStream(new FileOutputStream((continueRand?newGraphPrefix:graphPrefix)+".switch."+index));
					graphWriter.writeObject(switchGraph);
					if (from == 0)
						switchGraph.writeBalance(outputDir+File.separator+version+".switch.balance", index!=0);
				}
				
				// allow to interrupt this thread
				wait(1);
				
				if (massbalance) {
					// mass-balanced randomization
					MetabolicGraph mbGraph = null;
					if (continueRand) {
						ObjectInputStream reader = new ObjectInputStream(new FileInputStream(graphPrefix+".massbalance."+index));
						mbGraph = (MetabolicGraph)reader.readObject();
						reader.close();
					} else {
						mbGraph = graph.clone();
						if (index == from) {
							if (!mbGraph.compare(graph))
								throw new RuntimeException("Unidentical graph clone.");
						}
					}
					
					time = System.currentTimeMillis();
					// to print the substitutability differences include: substitutabilityFile+".differences"
					mbGraph.randomize(graph, classes, preservedV, p, !noniterative, depth, deltaG?outputDir+File.separator+version+".massbalance."+from+"-"+to+(continueRand?"["+newGraphIndex+"]":"")+".deltaG":null, null, index!=from, (index!=from||continueRand?null:randomizeLog));
//					if (modified != null)
//						modifiedRatio = modifiedRatio.add(modified);
					massBalanceTime += (System.currentTimeMillis()-time);
					// serialize graph
					ObjectOutputStream graphWriter = new ObjectOutputStream(new FileOutputStream((continueRand?newGraphPrefix:graphPrefix)+".massbalance."+index));
					graphWriter.writeObject(mbGraph);
//					mbGraph.write2(graphPrefix+".massbalance.network."+index, false);
				}
				
				// progress output
				if (runs>= 10 && (run*10/decimal) >= runs) {
					System.out.print((10*decimal) + "%");
					decimal++;
					fifth++;
					progress = true;
				} else if (runs>=50 && (run*50/fifth) >= runs) {
					System.out.print(".");
					fifth++;
					progress = true;
				}
			}
			
			if (progress)
				System.out.println();
			if (massbalance)
				System.out.println("Total mass-balanced randomization time: "+(float)massBalanceTime/60000f + " minutes.");
			if (switchRand)
				System.out.println("Total switch randomization time: "+(float)switchTime/60000f + " minutes.");
			
//			if (massbalance) {
//				BigDecimal averageModified = modifiedRatio.divide(new BigDecimal(String.valueOf(runs)), 6, RoundingMode.HALF_UP).multiply(new BigDecimal("100"));
//				System.out.println("Average of "+averageModified.floatValue()+"% reactions modified.");
//			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (OutOfMemoryError e) {
			e.printStackTrace();
			if (e.getMessage().contains("heap space"))
				System.err.println("Try increasing the JVM heap space by adding the java argument '-server' or '-Xmx1024m'.");
		} catch (InterruptedException e) {
			System.out.println("execution stopped.");
		}
		
		System.out.println("Total elapsed: " + ((float)(System.currentTimeMillis()-start)/60000f + " minutes."));
		System.out.println("See "+(!continueRand?graphLog+" and ":"")+randomizeLog+" for more information.");
		if (!massbalance && !switchRand)
			System.out.println("No randomized networks were generated. Add 'massbalance' argument to generate 1000 randomized networks with default parameters.");
		else
			System.out.println("Files were written to "+outputDir+File.separator+". Use "+version+" as network name for calculating graph properties.");
	}
}
