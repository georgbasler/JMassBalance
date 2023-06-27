/**
 * 
 */
package massbalance;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import org.jgrapht.Graphs;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.sbml.libsbml.ListOfReactions;
import org.sbml.libsbml.ListOfSpecies;
import org.sbml.libsbml.ListOfSpeciesReferences;
import org.sbml.libsbml.ListOfSpeciesTypes;
import org.sbml.libsbml.Model;
import org.sbml.libsbml.Reaction;
import org.sbml.libsbml.SBMLDocument;
import org.sbml.libsbml.SBMLError;
import org.sbml.libsbml.SBMLReader;
import org.sbml.libsbml.SBase;
import org.sbml.libsbml.Species;
import org.sbml.libsbml.SpeciesType;
import org.sbml.libsbml.XMLNode;

/**
 * <p>
 * This class contains static methods for parsing compound masses, writing graph information,
 * and merging property files.
 * 
 * <p>
 * The public methods can be called directly from the command line by executing
 * <code>java -jar massbalance.jar utilities [method] [args...]</code>,
 * where <code>method</code> is the method name and <code>args...</code>
 * are the method arguments.
 * 
 * <p>
 * The following method arguments are supported
 * and will be translated from the given string argument:
 * <p>
 * <code>MetabolicGraph</code>: pass the file name of the deserialized object.<br>
 * <code>int</code>: pass the integer argument.<br>
 * <code>String</code>: pass the String argument.<br>
 * <code>boolean</code>: pass 'true' or 'false'.<br>
 * <code>ArrayList&lt;String&gt;</code>: pass a list of strings as ['String1','String2',...].<br>
 *
 * <p>
 * An overview of the available methods is printed by executing the main class without arguments.
 * 
 * @author Georg Basler
 *
 */
public class Utilities {
	
	static final boolean DEBUG = false;
	static final String PARAMETERS = "Usage: <method> [arguments] ...";
	static final String CONFIG_FILE = "jmassbalance.config";
	static boolean configLoaded = false;
	
	public static boolean sbml = false;
	public static String sbmlMessage = "";
	public static Throwable sbmlThrowable;
	public static SBMLDocument sbmlDocument;
	
	/**
	 * Configuration variables, may be set from the config file.
	 */
	/** Chemical elements to be considered for the mass vectors. */
	static String[] ELEMENTS = new String[]{"C", "H", "N", "O", "P", "S"};;
	/** Hydrogen name. */
	static String HYDROGEN_NAME = "$Hydrogen";
	/** Mass vector of hydrogen. */
	static int[] HYDROGEN_MASS = new int[]{0, 1, 0, 0, 0, 0};
	/** Mass vectors of the four phosphate forms. */
	static int[][] PHOSPHATE_FORMS = new int[][]{{0, 0, 0, 4, 1, 0}, {0, 1, 0, 4, 1, 0}, {0, 2, 0, 4, 1, 0}, {0, 3, 0, 4, 1, 0}};
	/** Names of the phosphate forms. */
	static String[] PHOSPHATE_NAMES = new String[]{"$Phosphate", "$Hydrogen-Phosphate", "$Dihydrogen-Phosphate", "$Phosphoric-Acid"};
	/** Substrate/product-delimiter for reaction equations. */
	/** Side-delimiter for reversible reactions. */
	static String DELIMITER_PLUS = "$+$";
	static String DELIMITER_EQUALS = "$=$";
	/** Side-delimiter for irreversible (forward) reactions. */
	static String DELIMITER_FORWARD = "$>$";
	/** Begin-delimiter for stoichiometric coefficients. */
	static char DELIMITER_COEFFICIENT_START = '{';
	/** End-delimiter for stoichiometric coefficients. */
	static char DELIMITER_COEFFICIENT_END = '}';
	/** Delimiter for compound compartments. */
	static String DELIMITER_COMPOUND_COMPARTMENT = "$";
	/** Default compartment. */
	static String DEFAULT_COMPARTMENT = "Default_Compartment";
	/** Begin-delimiter for reaction compartments. */
	static char DELIMITER_REACTION_COMPARTMENT_START = '[';
	/** End-delimiter for reaction compartments. */
	static char DELIMITER_REACTION_COMPARTMENT_END = ']';
	/** Estimated ratio of the number of equivalence classes divided by the number of compounds. */
	static int RATIO_CLASSES_COMPOUNDS = 245; // 220: athaliana
	/** Estimated ration of the number of compounds divided by the number of reactions. */
	static float RATIO_COMPOUNDS_REACTIONS = 1;
	/** Estimated ratio of the number of equivalence classes divided by the number of reactions. */
	static float RATIO_CLASSES_REACTIONS = RATIO_CLASSES_COMPOUNDS * RATIO_COMPOUNDS_REACTIONS;
	/** Maximum factor to multiply a reaction stoichiometry **/
	static int MAX_STOICHIOMETRY_FACTOR = 7;
	/** Default number of compounds. */
	static int DEFAULT_NUMBER_OF_COMPOUNDS = 2700;
	/** Default number of reactions. */
	static int DEFAULT_NUMBER_OF_REACTIONS = 2900;
	/** Default number of classes. */
	static int DEFAULT_NUMBER_OF_CLASSES = 980000; // 410000;

	/**
	 * Static initializer for loading the configuration parameters.
	 */
	static {
		loadConfiguration();
	}
	
	/**
	 * <p>
 	 * The main method can be used to call the public methods of this class
 	 * directly from the command line by executing<br>
 	 * <code>java -jar massbalance.jar utilities [method] [args...]</code><br>
 	 * where <code>method</code> is the method name and <code>args...</code>
 	 * are the method arguments. The following method arguments are supported
 	 * and will be translated from the given string argument:
 	 * <p>
 	 * <code>MetabolicGraph</code>: pass the file name of the deserialized object.<br>
 	 * <code>int</code>: pass the integer argument.<br>
 	 * <code>String</code>: pass the String argument.<br>
 	 * <code>boolean</code>: pass 'true' or 'false'.<br>
	 * <code>ArrayList&lt;String&gt;</code>: pass a list of strings as ['String1','String2',...].<br>
	 * 
	 * <p>
	 * An overview of the available methods is printed by executing the main class without arguments.
	 */
	public static void main(String[] args) {
		ObjectInputStream graphReader = null;
		
		if (args.length < 1) {
			exit();
		}
		
		try {
			
			// parse the method string and arguments
			ArrayList<Class<?>> types = new ArrayList<Class<?>>(args.length);
			ArrayList<Object> params = new ArrayList<Object>(args.length);
			for (int i=0; i<args.length-1; i++) {
				try {
					// try to deserialize a graph
					graphReader = new ObjectInputStream(new FileInputStream(args[i+1]));
					MetabolicGraph graph = (MetabolicGraph)graphReader.readObject();
					graphReader.close();
					params.add(graph);
					types.add(MetabolicGraph.class);
				} catch (IOException ioe) {
						// not a graph
						if (graphReader != null)
							graphReader.close();
					try {
						// try a from-to int range
						if (args[i+1].contains("-") && Character.isDigit(args[i+1].charAt(0)) && Character.isDigit(args[i+1].charAt(args[i+1].length()-1))) {
							// try start and end indices
							params.add(Integer.parseInt(args[i+1].substring(0, args[i+1].indexOf("-"))));
							types.add(int.class);
							params.add(Integer.parseInt(args[i+1].substring(args[i+1].indexOf("-")+1, args[i+1].length())));
							types.add(int.class);
						} else {
							// try an integer
							try {
								params.add(Integer.parseInt(args[i+1]));
								types.add(int.class);
							} catch (NumberFormatException e) {
								// try a double
								params.add(Double.parseDouble(args[i+1]));
								types.add(double.class);
							}
						}
					} catch (NumberFormatException e) {
						// try a boolean
						if (args[i+1].equals("true") || args[i+1].equals("false")) {
							params.add(args[i+1].equals("true"));
							types.add(boolean.class);
						} else if (args[i+1].startsWith("[") && args[i+1].endsWith("]")) {
							// try a string list
							ArrayList<String> strings = new ArrayList<String>();
							int start = 1, end = 0;
							do {
								end = args[i+1].indexOf(",", start+1);
								strings.add(args[i+1].substring(start, (end != -1)?end:args[i+1].lastIndexOf("]")));
								start = end+1;
							} while (end != -1);
							params.add(strings);
							types.add(ArrayList.class);
						} else {
							// take argument as string
							params.add(args[i+1]);
							types.add(String.class);
						}
					}
				}
			}
			
			// call the method given as first argument with the remaining argument parameters
			Class<?> utilities = Class.forName("massbalance.Utilities");
			Class<?>[] classTypes = new Class[types.size()];
			types.toArray(classTypes);
			Method method = utilities.getDeclaredMethod(args[0], classTypes);
			if (!Modifier.isPublic(method.getModifiers()))
				exit();
			method.invoke(null, params.toArray());
			
			// draw a random number between 1 and args[0]
//			int max = Integer.parseInt(args[0]);
//			Random random = new Random();
//			System.out.println("Your throw [1,"+max+"]: "+(random.nextInt(max)+1));
			
//			printExampleReactions(args);
//			ObjectInputStream graphReader = new ObjectInputStream(new FileInputStream(args[0]+File.separator+"graphs"+File.separator+args[1]));
//			MetabolicGraph graph = (MetabolicGraph)graphReader.readObject();
//			graphReader.close();
//			countReactions(graph, true);
			
//			boolean biocyc = true;
//			checkNetwork(args, biocyc);
			
			// count unique strings
//			BufferedReader reader = new BufferedReader(new FileReader(args[0]));
//			String line;
//			HashSet<String> names = new HashSet<String>(); 
//			while ((line = reader.readLine()) != null)
//				names.add(line.trim());
//			System.out.println(names.size());
			
			// compound-masses mapping
//			compoundMassesMap(args[0], args[1], args[2]);
			
			// test graph files
//			testGraphs(args[0], args[1], !new String("switch").equals(args[2]));
			
//			maxSingleClass(args[0]);
//			writeCompartments(args[0]);
//			testSolveStoichiometrySingle();
//			testSolveStoichiometryDouble();
//			writeBioCycWeights(args[0], args[1]);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Sets the default configuration parameters, then tries to load the parameter values
	 * from the configuration file. The method is executed only once during initialization.
	 */
	protected static void loadConfiguration() {
		boolean debug = false;
		
		// try to set the configuration once per initialization
		if (configLoaded)
			return;
		else
			configLoaded = true;
		
		try {
			
			BufferedReader configReader = new BufferedReader(new FileReader(CONFIG_FILE));
			String line;
			while ((line = configReader.readLine()) != null) {
				
				if (line.startsWith("#") || line.trim().length() == 0)
					continue;
				
				String name = line.substring(0, line.indexOf("=")).trim();
				String value = line.substring(line.indexOf("=")+1).trim();
				
				try {
					// int types
					if (name.equals("RATIO_CLASSES_COMPOUNDS")) {
						RATIO_CLASSES_COMPOUNDS = Integer.parseInt(value);
					} else if (name.equals("RATIO_COMPOUNDS_REACTIONS")) {
						RATIO_COMPOUNDS_REACTIONS = Float.parseFloat(value);
					} else if (name.equals("MAX_STOICHIOMETRY_FACTOR")) {
						MAX_STOICHIOMETRY_FACTOR = Integer.parseInt(value);
					} else if (name.equals("DEFAULT_NUMBER_OF_COMPOUNDS")) {
						DEFAULT_NUMBER_OF_COMPOUNDS = Integer.parseInt(value);
					} else if (name.equals("DEFAULT_NUMBER_OF_REACTIONS")) {
						DEFAULT_NUMBER_OF_REACTIONS = Integer.parseInt(value);
					} else if (name.equals("DEFAULT_NUMBER_OF_CLASSES")) {
						DEFAULT_NUMBER_OF_CLASSES = Integer.parseInt(value);
					
					// char types
					} else if (name.equals("DELIMITER_REACTION_COMPARTMENT")) {
						if (value.length() == 2 && value.charAt(0) != value.charAt(1)) {
							DELIMITER_REACTION_COMPARTMENT_START = value.charAt(0);
							DELIMITER_REACTION_COMPARTMENT_END = value.charAt(1);
						} else {
							System.out.println("DELIMITER_REACTION_COMPARTMENT must consist of one opening and another closing character.");
						}
							
					} else if (name.equals("DELIMITER_COEFFICIENT")) {
						if (value.length() == 2 && value.charAt(0) != value.charAt(1)) {
							DELIMITER_COEFFICIENT_START = value.charAt(0);
							DELIMITER_COEFFICIENT_END = value.charAt(1);
						} else {
							System.out.println("DELIMITER_COEFFICIENT must consist of one opening and another closing character.");
						}
						
					// String types
					} else if (name.equals("HYDROGEN_NAME")) {
						HYDROGEN_NAME = value;
					} else if (name.equals("DELIMITER_PLUS")) {
						DELIMITER_PLUS = value;
					} else if (name.equals("DELIMITER_EQUALS")) {
						DELIMITER_EQUALS = value;
					} else if (name.equals("DELIMITER_FORWARD")) {
						DELIMITER_FORWARD = value;
					} else if (name.equals("DELIMITER_COMPOUND_COMPARTMENT")) {
						DELIMITER_COMPOUND_COMPARTMENT = value;
					} else if (name.equals("DEFAULT_COMPARTMENT")) {
						DEFAULT_COMPARTMENT = value;						

						
					// String[] types
					} else if (name.equals("ELEMENTS")) {
						String[] elements = value.split(",");
						// check if hydrogen and phosphate exists, necessary for fixing the balance
						boolean hasHydrogen = false, hasOxygen = false, hasPhosphate = false;
						for (int i=0; i<elements.length; i++) {
							if (elements[i].equals("H"))
								hasHydrogen = true;
							else if (elements[i].equals("O"))
								hasOxygen = true;
							else if (elements[i].equals("P"))
								hasPhosphate = true;
						}
						if (!hasHydrogen || !hasOxygen || !hasPhosphate)
							System.out.println("ELEMENTS must contain at least H, O, and P. Using default value.");
						else
							ELEMENTS = elements;
						
					} else if (name.equals("PHOSPHATE_NAMES")) {
						
						String[] elements = value.split(",");
						// check if the phosphate names are unique
						boolean unique = false;
						for (int i=0; i<elements.length-1; i++) {
							for (int j=i+1; j<elements.length; j++) {
								if (elements[i].equals(elements[j]))
									break;
								else if (i==elements.length-2 && j==elements.length-1)
									unique = true;
							}
						}
						if (!unique || elements.length != 4)
							System.out.println("PHOSPHATE_NAMES must be 4 unique strings. Using default values.");
						else
							PHOSPHATE_NAMES = elements;
						
					} else {
						System.out.println("Invalid parameter "+name+" in configuration file "+CONFIG_FILE);
					}
					
					// re-set RATIO_CLASSES_REACTIONS
					RATIO_CLASSES_REACTIONS = RATIO_CLASSES_COMPOUNDS * RATIO_COMPOUNDS_REACTIONS;
					// create the mass vectors for hydrogen, phosphate, hydrogen phosphate, dihydrogen phosphate, and phosphoric acid
					HYDROGEN_MASS = new int[ELEMENTS.length];
					PHOSPHATE_FORMS = new int[4][ELEMENTS.length];
					for (int i=0; i<ELEMENTS.length; i++) {
						if (ELEMENTS[i].equals("H")) {
							HYDROGEN_MASS[i] = 1;
							PHOSPHATE_FORMS[0][i] = 0;
							PHOSPHATE_FORMS[1][i] = 1;
							PHOSPHATE_FORMS[2][i] = 2;
							PHOSPHATE_FORMS[3][i] = 3;
						} else if (ELEMENTS[i].equals("O")) {
							PHOSPHATE_FORMS[0][i] = 4;
							PHOSPHATE_FORMS[1][i] = 4;
							PHOSPHATE_FORMS[2][i] = 4;
							PHOSPHATE_FORMS[3][i] = 4;
						} else if (ELEMENTS[i].equals("P")) {
							PHOSPHATE_FORMS[0][i] = 1;
							PHOSPHATE_FORMS[1][i] = 1;
							PHOSPHATE_FORMS[2][i] = 1;
							PHOSPHATE_FORMS[3][i] = 1;
						}
					}
					
				} catch (NumberFormatException e1) {
					System.out.println("Invalid value: "+value+" for parameter "+name+" in configuration file "+CONFIG_FILE+". Using default value.");
				}
			}
			
			configReader.close();
			
		} catch (FileNotFoundException e) {
			// continue with the default configuration
		} catch (IOException e) {
			System.err.println("I/O error loading the configuration file "+CONFIG_FILE);
			e.printStackTrace();
			System.exit(-1);
		}
		
		if (debug) {
			System.out.println("*** loadConfiguration ***");
			System.out.print("ELEMENTS = "); 
			for (int i=0; i<ELEMENTS.length; i++)
				System.out.print(ELEMENTS[i]+" ");
			System.out.println();
			System.out.println("HYDROGEN_NAME = "+HYDROGEN_NAME);
			System.out.print("PHOSPHATE_NAMES = ");
			for (int i=0; i<Utilities.PHOSPHATE_NAMES.length; i++)
				System.out.print(Utilities.PHOSPHATE_NAMES[i]+" ");
			System.out.println();
			System.out.print("PHOSPHATE_FORMS = ");
			for (int i=0; i<Utilities.PHOSPHATE_FORMS.length; i++) {
				for (int j=0; j<Utilities.ELEMENTS.length; j++) {
					System.out.print(Utilities.ELEMENTS[j]+Utilities.PHOSPHATE_FORMS[i][j]+" ");
					System.out.print(" ");
				}
				System.out.print("; ");
			}
			System.out.println();
			System.out.println("DELIMITER_PLUS = "+DELIMITER_PLUS);
			System.out.println("DELIMITER_EQUALS = "+DELIMITER_EQUALS);
			System.out.println("DELIMITER_FORWARD = "+DELIMITER_FORWARD);
			System.out.println("DELIMITER_COEFFICIENT = "+DELIMITER_COEFFICIENT_START+" - "+DELIMITER_COEFFICIENT_END);
			System.out.println("DELIMITER_COMPOUND_COMPARTMENT = "+DELIMITER_COMPOUND_COMPARTMENT);
			System.out.println("DELIMITER_REACTION_COMPARTMENT = "+DELIMITER_REACTION_COMPARTMENT_START+" - "+DELIMITER_REACTION_COMPARTMENT_END);
			System.out.println("RATIO_CLASSES_COMPOUNDS = "+RATIO_CLASSES_COMPOUNDS);
			System.out.println("RATIO_COMPOUNDS_REACTIONS = "+RATIO_COMPOUNDS_REACTIONS);
			System.out.println("RATIO_CLASSES_REACTIONS = "+RATIO_CLASSES_REACTIONS);
			System.out.println("MAX_STOICHIOMETRY_FACTOR = "+MAX_STOICHIOMETRY_FACTOR);
			System.out.println("DEFAULT_NUMBER_OF_COMPOUNDS = "+DEFAULT_NUMBER_OF_COMPOUNDS);
			System.out.println("DEFAULT_NUMBER_OF_REACTIONS = "+DEFAULT_NUMBER_OF_REACTIONS);
			System.out.println("DEFAULT_NUMBER_OF_CLASSES = "+DEFAULT_NUMBER_OF_CLASSES);
		}
	}
	
	/**
	 * Matrix multiplication of two-dimensional primitive double arrays. 
	 * @param left
	 * @param right
	 * @return
	 */
	protected static double[][] matrixMultiply(double[][] left, double[][] right) {
		if (left.length != right[0].length || left[0].length != right.length)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		
		int m = left.length;
		int n = right.length;
		double[][] result = new double[m][n];
		
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				// calculate the new i,j
				double value = 0;
				for (int k=0; k<m; k++)
					value += (left[i][k]*right[k][j]);
				result[i][j] = value;
			}
		}
		
		return result;
	}
	
	/**
	 * Performs matrix multiplication of a Fortran-style matrix array with a vector
	 * (used as external reverse communication driver for ARPACK function SNAUPD).
	 * The matrix is multiplied with the vector starting at index inOffset, the result
	 * is saved into the same vector starting at index outOffset. Therefore, the vector
	 * needs to be at least of size inOffset+outOffset.
	 * 
	 * @param matrix
	 * @param vector
	 * @return
	 */
	protected static void matrixVectorMultiplyF(float[] matrix, int m, float[] vector, int inOffset, int outOffset) {
		
		if (vector.length < inOffset+outOffset)
			throw new IllegalArgumentException("Vector has to be at least of size inOffset+outOffset="+(inOffset+outOffset)+".");
		
		if (m > vector.length/2)
			throw new IllegalArgumentException("Matrix dimension too large for input-output vector of size "+vector.length+".");
		
		for (int i=0; i<m; i++) {
			int iSum = 0;
			for (int j=0; j<m; j++)
				iSum += (matrix[m*j+i]*vector[inOffset+j]);
			vector[outOffset+i] = iSum;
		}
	}
	
	protected static void testMatrixVectorMultiplyF() {
//		int m = 3000;
//		float[] A = new float[m*m];
//		float[] V = new float[3*m];
//		float[] X = new float[m];
//		float[] Y = new float[m];
//		Random r = new Random();
//		for (int i=0; i<A.length; i++)
//			A[i] = r.nextFloat();
//		for (int i=0; i<V.length; i++)
//			V[i] = r.nextFloat();
//		for (int i=0; i<X.length; i++)
//			X[i] = r.nextFloat();
		float[] A = {1, 4, 7, 2, 5, 8, 3, 6, 9};
		int m = 3;
		float[] V = {0, 0, 0, 0, 0, 0, 2, 2, 2};
		float[] X = {1, 1, 1, 2, 2, 2, 3, 3, 3};
		float[] Y = new float[m];
		
		long time1 = System.currentTimeMillis();
		matrixVectorMultiplyF(A, m, V, 2*m, m);
		System.out.println("matrixVectorMultiplyF: "+((System.currentTimeMillis()-time1)/(float)1000)+"s.");
		long time2 = System.currentTimeMillis();
		// org.netlib.blas.Sgemv.sgemv("N", m, m, 1f, A_float, 0, m, X, 0, 1, 0f, Y, 0, 1);
		org.netlib.blas.Sgemv.sgemv("N", m, m, 1f, A, 0, m, X, 0, 3, 0f, Y, 0, 1);
		System.out.println("sgemv: "+((System.currentTimeMillis()-time2)/(float)1000)+"s.");
		
		System.out.println("V=");
		for (int i=0; i<V.length; i++)
			System.out.print(V[i]+"\t");
		System.out.println();

		System.out.println("X=");
		for (int i=0; i<X.length; i++)
			System.out.print(X[i]+"\t");
		System.out.println();
		
		System.out.println("Y=");
		for (int i=0; i<Y.length; i++)
			System.out.print(Y[i]+"\t");
		System.out.println();
	}
	
	/**
	 * Matrix addition of two-dimensional primitive double arrays. 
	 * @param left
	 * @param right
	 * @return
	 */
	protected static double[][] matrixSum(double[][] matrix1, double[][] matrix2) {
		if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length)
			throw new IllegalArgumentException("Matrix dimensions do not match.");
		
		int m = matrix1.length;
		int n = matrix1[0].length;
		double[][] result = new double[m][n];
		
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				// calculate the new i,j
				result[i][j] = matrix1[i][j]+matrix2[i][j];
			}
		}
		
		return result;
	}
	
	/**
	 * Calculates the edit distance of strings s1 and s2.
	 * @param s1
	 * @param s2
	 * @return
	 */
	protected static int editDistance(String s1, String s2) {
		int[][] matrix = new int[s1.length()+1][s2.length()+1];
		
		for (int i=0; i<matrix.length; i++) {
			for (int j=0; j<matrix[i].length; j++) {
				matrix[i][j] = (i==0) ? j : ((j==0) ? i : 0);
				if (i>0 && j>0) {
					if (s1.charAt(i-1) == s2.charAt(j-1))
						matrix[i][j] = matrix[i-1][j-1];
					else
						matrix[i][j] = Math.min(matrix[i][j-1]+1, Math.min(matrix[i-1][j-1]+1, matrix[i-1][j] + 1));
				}
			}
		}
		return matrix[s1.length()][s2.length()];
	}
	
	/**
	 * Performance test of different linear algebra packages for eigenvalue
	 * composition of the reaction-essentiality matrix.
	 * 
	 * Results on pittakos: 2 x Intel(R) Xeon(R) E5345 2.33GHz (8 cores, 16 GB memory, Linux Fedora 13):
	 * 
	 * UJMP:   47.994 seconds
	 * JAMA:   47.453 seconds.
	 * Colt:   51.008 seconds.
	 * JBLAS:  67.908 seconds.
	 * MTJ:    70.484 seconds.
	 * JLAPACK (float):        4.5 seconds.
	 * JLAPACK (double):       7.1 seconds.
	 *  
	 * @param graph
	 */
//	public static void performanceDecomposition(MetabolicGraph graph) {
//		Collection<Vertex> reactionsSet = new HashSet<Vertex>(graph.getReactions());
//		
//		System.out.print("Sorting...");
//		
//		// add the reversed reactions
//		for (Vertex reaction : graph.getReactions())
//			if (reaction.reversedReaction() != null)
//				reactionsSet.add(reaction.reversedReaction());
//		
//		// sort the reactions by name (and implicitly by reversedness)
//    	ArrayList<Vertex> reactionsList = new ArrayList<Vertex>(reactionsSet);
//    	Collections.sort(reactionsList, new Comparator<Vertex>() {
//    		public int compare(Vertex r1, Vertex r2) {
//    			return (r1.getName().compareTo(r2.getName()));
//    		}});
//    	
//    	System.out.println("Creating essentiality matrix...");
//    	
//		// construct the transposed reaction incidence matrix multiplied by the local essentialities and row-normalized
//    	int m = reactionsList.size();
//    	
//    	// initialize the matrices
//    	Jama.Matrix jamaMatrix = new Jama.Matrix(m, m);
//    	cern.colt.matrix.impl.SparseDoubleMatrix2D coltMatrix = new SparseDoubleMatrix2D(m, m);
////    	org.ujmp.core.Matrix ujmpMatrix = org.ujmp.core.MatrixFactory.sparse(m,m);
//    	no.uib.cipr.matrix.DenseMatrix mtjMatrix = new no.uib.cipr.matrix.DenseMatrix(m, m);
//    	org.jblas.DoubleMatrix jblasMatrix = new org.jblas.DoubleMatrix(m, m);
//    	double[][] jlapackMatrix = new double[m][m];
//    	float[][] jlapackMatrix_float = new float[m][m];
//    	shared.array.RealArray sharedMatrix = new shared.array.RealArray(m, m);
//    	
//		for (int i=0; i<m; i++) {
//			for (int j=0; j<m; j++) {
//				// set the local essentiality of reaction i for reaction j
//				double essentiality = graph.localEssentiality(reactionsList.get(i), reactionsList.get(j));
//				jamaMatrix.set(i, j, essentiality);
//				coltMatrix.set(i, j, essentiality);
////				ujmpMatrix.setAsDouble(essentiality, i, j);
//				mtjMatrix.set(i, j, essentiality);
//				jblasMatrix.put(i, j, essentiality);
//				jlapackMatrix[i][j] = essentiality;
//				jlapackMatrix_float[i][j] = (float)essentiality;
//				sharedMatrix.set(essentiality, i, j);
//			}
//		}
//		
//		// normalize colmuns and set the zero-columns to 1/m
//		for (int j=0; j<m; j++) {
//			double sum = 0;
//			for (int i=0; i<m; i++)
//				sum += jamaMatrix.get(i, j);
//			for (int i=0; i<m; i++) {
//				if (sum == 0) {
//					jamaMatrix.set(i, j, 1/m);
//					coltMatrix.set(i, j, 1/m);
////					ujmpMatrix.setAsDouble(1/m, i, j);
//					mtjMatrix.set(i, j, 1/m);
//					jblasMatrix.put(i, j, 1/m);
//					jlapackMatrix[i][j] = 1/m;
//					jlapackMatrix_float[i][j] = 1/m;
//					sharedMatrix.set(1/m, i, j);
//				} else {
//					jamaMatrix.set(i, j, jamaMatrix.get(i, j)/sum);
//					coltMatrix.set(i, j, coltMatrix.get(i, j)/sum);
////					ujmpMatrix.setAsDouble(ujmpMatrix.getAsDouble(i, j)/sum);
//					mtjMatrix.set(i, j, mtjMatrix.get(i, j)/sum);
//					jblasMatrix.put(i, j, jblasMatrix.get(i, j)/sum);
//					jlapackMatrix[i][j] = jlapackMatrix[i][j]/sum;
//					jlapackMatrix_float[i][j] = jlapackMatrix_float[i][j]/(float)sum;
//					sharedMatrix.set(sharedMatrix.get(i,j)/sum, i, j);
//				}
//			}
//		}
//		
//		System.out.println("Eigenvalue decomposition...");
//		
////		System.out.print("UJMP:\t");
////		long start = System.currentTimeMillis();
////		org.ujmp.core.Matrix[] ujmp = new org.ujmp.core.Matrix[m];
////		ujmpMatrix.eig();
////		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
////		System.gc();
//
//		System.out.print("Shared:\t");
//		long start = System.currentTimeMillis();
//		sharedMatrix.mEigs();
//		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
//		System.gc();
//		
//		System.out.print("JLAPACK (float):\t");
//		start = System.currentTimeMillis();
//		// NOTE: check args 1, 2, 4, 5, 
//		float[] TAU_float = new float[m];
//		float[] WORK_float = new float[m];
//		org.netlib.util.intW intW_float = new org.netlib.util.intW(0);
//		org.netlib.lapack.SGEHRD.SGEHRD(m, 1, m, jlapackMatrix_float, TAU_float, WORK_float, m, intW_float);
////		.dgehrd(m, 1, m, jlapackMatrix, m, TAU, WORK, m, arg8, arg9, arg10, arg11);
//		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
//		System.gc();
//		
//		System.out.print("JLAPACK (double):\t");
//		start = System.currentTimeMillis();
//		// NOTE: check args 1, 2, 4, 5, 
//		double[] TAU = new double[m];
//		double[] WORK = new double[m];
//		org.netlib.util.intW intW_double = new org.netlib.util.intW(0);
//		org.netlib.lapack.DGEHRD.DGEHRD(m, 1, m, jlapackMatrix, TAU, WORK, m, intW_double);
////		.dgehrd(m, 1, m, jlapackMatrix, m, TAU, WORK, m, arg8, arg9, arg10, arg11);
//		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
//		System.gc();
//		
////		System.out.print("JAMA:\t");
////		start = System.currentTimeMillis();
////		Jama.EigenvalueDecomposition jama = new Jama.EigenvalueDecomposition(jamaMatrix);
////		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
////		System.gc();
////		
////		System.out.print("Colt:\t");
////		start = System.currentTimeMillis();
////		cern.colt.matrix.linalg.EigenvalueDecomposition colt = new cern.colt.matrix.linalg.EigenvalueDecomposition(coltMatrix);
////		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
////		System.gc();
////		
////		System.out.print("JBLAS:\t");
////		start = System.currentTimeMillis();
////		org.jblas.Eigen.eigenvectors(jblasMatrix);
////		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
////		System.gc();
////		
////		System.out.print("MTJ:\t");
////		start = System.currentTimeMillis();
////		no.uib.cipr.matrix.EVD mtj = new no.uib.cipr.matrix.EVD(m, true, false);
////		try {
////			mtj.factor(mtjMatrix);
////		} catch (NotConvergedException e) {
////			e.printStackTrace();
////		}
////		System.out.println((float)(System.currentTimeMillis()-start)/(float)1000+" seconds.");
//		
//	}
	
	/**
	 * Print usage and exit.
	 */
	protected static void exit() {
		System.out.println(PARAMETERS);
		Class<?> utilities;
		try {
			utilities = Class.forName("massbalance.Utilities");
			System.out.println("Available methods:");
			Method[] methods = utilities.getMethods();
			for (int i=0; i<methods.length; i++) {
				if (methods[i].getDeclaringClass().equals(utilities) && Modifier.isPublic(methods[i].getModifiers()) && Modifier.isStatic(methods[i].getModifiers())) {
					boolean validTypes = true;
					// check for implemented parameter types
					Class<?>[] types = methods[i].getParameterTypes();
					for (int j=0; j<types.length; j++)
						if (!types[j].equals(int.class) && !types[j].equals(String.class) && !types[j].equals(boolean.class)
								&& !types[j].equals(MetabolicGraph.class) && !types[j].equals(ArrayList.class)) {
							validTypes = false;
							break;
						}
					if (validTypes) {
						String methodString = methods[i].toString();
						System.out.println("   "+methodString.substring(methodString.indexOf("Utilities")+10, methodString.indexOf(")")+1));
					}
				}
			}
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
		System.exit(-1);
	}
	
	/**
	 * Writes every tenth element from every line of the input file to the output file.
	 *  
	 * @param inputFile
	 * @param outputFile
	 * @throws IOException
	 */
	public static void select(String inputFile, String outputFile, boolean append) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, append));
		String line;
		int values = 0, lines = 0;
		
		while ((line = reader.readLine()) != null) {
			int start = 0, element = 0;
			int end = line.indexOf("\t", start);
			lines++;
			
			while (end != -1) {
				
				if (element % 10 == 0) {
					writer.write(line.substring(start, end)+"\t");
					values++;
				}
				start = end+1;
				end = line.indexOf("\t", start);
				element++;
				// last column
				if (end == -1 && element % 10 == 0) {
					writer.write(line.substring(start, line.length()));
					values++;
				}
			}
			writer.write("\n");
		}
		
		reader.close();
		writer.close();
		
		System.out.println("Copied "+values+" values in "+lines+" lines.");
	}
	
	/**
	 * Prints the lines in the given from-to range of the given inputFile.
	 * @param inputFile
	 * @param from
	 * @param to
	 * @throws IOException
	 */
	public static void writeLines(String inputFile, int from, int to) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(inputFile+".lines"+from+"-"+to, false));
		
		String line;
		int number = 1;
		while ((line = reader.readLine()) != null) {
			if (number >= from && number <= to)
				writer.write(line+"\n");
			number++;
		}
		
		reader.close();
		writer.close();
	}
	
	/**
	 * Appends the content from inputFile to the end of outputFile. Careful when applying
	 * to files with unfinished lines: remove the unfinished line before.
	 *  
	 * @param inputFile
	 * @param outputFile
	 * @throws IOException
	 * @return number of lines in the output file
	 */
	public static int appendFile(String inputFile, String outputFile) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, true));
		String line;
		int lines = 0;
		
		while ((line = reader.readLine()) != null) {
			writer.write(line+"\n");
			lines++;
		}
		
		reader.close();
		writer.close();
		return lines;
	}
	
	/**
	 * Appends the content from inputFile and the next numberOfFiles-1 inputFiles with
	 * successive index numbers to the end of outputFile. E.g., "in.1.file 3 outFile" appends
	 * the contents of "in.1.file", "in.2.file", and "in.3.file", to "outFile". Index ranges
	 * are also supported if the file name contains a dash, e.g. "in.0-9.file 3 outFile"
	 * appends the contents of "in.0-9.file", "in.10-19.file", and "in.20-29.file" to "outFile".
	 * 
	 * Careful when applying to files with unfinished lines: remove the unfinished line before.
	 *  
	 * @param firstInputFile
	 * @param numberOfFiles
	 * @param outputFile
	 * @throws IOException
	 */
	public static void appendFiles(String firstInputFile, int numberOfFiles, String outputFile) throws IOException {
		int from = Integer.MAX_VALUE, to = Integer.MAX_VALUE;
		String inputDir = new File(firstInputFile).getParent();
		firstInputFile = new File(firstInputFile).getName();
		int dash = firstInputFile.lastIndexOf("-");
		int lines = 0;
		
		if (dash == -1) {
			// parse index
			int index = 0;
			int end = firstInputFile.length()-1;
			while (!Character.isDigit(firstInputFile.charAt(end)))
				end--;
			int start = end;
			try {
				while (start >= 0)
					index = Integer.parseInt(firstInputFile.substring(start--, end+1));
			} catch (NumberFormatException e) {
				// parsed longest possible number
			}
			// append the files
			for (int i=0; i<numberOfFiles; i++) {
				System.out.println("Appending '"+firstInputFile+"'.");
				lines += appendFile(inputDir+File.separator+firstInputFile, outputFile);
				String indexString = String.valueOf(index);
				String newIndex = String.valueOf(++index);
				firstInputFile = firstInputFile.replace(indexString, newIndex);
			}
		} else {
			// parse from-to indices
			int start = dash-1;
			try {
				while (start >= 0)
					from = Integer.parseInt(firstInputFile.substring(start--, dash));
			} catch (NumberFormatException e) {
				// parsed longest possible number
			}
			// parse to-index
			int end = dash+2;
			try {
				while (end < firstInputFile.length())
					to = Integer.parseInt(firstInputFile.substring(dash+1, end++));
			} catch (NumberFormatException e) {
				// parsed longest possible number
			}
			
			// append the files
			int stepSize = to-from+1;
			for (int i=0; i<numberOfFiles; i++) {
				System.out.println("Appending '"+firstInputFile+"'.");
				lines += appendFile(inputDir+File.separator+firstInputFile, outputFile);
				String fromTo = ("."+from+"-"+to+".");
				from += stepSize;
				to += stepSize;
				String newFromTo = ("."+from+"-"+to+".");
				firstInputFile = firstInputFile.replace(fromTo, newFromTo);
			}
		}
		System.out.println("Output file: "+new File(outputFile).getName()+" ("+lines+" lines appended).");
	}
	
	/**
	 * Appends the tab-separated tokens from each line in inputFile2 to the corresponding
	 * lines in inputFile1 and prints the resulting lines to outputFile.
	 * 
	 * Lines starting with '#' are skipped. In addition, if offset<0, the inputFile1 is
	 * skipped by -offset lines, if offset>0, inputFile2 is skipped by offset lines. 
	 * 
	 * @param inputFile1
	 * @param inputFile2
	 * @param outputFile
	 * @param offset
	 * @throws IOException
	 */
	public static void appendColumns(String inputFile1, String inputFile2, String outputFile, int offset) throws IOException {
		BufferedReader reader1 = new BufferedReader(new FileReader(inputFile1));
		BufferedReader reader2 = new BufferedReader(new FileReader(inputFile2));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
		String line1 = null, line2 = null;
		boolean stop = false;
		
		if (offset < 0) {
			// advance inputFile1
			while ((line1 = reader1.readLine()) != null && offset < 0) {
				if (!line1.startsWith("#")) {
					offset++;
					writer.write(line1+"\n");
				}
			}
		} else if (offset > 0) {
			// advance inputFile2
			while ((line2 = reader2.readLine()) != null && offset > 0) {
				if (!line2.startsWith("#")) {
					offset--;
					writer.write(line2+"\n");
				}
			}
		} else {
			// read the first lines
			line1 = reader1.readLine();
			line2 = reader2.readLine();
		}
		
		do {
			// skip inputFile1 comments
			if (line1.startsWith("#"))
				continue;
			// skip inputFile2 comments
			while ((line2 = reader2.readLine()) != null && line2.startsWith("#"))
				;
			
			if (line2 == null) {
				stop = true;
				break;
			}
			
			writer.write(line1+"\t"+line2+"\n");
		} while ((line1 = reader1.readLine()) != null);
		
		if (stop)
			System.out.println(inputFile1+"\nhas more (unskipped) lines than\n"+inputFile2);
		else if ((line2 = reader2.readLine()) != null) {
			System.out.println(inputFile1+"\nhas less (unskipped) lines than\n"+inputFile2);
			System.out.println(line2+"\n");
		}
		
		reader1.close();
		reader2.close();
		writer.close();
	}
	
	/**
	 * Searches in inputFile2 for the first string of each line in inputFile1, and
	 * joins the remaining strings in outputFile.
	 * 
	 * @param inputFile1
	 * @param inputFile2
	 * @param outputFile
	 * @throws IOException
	 */
	public static void joinLines(String inputFile1, String inputFile2, String outputFile) throws IOException {
		BufferedReader reader1 = new BufferedReader(new FileReader(inputFile1));
		BufferedReader reader2 = new BufferedReader(new FileReader(inputFile2));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		String line;
		HashMap<String,String> map = new HashMap<String,String>(); 
		
		// create a hash map from the first key strings and the remaining values in inputFile2
		while ((line = reader2.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			map.put(line.substring(0, line.indexOf("\t")).trim(), line.substring(line.indexOf("\t")+1).trim());
		}
		
		// parse through inputFile1 and join the corresponding lines in outputFile
		while ((line = reader1.readLine()) != null) {
			String key = line.substring(0, line.indexOf("\t")).trim();
			String values1 = line.substring(line.indexOf("\t")+1).trim();
			String values2 = map.get(key);
			if (values2 != null) {
				// prefix the deltaG, and uncertainty
//				values2 = "D="+values2.replace("\t", "\tU=");
				writer.write(key+"\t"+values1+"\t"+values2+"\n");
			} else {
				System.out.println("Not found: "+key);
			}
		}
		
		reader1.close();
		reader2.close();
		writer.close();
	}
	
	/**
	 * Prints the largest single equivalence class.
	 * 
	 * @param classesFile
	 * @throws IOException
	 */
	public static void maxSingleClass(String classesFile) throws IOException, InterruptedException {
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = Equivalence.readClasses(classesFile);
		ArrayList<Integer> maxKey = classes.keySet().iterator().next();
		
		for (ArrayList<Integer> key : classes.keySet()) {
			if (classes.get(maxKey).getNames().size() < classes.get(key).getNames().size())
				maxKey = key;
		}
		System.out.println();
		EquivalenceClass eqClass = classes.get(maxKey);
		for (int i=0; i<eqClass.getNames().size(); i++) {
			System.out.print(eqClass.getName(i)+"\t");
			for (int j=0; j<ELEMENTS.length; j++)
				System.out.print(eqClass.getMass(i)[j]+"\t");
		}
		System.out.println();
		System.out.println(eqClass.getNames().size()+" compounds.");
	}
	
	/**
	 * Prints the largest pair equivalence class, excluding implicit pairs.
	 * 
	 * @param classesFile
	 * @throws IOException
	 */
	public static void maxExplicitDoubleClass(String classesFile) throws IOException, InterruptedException {
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = Equivalence.readClasses(classesFile);
		ArrayList<Integer> maxKey = classes.keySet().iterator().next();
		
		for (ArrayList<Integer> key : classes.keySet()) {
			if (classes.get(maxKey).getExplicitDoubleNames().size() < classes.get(key).getExplicitDoubleNames().size())
				maxKey = key;
		}
		EquivalenceClass eqClass = classes.get(maxKey);
		System.out.println();
		for (int i=0; i<eqClass.getExplicitDoubleNames().size(); i++) {
			System.out.print(eqClass.getDoubleName(i)+"\t");
			for (int j=0; j<ELEMENTS.length; j++)
				System.out.print(eqClass.getDoubleMass(i)[j]+"\t");
		}
		System.out.println();
		System.out.println(eqClass.getExplicitDoubleNames().size()/2+" pairs of compounds.");
	}
	
	/**
	 * Prints the largest pair equivalence class, including implicit pairs.
	 * 
	 * @param classesFile
	 * @throws IOException
	 */
	public static void maxDoubleClass(String classesFile) throws IOException, InterruptedException {
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = Equivalence.readClasses(classesFile);
		ArrayList<Integer> maxKey = classes.keySet().iterator().next();
		
		for (ArrayList<Integer> key : classes.keySet()) {
			if (classes.get(maxKey).getPairs().size() < classes.get(key).getPairs().size())
				maxKey = key;
		}
		EquivalenceClass eqClass = classes.get(maxKey);
		for (int i=0; i<eqClass.getExplicitDoubleNames().size(); i++) {
			System.out.print(eqClass.getDoubleName(i)+"\t");
			for (int j=0; j<ELEMENTS.length; j++)
				System.out.print(eqClass.getDoubleMass(i)[j]+"\t");
		}
		System.out.println();
		System.out.println(eqClass.getPairs().size()+" pairs of compounds.");
	}
	
	/**
	 * Prints the equivalence classes in which the given compound is contained.
	 * 
	 * @param compound
	 */
	public static void printEquivalenceClass(String classesFile, String compound) throws IOException, InterruptedException {
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = Equivalence.readClasses(classesFile);
		
		for (ArrayList<Integer> key : classes.keySet()) {
			boolean hit = false;
			for (String single : classes.get(key).getNames()) {
				if (single.equals(compound)) {
					hit = true;
					break;
				}
			}
			if (hit) {
				System.out.print("Single class: ");
				for (String single : classes.get(key).getNames()) {
					System.out.print(single+"; ");
				}
				System.out.println();
				hit = false;
			}
			
			for (String[] pair : classes.get(key).getPairs().keySet()) {
				if (pair[0].equals(compound) || pair[1].equals(compound)) {
					hit = true;
					break;
				}
			}
			if (hit) {
				System.out.print("Pairs class: ");
				for (String[] pair : classes.get(key).getPairs().keySet()) {
					System.out.print(pair[0]+" + "+pair[1]+"; ");
				}
				System.out.println();
			}
		}
		
	}
	
	/**
	 * Prints the number of single and double mass equivalence classes and their average sizes.
	 * 
	 * @param classesFile
	 * @throws IOException
	 */
	public static void avgClassSizes(String classesFile) throws IOException, InterruptedException {
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = Equivalence.readClasses(classesFile);
		double singleSizesSum = 0, doubleSizesSum = 0;
		int singleClasses = 0, doubleClasses = 0;
		
		for (ArrayList<Integer> key : classes.keySet()) {
			if (classes.get(key).getNames().size() > 0) {
				singleClasses++;
				singleSizesSum += classes.get(key).getNames().size();
			}
			if (classes.get(key).getPairs().size() > 0) {
				doubleClasses++;
				doubleSizesSum += classes.get(key).getPairs().size();
			}
		}
		System.out.println(singleClasses+" single classes, average size: "+(singleSizesSum/singleClasses));
		System.out.println(doubleClasses+" double classes, average size: "+(doubleSizesSum/doubleClasses));
	}
	
	/**
	 * Checks readClasses() by comparing the created classes to the written and read classes.
	 * 
	 * @param inputDir
	 * @param outputDir
	 * @return
	 * @throws IOException
	 */
	protected boolean checkReadClasses(String inputFile, String outputDir) throws IOException, InterruptedException {
		
		String version = Utilities.getVersion(inputFile);
		String classesFile = outputDir+File.separator+version+".classes";
		
		Equivalence.main(new String[]{inputFile, outputDir});
		HashMap<ArrayList<Integer>, EquivalenceClass> classes1 = Equivalence.classes;
		
		// read the classes
		HashMap<ArrayList<Integer>, EquivalenceClass> classes2 = Equivalence.readClasses(classesFile);
		// compare
		return Equivalence.compareClasses(classes1, classes2);
	}
	
	/**
	 * Prints all compartments in the given network.
	 * 
	 * @param reactionsFile
	 * @throws IOException
	 */
	public static void writeCompartments(String reactionsFile) throws IOException {
		BufferedReader reactionsReader = new BufferedReader(new FileReader(reactionsFile));
		String line;
		HashSet<String> compounds = new HashSet<String>();
		HashSet<String> compartments = new HashSet<String>();
		
		while ((line = reactionsReader.readLine()) != null) {
			
			if (line.startsWith("#"))
				continue;
			// parse the reaction equation
			String equation = line.substring(line.indexOf("\t")+1, line.length());
			boolean newToken = false;
			
			for (int i=0, j=0; i<equation.length(); i++) {
				boolean lastToken = false;
				int increment = 0;
				
				if (i == equation.length()-1)
					lastToken = true;
				
				if (equation.charAt(i) == DELIMITER_COEFFICIENT_START) {
					// skip the stoichiometric coefficient
					i = equation.indexOf(DELIMITER_COEFFICIENT_END, i);
					j = i+1;
					
				} else if (lastToken && newToken || (increment = Utilities.isDelimiter(equation, i)) != -1) {
					String compartment = null;
					// add the previous token
					String compound = equation.substring(j, lastToken ? i+1 : i).trim();
					// divide into name and compartment
					String name = compound.substring(0, compound.indexOf(DELIMITER_COMPOUND_COMPARTMENT));
					compartment = compound.substring(compound.indexOf(DELIMITER_COMPOUND_COMPARTMENT)+1, compound.length());
					
					compounds.add(name);
					compartments.add(compartment);
					System.out.println("Compound: "+name+", compartment: "+compartment);
					
					newToken = false;
					i += (increment-1);
					
				} else if (newToken == false) {
					// beginning of a new token
					newToken = true;
					j = i;
				}
			}
		}
		
		System.out.print("Network has "+compartments.size()+" compartments: ");
		for (String compartment : compartments)
			System.out.print(compartment+" ");
		System.out.println();
		
	}
	
	/**
	 * Tests the serialized graphs in the given directory of the given species version. If
	 * massbalance==true massbalance randomized files are tested, otherwise switch randomized files.
	 * 
	 * @param inputDir
	 * @param version
	 * @param massbalance
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void testGraphs(String inputDir, String version, boolean massbalance, int from, int to) throws IOException, ClassNotFoundException {
		int success = 0, errors = 0, edges;
		String randomization, file = "";
		
		if (!new File(inputDir).isDirectory())
			throw new IllegalArgumentException("Not a directory: "+inputDir);
		
		if (massbalance) {
			randomization = "massbalance";
		} else {
			randomization = "switch";
		}
		
		file = inputDir+File.separator+version;
		ObjectInputStream reader = new ObjectInputStream(new FileInputStream(file));
		MetabolicGraph graph = (MetabolicGraph)reader.readObject();
		reader.close();
		edges = graph.edgeSet().size();
		
		for (int i=from; i<=to; i++) {
			File randomizedFile = new File(inputDir+File.separator+version+"."+randomization+"."+i);
			try {
				ObjectInputStream randomizedReader = new ObjectInputStream(new FileInputStream(randomizedFile));
				MetabolicGraph randomizedGraph = (MetabolicGraph)randomizedReader.readObject();
				randomizedReader.close();
				if (edges != randomizedGraph.edgeSet().size()) {
					errors++;
					System.out.println("Different number of edges: "+edges+" in "+version+"."+randomization+"."+(i-1)+", "+randomizedGraph.edgeSet().size()+" in "+version+"."+randomization+"."+i);
				} else {
					success++;
				}
			} catch (IOException e) {
				errors++;
				System.out.println("Error reading "+randomizedFile.getName()+": "+e); 
			}
		}
		System.out.println("Read "+success+" files successfully, "+errors+" error(s).");
	}
	
	/**
	 * Creates a graph from the biocyc or flat file input directory and deserializes it to the output file name. 
	 * 
	 * @param networkFile
	 * @param outputFile
	 * @param reversible
	 * @param compartments
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void createGraph(String networkFile, String outputFile, boolean reversible, boolean compartments) throws IOException, ClassNotFoundException {
		String compoundsFile, pathwayFile;
		HashMap<String, int[]> massesMap;
		
		String infoFile = networkFile+".checkNetwork.log";
		MetabolicGraph graph;
		
		if (networkFile.endsWith(".xml")) {
			// SBML format
			SBMLDocument sbmlDoc = Utilities.initSBML(networkFile);
			Model sbmlModel = sbmlDoc.getModel();
			massesMap = Utilities.createMasses(sbmlModel, infoFile);
			graph = Properties.createSBMLGraph(networkFile, infoFile, massesMap, reversible, compartments, true);
			
		} else if (networkFile.endsWith(".dat")) {
			// biocyc format
			compoundsFile = networkFile.substring(0, networkFile.lastIndexOf(File.separator)+1)+"compounds.dat";
			pathwayFile = networkFile.substring(0, networkFile.lastIndexOf(File.separator)+1)+"pathways.dat";
			if (!new File(compoundsFile).exists()) {
				System.err.println("Missing compounds file: "+compoundsFile);
				System.exit(-1);
			}
			massesMap = Utilities.createMassesBioCyc(compoundsFile, networkFile, infoFile);
			graph = Properties.createBiocycGraph(networkFile, pathwayFile, infoFile, massesMap, reversible, compartments, true);
			
		} else {
			// text format
			compoundsFile = networkFile+".compounds";
			if (!new File(compoundsFile).exists()) {
				System.err.println("Missing compounds file: "+compoundsFile);
				System.exit(-1);
			}
			massesMap = Utilities.createMasses(networkFile, compoundsFile, infoFile);
			graph = Properties.createGraph(networkFile, compoundsFile, infoFile, massesMap, reversible, compartments, true);
		}
		
		ObjectOutputStream writer = new ObjectOutputStream(new FileOutputStream(outputFile));
		writer.writeObject(graph);
		writer.close();
	}
	
	protected static void checkNetwork(String networkFile, String outputDir, boolean reversible, boolean compartments) throws IOException, ClassNotFoundException {
		String pathwayFile = getPathwaysFile(networkFile);
		String version = getVersion(networkFile);
		String infoFile = outputDir+File.separator+version+".checkNetwork.log";
		
		HashMap<String, int[]> massesMap;
		MetabolicGraph graph;
		
		if (networkFile.endsWith(".xml")) {
			// SBML format
			SBMLDocument sbmlDoc = Utilities.initSBML(networkFile);
			Model sbmlModel = sbmlDoc.getModel();
			massesMap = Utilities.createMasses(sbmlModel, infoFile);
			graph = Properties.createSBMLGraph(networkFile, infoFile, massesMap, reversible, compartments, true);
			
		} else if (networkFile.endsWith(".dat")) {
			// biocyc format
			massesMap = Utilities.createMassesBioCyc(Utilities.getCompoundsFile(networkFile), networkFile, infoFile);
			graph = Properties.createBiocycGraph(networkFile, pathwayFile, infoFile, massesMap, false, compartments, true);
			
		} else {
			// text format
			massesMap = Utilities.createMasses(networkFile, Utilities.getCompoundsFile(networkFile), infoFile);
			graph = Properties.createGraph(networkFile, Utilities.getCompoundsFile(networkFile), infoFile, massesMap, reversible, compartments, true);
			
		}
		
		// read randomized graphs
		String graphPrefix = outputDir+File.separator+version+"-graphs"+File.separator+version;
		ObjectInputStream mbReader = new ObjectInputStream(new FileInputStream(graphPrefix+".massbalance.0"));
		MetabolicGraph mbGraph = (MetabolicGraph)mbReader.readObject();
		mbReader.close();
		ObjectInputStream sReader = new ObjectInputStream(new FileInputStream(graphPrefix+".switch.0"));
		MetabolicGraph sGraph = (MetabolicGraph)sReader.readObject();
		sReader.close();
		
		System.out.println("ORIGINAL: ");
		countReactions(graph, reversible);
		System.out.println();
		System.out.println("MASSBALANCED: ");
		countReactions(mbGraph, reversible);
		System.out.println();
		System.out.println("SWITCH: ");
		countReactions(sGraph, reversible);
	}
	
	/**
	 * Checks certain properties of the original network and the given range of randomized networks,
	 * e.g. the number of isolated compound nodes or the number of reactions containing a cycle.
	 * 
	 * @param networkName
	 * @param inputDir
	 * @param randomization
	 * @param from
	 * @param to
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void checkNetworks(String networkName, String inputDir, String randomization, int from, int to, String outputFile) throws IOException, ClassNotFoundException {
		String originalFile = inputDir+File.separator+networkName+"-graphs"+File.separator+networkName;
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
		ObjectInputStream graphReader = new ObjectInputStream(new FileInputStream(originalFile));
		MetabolicGraph originalGraph = (MetabolicGraph)graphReader.readObject();
		graphReader.close();
		
		ArrayList<Integer> originalComponents = originalGraph.connectedComponents(true);
		writer.write(String.valueOf(originalComponents.size())+"\n");
//		System.out.println("Original cycles: ");
//		countReactionsWithCycles(originalGraph);
		
		for (int i=from; i<=to; i++) {
			String file = originalFile+"."+randomization+"."+String.valueOf(i);
			ObjectInputStream reader = new ObjectInputStream(new FileInputStream(file));
			MetabolicGraph graph = (MetabolicGraph)reader.readObject();
			reader.close();
			writer.write(graph.connectedComponents(false).size()+"\n");
//			System.out.print("cycles: ");
//			countReactionsWithCycles(graph);
		}
		
		writer.close();
	}
	
	/**
	 * Counts the number of isolated nodes in the graph, i.e. nodes with an in-degree and out-degree of 0.
	 * @param graph
	 */
	public static int countIsolatedNodes(MetabolicGraph graph) {
		int isolated = 0;
		for (Vertex v : graph.vertexSet()) {
			if (graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) == 0)
				isolated++;
		}
		return isolated;
	}
	
	/**
	 * Counts the number of reactions which contain a cycle, i.e. a substrate which is also a product.
	 * @param graph
	 */
	public static void countReactionsWithCycles(MetabolicGraph graph) {
		int cycles = 0;
		for (Vertex reaction : graph.getReactions()) {
			if (graph.hasCycle(reaction))
				cycles++;
		}
		System.out.println(cycles);
	}
	
	/**
	 * Determines the mass-balance of the reactions in the given graph and prints
	 * unbalanced reactions.
	 * 
	 * @param graph
	 */
	public static void checkBalance(MetabolicGraph graph) {
		System.out.println("Checking atom balance...");
		
		for (Vertex reaction : graph.getReactions()) {
			boolean balanced = true, hydrogenBalanced = true;
			int[] balance = graph.balance(reaction);
			if (balance == null)
				continue;
			String message = "";
			for (int i=0; i<balance.length; i++) {
				if (balance[i] != 0) {
					if (ELEMENTS[i].equals("H"))
						hydrogenBalanced = false;
					balanced = false;
					message += (Utilities.ELEMENTS[i]+"["+balance[i]+"] ");
				}
			}
			if (!balanced) {
				System.out.println("Reaction "+reaction.getName()+" is unbalanced: "+message);
				for (Vertex p : Graphs.predecessorListOf(graph, reaction))
					System.out.print(DELIMITER_COEFFICIENT_START+graph.getEdgeWeight(graph.getEdge(p,reaction))+DELIMITER_COEFFICIENT_END+" "+p.getName()+" ");
				System.out.print(" != ");
				for (Vertex s : Graphs.successorListOf(graph, reaction))
					System.out.print(DELIMITER_COEFFICIENT_START+graph.getEdgeWeight(graph.getEdge(s,reaction))+DELIMITER_COEFFICIENT_END+" "+s.getName()+" ");
				System.out.println();
			} else if (!hydrogenBalanced)
				throw new RuntimeException("Reaction "+reaction.getName()+" is hydrogen unbalanced (ignore if nofix).");
		}
		
		System.out.println("No hydrogen unbalanced reactions.");
	}
	
	/**
	 * Counts the number of reactions, reversible, cyclic, zero-degree reactions, and the average
	 * and maximum reaction degrees.
	 * 
	 * @param graph
	 * @param print
	 * @throws IOException
	 */
	public static void countReactions(MetabolicGraph graph, boolean print) throws IOException {
		// count reversible/irreversible and transport reactions
		// original graph
		int reversibleReactions = 0, cyclicReactions = 0, zeroDegreeReactions = 0;
		for (Vertex reaction : graph.getReactions()) {
			if (reaction.reversedReaction() != null)
				reversibleReactions++;
			if (graph.inDegreeOf(reaction) == 0 || graph.outDegreeOf(reaction) == 0) {
				if (print && graph.hasCompartments()) 
					System.out.println("Reaction "+reaction.getName()+" has in-degree "+graph.inDegreeOf(reaction)+", out-degree "+graph.outDegreeOf(reaction));
				zeroDegreeReactions++;
			}
			if (graph.isCyclic(reaction, print))
				cyclicReactions++;
		}
		System.out.println("Total reactions: "+graph.getReactions().size()+", reversible: "+reversibleReactions+", cyclic: "+cyclicReactions+", zero-degree reactions: "+zeroDegreeReactions+", average reaction degree: "+graph.averageReactionDegree()+", maximum reaction degree: "+graph.maxReactionDegree());
	}
	
	/**
	 * Prints the indices in the sorted lists of reactions from the given reaction list.
	 * 
	 * @param graph
	 * @param reactions
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void printReactionIndices(MetabolicGraph graph, ArrayList<String> reactions) throws IOException, ClassNotFoundException {
		
		Collection<Vertex> reactionsSet = graph.getReactions();
    	
		// sort the reactions by name (and implicitly by reversedness)
    	ArrayList<Vertex> reactionsList = new ArrayList<Vertex>(reactionsSet);
    	Collections.sort(reactionsList, new VertexComparator());
    	
    	// print the given compound's indices
    	for (String v : reactions) {
			// count reaction index by including reversed reactions
			System.out.print(v+"\t");
			Vertex reaction = graph.getReaction(v);
			if (reaction == null) {
				System.out.println("Reaction not found: "+v);
			} else {
				int i=0;
				for (Vertex r : reactionsList) {
					if (reaction.equals(r)) {
						System.out.print(i);
						break;
					}
					if (r.reversedReaction() != null)
						i++;
					i++;
				}
				System.out.println();
			}
    	}
	}
	
	/**
	 * Prints the indices in the sorted lists of compounds from the given compounds list.
	 * 
	 * @param graph
	 * @param compounds
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void printCompoundIndices(MetabolicGraph graph, ArrayList<String> compounds) throws IOException, ClassNotFoundException {
		
		Collection<Vertex> compoundsSet = graph.getCompounds();
		
		// sort the compounds by name and compartment
		ArrayList<Vertex> compoundsList = new ArrayList<Vertex>(compoundsSet);
    	Collections.sort(compoundsList, new VertexComparator());
    	
    	// print the given compound's indices
    	for (String v : compounds) {
    		
    		if (v.contains(DELIMITER_COMPOUND_COMPARTMENT) || (!graph.hasCompartments() && graph.getCompound(v, null) != null)) {
    			// print compound index
				if (graph.hasCompartments()) {
					if (!v.contains(DELIMITER_COMPOUND_COMPARTMENT))
						throw new IllegalArgumentException("No compartment specified for compound "+v);
					
					String name = v.substring(0, v.indexOf(DELIMITER_COMPOUND_COMPARTMENT));
					String compartment = v.substring(v.indexOf(DELIMITER_COMPOUND_COMPARTMENT)+1);
					Vertex compound = graph.getCompound(name, compartment);
					System.out.println(compoundsList.indexOf(compound));
					
				} else {
					Vertex compound = graph.getCompound(v, null);
					System.out.println(compound+": "+compoundsList.indexOf(compound));
				}
    		} else {
    			System.out.println("Compound not found: "+v);
    		}
    	}
	}
	
	/**
	 * Prints all possible substitutions of the given reaction.
	 * 
	 * @param graph
	 * @param reaction
	 * @param classes
	 */
	protected static void printSubstitutions(MetabolicGraph graph, Vertex reaction, HashMap<ArrayList<Integer>, EquivalenceClass> classes) {
		
		// print original reaction
		printReaction(graph, reaction, true, false);
		
		// print substrate substitutions
		int inDegree = graph.inDegreeOf(reaction);
		DefaultWeightedEdge[] substrateEdges = new DefaultWeightedEdge[inDegree];
		graph.incomingEdgesOf(reaction).toArray(substrateEdges);
		MetabolicGraph.Substitutions substrateSubstitutions = graph.getSubstitutions(reaction, new HashSet<Vertex>(), inDegree, 0, substrateEdges, classes, false, graph.hasCycle(reaction));
		for (int i=0; i<substrateSubstitutions.edgeIndices.size(); i++) {
			System.out.println(graph.getEdgeSource(substrateEdges[substrateSubstitutions.edgeIndices.get(i)]).getName() + " ==> " + substrateSubstitutions.substitutes.get(i).getName() + " [edge index " +substrateSubstitutions.edgeIndices.get(i)+"]");
		}
		for (int i=0; i<substrateSubstitutions.edgePairIndices.size()-1; i+=2) {
			System.out.println(graph.getEdgeSource(substrateEdges[substrateSubstitutions.edgePairIndices.get(i)]).getName() + " + " + graph.getEdgeSource(substrateEdges[substrateSubstitutions.edgePairIndices.get(i+1)]).getName()
			+ " ==> " + substrateSubstitutions.substitutePairs.get(i).getName() + " + " + substrateSubstitutions.substitutePairs.get(i+1).getName() + " [edge indices " +substrateSubstitutions.edgePairIndices.get(i)+","+substrateSubstitutions.edgePairIndices.get(i+1)+"]");
		}
		
		// print product substitutions
		int outDegree = graph.outDegreeOf(reaction);
		DefaultWeightedEdge[] productEdges = new DefaultWeightedEdge[outDegree];
		graph.outgoingEdgesOf(reaction).toArray(productEdges);
		MetabolicGraph.Substitutions productSubstitutions = graph.getSubstitutions(reaction, new HashSet<Vertex>(), outDegree, 1, productEdges, classes, false, graph.hasCycle(reaction));
		for (int i=0; i<productSubstitutions.edgeIndices.size(); i++) {
			System.out.println(graph.getEdgeTarget(productEdges[productSubstitutions.edgeIndices.get(i)]).getName() + " ==> " + productSubstitutions.substitutes.get(i).getName() + " [edge index " +productSubstitutions.edgeIndices.get(i)+"]");
		}
		for (int i=0; i<productSubstitutions.edgePairIndices.size()-1; i+=2) {
			System.out.println(graph.getEdgeTarget(productEdges[productSubstitutions.edgePairIndices.get(i)]).getName() + " + " + graph.getEdgeTarget(productEdges[productSubstitutions.edgePairIndices.get(i+1)]).getName()
			+ " ==> " + productSubstitutions.substitutePairs.get(i).getName() + " + " + productSubstitutions.substitutePairs.get(i+1).getName() + " [edge indices " +productSubstitutions.edgePairIndices.get(i)+","+productSubstitutions.edgePairIndices.get(i+1)+"]");
		}
	}
	
	/**
	 * Finds and prints all transport reactions in which the given compounds are involved.
	 * If exclusive==true, than only transporters which exclusively carry the given compound
	 * (in- and out-degree is 1) are printed.
	 * 
	 * @param graph
	 * @param compoundNames
	 * @param exclusive
	 */
	public static void findTransporters(MetabolicGraph graph, ArrayList<String> compoundNames, boolean exclusive) throws IOException, ClassNotFoundException {
		
		Collection<Vertex> reactionsSet = graph.getReactions();
    	
		// sort the reactions by name (and implicitly by reversedness)
    	ArrayList<Vertex> reactionsList = new ArrayList<Vertex>(reactionsSet);
    	Collections.sort(reactionsList, new VertexComparator());
    	
		// find all given compounds
		for (String name : compoundNames) {
			
			System.out.println("Finding "+(exclusive?"exclusive ":"")+"transporters for: "+name);
			
			// search all pairs of compartments
			Set<String> compartments = graph.getCompartments().keySet();
			ArrayList<String> compartmentsList = new ArrayList<String>(compartments);
			Collections.sort(compartmentsList);
			
			for (int i=0; i<compartmentsList.size()-1; i++) {
				String compartment = compartmentsList.get(i); 
					
				for (int j=i+1; j<compartmentsList.size(); j++) {
					String compartment2 = compartmentsList.get(j);
					
					Vertex compound = graph.getCompound(name, compartment);
					if (compound != null) {
						HashSet<Vertex> visited = new HashSet<Vertex>();
						// search the predecessor reactions
						for (Vertex reaction : Graphs.predecessorListOf(graph, compound)) {
							if (!reaction.isReversed() && !visited.contains(reaction) && (!exclusive || graph.inDegreeOf(reaction)==1 && graph.outDegreeOf(reaction)==1)) {
								visited.add(reaction);
								// check for transporter: any other compartment must contain the same compound name 
								if (Graphs.predecessorListOf(graph, reaction).contains(graph.getCompound(name, compartment2))) {
									printReaction(graph, reaction, true, true);
									ArrayList<String> reactionList = new ArrayList<String>(1);
									reactionList.add(reaction.getName());
									printReactionIndices(graph, reactionList);
								}
							}
						}
						// search the successor reactions
						for (Vertex reaction : Graphs.successorListOf(graph, compound)) {
							if (!reaction.isReversed() && !visited.contains(reaction) && (!exclusive || graph.inDegreeOf(reaction)==1 && graph.outDegreeOf(reaction)==1)) {
								visited.add(reaction);
								if (Graphs.successorListOf(graph, reaction).contains(graph.getCompound(name, compartment2))) {
									printReaction(graph, reaction, true, true);
									ArrayList<String> reactionList = new ArrayList<String>(1);
									reactionList.add(reaction.getName());
									printReactionIndices(graph, reactionList);
								}
							}
						}
					}
				}
			}
		}
		
	}
	
	/***
	 * Prints the sorted stoichiometric matrix indices of all compounds in a compartment
	 * starting with "ext", "ex", or "e".
	 * 
	 * @param graph
	 */
	public static void printExternalCompoundIndices(MetabolicGraph graph) {
		String external = null;
		
		if (!graph.hasCompartments())
			throw new IllegalArgumentException("Graph has no compartments.");
		
		Collection<Vertex> compoundsSet = graph.getCompounds();
		
		// sort the compounds by name and compartment
		ArrayList<Vertex> compoundsList = new ArrayList<Vertex>(compoundsSet);
    	Collections.sort(compoundsList, new VertexComparator());
		
		// find the external compartment string
		for (String compartment : graph.getCompartments().keySet()) {
			if (compartment.startsWith("ext"))
				external = compartment;
			else if (compartment.startsWith("ex"))
				external = compartment;
			else if (compartment.startsWith("e"))
				external = compartment;
		}
		
		System.out.println("Compounds in compartment "+external+":");
		int i=0;
		for (Vertex compound : compoundsList) {
			if (compound.getCompartment().equals(external))
				System.out.println(i);
			i++;
		}
	}
			
	/**
	 * Print the summary statistics for the given graph:
	 * 
	 * - number of compounds and those containing only CHNOPS
	 * - number of reactions, unbalanced and unannotated reactions
	 * - reversibility of the graph, number of vertices and edges
	 * - sizes of connected and strongly connected components
	 * 
	 * @param graph
	 */
	public static void summary(MetabolicGraph graph) {
		int chnops = 0, unbalanced = 0, unannotated = 0;
		
		for (Vertex compound : graph.getCompounds()) {
			if (compound.getMass() != null)
				chnops++;
		}
		for (Vertex reaction : graph.getReactions()) {
			int[] balance = graph.balance(reaction);
			if (balance == null)
				unannotated++;
			else {
				for (int i=0; i<balance.length; i++) {
					if (balance[i] != 0) {
						unbalanced++;
						break;
					}
				}
			}
		}
		
		System.out.println("Compounds: "+graph.getCompounds().size()+", CHNOPS: "+chnops);
		System.out.println("Reactions: "+graph.getReactions().size()+", unbalanced: "+unbalanced+", unannotated: "+unannotated);
		System.out.println("Graph is "+(graph.isReversible()?"":"not ")+"reversible, "+graph.vertexSet().size()+" vertices, "+graph.edgeSet().size()+" edges.");
		
		// print the (strongly) connected components
		graph.connectedComponents(true);
		StrongConnectivityInspector<Vertex, DefaultWeightedEdge> connectivity = new StrongConnectivityInspector<Vertex, DefaultWeightedEdge>(graph);
		List<Set<Vertex>> components = connectivity.stronglyConnectedSets();
		graph.cycleCount(components, 2, false, true);
	}
	
	/**
	 * Prints some example reactions of the graph with in- and out-degree smaller than 4.
	 * 
	 * @param networkFile
	 * @throws IOException
	 */
	public static void printExampleReactions(String networkFile) throws IOException {
		// print some simple example reactions
		String version = getVersion(networkFile);
		String compoundsFile = getCompoundsFile(networkFile);
		String pathwayFile = getPathwaysFile(networkFile);
		String infoFile = networkFile+".printExampleReactions.log";
		boolean reversible = (version.equals("athaliana") || version.equals("ymn2_0"));
		boolean compartments = version.equals("barley");
		
		HashMap<String, int[]> massesMap;
		
		
		createMassesBioCyc(compoundsFile, networkFile, infoFile);
		MetabolicGraph graph;

		if (networkFile.endsWith(".xml")) {
			// SBML format
			SBMLDocument sbmlDoc = Utilities.initSBML(networkFile);
			Model sbmlModel = sbmlDoc.getModel();
			massesMap = Utilities.createMasses(sbmlModel, infoFile);
			graph = Properties.createSBMLGraph(networkFile, infoFile, massesMap, reversible, compartments, true);
			
		} else if (networkFile.endsWith(".dat")) {
			// biocyc format
			massesMap = Utilities.createMassesBioCyc(Utilities.getCompoundsFile(networkFile), networkFile, infoFile);
			graph = Properties.createBiocycGraph(networkFile, pathwayFile, infoFile, massesMap, false, compartments, true);
			
		} else {
			// text format
			massesMap = Utilities.createMasses(networkFile, Utilities.getCompoundsFile(networkFile), infoFile);
			graph = Properties.createGraph(networkFile, Utilities.getCompoundsFile(networkFile), infoFile, massesMap, reversible, compartments, true);
			
		}
		
		for (Vertex reaction : graph.getReactions()) {
			int atoms = 0;
			reactions:
			if (graph.inDegreeOf(reaction)+graph.outDegreeOf(reaction) <= 3) {
				for (Vertex compound : Graphs.predecessorListOf(graph, reaction)) {
					if (compound.getMass() == null)
						break reactions;
					for (int i=0; i<ELEMENTS.length; i++)
						atoms += compound.getMass()[i];
				}
				for (Vertex compound : Graphs.successorListOf(graph, reaction)) {
					if (compound.getMass() == null)
						break reactions; 
					for (int i=0; i<ELEMENTS.length; i++)
						atoms += compound.getMass()[i];
				}
			}
			if (atoms >0 && atoms <= 26)
				System.out.println(reaction.getName()+" ("+atoms+" atoms)");
			
		}
	}
	
	/**
	 * Prints the stoichiometric coefficients of the given graph,
	 * i.e., its edge weights.
	 * 
	 * @param graph
	 */
	public static void printCoefficients(MetabolicGraph graph) {
		for (DefaultWeightedEdge edge : graph.edgeSet())
			System.out.println(graph.getEdgeWeight(edge));
	}
	
	/**
	 * Print all compounds with masses and the largest integer factorization.
	 * 
	 * @param compoundsFile
	 * @param networkFile
	 * @param infoFile
	 * @throws IOException
	 */
	public static void printMasses(String compoundsFile, String networkFile, String infoFile) throws IOException {
		int[] primes = {2, 3, 5, 7, 11, 13};
		BigDecimal max = BigDecimal.ZERO;
		String maxFactorizationCompound = "", maxAtomsCompound = "";
		int maxAtoms = 0;
		
		HashMap<String, int[]> massesMap = createMassesBioCyc(compoundsFile, networkFile, infoFile);
		
		for (String compound : massesMap.keySet()) {
			System.out.print(compound+"\t");
			int[] mass = massesMap.get(compound);
			if (mass != null) {
				BigDecimal factorization = BigDecimal.ONE; 
				for (int i=0; i<ELEMENTS.length; i++) {
					System.out.print(ELEMENTS[i]+mass[i]+" ");
					factorization = factorization.multiply(new BigDecimal(Math.pow(primes[i], mass[i])));
					if (mass[i] > maxAtoms) {
						maxAtoms = mass[i];
						maxAtomsCompound = compound + " [" + ELEMENTS[i] + "]";
					}
				}
				if (factorization.compareTo(max) == 1) {
					max = factorization;
					maxFactorizationCompound = compound;
				}
			}
			System.out.println();
		}
		
		System.out.println("Largest factorization: "+maxFactorizationCompound+": "+max.stripTrailingZeros());
		System.out.println("Largest number of individual atoms: "+maxAtomsCompound+": "+maxAtoms);
	}
	
	/**
	 * Print all reactions in the given graph that have a coefficient at least as large as min.
	 * Reactions with more than one large coefficient are repeatedly printed.
	 * 
	 * @param graph
	 * @param min
	 */
	public static void printLargeCoefficients(MetabolicGraph graph, int min) {
		for (DefaultWeightedEdge edge : graph.edgeSet()) {
			double weight = graph.getEdgeWeight(edge);
			if (weight >= min) {
				Vertex reaction = (graph.getEdgeSource(edge).getType() == Vertex.REACTION ? graph.getEdgeSource(edge) : graph.getEdgeTarget(edge));
				printReaction(graph, reaction, true, true);
			}
		}
	}
	
	/**
	 * Searches for the largest coefficient in the network and prints the 
	 * corresponding reaction.
	 * 
	 * @param graph
	 */
	public static void printMaxCoefficient(MetabolicGraph graph) {
		DefaultWeightedEdge maxEdge = null;
		double maxWeight = -1;
		
		for (DefaultWeightedEdge edge : graph.edgeSet()) {
			double weight = graph.getEdgeWeight(edge);
			if (maxWeight < weight) {
				maxWeight = weight;
				maxEdge = edge;
			}
		}
		
		printReaction(graph, (graph.getEdgeSource(maxEdge).getType()==Vertex.REACTION ? graph.getEdgeSource(maxEdge) : graph.getEdgeTarget(maxEdge)), true, true);
	}
	
	/**
	 * Searches for the largest deltaGr in the network and prints the 
	 * corresponding reaction. If reversible==true, then the deltaGr of
	 * reversible reactions is 0, otherwise, the maximum is determined
	 * for each direction independently. If excludeTransporters, transport
	 * reactions as well as import/export reactions are skipped.
	 * 
	 * @param graph
	 */
	public static void printMaxDeltaGr(MetabolicGraph graph, boolean reversible, boolean excludeTransprorters) {
		double maxDeltaGr = -Double.MAX_VALUE;
		Vertex maxReaction = null;
		
		for (Vertex reaction : graph.getReactions()) {
			
			// skip transporters and import/export reactions
			if (excludeTransprorters &&
					(graph.inDegreeOf(reaction) == 0 || graph.outDegreeOf(reaction) == 0 
					 || graph.getCompartments(reaction).keySet().size() > 1)) {
				continue;
			}
			
			Double deltaGr = graph.getDeltaGr(reaction, reversible);
			if (deltaGr != null && deltaGr > maxDeltaGr) {
				maxDeltaGr = deltaGr;
				maxReaction = reaction;
			}
			// if reversible==false, check the reversed direction independently
			if (!reversible && reaction.reversedReaction() != null) {
				deltaGr = graph.getDeltaGr(reaction.reversedReaction(), false);
				if (deltaGr != null && deltaGr > maxDeltaGr) {
					maxDeltaGr = deltaGr;
					maxReaction = reaction.reversedReaction();
				}
			}

		}
		
		System.out.println("Largest deltaGr: "+maxDeltaGr);
		printReaction(graph, maxReaction, true, true);
	}
	
	public static void printMass(MetabolicGraph graph, String compoundString) {
		Vertex compound;
		if (graph.hasCompartments()) {
			if (!compoundString.contains(DELIMITER_COMPOUND_COMPARTMENT)) {
				System.out.println("No compartment specified.");
				return;
			}
			String name = compoundString.substring(0, compoundString.indexOf(DELIMITER_COMPOUND_COMPARTMENT));
			String compartment = compoundString.substring(compoundString.indexOf(DELIMITER_COMPOUND_COMPARTMENT)+1);
			compound = graph.getCompound(name, compartment);
		} else {
			compound = graph.getCompound(compoundString, null);
		}
		
		if (compound == null) {
			System.out.println("No such compound.");
			return;
		}
		
		System.out.print(compound.getName()+(graph.hasCompartments()?DELIMITER_COMPOUND_COMPARTMENT+compound.getCompartment():"")+" (");
		int[] mass = compound.getMass();
		if (mass != null) {
			for (int i=0; i<ELEMENTS.length; i++)
				if (mass[i] != 0)
					System.out.print(ELEMENTS[i]+mass[i]);
		} else {
			System.out.print("null");
		}
		System.out.print(")");
		
	}
	
	/**
	 * Prints all reaction equations of the given graph.
	 * 
	 * @param graph
	 */
	public static void printReactions(MetabolicGraph graph) {
		for (Vertex reaction : graph.getReactions())
			printReaction(graph, reaction, false, false);
	}
	
	/**
	 * Prints the equation of the given reaction.
	 *  
	 * @param graph
	 * @param reaction
	 */
	public static void printReaction(MetabolicGraph graph, Vertex reaction, boolean printMass, boolean printDeltaG) {
		
		if (reaction == null) {
			System.out.println("No such reaction.");
			return;
		}
		
		System.out.print(reaction.getName()+(printDeltaG ? " ("+graph.getDeltaGr(reaction, false)+")" : "")+"\t");
		for (Iterator<DefaultWeightedEdge> edgeIt = graph.incomingEdgesOf(reaction).iterator(); edgeIt.hasNext();) {
			// print the substrate
			DefaultWeightedEdge edge = edgeIt.next();
			Vertex substrate = graph.getEdgeSource(edge);
			double edgeWeight = graph.getEdgeWeight(edge);
			System.out.print((edgeWeight != 1 ? DELIMITER_COEFFICIENT_START+String.valueOf(edgeWeight)+DELIMITER_COEFFICIENT_END+" " : "")+substrate.getName()+(graph.hasCompartments() ? DELIMITER_COMPOUND_COMPARTMENT+substrate.getCompartment() : "")+(printDeltaG ? " ["+String.valueOf(substrate.getDeltaGf())+"] " : ""));
			// print the mass
			if (printMass) {
				System.out.print("(");
				int[] mass = substrate.getMass();
				if (mass != null) {
					for (int i=0; i<ELEMENTS.length; i++)
						if (mass[i] != 0)
							System.out.print(ELEMENTS[i]+mass[i]);
				} else {
					System.out.print("null");
				}
				System.out.print(")");
			}
			
			if (edgeIt.hasNext())
				System.out.print(" + ");
		}
		if (reaction.reversedReaction() != null)
			System.out.print(" = ");
		else
			System.out.print(" > ");
		for (Iterator<DefaultWeightedEdge> edgeIt = graph.outgoingEdgesOf(reaction).iterator(); edgeIt.hasNext();) {
			// print the product
			DefaultWeightedEdge edge = edgeIt.next();
			Vertex product = graph.getEdgeTarget(edge);
			double edgeWeight = graph.getEdgeWeight(edge);
			System.out.print((edgeWeight != 1 ? DELIMITER_COEFFICIENT_START+String.valueOf(edgeWeight)+DELIMITER_COEFFICIENT_END+" " : "")+product.getName()+(graph.hasCompartments() ? DELIMITER_COMPOUND_COMPARTMENT+product.getCompartment() : "")+(printDeltaG ? " ["+String.valueOf(product.getDeltaGf())+"] " : ""));
			// print the mass
			if (printMass) {
				System.out.print("(");
				int[] mass = product.getMass();
				if (mass != null) {
					for (int i=0; i<ELEMENTS.length; i++)
						if (mass[i] != 0)
							System.out.print(ELEMENTS[i]+mass[i]);
				} else {
					System.out.print("null");
				}
				System.out.print(")");
			}
			
			if (edgeIt.hasNext())
				System.out.print(" + ");
		}
		System.out.println();
	}
	
	/**
	 * Prints the reaction equation from the given graph.
	 * 
	 * @param graph
	 * @param reaction
	 */
	public static void printReaction(MetabolicGraph graph, String reaction, boolean printMass, boolean printDeltaG) {
		printReaction(graph, graph.getReaction(reaction), printMass, printDeltaG);
	}
	
	/**
	 * Prints the reaction equations from the given graph.
	 * 
	 * @param graph
	 * @param reactions
	 */
	public static void printReactions(MetabolicGraph graph, ArrayList<String> reactions) {
		for (String reaction : reactions)
			printReaction(graph, graph.getReaction(reaction), true, true);
	}
	
	/**
	 * Prints all reactions with in-degree+out-degree < 3 in the graph.
	 * 
	 * @param graph
	 */
	public static void printSmallDegreeReactions(MetabolicGraph graph) {
		for (Vertex reaction : graph.getReactions()) {
			if (graph.inDegreeOf(reaction)+graph.outDegreeOf(reaction) < 3)
				printReaction(graph, reaction, true, true);
		}
	}
	
	/**
	 * Prints the number of compounds given by number with largest degree.
	 * 
	 * @param graph
	 * @param threshold
	 */
	public static void printHubCompounds(final MetabolicGraph graph, int number) {
		ArrayList<Vertex> compounds = new ArrayList<Vertex>(graph.getCompounds());
    	Collections.sort(compounds, new Comparator<Vertex>() {
    		public int compare(Vertex c1, Vertex c2) {
    			return (new Integer(graph.inDegreeOf(c2)+graph.outDegreeOf(c2)).compareTo(new Integer(graph.inDegreeOf(c1)+graph.outDegreeOf(c1))));
    		}});
    	
    	for (int i=0; i<number; i++) {
    		Vertex hub = compounds.get(i);
    		System.out.println(hub.getName()+(graph.hasCompartments()?Utilities.DELIMITER_COMPOUND_COMPARTMENT+hub.getCompartment():"")+"\t"+String.valueOf(graph.inDegreeOf(hub)+graph.outDegreeOf(hub)));
    	}
	}
	
	
	/**
	 * Prints a particular type of pathway, if present:
	 * 
	 * A -> B+D
	 * B+D -> C
	 * 
	 * where D is a cofactor (has a large degree), and B not.
	 * 
	 * @param graph
	 */
	public static void printPathway(MetabolicGraph graph) {
		HashSet<Vertex> hubs = new HashSet<Vertex>();
		
		// start from the hubs
		for (Vertex d : graph.getCompounds()) {
			if (graph.inDegreeOf(d)+graph.outDegreeOf(d) >= 500)
				hubs.add(d);
		}
		
		// from each hub, search for the pathway
		for (Vertex d : hubs) {
			for (Vertex r1 : Graphs.predecessorListOf(graph, d)) {
				
				int degree1 = graph.inDegreeOf(r1)+graph.outDegreeOf(r1);
//				if (graph.inDegreeOf(r1)==1 && graph.outDegreeOf(r1)==2) {
				if (degree1<=4) {
					
					for (Vertex b : Graphs.successorListOf(graph, r1)) {
						if (!d.equals(b) && graph.inDegreeOf(b)+graph.outDegreeOf(b) >= 30 && graph.inDegreeOf(b)+graph.outDegreeOf(b) <= 400) {
							for (Vertex r2 : Graphs.successorListOf(graph, b)) {
								if (r1.reversedReaction()==null || r2.reversedReaction()==null || !r1.reversedReaction().equals(r2)) {
									int degree2 = graph.inDegreeOf(r2)+graph.outDegreeOf(r2);
									if (Graphs.predecessorListOf(graph, r2).contains(d) && degree1+degree2 <= 7) {//graph.inDegreeOf(r2)==2 && graph.outDegreeOf(r2) == 1) {
										printReaction(graph, r1, false, false);
										printReaction(graph, r2, false, false);
										System.out.println("----------------------------------------------");
										
									}
								}
							}
						}
					}
				}
			}
			
		}
	}
	
	/**
	 * Prints reactions containing TCA cycle intermediaries from ecoli, satisfying
	 * additional conditions.
	 *  
	 * @param graph
	 */
	public static void printTCACycleReactions(MetabolicGraph graph) {
		String[] tca = new String[]{"accoa", "cit", "icit", "akg", "succoa", "succ", "fum", "mal-L", "oaa"};
//		String[] tca = new String[]{"cit", "icit", "mal-L"};
		ArrayList<String> tcaIntermediaries = new ArrayList<String>(Arrays.asList(tca));
		ArrayList<Vertex> reactions = new ArrayList<Vertex>();
		boolean hasUnary = false;
		boolean hasHighDeltaGr = false;
		boolean hasLowDeltaGr = false;
		
		// collect reactions of interest
		for (Vertex reaction : graph.getReactions()) {
			boolean valid = false;
			ArrayList<String> substrates = new ArrayList<String>(graph.inDegreeOf(reaction));
			
			// skip reactions with in- or out-degree bigger than 3
			if (graph.inDegreeOf(reaction) > 3 || graph.outDegreeOf(reaction) > 3 )
				continue;
			
			// check for substrate intermediaries
			for (Vertex substrate : Graphs.predecessorListOf(graph, reaction)) {
				if (tcaIntermediaries.contains(substrate.getName())) {
					substrates.add(substrate.getName());
					valid = true;
					break;
				}
			}
			// check for product intermediaries			
//			if (valid) {
//				valid = false;
				for (Vertex product : Graphs.successorListOf(graph, reaction)) {
					// skip reactions with circular intermediaries
					if (substrates.contains(product.getName())) {
						valid = false;
						break;
					}
					if (tcaIntermediaries.contains(product.getName()))
						valid = true;
				}
//			}
			
			if (valid) {
				if (graph.inDegreeOf(reaction) == 1 || graph.outDegreeOf(reaction) == 1)
					hasUnary = true;
				Double deltaGr = graph.getDeltaGr(reaction, true);
				if (deltaGr != null && deltaGr < -100)
					hasLowDeltaGr = true;
				if (deltaGr != null && deltaGr > 50)
					hasHighDeltaGr = true;
				reactions.add(reaction);
			}
		}
		
		// print interesting reactions
		if (hasUnary && hasLowDeltaGr && hasHighDeltaGr) {
			for (Vertex reaction : reactions)
				printReaction(graph, reaction, true, true);
		}
	}
	
	//////////////////////////
	// Mass parsing methods.//
	//////////////////////////
	
	protected static HashMap<String, int[]> createMasses(String networkFile, String compoundsFile, String infoFile) throws IOException {
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader networkReader = new BufferedReader(new FileReader(networkFile));
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(infoFile));
		HashMap<String, int[]> massesMap = new HashMap<String, int[]>((int)(DEFAULT_NUMBER_OF_COMPOUNDS/0.75f+1), 0.75f);
		boolean hasHydrogen = false;
		boolean[] hasPhosphate = {false, false, false, false};
		String line;
		int annotated = 0;
		
		String networkVersion = getVersion(networkFile);
		infoWriter.write("# "+networkVersion+"\n");
		
		// collect the compounds actually in the network
		HashSet<String> existing = new HashSet<String>((int)(Utilities.DEFAULT_NUMBER_OF_COMPOUNDS/0.75f+1));
		while ((line = networkReader.readLine()) != null) {
			
			if (line.startsWith("#"))
				continue;
			
			// extract reaction equation string
			int start = line.indexOf("\t")+1;
			int end = (line.indexOf("\t", start) >= 0 ? line.indexOf("\t", start) : line.length());
			String equation = line.substring(start, end);
			
			// collect the tokens from the equation
			existing.addAll(parseCompounds(equation));
		}
		
		// collect the masses from the compounds file
		while ((line = compoundsReader.readLine()) != null) {
			
			if (line.startsWith("#"))
				continue;
			
			if (line.indexOf("\t") < 1)
				throw new RuntimeException("Invalid line in "+compoundsFile+" (no TAB space): "+line);
			
			// extract compound name and mass
			String compound = line.substring(0, line.indexOf("\t")).trim();
			
			if (existing.contains(compound)) {
				int end = line.indexOf("\t", line.indexOf("\t")+1) > -1 ? line.indexOf("\t", line.indexOf("\t")+1) : line.length();
				int[] mass;
				if (line.substring(line.indexOf("\t")+1, end).contains("InChI="))
					mass = parseInChIString(compound, line.substring(line.indexOf("\t")+1, end), infoWriter);
				else
					mass = parseMass(line.substring(line.indexOf("\t")+1, end), infoWriter);
				if (mass != null) {
					// determine hydrogen and phosphate forms
					if (Properties.isHydrogen(mass))
						hasHydrogen = true;
					int type = Properties.isPhosphate(mass);
					if (type >= 0 && type <= 3)
						hasPhosphate[type] = true;
					if (massesMap.put(compound, mass) == null)
						annotated++;
				} else {
					infoWriter.write("Invalid mass for compound "+compound+"\n");
				}
			}
		}
		
		// add missing hydrogen and phosphate forms
		if (!hasHydrogen)
			massesMap.put(Utilities.HYDROGEN_NAME, Utilities.HYDROGEN_MASS);
		for (int i=0; i<hasPhosphate.length; i++) {
			if (!hasPhosphate[i])
				massesMap.put(Utilities.PHOSPHATE_NAMES[i], Utilities.PHOSPHATE_FORMS[i]);
		}
		
		if (DEFAULT_NUMBER_OF_COMPOUNDS < massesMap.size())
			System.out.println("Initial masses hash capacity "+DEFAULT_NUMBER_OF_COMPOUNDS+" smaller than final map size ("+massesMap.size()+"). Increase DEFAULT_NUMBER_OF_COMPOUNDS to at least "+massesMap.size()+" for better performance.");
		
		float ratioAnnotated = (float)annotated*100/existing.size();
		String message = existing.size()+" compounds in the original network, "+annotated+" have annotated mass with ";
		for (int i=0; i<ELEMENTS.length; i++)
			message += ELEMENTS[i]+(i<ELEMENTS.length-1?",":"");
		message += " ("+ratioAnnotated+"%).";
		System.out.println(message);
		infoWriter.write(message+"\n");
		
		compoundsReader.close();
		networkReader.close();
		infoWriter.close();
		
		return massesMap;
	}
	
	protected static ArrayList<String> parseCompounds(String equation) {
		boolean debug = false;
		int sideDelimiter, delimiterLength;
		
		if (!equation.contains(DELIMITER_EQUALS) && !equation.contains(DELIMITER_FORWARD)) {
			System.out.println("Missing or incorrect side delimiter in reaction equation: "+equation);
			return new ArrayList<String>(0);
		}
		
		ArrayList<String> compounds = new ArrayList<String>();
		
		// remove the stoichiometric coefficients
		equation = equation.replaceAll("\\"+DELIMITER_COEFFICIENT_START+"[^\\"+DELIMITER_COEFFICIENT_END+"]*"+"\\"+DELIMITER_COEFFICIENT_END, "");
		// remove the reaction compartment delimiter
		equation = equation.replaceAll("\\"+DELIMITER_REACTION_COMPARTMENT_START+"[^\\"+DELIMITER_REACTION_COMPARTMENT_END+"]*"+"\\"+DELIMITER_REACTION_COMPARTMENT_END, "");
		
		// extract the position and length of the reaction side delimiter
		if (equation.contains(DELIMITER_EQUALS)) {
			sideDelimiter = equation.indexOf(DELIMITER_EQUALS);
			delimiterLength = DELIMITER_EQUALS.length();
		} else {
			sideDelimiter = equation.indexOf(DELIMITER_FORWARD);
			delimiterLength = DELIMITER_FORWARD.length();
		}
		
		// split the equation string by the side delimiter
		String leftSide = equation.substring(0, sideDelimiter).trim();
		String rightSide = equation.substring(sideDelimiter+delimiterLength).trim();
		
		// extract the compound names
		String[] substrates = leftSide.split("\\Q"+DELIMITER_PLUS+"\\E");
		String[] products = rightSide.split("\\Q"+DELIMITER_PLUS+"\\E");
		
		// add the trimmed compound names without compartment delimiters
		for (int i=0; i<substrates.length; i++)
			compounds.add(substrates[i].trim().replaceFirst("\\"+DELIMITER_COMPOUND_COMPARTMENT+".*", ""));
		for (int i=0; i<products.length; i++)
			compounds.add(products[i].trim().replaceFirst("\\"+DELIMITER_COMPOUND_COMPARTMENT+".*", ""));
		
		if (debug) {
			System.out.print(equation+":\t");
			for (String compound : compounds)
				System.out.print(compound+", ");
			System.out.println();
		}
		
		return compounds;
	}
	
	/**
	 * Reads the mass vectors from the given SBML file into a hash map containing the
	 * species id as String key and the mass vector as int[] value. The SBML Level 2
	 * recommendation is to annotate the mass as InChI elements of the second level
	 * child of an rdf:Description element (see <a href="http://sbml.org/Community/Wiki/About_annotations_in_Level_2" target="_blank">).
	 * However, Herrgard et al., Nat Biotechnol 26(10), 1155-1160 Oct (2008) use an
	 * InChI element inside the annotation element of the corresponding SpeciesType.
	 * In version 4 of YeastNet (http://www.comp-sys-bio.org/yeastnet/), they add the
	 * mass formula as part of the notes body. All three types of annotation, as well
	 * as combinations thereof, are supported by this method. 
	 * 
	 * @param sbmModel
	 * @param infoFile
	 * @return
	 * @throws IOException
	 */
	protected static HashMap<String, int[]> createMasses(Model sbmlModel, String infoFile) throws IOException {
		boolean useTypes = false;
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(infoFile));
		HashMap<String, int[]> massesMap = new HashMap<String, int[]>();
		boolean hasHydrogen = false;
		boolean[] hasPhosphate = {false, false, false, false};
		int annotated = 0;
		ListOfSpeciesTypes speciesTypes = null;
		
		HashSet<SBase> existing = new HashSet<SBase>(Utilities.DEFAULT_NUMBER_OF_COMPOUNDS);
		ListOfReactions reactions = sbmlModel.getListOfReactions();
		ListOfSpecies species = sbmlModel.getListOfSpecies();
		
		// The SpeciesType object class is only available in SBML Level 2 Versions 2-4. It is not available in Level 1 nor Level 3.
		if (sbmlModel.getLevel() == 2 && sbmlModel.getVersion() >= 2 && sbmlModel.getVersion() <= 4) {
			speciesTypes = sbmlModel.getListOfSpeciesTypes();
			// check if we have SpeciesTypes as in YeastNet (Herrgard et al., Nat Biotechnol 26(10), 1155-1160 Oct (2008))
			if (speciesTypes.size() > 0)
				useTypes = true;
		}
		
		// collect the compounds actually in the network
		for (int i=0; i<reactions.size(); i++) {
			Reaction reaction = reactions.get(i);
			// collect the substrate and product names
			ListOfSpeciesReferences substrates = reaction.getListOfReactants();
			ListOfSpeciesReferences products = reaction.getListOfProducts();
			
			for (int j=0; j<substrates.size()+products.size(); j++) {
				// add the compound and type to the set of existing compounds/types
				String compoundSpecies = (j<substrates.size() ? substrates.get(j).getSpecies() : products.get(j-substrates.size()).getSpecies());
				Species compound = species.get(compoundSpecies);
				if (compound != null) {
					if (useTypes) {
						SpeciesType compoundType = speciesTypes.get(compound.getSpeciesType());
						if (compoundType != null)
							existing.add(compoundType);
					} else {
						existing.add(compound);
					}
				} else {
					infoWriter.write("Reaction "+reaction.getName()+" contains unannotated species "+compoundSpecies+".\n");
				}
			}
		}
		
		// get the sum formulas of the existing compound species
		for (SBase sBase : existing) {
			
			int[] mass = null;
			String InChI = null, compoundName = null, compoundId = null;
			
			if (useTypes) {
				// YeastNet format: search the corresponding SpeciesType for an "inchi" annotation
				SpeciesType compoundType = (SpeciesType)sBase;
				
				if (compoundType != null && existing.contains(compoundType)) {
					XMLNode notes = compoundType.getNotes();
					XMLNode annotation;
					
					// parse inchi from the notes body subelement (YeastNet 4.0) 
					if (notes != null) {
						for (int j=0; j<notes.getNumChildren(); j++) {
							XMLNode child = notes.getChild(j);
							for (int k=0; k<child.getNumChildren(); k++) {
								XMLNode grandchild = child.getChild(k);
								String formula = grandchild.getChild(0).getCharacters();
								if (formula.startsWith("formula:")) {
									InChI = formula.substring(8,formula.length());
									break;
								}
							}
							if (InChI != null)
								break;
						}
						
					} else if ((annotation = compoundType.getAnnotation()) != null) {
						// parse inchi from the annotation subelement (YeastNet 2.0)
						for (int j=0; j<annotation.getNumChildren(); j++) {
							XMLNode child = annotation.getChild(j);
							if (child.getName().equals("inchi")) {
								InChI = child.getChild(0).getCharacters();
								break;
							}
						}
					}
					// store the name and id of the SpeciesType
					compoundId = compoundType.getId();
					compoundName = compoundType.getName();
				}
				
			} else {
				
				XMLNode annotation, notes;
				Species compound = (Species)sBase;
				
				if ((notes = compound.getNotes()) != null) {
					// parse mass from "html" subelements of notes (some BiGG SBML versions)
					for (int j=0; j<notes.getNumChildren(); j++) {
						XMLNode child = notes.getChild(j);
						if (child.getName().equals("html") || child.getName().equals("body")) {
							InChI = child.getChild(0).getChild(0).getCharacters();
							if (InChI.startsWith("FORMULA: "))
								InChI = InChI.substring(9);
							break;
						}
					}
					// store the name and id of the Species
					compoundId = compound.getId();
					compoundName = compound.getName();
					
				} else if ((annotation = compound.getAnnotation()) != null) {
					// parse mass from "inchi" subelements of annotation (other SBML versions)
					for (int j=0; j<annotation.getNumChildren(); j++) {
						XMLNode child = annotation.getChild(j);
						if (child.getName().equals("inchi")) {
							InChI = child.getChild(0).getCharacters();
							break;
						}
						// search rdf:Description grandchildren for "inchi" element (SBML Level 2 recommendation)
						for (int k=0; k<child.getNumChildren(); k++) {
							XMLNode grandchild = child.getChild(k);
							XMLNode greatgrandchild = grandchild.getChild(0);
							if (greatgrandchild.getName().equals("inchi")) {
								InChI = greatgrandchild.getChild(0).getCharacters();
								break;
							}
						}
						if (InChI != null)
							break;
					}
					// store the name and id of the Species
					compoundId = compound.getId();
					compoundName = compound.getName();
				}
			}
			
			if (InChI != null) {
				mass = parseInChIString("", InChI,infoWriter);
				if (mass != null) {
					// determine hydrogen and phosphate forms
					if (Properties.isHydrogen(mass))
						hasHydrogen = true;
					int type = Properties.isPhosphate(mass);
					if (type >= 0 && type <= 3)
						hasPhosphate[type] = true;
					if (massesMap.put(compoundId, mass) == null)
						annotated++;
				} else {
					infoWriter.write("Invalid element in compound "+compoundName+" (id "+compoundId+"), InChI: "+InChI+"\n");
				}
				
			} else {
				infoWriter.write(compoundName+" (id "+compoundId+") has no annotated mass.\n");
			}
		}
		
		// add missing hydrogen and phosphate forms
		if (!hasHydrogen)
			massesMap.put(Utilities.HYDROGEN_NAME, Utilities.HYDROGEN_MASS);
		for (int i=0; i<hasPhosphate.length; i++) {
			if (!hasPhosphate[i])
				massesMap.put(Utilities.PHOSPHATE_NAMES[i], Utilities.PHOSPHATE_FORMS[i]);
		}
		
		if (DEFAULT_NUMBER_OF_COMPOUNDS < massesMap.size())
			System.out.println("Initial masses hash capacity "+DEFAULT_NUMBER_OF_COMPOUNDS+" smaller than final map size ("+massesMap.size()+"). Increase DEFAULT_NUMBER_OF_COMPOUNDS to at least "+massesMap.size()+" for better performance.");
		
		float ratioAnnotated = (float)annotated*100/existing.size();
		String message = existing.size()+" compounds in the original network, "+annotated+" have annotated mass with ";
		for (int i=0; i<ELEMENTS.length; i++)
			message += ELEMENTS[i]+(i<ELEMENTS.length-1?",":"");
		message += " ("+ratioAnnotated+"%).";
		System.out.println(message);
		infoWriter.write(message+"\n");
		infoWriter.close();
		
		return massesMap;
	}
	
	/**
		 * Reads the mass vectors from the given BioCyc compounds file into a hash map containing the
		 * UNIQUE-ID as String key and the mass vector as int[] value.
		 * 
		 * @param compoundsFile
		 * @return
		 * @throws IOException
		 */
		protected static HashMap<String, int[]> createMassesBioCyc(String compoundsFile, String reactionsFile, String infoFile) throws IOException {
			BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
			BufferedReader reactionsReader = new BufferedReader(new FileReader(reactionsFile));
			BufferedWriter infoWriter = new BufferedWriter(new FileWriter(infoFile));
			HashMap<String, int[]> massesMap = new HashMap<String, int[]>((int)(DEFAULT_NUMBER_OF_COMPOUNDS/0.75f+1), 0.75f);
			String line, name = "";
			int numAnnotated = 0; // , compounds = 0, notInNetwork = 0, invalid = 0
			int[] mass = new int[ELEMENTS.length];
			boolean exists = false, valid = false, annotated = false;
			boolean[] hasPhosphate = {false, false, false, false};
			boolean hasHydrogen = false;
			String invalidElement = "";
			TreeMap<String, Integer> invalidElements = new TreeMap<String, Integer>();
			
			infoWriter.write("# "+getVersion(compoundsFile)+"\n");
			
			// collect the compounds actually in the network
			HashSet<String> existing = new HashSet<String>((int)(Utilities.DEFAULT_NUMBER_OF_COMPOUNDS/0.75f+1));
			boolean left = false, right = false, smallMoleculeReaction = false;
			while ((line = reactionsReader.readLine()) != null) {
				
				if (line.length() >= 32 && line.substring(0, 32).equals("TYPES - Small-Molecule-Reactions"))
					smallMoleculeReaction = true;
				if (line.length() >= 12 && line.substring(0, 12).equals("UNIQUE-ID - "))
					smallMoleculeReaction = false;
				
				// only consider small molecule reactions
				if (smallMoleculeReaction) {
					if (line.length() >= 7 && line.substring(0, 7).equals("LEFT - "))
						left = true;
					else if (line.length() >= 8 && line.substring(0, 8).equals("RIGHT - "))
						right = true;
					
					if (left || right) {
						String compound = (left? line.substring(7) : line.substring(8));
						// remove the pipes				
						if (compound.charAt(0) == '|' && compound.charAt(compound.length()-1) == '|')
							compound = compound.substring(1, compound.length()-1);
						existing.add(compound);
						left = false;
						right = false;
					}
				}
			}
			
			while ((line = compoundsReader.readLine()) != null) {
				
				if (line.length() >= 2 && line.substring(0, 2).equals("//")) {
					// compound parsing finished: add name and mass to the equivalence class
					
					if (!exists) {
//						notInNetwork++;
						infoWriter.write("Compound "+name+" in compounds file, but not in network.\n");
					} else if (!annotated) {
//						compounds++;
						infoWriter.write("Compound "+name+" has no annotated mass.\n");
					} else if (!valid) {
//						invalid++;
//						compounds++;
						infoWriter.write("Compound "+name+" contains chemical element(s) other than ");
						for (int i=0; i<ELEMENTS.length; i++)
							infoWriter.write(ELEMENTS[i]);
						infoWriter.write(" ("+invalidElement+")\n");
					} else {
//						compounds++;
						numAnnotated++;
						// determine hydrogen and phosphate forms
						int type = Properties.isPhosphate(mass);
						if (type >= 0 && type <= 3)
							hasPhosphate[type] = true;
						if (Properties.isHydrogen(mass))
							hasHydrogen = true;
						// remove mass vectors containing only zeroes
						if (Arrays.equals(mass, new int[ELEMENTS.length]))
							massesMap.put(name, null);
						else
							massesMap.put(name, mass);
					}
					
					// reset name and formula
					name = "";
					mass = new int[ELEMENTS.length];
					
				} else if (line.length() >= 12 && line.substring(0, 12).equals("UNIQUE-ID - ")) {
					// extract compound name and enable formula parsing
					name = line.substring(12, line.length());
					exists = existing.contains(name);
					valid = true;
					annotated = false;
	
				} else if (exists && valid && line.length() >= 19 && line.substring(0, 19).equals("CHEMICAL-FORMULA - ")) {
					// set element index to current element, valid remains false if the formula results invalid
					annotated = true;
					int element = -1;
					valid = false;
					for (int i=0; i<ELEMENTS.length; i++) {
						if (line.substring(20, line.indexOf(" ", 20)).equalsIgnoreCase(ELEMENTS[i]) && line.charAt(20+ELEMENTS[i].length()) == ' ') {
							element = i;
							valid = true;
							break;
						}
					}
					if (!valid) {
						invalidElement = line.substring(20, line.indexOf(' ', 20));
						if (invalidElements.containsKey(invalidElement))
							invalidElements.put(invalidElement, invalidElements.get(invalidElement)+1);
						else
							invalidElements.put(invalidElement, 1);
						continue;
					}
					// extract the element
					mass[element] = Integer.parseInt(line.substring(20+ELEMENTS[element].length()+1, line.lastIndexOf(')')));
				}
			}
			
			// add missing hydrogen and phosphate forms
			if (!hasHydrogen)
				massesMap.put(Utilities.HYDROGEN_NAME, Utilities.HYDROGEN_MASS);
			for (int i=0; i<hasPhosphate.length; i++) {
				if (!hasPhosphate[i])
					massesMap.put(Utilities.PHOSPHATE_NAMES[i], Utilities.PHOSPHATE_FORMS[i]);
			}
			
			// determine the most frequent undefined element; may be removed for better performance
	//		String topElement = "", secondElement = "";
			int count = 0, count2 = 0;
			infoWriter.write("Non-");
			for (int i=0; i<ELEMENTS.length; i++)
				infoWriter.write(ELEMENTS[i]);
			infoWriter.write(" elements: ");
			for (Iterator<String> elements = invalidElements.keySet().iterator(); elements.hasNext();) {
				String key = elements.next();
				infoWriter.write(key);
				if (elements.hasNext())
					infoWriter.write(", ");
				else
					infoWriter.write("\n");
				Integer value = invalidElements.get(key);
				if (value.intValue() > count) {
	//				topElement = key;
					count = value.intValue();
				} else if (value.intValue() > count2) {
	//				secondElement = key;
					count2 = value.intValue();
				}
			}
			
			if (DEFAULT_NUMBER_OF_COMPOUNDS < massesMap.size())
				System.out.println("Initial masses hash capacity "+DEFAULT_NUMBER_OF_COMPOUNDS+" smaller than final map size ("+massesMap.size()+"). Increase DEFAULT_NUMBER_OF_COMPOUNDS to at least "+massesMap.size()+" for better performance.");
			
			float ratioAnnotated = (float)numAnnotated*100/existing.size();		
			String message = existing.size()+" compounds in the original network, "+numAnnotated+" have annotated mass with ";
			for (int i=0; i<ELEMENTS.length; i++)
				message += ELEMENTS[i]+(i<ELEMENTS.length-1?",":"");
			message += " ("+ratioAnnotated+").";
	//		("+ratioAnnotated+"%), "+numValid+" contain only "+new String(ELEMENTS)+" ("+ratioValid+"%). Most frequent undefined elements are "+topElement+" ("+count+") and "+secondElement+" ("+count2+").\n"; 
			System.out.println(message);
			infoWriter.write(message+"\n");
			
			compoundsReader.close();
			infoWriter.close();
			return massesMap;
		}

	/**
	 * Parses an InChI string and returns an int[ELEMENTS.length] containing the mass vector. 
	 * 
	 * @param InChI
	 * @return Chemical formula string.
	 */
	public static int[] parseInChIString(String name, String InChI, BufferedWriter infoWriter) throws IOException {
		boolean debug = false;
		int[] mass = new int[ELEMENTS.length];
		String massString;
		int start = InChI.indexOf("/")+1;
		int end = start;
		
		while (end < InChI.length()
				&& (Character.isLetter(InChI.charAt(end)) || Character.isDigit(InChI.charAt(end)) || InChI.charAt(end) == '.'))
			end++;
		massString = InChI.substring(start, end);
		
		if (massString.length() == 0)
			return null;
		
		if (massString.equals("p") && InChI.charAt(end) == '+') {
			// hydrogen is represented as a positive charge p+x: replace by Hx
			start = start+2;
			end++;
			while (end < InChI.length() && Character.isDigit(InChI.charAt(end)))
				end++;
			int charge = Integer.parseInt(InChI.substring(start, end));
			massString = "H"+charge;
		}
		
		if (debug)
			System.out.print("Compound: "+name+", massString: "+massString+", ");
		
		for (int i=0; i<massString.length(); i++) {
			
			// skip digits, dots, and other stuff
			if (Character.isLetter(massString.charAt(i))) {
				boolean valid = false;
				// check if the element is valid
				int element = 0;
				// match any defined ELEMENT
				for (; element < ELEMENTS.length; element++) {
					String e = ELEMENTS[element];
					// skip if element mismatches, element length is larger than the remaining massString,
					// or if the next character is a lower case character
					if (massString.length() >= i+e.length() && e.equals(massString.substring(i, i+e.length()))
							&& (i+e.length() == massString.length() || !Character.isLowerCase(massString.charAt(i+e.length())))) {
						valid = true;
						i += e.length()-1;
						break;
					}
				}
				// parse the mass
				if (valid) {
					end = i+1;
					while (end < massString.length() && Character.isDigit(massString.charAt(end)))
						end++;
					if (end > i+1)
						mass[element] += Integer.parseInt(massString.substring(i+1, end));
					else 
						mass[element] += 1;
				} else {
					
					if (debug || infoWriter != null) {
						// print the missing elements
						int j=i+1;
						for (; j<massString.length(); j++) {
							try {
								Integer.parseInt(massString.substring(j, j+1));
								break;
							} catch (NumberFormatException e) {
								if (Character.isUpperCase(massString.charAt(j)))
									break;
							}
						}
						if (debug)
							System.out.println("Invalid element: "+massString.substring(i, j)+"\n");
						if (infoWriter != null)
							infoWriter.write("Invalid element: "+massString.substring(i, j)+"\n");
					}
					
					return null;
				}
			}
		}
		
		// remove mass vectors containing only zeroes
		if (Arrays.equals(mass, new int[ELEMENTS.length]))
			return null;
		
		if (debug) {
			System.out.print("parsed mass: ");
			for (int i=0; i<ELEMENTS.length; i++)
				if (mass[i] != 0)
					System.out.print(ELEMENTS[i]+String.valueOf(mass[i]));
			System.out.println();
		}
		
		return mass;
	}
		
	/**
	 * Parses a string of format CXHYNZ..., where C,H,N are chemical elements followed by the
	 * number of atoms X,Y,Z. The elements need to be in a specific order. Dots '.' are treated as 1.
	 * Returns an int[ELEMENTS.length] containing the mass vector.
	 * 
	 * @param massString
	 * @return The mass vector.
	 */
	public static int[] parseMass(String massString, BufferedWriter infoWriter) throws IOException {
		boolean debug = false;
		int[] mass = new int[ELEMENTS.length];
		
		if (massString.length() == 0)
			return null;
		
		if (debug)
			System.out.print("Parsing "+massString+": ");
		
		// treat sum formulas enclosed by 'n' as with n=1 (same equivalence class)
		if (massString.charAt(0) == '(' && massString.charAt(massString.length()-2) == ')' && massString.charAt(massString.length()-1) == 'n')
			massString = massString.substring(1, massString.length()-2);
		
		for (int i=0; i<massString.length(); i++) {
			
			// skip digits and dots
			if (!Character.isDigit(massString.charAt(i)) && massString.charAt(i) != '.') {
				boolean valid = false;
				// check if the element is valid
				int element = 0;
				// match any defined ELEMENT
				for (; element < ELEMENTS.length; element++) {
					String e = ELEMENTS[element];
					// skip if element mismatches, element length is larger than the remaining massString,
					// or if the next character is a lower case character
					if (massString.length() >= i+e.length() && e.equals(massString.substring(i, i+e.length()))
							&& (i+e.length() == massString.length() || !Character.isLowerCase(massString.charAt(i+e.length())))) {
						valid = true;
						i += e.length()-1;
						break;
					}
				}
				// parse the mass
				if (valid && mass[element] == 0) {
					int end = i+1;
					while (end < massString.length() && Character.isDigit(massString.charAt(end)))
						end++;
					if (end > i+1)
						mass[element] = Integer.parseInt(massString.substring(i+1, end));
					else 
						mass[element] = 1;
				} else {
					
					if (debug || infoWriter != null) {
						// print the missing elements
						int j=i+1;
						for (; j<massString.length(); j++) {
							try {
								Integer.parseInt(massString.substring(j, j+1));
								break;
							} catch (NumberFormatException e) {
								if (Character.isUpperCase(massString.charAt(j)))
									break;
							}
						}
						if (debug)
							System.out.println("Invalid element: "+massString.substring(i, j)+"\n");
						if (infoWriter != null)
							infoWriter.write("Invalid element: "+massString.substring(i, j)+"\n");
					}
					
					return null;
				}
			}
		}
		
		// remove mass vectors containing only zeroes
		if (Arrays.equals(mass, new int[ELEMENTS.length]))
			return null;
		
		if (debug) {
			for (int i=0; i<ELEMENTS.length; i++)
				System.out.print(mass[i]+" ");
			System.out.println();
		}
		
		return mass;
	}
	
	
	/**
	 * Tests whether an equation side delimiter starts at the given index i
	 * in the String s. Returns the length of the found delimiter or -1, if no
	 * delimiter was found at the given position.
	 * 
	 * @param s
	 * @param i
	 * @return
	 */
	protected static int isSideDelimiter(String s, int i) {
		int equalsLength = DELIMITER_EQUALS.length();
		int forwardLength = DELIMITER_FORWARD.length();
				
		if (i <= s.length()-equalsLength && DELIMITER_EQUALS.equals(s.substring(i, i+equalsLength)))
			return equalsLength;
		else if (i <= s.length()-forwardLength && DELIMITER_FORWARD.equals(s.substring(i, i+forwardLength)))
			return forwardLength;
		else
			return -1;
	}
			
	/**
	 * Tests whether an equation side delimiter or plus starts at the given index i
	 * in the String s. Returns the length of the found delimiter or -1, if no
	 * delimiter was found at the given position.
	 * 
	 * @param s
	 * @param i
	 * @return
	 */
	protected static int isDelimiter(String s, int i) {
		int plusLength = DELIMITER_PLUS.length();
		
		if (i <= s.length()-plusLength && DELIMITER_PLUS.equals(s.substring(i, i+plusLength)))
			return plusLength;
		else
			return isSideDelimiter(s, i);
	}
	
	/**
	 * Parses a string and returns the stoichiometric coefficient, or -1 if
	 * the string is not a stoichiometric coefficient, i.e., a double number
	 * enclosed by DELIMITER_COEFFICIENT_START and DELIMITER_COEFFICIENT_END.
	 * 
	 * @param s
	 * @return The stoichiometric coefficient, or -1 if the string is not a
	 * stoichiometric coefficient.
	 */
	public static double parseCoefficient(String s) {
		double coefficient = -1;
		if (s.length() > 2 && s.charAt(0) == DELIMITER_COEFFICIENT_START && s.charAt(s.length()-1) == DELIMITER_COEFFICIENT_END) {
			 try {
				 coefficient = Double.parseDouble(s.substring(1, s.length()-1));
			 } catch (NumberFormatException e) {
				 System.out.println("invalid coefficient: "+s);
				 return -1;
			 }
		}
		if (coefficient == -1)
			System.out.println("invalid coefficient: "+s);
		return coefficient;
	}
	
	/**
	 * Calculates the greatest common divisor the given values.
	 * 
	 * @param values
	 * @return The greatest common divisor.
	 */
	public static long gcd(long[] values) {
		long m = values[0];
		// calculate m, the greatest common divisor of values
		for (int i=0; i<values.length-1; i++) {
			long n = values[i+1];
			long r;
			while (n != 0) {
				r = m % n;
				m = n;
				n = r;
			}
		}
		return m;
	}
	
	/**
	 * Returns the least common integer multiple for x and y, determined
	 * by dividing x*y by the greatest common divisors.
	 * 
	 * @param x
	 * @param y
	 * @return The least common integer multiple of x and y.
	 */
	public static long lcm(int x, int y) {
		long xy = x*y;
		long gcd = gcd(new long[]{(long)x,(long)y});
		long lcm = xy/gcd;
		
//		System.out.println("x="+x+", y="+y+", xy="+xy+", gcd="+gcd+", lcm="+lcm);
		
		// check for numerical imprecision
		if (lcm < 0 || lcm%x != 0 || lcm%y != 0)
			return -1;
		
		return lcm;
	}
	
	/**
	 * Calculates the binomial coefficient (n over k).
	 * 
	 * @param n
	 * @param k
	 * @return The binomial coefficient.
	 */
	public static long binomialCoefficient(long n, long k) {
		long a = 1, b = 1;
		
		if (k<0)
			throw new IllegalArgumentException("n must be positive.");
		if (n < k)
			throw new IllegalArgumentException("n="+n+" must be larger than k="+k+".");
		if (k==0 || k==n)
			return 1;
		if (k==1 || k==n-1)
			return n;
		
		for (long i=n-k+1; i<=n; i++) {
			a *= i;
			if (a < 1 || a == Long.MAX_VALUE)
				throw new ArithmeticException("long overflow. n="+n+", k="+k+", a="+a);
		}
		for (long i=1; i<=k; i++) {
			b *= i;
			if (b < 1 || b == Long.MAX_VALUE)
				throw new ArithmeticException("long overflow. n="+n+", k="+k+", b="+b);
		}
		
		if ((a % b) != 0)
			throw new ArithmeticException("Imprecise division: a="+a+", b="+b+" (n="+n+", k="+k+")");
		
		return a/b;
	}
	
	/**
	 * Tests a and b for linear dependence.
	 * 
	 * @param a
	 * @param b
	 * @return The factor f=b/a, or 0, if the vectors are linearly independent.
	 */
	public static double linearDependent(double[] a, double[] b) {
		double f=0;
		if (a.length != b.length)
			throw new IllegalArgumentException("Vectors must have equal length.");
		
		for (int i=0; i<a.length; i++) {
			
			if (a[i] == 0 && b[i]== 0)
				continue;
			if (a[i] == 0 ^ b[i] == 0)
				return 0;
			if (f == 0) {
				f = b[i]/a[i];
			} else {
				if (Math.abs(f-(b[i]/a[i])) > 0.0001)
					return 0;
			}
		}
		
		return f;
	}
	
	/**
	 * Maps the compound names given in indexFile to the pair of corresponding names given in mappingFile,
	 * and prints the indices of the mapped compounds in the sorted stoichiometric matrix.
	 * 
	 * @param graph
	 * @param mappingFile
	 * @param indexFile
	 * @throws IOException
	 */
	public static void printMappedIndices(MetabolicGraph graph, String mappingFile, String indexFile) throws IOException {
		BufferedReader mappingReader = new BufferedReader(new FileReader(mappingFile));
		BufferedReader indexReader = new BufferedReader(new FileReader(indexFile));
		String line;
		HashMap<String, String> mapping = new HashMap<String,String>();
		Collection<Vertex> compoundsSet = graph.getCompounds();
		
		// sort the compounds by name and compartment
		ArrayList<Vertex> compoundsList = new ArrayList<Vertex>(compoundsSet);
    	Collections.sort(compoundsList, new VertexComparator());

		// map the keys and values
		while ((line = mappingReader.readLine()) != null) {
//			System.out.println("key: "+line.substring(line.indexOf("\t")).trim()+", value: "+line.substring(0,line.indexOf("\t")).trim());
			mapping.put(line.substring(line.indexOf("\t")).trim(), line.substring(0,line.indexOf("\t")).trim());
		}

		// find the mapped values
		while ((line = indexReader.readLine()) != null) {
			String mappedName = mapping.get(line.trim());
			Vertex v = graph.getCompound(mappedName, "c");
			if (compoundsList.indexOf(v) != -1) {
				System.out.println(compoundsList.indexOf(v));
			} else {
				// try 
				
				System.out.println("Not found: "+line.trim()+"."+" Mapping: "+mappedName);
			}
		}
	}
	
	/**
	 * Convenience method to map compound strings with their masses from 2 files: one contains the compound strings and
	 * aliases, the other maps the aliases to sum formulas. Writes a table of compound names and corresponding mass strings. 
	 * 
	 * @param compoundsFile
	 * @param massesFile
	 * @param outputFile
	 * @throws IOException
	 */
	public static void compoundAliasMap(String compoundsFile, String massesFile, String outputFile) throws IOException {
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader massesReader = new BufferedReader(new FileReader(massesFile));
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
		
		String line;
		int compounds = 0, masses = 0;
		HashMap<String,String> compoundMasses = new HashMap<String,String>(); 
		
		// read masses and their identifiers into a hash map
		while ((line = massesReader.readLine()) != null) {
			
			if (line.indexOf("\t") >= 0 && line.indexOf("\t") != line.length()-1) { 
				String name = line.substring(0, line.indexOf("\t"));
				String mass = line.substring(line.indexOf("\t")+1);
				compoundMasses.put(name, mass);				
			} else {
				System.out.println("No mass for "+line);
			}
		}
		
		// read compound names and print the mapped names with masses
		while ((line = compoundsReader.readLine()) != null) {
			compounds++;
			String name, alias1, alias2 = "";
			
			int tab = line.indexOf("\t");
			name = line.substring(0, tab);
			
			if (line.indexOf("\t", tab+1) < 0) {
				alias1 = line.substring(tab+1, line.length());
			} else {
				tab = line.indexOf("\t", tab+1);
				alias1 = line.substring(line.indexOf("\t")+1, tab);
				alias2 = line.substring(tab+1, line.length());
			}
			
			// remove the pipes				
			if (alias1.charAt(0) == '|' && alias1.charAt(alias1.length()-1) == '|')
				alias1 = alias1.substring(1, alias1.length()-1);
			if (alias2.charAt(0) == '|' && alias2.charAt(alias2.length()-1) == '|')
				alias2 = alias2.substring(1, alias2.length()-1);
			
			String mass = compoundMasses.get(alias1);
			if (mass == null || mass.equals(""))
				// try second alias
				mass = compoundMasses.get(alias2);
			if (mass == null || mass.equals(""))
				System.out.println("No mass found for: "+name+", alias1: "+alias1+", alias2: "+alias2);
			else {		
				outputWriter.write(name+"\t"+mass+"\n");
				masses++;
			}
		}
		
		float ratioMasses = (float)masses*100/compounds;
		System.out.println("Found "+masses+" masses for "+compounds+" compounds ("+ratioMasses+").");
		compoundsReader.close();
		massesReader.close();
		outputWriter.close();
		
	}
	
	/**
	 * Searches the compounds from compoundsFile in the databaseFile by case-insensitive name and sum formula,
	 * and prints the compounds with deltaG to outputFile. A difference in hydrogen atoms is tolerated as
	 * specified by hydrogenTolerance. Used to append the deltaG values from the KEGG compounds to the compounds
	 * file of a network.
	 * 
	 * @param compoundsFile
	 * @param databaseFile
	 * @param outputFile
	 * @param hydrogenTolerance
	 */
	public static void appendDeltaG(String compoundsFile, String databaseFile, String outputFile, int hydrogenTolerance) throws IOException {
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		int counter = 0;
		String line;
		HashMap<String,String> formulas = new HashMap<String,String>();
		HashMap<String,String> deltaGs = new HashMap<String,String>();
		HashMap<String,String> uncertainties = new HashMap<String,String>();
		
		// store the lower-case compound synonyms, formulas, deltaG, and uncertainties
		while ((line = databaseReader.readLine()) != null) {
			// get tokens (don't use lower-case here, as we need case-sensitive sum formulas)
			ArrayList<String> tokens = parseTokens(line, "\t", false);
			ArrayList<String> names = new ArrayList<String>(tokens.size());
			
			// parse the names
			for (String token : tokens) {
				if (token.length() > 1 && token.charAt(1) != '=')
					names.add(token.toLowerCase());
			}
			// parse the formula, deltaG, and uncertainty
			for (String token : tokens) {
				try {
					// add a valid formula, deltaG, or uncertainty for every synonym
					if (token.startsWith("F=")) {
						if (parseMass(token.substring(2), null) == null)
							continue;
						for (String name : names)
							formulas.put(name, token);
					} else if (token.startsWith("D=")) {
						Double.parseDouble(token.substring(2));
						for (String name : names)
							deltaGs.put(name, token);
					} else if (token.startsWith("U=")) {
						Double.parseDouble(token.substring(2));
						for (String name : names)
							uncertainties.put(name, token);
					}
				} catch (NumberFormatException e) {
					// not a valid deltaG/uncertainty
				}
			}
		}
		
		// print the mappings
//		for (String name : formulas.keySet())
//			System.out.println(name+"\tF="+formulas.get(name)+"\tD="+deltaG.get(name)+"\tU="+uncertainties.get(name));
		
		// iterate over the compounds file, compare the sum formula and write with deltaG to outputfile
		while ((line = compoundsReader.readLine()) != null) {
			
			if (line.startsWith("#")) {
				writer.write(line+"\n");
				continue;
			}
			
			// get tokens (don't use lower-case here, as we need case-sensitive sum formulas)
			ArrayList<String> tokens = parseTokens(line, "\t", false);
			String compound = tokens.get(0).toLowerCase();
			// parse the sum formula of the original compound
			int[] formula = null;
			for (String token : tokens) {
				if (token.length() > 1 && token.startsWith("F=")) {
					formula = parseMass(token.substring(2), null);
					break;
				}
			}
			// get the sum formula of the database compound
			String formulaDBString = formulas.get(compound);
			int[] formulaDB = null;
			if (formulaDBString != null && formulaDBString.length() > 1)
				formulaDB = parseMass(formulaDBString.substring(2), null);
			// try to find the deltaG and uncertainty of the compound
			String deltaGDB = deltaGs.get(compound);
			String uncertaintyDB = uncertainties.get(compound);
			
			if (deltaGDB == null) {
				System.out.println("Not found: "+compound);
			} else if (formula == null || formulaDB == null) {
				System.out.println("No formula in "+(formula == null?"compounds file":"database file")+": "+compound);
			} else {
				// compare the sum formula, allow the specified hydrogenTolerance
				boolean equal = true;
				for (int i=0; i<ELEMENTS.length; i++) {
					if ((!ELEMENTS[i].equals("H") && formula[i] != formulaDB[i]) || (ELEMENTS[i].equals("H") && Math.abs(formula[i]-formulaDB[i]) > hydrogenTolerance)) {
						equal = false;
						break;
					}
				}
				if (equal) {
					// add the deltaG
					line = line+"\t"+deltaGDB+"\t"+uncertaintyDB;
					counter++;
				} else {
					System.out.print("Unequal sum formula for "+compound+": ");
					for (int i=0; i<ELEMENTS.length; i++)
						System.out.print(ELEMENTS[i]+formula[i]);
					System.out.print(" != ");
					for (int i=0; i<ELEMENTS.length; i++)
						System.out.print(formulaDB[i]);
					System.out.println();
				}
			}
			writer.write(line+"\n");
		}
		
		System.out.println("Added "+counter+" deltaG values.");
		
		compoundsReader.close();
		databaseReader.close();
		writer.close();
	}
	
	/**
	 * Searches a KEGG ID in compoundsFile, matches the compound in databaseFile and appends
	 * the deltaG and uncertainties to the compounds in outputFile. Used to append the deltaG values
	 * from the KEGG compounds to the compounds file of a network, matching by ID.
	 * 
	 * @param compoundsFile
	 * @param databaseFile
	 * @param outputFile
	 */
	public static void appendDeltaGFromKEGG(String compoundsFile, String databaseFile, String outputFile) throws IOException {
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		String line;
		HashMap<String,String> deltaGs = new HashMap<String,String>();
		HashMap<String,String> uncertainties = new HashMap<String,String>();
		int added = 0;
		
		// store the KEGG IDs, deltaG, and uncertainties
		while ((line = databaseReader.readLine()) != null) {
				
			String id = line.substring(0, line.indexOf("\t")).trim();
			int start = line.indexOf("\t")+1, end = 0;
			
			// parse the tab-separated tokens
			do {
				end = line.indexOf("\t", start);
				if (end == -1)
					end = line.length();
				
				String token = line.substring(start, end);
				start = end+1;
				// extract the value and add it to the corresponding map
				if (token.startsWith("D="))
					deltaGs.put(id, token);
				else if (token.startsWith("U="))
					uncertainties.put(id, token);
				
			} while (end < line.length());
		}
		
		System.out.println("Database file has "+deltaGs.size()+" DeltaG values.");
		
		// iterate over the compounds file, search for the id and write with deltaG to outputfile
		while ((line = compoundsReader.readLine()) != null) {
			int start = 0, end = 0;
			
			if (line.startsWith("#")) {
				writer.write(line+"\n");
				continue;
			}
			
			// parse the tab-separated tokens
			do {
				end = line.indexOf("\t", start);
				if (end == -1)
					end = line.length();
				
				String token = line.substring(start, end);
				start = end+1;
				// it is a KEGG id, then append the deltaG and uncertainty 
				if (token.matches("C\\d\\d\\d\\d\\d")) {
					String deltaG = deltaGs.get(token);
					String uncertainty = uncertainties.get(token);
					if (deltaG != null) {
						line = line + "\t" + deltaG + "\t" + uncertainty;
						added++;
					}
					break;
				}
				
			} while (end < line.length());
			
			writer.write(line+"\n");
		}
		
		System.out.println("Appended "+added+" values in "+outputFile);
		
		compoundsReader.close();
		databaseReader.close();
		writer.close();
	}
	
	/**
	 * Parses the sum formulas and deltaG values from the databaseFile and determines
	 * the differences in deltaGf between compounds which have an equal sum formula
	 * up to the specified hydrogenTolerance in the number of hydrogen atoms.
	 * 
	 * @param databaseFile
	 * @param hydrogenTolerance
	 * @throws IOException
	 */
	public static void deltaGDifferences(String databaseFile, int hydrogenTolerance) throws IOException {
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		String line;
		
		class Compound {
			int[] mass;
			Double deltaG;
		};
		
		// parse the compound masses and deltaGf
		ArrayList<Compound> compounds = new ArrayList<Compound>();
		while ((line = databaseReader.readLine()) != null) {
			
			Compound compound = new Compound();
			ArrayList<String> tokens = parseTokens(line, "\t", false);
			for (String token : tokens) {
				if (token.length() > 1) {
					if (token.startsWith("F=")) {
						int[] mass = parseMass(token.substring(2), null);
						if (mass != null)
							compound.mass = mass;
					} else if (token.startsWith("D=")) {
						try {
							compound.deltaG = Double.parseDouble(token.substring(2));
						} catch (NumberFormatException e) {
							// invalid deltaG
						}
					}
				}
			}
			if (compound.mass != null && compound.deltaG != null)
				compounds.add(compound);
		}
		databaseReader.close();
		
		// identify the differences in deltaGf between compounds with equal/similar sum formula
		ArrayList<Double> differences = new ArrayList<Double>((compounds.size()*(compounds.size()-1))/2);
		double maxDifference = 0, sumDifference = 0;
		int numCompounds = compounds.size();
		int counter = 0;
		
		for (int i=0; i<numCompounds-1; i++) {
			for (int j=i+1; j<numCompounds; j++) {
				Compound compound1 = compounds.get(i);
				Compound compound2 = compounds.get(j);
				// compare the sum formulas
				boolean equal = true;
				for (int k=0; k<ELEMENTS.length; k++) {
					if ((!ELEMENTS[k].equals("H") && compound1.mass[k] != compound2.mass[k])
							|| (ELEMENTS[k].equals("H") && Math.abs(compound1.mass[k]-compound2.mass[k]) > hydrogenTolerance)) {
						equal = false;
						break;
					}
				}
				// store the deltaG difference
				if (equal) {
					double difference = Math.abs(compound1.deltaG-compound2.deltaG);
					differences.add(difference);
					if (difference > maxDifference)
						maxDifference = difference;
					sumDifference += difference;
					counter++;
					if (difference >= 30) {
						System.out.println("Difference "+difference+" = ("+compound1.deltaG+")-("+compound2.deltaG+")");
						System.out.print("mass1: ");
						for (int k=0; k<ELEMENTS.length; k++)
							if (compound1.mass[k] != 0)
								System.out.print(ELEMENTS[k]+compound1.mass[k]);
						System.out.print(", mass2: ");
						for (int k=0; k<ELEMENTS.length; k++)
							if (compound2.mass[k] != 0)
								System.out.print(ELEMENTS[k]+compound2.mass[k]);
						System.out.println();
					}
				}
			}
		}
		
		for (Double difference : differences)
			System.out.println(String.valueOf(difference));
		
		System.out.println("Average difference: "+String.valueOf((float)(sumDifference/counter))+", maximum difference: "+maxDifference);
	}
	
	/**
	 * Parses the input string and returns a list of tokens separated by delimiter.
	 * If lowercase==true, all characters are converted to lower-case.
	 * 
	 * @param input
	 * @param delimiter
	 * @param lowercase
	 * @return ArrayList containing the parsed tokens.
	 */
	public static ArrayList<String> parseTokens(String input, String delimiter, boolean lowercase) {
		ArrayList<String> tokens = new ArrayList<String>();
		int start = 0, end = 0;
		// parse the delimiter-separated tokens
		do {
			end = input.indexOf(delimiter, start);
			if (end == -1)
				end = input.length();
			
			String token = input.substring(start, end);
			if (lowercase)
				token = token.toLowerCase();
			start = end+1;
			tokens.add(token);
			
		} while (end < input.length());
		
		return tokens;
	}
	
	/**
	 * Parses the inputMap file for synonyms and searches all synonyms in databaseFile.
	 * If a synonym matches, the database line is appended by the name from inputMap
	 * and written to outputFile. Used to add new compound synonyms from a mapping file
	 * of a network to the database, which can then be used by appendDeltaG to add
	 * the deltaG values to the compounds list.
	 * 
	 * @param inputMap
	 * @param databaseFile
	 * @param outputFile
	 */
	public static void joinSynonyms(String inputMap, String databaseFile, String outputFile, boolean caseSensitive, int minMatch, int minAdd) throws IOException {
		BufferedReader mapReader = new BufferedReader(new FileReader(inputMap));
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		HashSet<HashSet<String>> newSynonyms = new HashSet<HashSet<String>>();
		HashSet<String> added = new HashSet<String>();
		String line;
		int counter = 0;
		
		// create a 2-dimensional array of the new compound synonyms
		while ((line = mapReader.readLine()) != null) {
			ArrayList<String> tokens = parseTokens(line, "\t", !caseSensitive);
			HashSet<String> names = new HashSet<String>(tokens.size());
			for (String token : tokens) {
				if (token.length() > 1 && token.charAt(1) != '=')
					names.add(token);
			}
			newSynonyms.add(names);
		}
		
		// find the names from the inputMap and write the joined line
		while ((line = databaseReader.readLine()) != null) {
			ArrayList<String> dbTokens = parseTokens(line, "\t", false);
			ArrayList<String> dbTokensLowerCase = parseTokens(line, "\t", !caseSensitive);
			HashSet<String> found = null;
			String foundToken = null;
			
			// search all new synonyms for a hit
			for (HashSet<String> newSet : newSynonyms) {
				for (String newToken : newSet) {
					if (dbTokensLowerCase.contains(newToken)) {
						found = newSet;
						foundToken = newToken;
						break;
					}
				}
				if (found != null)
					break;
			}
			
			// add the synonyms not already present
			if (found != null) {
				for (String f : found) {
					if (f.length() > 1 && f.charAt(1)!='=' && !dbTokensLowerCase.contains(f)) {
						if (added.contains(f))
							System.out.println("Multiple match: '"+f+"' by '"+foundToken+"'.");
						else if (f.length() < minAdd)
							System.out.println("Name too short: '"+f+"', matched by '"+foundToken+"'.");
						else if (foundToken.length() < minMatch)	
							System.out.println("Match too short: '"+foundToken+"' for '"+f+"'.");
						else {
							dbTokens.add(1, f);
							added.add(f);
							System.out.println("Added: '"+f+"', matched by '"+foundToken+"'.");
							counter++;							
						}
					}
				}
			}
			
			// write the new line
			for (String token : dbTokens)
				writer.write(token+"\t");
			writer.write("\n");
		}
		
		System.out.println("Added "+counter+" synonyms.");
		
		writer.close();
		mapReader.close();
		databaseReader.close();
		
	}
	
	/**
	 * Parses the compounds file for synonyms (more than one name for a compound) and searches 
	 * similar synonyms with similar sum formula in databaseFile.
	 * 
	 * The names and sum formulas of the compounds file are collected, and then compared with
	 * the names and sum formulas in the database file. The distance between names from the
	 * compounds file and names from the database file is determined by the edit distance
	 * (Levenshtein distance) divided by the length of the larger name. If the distance is
	 * below maxDistance, and the sum formulas equal
	 * up to a difference of the number of hydrogen atoms specified by hydrogenTolerance, then
	 * all synonyms of the compound are added to the corresponding line in the database file.
	 * The output is written to outputFile.
	 *  
	 * Used to add new compound synonyms from a compounds file of a network to a KEGG database
	 * file with synonyms, which can then be used by appendDeltaG to add the deltaG values to
	 * the compounds.
	 * 
	 * @param compoundsFile
	 * @param databaseFile
	 * @param outputFile
	 * @param maxDistance
	 * @param hydrogenTolerance
	 */
	public static void joinSimilarSynonyms(String compoundsFile, String databaseFile, String outputFile, double maxDistance, int hydrogenTolerance) throws IOException {
		boolean debug = false;
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		HashSet<HashSet<String>> compounds = new HashSet<HashSet<String>>();
		HashSet<String> added = new HashSet<String>();
		String line;
		int counter = 0;
		
		// collect the compound names and formulas from compoundsFile
		while ((line = compoundsReader.readLine()) != null) {
			ArrayList<String> tokenArray = parseTokens(line, "\t", false);
			HashSet<String> tokens = new HashSet<String>(tokenArray);
			compounds.add(tokens);
		}
		
		// parse database file
		while ((line = databaseReader.readLine()) != null) {
			
			int[] newFormula = new int[ELEMENTS.length], dbFormula = new int[ELEMENTS.length];
			ArrayList<String> dbTokens = parseTokens(line, "\t", false);
			ArrayList<String> dbNames = new ArrayList<String>(dbTokens.size());
			
			// collect the lower-case names and sum formula from the database
			for (String dbToken : dbTokens) {
				if (dbToken.length() > 1) {
					if (dbToken.startsWith("F="))
						dbFormula = parseMass(dbToken.substring(2), null);
					else if (dbToken.charAt(1) != '=')
						dbNames.add(dbToken.toLowerCase());
				}
			}
			
			// search all compounds for an approximate name and formula
			for (HashSet<String> tokens : compounds) {
				
				boolean formulaHit = false;
				double minDistance = 1;
				String matchingDBToken = null, newToken = null;
				
				for (String token : tokens) {
					if (token.startsWith("F=")) {
						// compare the sum formulas
						newFormula = parseMass(token.substring(2), null);
						if (dbFormula != null && newFormula != null) {
							formulaHit = true;
							for (int i=0; i<ELEMENTS.length; i++) {
								if ((!ELEMENTS[i].equals("H") && dbFormula[i] != newFormula[i]) || (ELEMENTS[i].equals("H") && Math.abs(dbFormula[i]-newFormula[i]) > hydrogenTolerance)) {
									formulaHit = false;
									break;
								}
							}
						}
						if (dbFormula == null || newFormula == null || !formulaHit)
							break;
						
					} else if (token.charAt(1) != '=') {

						// compare the case-insensitive names
						for (String dbName : dbNames) {
							
							// replace all special characters by $ to improve the matching distance
							for (int i=0; i<dbName.length(); i++) {
								char c = dbName.charAt(i);
								if (!Character.isLetterOrDigit(c))
									dbName = dbName.replace(c, '$');
							}
							for (int i=0; i<token.length(); i++) {
								char c = token.charAt(i);
								if (!Character.isLetterOrDigit(c))
									token = token.replace(c, '$');
							}
							
							// replace multiple $ by one 
							dbName = dbName.replaceAll("(\\$)+", "\\$");
							token = token.replaceAll("(\\$)+", "\\$");
							
							// remove leading and trailing $
							if (dbName.charAt(0)=='$')
								dbName = dbName.substring(1);
							if (dbName.charAt(dbName.length()-1)=='$')
								dbName = dbName.substring(0,dbName.length()-1);
							if (token.charAt(0)=='$')
								token = token.substring(1);
							if (token.charAt(token.length()-1)=='$')
								token = token.substring(0,token.length()-1);
							
							if (token.charAt(token.length()-1)=='$' || (dbName.charAt(dbName.length()-1)=='$'))
								throw new RuntimeException("shit");
							
							// normalize the edit distance by the longer string, which corresponds to the maximum edit distance
							double distance = (double)editDistance(dbName.toLowerCase(),token.toLowerCase()) / Math.max(dbName.length(), token.length());
							if (distance < minDistance) {
								minDistance = distance;
								matchingDBToken = dbName;
								newToken = token;
							}
						}
					}
				}
				
				if (formulaHit && minDistance <= maxDistance) {
					
					if (debug) {
						// print the match
						System.out.print("Matching name: "+newToken+" (");
						for (int i=0; i<newFormula.length; i++)
							if (newFormula[i] != 0)
								System.out.print(ELEMENTS[i]+String.valueOf(newFormula[i]));
						System.out.print(") ~ "+matchingDBToken+" (");
						for (int i=0; i<dbFormula.length; i++)
							if (dbFormula[i] != 0)
								System.out.print(ELEMENTS[i]+String.valueOf(dbFormula[i]));					
						System.out.print(" (distance "+minDistance+").");
						System.out.println();
					}
							
					// add the synonyms which are not already in the database
					for (String token : tokens) {
						if (token.length() > 1 && token.charAt(1) != '=' && !dbTokens.contains(token)) {
							// allow for multiple matches (ideally we would add the new tokens only to the database names with the best match)
							if (added.contains(token))
								System.out.println("Multiple match: '"+token+"' by '"+matchingDBToken+"'.");
//							else if (token.length() < minLength) {
//								if (debug)
//									System.out.println("Skipping short compound: '"+token+"', matched by '"+matchingDBToken+"'.");
							dbTokens.add(1, token);
							added.add(token);
							counter++;
							if (debug)
								System.out.println("Added: '"+token+"', matched by '"+newToken+" (synonym) with "+matchingDBToken+" (database)'.");
						}
					}
				}
				
			}
			// write the new line
			for (String token : dbTokens)
				writer.write(token+"\t");
			writer.write("\n");
		}
		
		System.out.println("Added "+counter+" new synonyms.");
		
		writer.close();
		compoundsReader.close();
		databaseReader.close();
		
	}
	
	public static void matchStrings(String templateFile, String databaseFile, String outputFile, double maxDistance) throws IOException {
		BufferedReader templateReader = new BufferedReader(new FileReader(templateFile));
		BufferedReader databaseReader = new BufferedReader(new FileReader(databaseFile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		ArrayList<String> templates = new ArrayList<String>();
		HashMap<String,String> matches = new HashMap<String,String>();
		HashMap<String,Double> distances = new HashMap<String,Double>();
		String line;
		int hits = 0, exactHits = 0;
		
		// collect the template names
		while ((line = templateReader.readLine()) != null) {
			if (line.startsWith("#") || line.length() == 0)
				continue;
			int end = line.indexOf("\t");
			if (end == -1)
				end = line.length();
			templates.add(line.substring(0, end).trim());
		}
		
		// find templates in the databse
		while ((line = databaseReader.readLine()) != null) {
			if (line.startsWith("#") || line.length() == 0)
				continue;
			
			// extract the database string
			int end = line.indexOf("\t");
			if (end == -1)
				end = line.length();
			String db = line.substring(0, end).trim();
			
			for (String template : templates) {
				double distance = (double)editDistance(db.toLowerCase(),template.toLowerCase()) / Math.max(db.length(), template.length());
				if (distance <= maxDistance) {
					// add the match if it is the first one or better than the existing match
					if (!distances.containsKey(template) || distances.get(template) > distance) {
						matches.put(template, db);
						distances.put(template, distance);
					}
				}
			}
		}
		
		// print the best matches
		for (String match : matches.keySet()) {
			writer.write(match+"\t"+matches.get(match)+"\t"+String.valueOf(distances.get(match))+"\n");
			if (distances.get(match) == 0)
				exactHits++;
			else
				hits++;
		}
		
		System.out.println(exactHits +" exact matches, "+hits+" approximate matches with distance <= "+maxDistance);
		
		writer.close();
		templateReader.close();
		databaseReader.close();
		
	}
	
	/**
	 * Maps a list of compounds from compoundsFile to a list containing compound synonyms and masses in massesFile.
	 * There may be several tab-separated synonyms before the mass.
	 * 
	 * The outputFile contains a list of the compounds with found masses. 
	 * 
	 * @param compoundsFile
	 * @param massesFile
	 * @param outputFile
	 * @throws IOException
	 */
	public static void compoundMassesMap(String compoundsFile, String massesFile, String outputFile) throws IOException {
		boolean debug = false;
		BufferedReader compoundsReader = new BufferedReader(new FileReader(compoundsFile));
		BufferedReader massesReader = new BufferedReader(new FileReader(massesFile));
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
		
		String line;
		int compounds = 0, masses = 0;
		HashMap<String,String> compoundMasses = new HashMap<String,String>(); 
		
		// read masses and their identifiers into a hash map
		while ((line = massesReader.readLine()) != null) {
			
			if (line.indexOf("\t") >= 0 && line.indexOf("\t") != line.length()-1) {
				int start = 0;
				int tab = line.indexOf("\t");
				String name = line.substring(start, tab).toLowerCase();
				String mass = line.substring(line.lastIndexOf("\t")+1);
				if (debug)
					System.out.println("Putting compound: "+name+", mass: "+mass);
				compoundMasses.put(name, mass);
				// add synonyms
				start = tab+1;
				while ((tab = line.indexOf("\t", start)) != line.lastIndexOf("\t")) {
					name = line.substring(start, tab).toLowerCase();
					mass = line.substring(line.lastIndexOf("\t")+1);
					if (debug)
						System.out.println("Putting compound: "+name+", mass: "+mass);
					compoundMasses.put(name, mass);
					start = tab+1;
				}
				
			} else {
				System.out.println("No mass for "+line);
			}
		}
		
		while ((line = compoundsReader.readLine()) != null) {
			compounds++;
			
			String mass = compoundMasses.get(line.trim().toLowerCase());
			if (mass == null || mass.equals(""))
				System.out.println("No mass found for: "+line);
			else {
				outputWriter.write(line.trim()+"\t"+mass+"\n");
				masses++;
			}
		}
		
		float ratioMasses = (float)masses*100/compounds;
		System.out.println("Found "+masses+" masses for "+compounds+" compounds ("+ratioMasses+").");
		compoundsReader.close();
		massesReader.close();
		outputWriter.close();
	}
	
	public static void reachable(ArrayList<String> substrateNames, ArrayList<String> productNames, MetabolicGraph graph) throws IOException {
		Set<Vertex> substrates = new HashSet<Vertex>();
		Set<Vertex> products = new HashSet<Vertex>();
		
		// get the substrate vertices
		for (String substrate : substrateNames) {
			String compartment = null;
			if (substrate.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
				compartment = graph.hasCompartments() ? substrate.substring(substrate.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1, substrate.length()) : null;
				substrate = substrate.substring(0, substrate.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
			}
			Vertex compound = graph.getCompound(substrate, compartment);
			if (compound == null) {
				// not reachable: compound from the reaction is not in the network
				System.out.println("No such compound: "+substrate+(compartment == null ? "":compartment)+" in the network.");
				return;
			}
			substrates.add(compound);
		}
		
		// get the product vertices
		for (String product : productNames) {
			String compartment = null;
			if (product.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
				compartment = graph.hasCompartments() ? product.substring(product.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1, product.length()) : null;
				product = product.substring(0, product.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
			}
			Vertex compound = graph.getCompound(product, compartment);
			if (compound == null) {
				// not reachable: compound from the reaction is not in the network
				System.out.println("No such compound: "+product+(compartment == null ? "":compartment)+" in the network.");
				return;
			}
			products.add(compound);
		}
		
		// calculate the scope of the substrates until all products are included or no more compounds are added
		boolean reachable = false, added = true; 
		int seedSize = substrates.size();
		HashSet<Vertex> scope = new HashSet<Vertex>(5*seedSize, 0.75f);
		HashSet<Vertex> visited = new HashSet<Vertex>(5*seedSize, 0.75f);
		HashSet<Vertex> tempSeed = new HashSet<Vertex>(3*seedSize, 0.75f);
		// copying the set gives a speedup, due to optimized initial size?
		substrates = new HashSet<Vertex>(substrates);
		
		// repeat until no more compounds are added
		while (added) {
			added = false;
			
			// iterate over the seed compounds
			for (Iterator<Vertex> seedIterator = substrates.iterator(); seedIterator.hasNext();) {
				// iterate over the reactions reachable from the seed
				for (Iterator<DefaultWeightedEdge> reachableIterator = graph.outgoingEdgesOf(seedIterator.next()).iterator(); reachableIterator.hasNext();) {
					Vertex reaction = graph.getEdgeTarget(reachableIterator.next());
					if (!visited.contains(reaction) && !scope.contains(reaction)) {
						boolean producible = true;
						// determine whether the reachable reaction is producible by the seed
						for (Iterator<DefaultWeightedEdge> producibleIterator = graph.incomingEdgesOf(reaction).iterator(); producibleIterator.hasNext();) {
							if (!substrates.contains(graph.getEdgeSource(producibleIterator.next()))) {
								producible = false;
								break;
							}
						}
						// add the producible compounds to the seed
						if (producible) {
							for (Iterator<DefaultWeightedEdge> producedIterator = graph.outgoingEdgesOf(reaction).iterator(); producedIterator.hasNext();)
								tempSeed.add(graph.getEdgeTarget(producedIterator.next()));
							scope.add(reaction);
							// test if all products are in the current scope = termpSeed + current substrates
							Set<Vertex> union = new HashSet<Vertex>(tempSeed);
							union.addAll(substrates);
							if (union.containsAll(products)) {
								reachable = true;
								break;
							}
							
						} else {
							visited.add(reaction);
						}
					}
				}
				if (reachable)
					break;
			}
			if (reachable)
				break;
			added = substrates.addAll(tempSeed);
			visited.clear();
			tempSeed.clear();
		}
		
		if (reachable) {
			System.out.println("The products are reachable from the substrates using reactions:");
			for (Vertex reaction : scope) {
				System.out.println(reaction.getName());
			}
		} else {
			System.out.println("The products are not reachable from the substrates.");
		}
	}
	
	/**
	 * Parses a reaction equation string and collects the substrates, products, and stoichiometric coefficient.
	 * Currently not used, but could be called from the corresponding parsers. Add return construct.
	 * 
	 * @param reaction
	 */
	protected void parseReactionEquation(String reaction, MetabolicGraph graph) {
		boolean debug = false;
		boolean newToken = false, left = true;
		HashMap<Vertex,Integer> substrates = new HashMap<Vertex,Integer>();
		HashMap<Vertex,Integer> products = new HashMap<Vertex,Integer>();
		
		// parse the reaction equation
		for (int i=0, j=0; i<reaction.length(); i++) {
			boolean lastToken = false;
			int coefficient = 1;
			String compartment = null;
			int increment = 0;
			
			if (i == reaction.length()-1)
				lastToken = true;
			
			if (reaction.charAt(i) == Utilities.DELIMITER_COEFFICIENT_START) {
				// parse the stoichiometric coefficient
				coefficient = (int)Utilities.parseCoefficient(reaction.substring(i, reaction.indexOf(Utilities.DELIMITER_COEFFICIENT_END, i+1)+1));
				if (coefficient <= 0)
					coefficient = 1;
				i = reaction.indexOf(Utilities.DELIMITER_COEFFICIENT_END, i);
				j = i+1;
				
			} else if (lastToken && newToken || (increment = Utilities.isDelimiter(reaction, i)) != -1) {
				// add the previous token
				String name = reaction.substring(j, lastToken ? i+1 : i).trim();
				
				if (name.length() > 0) {
				
					if (name.contains(Utilities.DELIMITER_COMPOUND_COMPARTMENT)) {
						compartment = graph.hasCompartments() ? name.substring(name.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT)+1, name.length()) : null;
						name = name.substring(0, name.indexOf(Utilities.DELIMITER_COMPOUND_COMPARTMENT));
					}
					Vertex compound = graph.getCompound(name, compartment);
					
					if (compound == null) {
						// no decomposition: compound from the reaction is not in the database/network
						System.out.println("No decomposition: compound "+name+(compartment == null ? "":compartment)+" not in database/network.");
					}
					// add the compound to the equation
					if (left) {
						substrates.put(compound, coefficient);
					} else {
						products.put(compound, coefficient);
					}
				}
				
				// parse the equation side delimiter
				if (Utilities.isSideDelimiter(reaction, i) != -1)
					left = false; 
				
				newToken = false;
				
				if (increment > 0)
					i += (increment-1);
				
			} else if (newToken == false) {
				// beginning of a new token
				newToken = true;
				j = i;
			}
		}
		
		if (debug) {
			// print the parsed reaction
			for (Iterator<Vertex> it = substrates.keySet().iterator(); it.hasNext(); ) {
				Vertex substrate = it.next();
				System.out.print(String.valueOf(DELIMITER_COEFFICIENT_START)+(int)substrates.get(substrate)+String.valueOf(DELIMITER_COEFFICIENT_END)+" "+substrate.getName()+(substrate.getCompartment()==null?"":substrate.getCompartment())+" ");
				if (it.hasNext())
					System.out.print("+ ");
				else 
					System.out.print("--> ");
			}
			for (Iterator<Vertex> it = products.keySet().iterator(); it.hasNext(); ) {
				Vertex product = it.next();
				System.out.print(String.valueOf(DELIMITER_COEFFICIENT_START)+(int)products.get(product)+String.valueOf(DELIMITER_COEFFICIENT_END)+" "+product.getName()+(product.getCompartment()==null?"":product.getCompartment())+" ");
				if (it.hasNext())
					System.out.print("+ ");
				else 
					System.out.println();
			}
		}
	}
	
	/**
	 * Checks if there is SBML support, then validates and returns
	 * the SBMLDocument. Terminates if there is no SBML support,
	 * or if the document contains fatal errors
	 * (see {@link Utilities#validateSBML(SBMLDocument)}). 
	 * 
	 * For some reason, the SBML Model is not stored properly as
	 * static variable. Therefore, we store the SBMLDocument, which
	 * avoids multiple validations.
	 * 
	 * @param sbmlFile
	 * @return The SBMLDocument.
	 */
	public static SBMLDocument initSBML(String sbmlFile) {
		// check if we have SBML support
		if (!Utilities.sbml) {
			if (Utilities.sbmlThrowable != null)
				System.err.println(Utilities.sbmlThrowable.getMessage());
			System.err.println(Utilities.sbmlMessage);
			System.err.println("For SBML support, libSBML needs to be properly installed. See http://sbml.org/Software/libSBML/docs/cpp-api/libsbml-installation.html.");
			System.exit(-1);
		}
		
		// return the sbml doc
		if (sbmlDocument != null)
			return sbmlDocument;
		else {
			try {
				// create and validate the sbml doc
				SBMLReader reader = new SBMLReader();
				// store the SBMLDocument, as the model is not properly stored
				sbmlDocument = reader.readSBML(sbmlFile);
				Utilities.validateSBML(sbmlDocument);
				return sbmlDocument;
			} catch (UnsatisfiedLinkError e) {
				System.err.println("For SBML support, libSBML needs to be properly installed. See http://sbml.org/Software/libSBML/docs/cpp-api/libsbml-installation.html.");
				throw e;
			}
		}
	}
	
	/**
	 * Performs an internal consistency check of the SBML Document.
	 * Terminates if a fatal error exists, prints out any errors and
	 * the number of warnings.
	 * 
	 * @param sbmlDoc
	 */
	public static void validateSBML(SBMLDocument sbmlDoc) {
		
		long totalIssues = sbmlDoc.checkInternalConsistency();
		long warnings = 0, errors = 0;
		HashMap<String, ArrayList<Long>> errorMessages = new HashMap<String, ArrayList<Long>>();
		
		for (int i=0; i<sbmlDoc.getNumErrors(); i++) {
			SBMLError error = sbmlDoc.getError(i);
			if (error.isFatal()) {
				// exit on fatal error
				System.err.println("Invalid SBML file. Consistency check (org.sbml.libsbml.SBMLDocument.checkInternalConsistency()) had "+totalIssues+" error(s).");
				System.err.println("Fatal error in column "+error.getColumn()+", line "+error.getLine()+": ");
				System.err.println(error.getMessage());
				System.exit(-1);
			} else if (error.isError()) {
				// print messages of errors
				ArrayList<Long> lines = errorMessages.get(error.getMessage());
				if (lines == null)
					lines = new ArrayList<Long>();
				lines.add(error.getLine());
				errorMessages.put(error.getMessage(), lines);
				errors++;
			} else if (error.isWarning()) {
				warnings++;
			}
		}
		
		for (String message : errorMessages.keySet()) {
			System.out.print("SBML error in line(s) ");
			int lines = 1;
			for (Iterator<Long> it = errorMessages.get(message).iterator(); it.hasNext() && lines <= 3; lines++) {
				long line = it.next();
				System.out.print(line);
				if (it.hasNext() && lines<3)
					System.out.print(", ");
				else if (lines==3)
					System.out.println("... ("+(errorMessages.get(message).size()-3)+" more).");
				else
					System.out.println();
			}
			System.out.print(message);
			if (message.contains("<model> definition"))
				System.out.println("A possible cause is a read error of a mounted network drive. Try copying the XML file to the local file system.");
		}
		
		if (sbmlDoc.getModel() == null) {
			System.err.println("Invalid SBML file (no model).");
			System.exit(-1);
		}
		
		if (errors != 0 || warnings != 0)
			System.out.println("SBML file had "+errors+" errors and "+warnings+" warnings. Generated network may contain errors.");
	}
	
	/**
	 * Prints the molecular weights of all compounds in the given biocyc directory.
	 *  
	 * @param compoundsDir
	 * @param outputPath
	 * @throws IOException
	 */
	public static void writeBioCycWeights(String compoundsDir, String outputPath) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(compoundsDir+File.separator+"compounds.dat"));
		String line;
		String version = getVersion(compoundsDir+File.separator+"compounds.dat");
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath+File.separator+version+".molweights", false));
		int weights = 0;
		
		while ((line = reader.readLine()) != null) {
			
			if (line.length() >= 19 && line.substring(0, 19).equals("MOLECULAR-WEIGHT - ")) {
				writer.write(line.substring(19, line.length()).trim()+"\n");
				weights++;
			}
		}
		
		reader.close();
		writer.close();
		
		System.out.println(weights+" molecular weights written.");
	}
	
	/**
	 * Writes the hash map of compound names and mass vectors into a file.
	 * 
	 * NOTE: This method is unnecessary if the masses are read from the compounds.dat both for creating equivalence classes
	 * and creating the graph.
	 */
	protected static void writeMasses(HashMap<String, int[]> massesMap, String version, String massesFile) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(massesFile));
		
		writer.write(version+"\n");
		for (Iterator<String> compounds = massesMap.keySet().iterator(); compounds.hasNext();) {
			String compound = compounds.next();
			writer.write(compound);
			for (int i=0; i<ELEMENTS.length; i++)
				writer.write("\t"+massesMap.get(compound)[i]);
			writer.write("\n");
		}
	}
	
	protected static String getVersion(String file) throws IOException {
		String version = null;
		
		if (file.endsWith(".xml")) {
			// sbml file: extract the version from the model, if possible
			SBMLDocument sbmlDoc = Utilities.initSBML(file);
			Model sbmlModel = sbmlDoc.getModel();
			if (sbmlModel.getId() == null || sbmlModel.getId().length() == 0)
				version = file.substring(file.lastIndexOf(File.separator)+1, file.lastIndexOf("."));
			else
				version = sbmlModel.getId();
			
		} else if (file.endsWith(".dat")) {
			// biocyc file: extract the version from concatenating the Species, Database, and Version strings. 
			String line;
			BufferedReader reader = new BufferedReader(new FileReader(file));
			while ((line = reader.readLine()) != null) {
				if (line.length() >= 10 && line.substring(0,10).equals("# Species:")) {
					version = line+"\n"+reader.readLine()+"\n"+reader.readLine();
					if (version == null || !version.contains("# Database:") || !version.contains("# Version:")) {
						reader.close();
						throw new RuntimeException("Invalid Species/Database/Version string in "+file);
					}
					version = version.replace("# Species: ", "").replace("# Database: ", "").replace("# Version: ", "").replace("\n", "").replace("#", "").replace(" ", "");
					break;
				}
			}
			reader.close();
			
		} else {
			// text file: use the file prefix as version
			version = file.substring(file.lastIndexOf(File.separator)+1, file.length());
//			BufferedReader reader = new BufferedReader(new FileReader(file));
//			String line = reader.readLine();
//			reader.close();
//			if (line.contains("#"))
//				version = line.substring(line.indexOf("#")+1, line.length()).trim();
//			else
//				throw new RuntimeException("No species/version string in the first line of "+file);
		}
		
		return version;
	}
	
	/**
	 * Returns the compounds file name corresponding to the given network file name.
	 * 
	 * @param networkFile
	 * @return
	 */
	protected static String getCompoundsFile(String networkFile) {
		String compoundsFile = null;
		
		if (networkFile.endsWith(".xml")) {
			throw new IllegalArgumentException("No compounds file required for SBML files.");
			
		} else if (networkFile.endsWith(".dat")) {
			// biocyc format
			compoundsFile = networkFile.substring(0, networkFile.lastIndexOf(File.separator)+1)+"compounds.dat";
			
		} else {
			// text format
			compoundsFile = networkFile+".compounds";
		}
		
		// a compounds file is mandatory
		if (!new File(compoundsFile).exists()) {
			System.err.println("Missing compounds file: "+compoundsFile);
			System.exit(-1);
		}
		
		return compoundsFile;
	}
	
	protected static String getPathwaysFile(String networkFile) {
		String pathwaysFile = null;
		
		if (networkFile.endsWith(".xml")) {
			throw new IllegalArgumentException("No compounds file required for SBML files.");
			
		} else if (networkFile.endsWith(".dat")) {
			// biocyc format
			pathwaysFile = networkFile.substring(0, networkFile.lastIndexOf(File.separator)+1)+"pathways.dat";
			
		} else {
			// text format
			throw new IllegalArgumentException("Pathways file is only required for biocyc files (.dat).");
		}
		
		// a pathways file is optional
		if (!new File(pathwaysFile).exists()) {
			System.out.println("Missing pathways file: "+pathwaysFile+" (optional). Used for inferring unannotated reaction directions.");
			return null;
		}
		
		return pathwaysFile;
	}
	
	protected static void printEquationSystem(int[] a1, int[] a2, int coeff1, int coeff2, int[] b1, int[] b2, boolean printTrivial) {
		boolean trivial = false;
		if (!printTrivial) {
			for (int l=0; l<Utilities.ELEMENTS.length; l++) {
				if (a1[l] == 0 && a2[l] == 0 && b1[l] == 0 && b2[l] == 0)
					continue;
				// do not print trivial unsolvable systems
				if ((a1[l]+a2[l] > coeff1*b1[l]+coeff2*b2[l]) || ((coeff1*b1[l]+coeff2*b2[l]) % (a1[l]+a2[l]) != 0)) {
					trivial = true;
					break;
				}
			}
		}
		if (printTrivial || !trivial) {
			System.out.println("a[i,0]\ta[i,1]\tb[i]");
			for (int l=0; l<Utilities.ELEMENTS.length; l++)
				System.out.println(a1[l]+"\t"+a2[l]+"\t"+((coeff1*b1[l])+(coeff2*b2[l]))+" ("+coeff1+"*"+(b1[l])+" + "+coeff2+"*"+b2[l]+")");
		}
	}
	
	/**
	 * Test routines for regression tests of MetabolicGraph.solveStoichiometry(int[], int, int[]).
	 * 
	 * @return True, if all tests completed successfully, false otherwise.
	 */
	static boolean testSolveStoichiometrySingle() {
		int test = 0, randomTests = 0;
		MetabolicGraph graph = new MetabolicGraph(false, false);
		System.out.println("Testing solveStoichiometry(int[], int, int[])..."); 
		
		{
			// random tests. target mass (b) is always equal or a multiple of the source mass (a).
			Random random = new Random();
			int max = 202;
			int maxCoefficient = 11;
			for (int i=0; i<10000; i++) {
				randomTests++;
				int[] a = {random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max)};
				int coeff = random.nextInt(maxCoefficient)+1;
				int factor = random.nextInt(maxCoefficient)+1;
				int[] b = {factor*a[0], factor*a[1], factor*a[2], factor*a[3], factor*a[4], factor*a[5]};
				int solution = graph.solveStoichiometry(a, coeff, b)[0];
				if (solution != factor*coeff) {
					System.out.println("Failed random test "+randomTests);
					System.out.println("a\tb");
					for (int j=0; j<ELEMENTS.length; j++)
						System.out.println(a[i]+"\t"+b[i]);
					System.out.println(""+"solution="+solution);
					return false;
				}
			}
		}
		
		{
			// test 1: single mass vector
			test++;
			System.out.println("test "+test);
			int[] substitute = {1, 7, 0, 0, 0, 0};
			int coeff = 3;
			int[] original = {3, 21, 0, 0, 0, 0};
			int solution = graph.solveStoichiometry(substitute, coeff, original)[0];
			if (solution != 9) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		System.out.println(test+" tests and "+randomTests+" random tests completed successfully.");
		return true;
	}
	
	/**
	 * Test routines for regression tests of MetabolicGraph.solveStoichiometry(int[], int[], int, int, int[], int[]).
	 * 
	 * @return True, if all tests completed successfully, false otherwise.
	 */
	static boolean testSolveStoichiometryDouble() {
		int test = 0, randomTests = 0;
		long start = System.nanoTime();
		MetabolicGraph graph = new MetabolicGraph(false, false);
		System.out.println("Testing solveStoichiometry(int[], int[], int, int, int[], int[])"); 
		
		{
			// test 1: all ones
			test++;
			System.out.println("test "+test);
			int[] a1 = {1, 1, 1, 1, 1, 1};
			int[] a2 = {1, 1, 1, 1, 1, 1};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {1, 1, 1, 1, 1, 1};
			int[] b2 = {1, 1, 1, 1, 1, 1};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 1) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 2: regular matrix with invalid local solution
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 6, 0, 2, 0, 0};
			int[] a2 = {1, 2, 0, 1, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {0, 1, 0, 0, 0, 0};
			int[] b2 = {8, 15, 0, 6, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 2 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 3: division by zero? (AKBLIG-RXN from EcoCyc13.5)
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 1, 0, 1, 0, 0};
			int[] a2 = {8, 13, 2, 5, 0, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {0, 1, 0, 0, 0, 0};
			int[] b2 = {4, 6, 1, 3, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 4: non-zero row and 1-zero row:
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 2, 0, 0, 0, 0};
			int[] a2 = {1, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {1, 1, 0, 0, 0, 0};
			int[] b2 = {1, 0, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 5: non-zero row and 1-zero row:
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 2, 0, 0, 0, 0};
			int[] a2 = {0, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {2, 1, 0, 0, 0, 0};
			int[] b2 = {1, 0, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 6: Linearly dependent rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {2, 1, 1, 0, 0, 0};
			int[] a2 = {2, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {6, 2, 1, 0, 0, 0};
			int[] b2 = {0, 1, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 7: Linearly dependent rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {2, 1, 1, 0, 0, 0};
			int[] a2 = {2, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {8, 2, 2, 0, 0, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 2 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 8: Linearly dependent rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {2, 1, 0, 0, 0, 0};
			int[] a2 = {2, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {8, 2, 0, 0, 0, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 3) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 9: Linearly dependent rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {1, 1, 0, 0, 0, 0};
			int[] a2 = {1, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {1, 1, 0, 0, 0, 0};
			int[] b2 = {1, 1, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 1) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 10: 3 linearly dependent rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {10, 15, 5, 13, 2, 1};
			int[] a2 = {10, 15, 5, 10, 2, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {10, 14, 5, 10 ,1, 1};
			int[] b2 = {10, 16, 5, 13, 3, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 1) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 11: all linearly dependent rows
			test++;
			System.out.println("test "+test);
			int[] a1 = {10, 15, 5, 11, 2, 1};
			int[] a2 = {10, 15, 5, 11, 2, 1};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {30, 24, 15, 32 ,5, 3};
			int[] b2 = {10, 36, 5, 12, 3, 1};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			// NOTE: this solution is the first found by exhaustive search
			if (solution[0] != 1 || solution[1] != 3) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 12: all linearly dependent or zero rows
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 0, 5, 11, 2, 1};
			int[] a2 = {0, 2, 0, 11, 2, 1};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {0, 0, 15, 33 ,6, 3};
			int[] b2 = {0, 4, 0, 22, 4, 2};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			// NOTE: this solution is the first found by exhaustive search
			if (solution[0] != 3 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				System.out.println("solution[0]="+solution[0]+", solution[1]="+solution[1]);
				return false;
			}
		}
		
		{
			// test 13: last row linearly dependent
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 0, 0, 0, 0, 2};
			int[] a2 = {0, 0, 0, 0, 0, 2};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {0, 0, 0, 0, 0, 3};
			int[] b2 = {0, 0, 0, 0, 0, 5};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			// NOTE: this solution is the first found by exhaustive search
			if (solution[0] != 1 || solution[1] != 3) {
				System.out.println("Failed test "+test);
				System.out.println("solution[0]="+solution[0]+", solution[1]="+solution[1]);
				return false;
			}
		}
		
		{
			// test 14: 2 1-zero rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {2, 0, 0, 0, 0, 0};
			int[] a2 = {0, 2, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {2, 0, 0, 0, 0, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 1) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 15: 2 1-zero rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 0, 0, 0, 4, 0};
			int[] a2 = {0, 0, 0, 0, 0, 6};
			int coeff1 = 6, coeff2 = 3;
			int[] b1 = {0, 0, 0, 0, 1, 0};
			int[] b2 = {0, 0, 0, 0, 0, 3};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 3 || solution[1] != 3) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 16: 1-zero row and non-zero row:
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 0, 5, 0, 0, 0};
			int[] a2 = {3, 0, 7, 0, 0, 0};
			int coeff1 = 3, coeff2 = 4;
			int[] b1 = {6, 0, 12, 0, 0, 0};
			int[] b2 = {0, 0, 9, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 6 || solution[1] != 6) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 17: non-zero row and 1-zero row:
			test++;
			System.out.println("test "+test);
			int[] a1 = {5, 0, 0, 0, 0, 0};
			int[] a2 = {7, 0, 0, 0, 0, 3};
			int coeff1 = 3, coeff2 = 4;
			int[] b1 = {12, 0, 0, 0, 0, 6};
			int[] b2 = {9, 0, 0, 0, 0, 0};	
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 6 || solution[1] != 6) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 18: 2 non-zero rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 1, 0, 0, 2, 0};
			int[] a2 = {2, 1, 0, 0, 0, 1};
			int coeff1 = 2, coeff2 = 1;
			int[] b1 = {6, 1, 0, 0, 1, 0};
			int[] b2 = {1, 4, 0, 0, 0, 1};	
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 19: 2 non-zero rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 0, 0, 0, 0, 0};
			int[] a2 = {0, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {6, 0, 0, 0, 0, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};	
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 2 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 20: 2 non-zero rows:
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 1, 0, 0, 0, 0};
			int[] a2 = {0, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {6, 0, 0, 0, 0, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 21
			test++;
			System.out.println("test "+test);
			int[] a1 = {3, 2, 0, 0, 0, 0};
			int[] a2 = {0, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {2, 1, 0, 0, 0, 0};
			int[] b2 = {1, 0, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 22: only one non-zero row
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 2, 0, 0, 0, 0};
			int[] a2 = {0, 1, 0, 0, 0, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {0, 2, 0, 0, 0, 0};
			int[] b2 = {0, 1, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 23: high numbers cause precision problems (RXN0-2001 from EcoCyc13.5)
			test++;
			System.out.println("test "+test);
			int[] a1 = {110, 196, 2, 39, 2, 0};
			int[] a2 = {60, 100, 1, 7, 1, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {110, 196, 2, 39, 2, 0};
			int[] b2 = {60, 100, 1, 7, 1, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 1 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 24: non-zero row enclosed by linearly dependent rows
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 2, 1, 2, 0, 0};
			int[] a2 = {0, 2, 0, 2, 0, 0};
			int coeff1 = 1, coeff2 = 1;
			int[] b1 = {0, 4, 2, 4, 0, 0};
			int[] b2 = {0, 4, 0, 4, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution[0] != 2 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				System.out.println("solution[0]="+solution[0]+", solution[1]="+solution[1]);
				return false;
			}
		}
		
		{
			// test 25: unsolvable even with non-integers
			test++;
			System.out.println("test "+test);
			int[] a1 = {9, 13, 3, 4, 0, 0};
			int[] a2 = {6, 10, 0, 12, 2, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {15, 21, 3, 16, 2, 0};
			int[] b2 = {0, 2, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 26: unsolvable even with non-integers
			test++;
			System.out.println("test "+test);
			int[] a1 = {6, 10, 0, 6, 0, 0};
			int[] a2 = {2, 6, 0, 2, 0, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {1, 1, 0, 1, 0, 0};
			int[] b2 = {0, 1, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 27: unsolvable even with non-integers
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 0, 0, 4, 1, 0};
			int[] a2 = {6, 12, 0, 5, 0, 0};
			int coeff1 = 1, coeff2 = 2;
			int[] b1 = {6, 11, 0, 9, 1, 0};
			int[] b2 = {0, 1, 0, 0, 0, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution != null) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 28: reversible substitution - 1
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 1, 0, 0, 0, 0};
			int[] a2 = {21, 26, 7, 14, 2, 0};
			int coeff1 = 3, coeff2 = 2;
			int[] b1 = {0, 1, 0, 0, 0, 0};
			int[] b2 = {21, 27, 7, 14, 2, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution == null || solution[0] != 5 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// test 29: reversible substitution - 2
			test++;
			System.out.println("test "+test);
			int[] a1 = {0, 1, 0, 0, 0, 0};
			int[] a2 = {21, 27, 7, 14, 2, 0};
			int coeff1 = 5, coeff2 = 2;
			int[] b1 = {0, 1, 0, 0, 0, 0};
			int[] b2 = {21, 26, 7, 14, 2, 0};
			int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
			if (solution == null || solution[0] != 3 || solution[1] != 2) {
				System.out.println("Failed test "+test);
				return false;
			}
		}
		
		{
			// random tests.
			Random random = new Random();
			int max = 202;
			int maxCoefficient = 11;
			for (int i=0; i<10000; i++) {
				randomTests++;
				int[] a1 = {random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max)};
				int[] a2 = {random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max), random.nextInt(max)};
				Integer[] a1norm = new Integer[ELEMENTS.length], a2norm = new Integer[ELEMENTS.length];
				Equivalence.normalize(a1).toArray(a1norm);
				Equivalence.normalize(a2).toArray(a2norm);
				int coeff1 = random.nextInt(maxCoefficient)+1, coeff2 = random.nextInt(maxCoefficient)+1;
				int factor1 = random.nextInt(maxCoefficient)+1, factor2 = random.nextInt(maxCoefficient)+1;
				int[] b1 = {factor1*a1norm[0], factor1*a1norm[1], factor1*a1norm[2], factor1*a1norm[3], factor1*a1norm[4], factor1*a1norm[5]};
				int[] b2 = {factor2*a2norm[0], factor2*a2norm[1], factor2*a2norm[2], factor2*a2norm[3], factor2*a2norm[4], factor2*a2norm[5]};
				int[] solution = graph.solveStoichiometry(a1, a2, coeff1, coeff2, b1, b2, false);
				for (int j=0; j<ELEMENTS.length; j++) {
					if (solution[0]*a1[j] != solution[2]*coeff1*b1[j] || solution[1]*a2[j] != solution[2]*coeff2*b2[j]) {
						System.out.println("Failed random test "+randomTests);
						printEquationSystem(a1, a2, coeff1, coeff2, b1, b2, true);
						System.out.println("solution[0]="+solution[0]+", solution[1]="+solution[1]+", factor1="+factor1+", factor2="+factor2);
						return false;
					}
				}
			}
		}
		
		System.out.println(test+" tests and "+randomTests+" random tests completed successfully.");
		System.out.println("Runtime: "+(System.nanoTime()-start)/1000000+" ms.");
		return true;
	}
	
	/**
	 * The following static block is needed in order to load the
	 * libSBML Java interface library when the application starts.
	 * If it fails, SMBL is set to false. The error messages
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
