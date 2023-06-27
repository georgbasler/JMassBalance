/**
 * 
 */
package massbalance;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.StringTokenizer;

import org.sbml.libsbml.Model;
import org.sbml.libsbml.SBMLDocument;

/**
 * Class for generating, reading, and writing the mass equivalence classes required for mass balanced
 * randomization.
 * 
 * 
 * @author Georg Basler
 * 
 */
public class Equivalence {

	static final String PARAMETERS = "Parameters required: <networkFile> <outputDir> [-check[-exhaustive]]";
	
	static final boolean DOUBLE = true;
	
	// global classes and masses are needed for calling from Randomize.java
	static HashMap<ArrayList<Integer>, EquivalenceClass> classes;
	static HashMap<String, int[]> massesMap;

	/**
	 * Static initializer for loading the configuration parameters.
	 */
	static {
		Utilities.loadConfiguration();
	}
	
	/**
	 * <p>
	 * Parses the mass vectors and generates the mass equivalence classes
	 * from the given files. Is called implicitly by {@link Randomize#main(String[])}.
	 * 
	 * <p>
     * <strong>Usage:</strong><br>
     * <code>networkFile outputDir</code><br>
	 * <br>
	 * 
	 * @param args
	 */
	public static void main(String[] args) throws InterruptedException {
	      
		if (args.length < 2) {
			System.out.println(PARAMETERS);
			System.exit(-1);
		}

		String networkFile = args[0];
		String outputDir = args[1];
		
		try {
			if (!new File(networkFile).exists() || new File(networkFile).isDirectory()) {
				System.err.println("File not found: "+networkFile);
				System.out.println(PARAMETERS);
				System.exit(-1);
			}
			if (!new File(outputDir).isDirectory())
				new File(outputDir).mkdir();
			
			String version = Utilities.getVersion(networkFile);
			String infoFile = outputDir+File.separator+version+".equivalence.log";
			new File(infoFile).delete();
			
			// check the format of the network file
			if (networkFile.endsWith(".xml")) {
				// SBML format
				SBMLDocument sbmlDoc = Utilities.initSBML(networkFile);
				Model sbmlModel = sbmlDoc.getModel();
				massesMap = Utilities.createMasses(sbmlModel, infoFile);
				
			} else if (networkFile.endsWith(".dat")) {
				// biocyc format
				massesMap = Utilities.createMassesBioCyc(Utilities.getCompoundsFile(networkFile), networkFile, infoFile);
				
			} else {
				// text format
				massesMap = Utilities.createMasses(networkFile, Utilities.getCompoundsFile(networkFile), infoFile);
			}
			
			// create the equivalence classes
			String classesFile = outputDir+File.separator+version+".classes";
			System.out.print("Creating equivalence classes: "+classesFile+". ");
			long time = System.currentTimeMillis();
			classes = createClasses(massesMap);
			
			// compare to other equivalence classes
//			HashMap<ArrayList<Integer>, EquivalenceClass> classes2 = Utilities.readClasses(classesFile+".old", numberOfCompounds);
//			System.out.println(Utilities.compareClasses(classes, classes2));
			
			if (args.length > 2 && args[2].startsWith("-check"))
				consistency(classes, massesMap.size(), args[2].equals("-check-exhaustive"));
			
			// save the classes to file
			writeClasses(classes, version, classesFile, infoFile);
			
			long elapsed = System.currentTimeMillis()-time;
			System.out.println(" ("+((float)elapsed/1000f) + " seconds).");
			System.out.println("See "+infoFile+" for more information.");
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Generates the mass equivalence classes from a given map of compound names and mass vectors. 
	 * 
	 * @param massesMap HashMap of compound names and their mass vectors.
	 * @return HashMap containing the basis vectors of equivalence classes as keys and EquivalenceClass as values.
	 */
	public static HashMap<ArrayList<Integer>, EquivalenceClass> createClasses(HashMap<String, int[]> massesMap) throws InterruptedException {
		// Choosing an initial hash map capacity such that (initialCapacity > numberOfClasses * loadFactor) avoids all rehash operations
		int initialCapacity = (int)(1.1*Utilities.RATIO_CLASSES_COMPOUNDS*massesMap.size());
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = new HashMap<ArrayList<Integer>, EquivalenceClass>((int)(initialCapacity/0.75f+1), 0.75f);
		long counter = 0;
		
		for (Iterator<String> compounds = massesMap.keySet().iterator(); compounds.hasNext();) {
			String compound = compounds.next();
			int[] mass = massesMap.get(compound);
			ArrayList<Integer> newKey = normalize(mass);
			EquivalenceClass eqClass = classes.get(newKey);
			
			if (eqClass == null) {
				eqClass = new EquivalenceClass(compound, mass);
				classes.put(newKey, eqClass);
			} else {
				eqClass.addCompound(compound, mass);
			}
	
			HashMap<ArrayList<Integer>, EquivalenceClass> temp = new HashMap<ArrayList<Integer>, EquivalenceClass>((int)(initialCapacity/0.75f+1), 0.75f);
			Set<ArrayList<Integer>> keys = classes.keySet();
			for (Iterator<ArrayList<Integer>> keysIterator = keys.iterator(); keysIterator.hasNext(); ) {
				
				ArrayList<Integer> key = keysIterator.next();
				
				// calculate sum mass with every compound from a different equivalence class
				// (pairs with both compounds from the same class are considered implicitly) 
				if (!key.equals(newKey)) {
					eqClass = classes.get(key);
					ArrayList<int[]> masses = eqClass.getMasses();
					
					int i=0;
					for (Iterator<int[]> massesIterator = masses.iterator(); massesIterator.hasNext(); i++) {
						int[] sum = sumMasses(mass, massesIterator.next());
						ArrayList<Integer> sumKey = normalize(sum);
						EquivalenceClass doubleClass = classes.get(sumKey);
						
						if (doubleClass == null) {
							// new class: check if the temp hash map contains it
							doubleClass = temp.get(sumKey);
							if (doubleClass == null) {
								// create new class and add it to temp hash map.
								// the temp hash map is necessary to avoid ConcurrentModificationException when iterating and adding. 
								doubleClass = new EquivalenceClass(compound, mass, eqClass.getName(i), eqClass.getMass(i));
								temp.put(sumKey, doubleClass);
							} else {
								doubleClass.addPair(compound, mass, eqClass.getName(i), eqClass.getMass(i));
							}
							// classes.put(sumKey, doubleClass); ConcurrentModificationException due to modification while iteration
						} else {
							doubleClass.addPair(compound, mass, eqClass.getName(i), eqClass.getMass(i));
						}
					}
				}
				
				if (counter % 100 == 0) {
					// allow to interrupt this thread
					Thread t = Thread.currentThread();
					synchronized(t) {
						t.wait(1);
					}
				}
				counter++;
			}
			
			// add the new, temporarily stored equivalence classes
			classes.putAll(temp);
		}
		
		if (initialCapacity < classes.size())
			System.out.print("Initial classes hash capacity "+initialCapacity+" smaller than final map size ("+classes.size()+"). Increase RATIO_CLASSES_COMPOUNDS to at least "+(int)(classes.size()/(1.1*massesMap.size())+1)+" for better performance. ");
		
		return classes;
	}
	
	
	/**
	 * Reads the equivalence classes from the classesFile and creates a hash map of the normalized mass vector as key
	 * and the EquivalenceClass as value. Single compounds of a class are expected in odd line numbers, pairs of compounds in
	 * even line numbers (lines prefixed with # are regarded comments).
	 * 
	 * @param classesFile
	 * @return HashMap containing the basis vectors of equivalence classes as keys and EquivalenceClass as values.
	 * @throws IOException
	 */
	public static HashMap<ArrayList<Integer>, EquivalenceClass> readClasses(String classesFile) throws IOException, InterruptedException {
		BufferedReader reader = new BufferedReader(new FileReader(classesFile));
		HashMap<ArrayList<Integer>, EquivalenceClass> classes = new HashMap<ArrayList<Integer>, EquivalenceClass>((int)(Utilities.DEFAULT_NUMBER_OF_CLASSES/0.75f), 0.75f);
		int counter = 0, lineNumber = 0;
		String line;
		boolean single = true; // single or double compound row  
		boolean odd = false; // first or second of a pair of compounds
		ArrayList<Integer> key = null;
		EquivalenceClass eqClass = null;
		
		// give some time to interrupt this thread
		Thread t = Thread.currentThread();
		synchronized(t) {
			t.wait(1);
		}
		
		while ((line = reader.readLine()) != null) {
			lineNumber++;
			
			if (line.startsWith("#"))
				continue;
			
			int pos = 0;
			// number of considered elements +1
			int elements = Utilities.ELEMENTS.length+1;
			int[] mass = null;
			ArrayList<String> names = new ArrayList<String>();
			ArrayList<int[]> masses = new ArrayList<int[]>();
			StringTokenizer tokenizer = new StringTokenizer(line, "\t");
			
			if (single) {
				key = null;
				eqClass = new EquivalenceClass();
			}
			
			while (tokenizer.hasMoreTokens()) {
					
				if (pos % elements == 0) {
					// every n+1st token: compound name
					names.add(tokenizer.nextToken());
					if (!single)
						odd = !odd;
					mass = new int[Utilities.ELEMENTS.length];
				} else {
					// other tokens: masses
					try {
						mass[(pos%elements)-1] = Integer.parseInt(tokenizer.nextToken());
					} catch (NumberFormatException e) {
						System.err.println("Number of tokens in classes file (line "+lineNumber+") does not match the specified number of atoms ("+Utilities.ELEMENTS.length+"). Make sure to use the same elements in the config file that were used for generating the equivalence classes.");
						System.exit(-1);
					}
				}
				
				if (pos % elements == Utilities.ELEMENTS.length) {
					// finished a compound: add mass
					if (single)
						counter++;
					masses.add(mass);
					
					if (key == null) {
						// determine key for single or double compounds
						if (single)
							key = Equivalence.normalize(mass);
						else if (!odd)
							key = Equivalence.normalize(sumMasses(masses.get(0), masses.get(1)));
					}
				}
				pos++;
			}
			
			// free some memory
			tokenizer = null;
			mass = null;
			
			if (pos % elements != 0)
				throw new RuntimeException("Number of tokens in classes file (line "+lineNumber+") does not match the specified number of atoms ("+Utilities.ELEMENTS.length+"). Make sure to use the same elements in the config file that were used for generating the equivalence classes.");
			
			// add compounds to class
			if (pos != 0) {
//				names.trimToSize();
//				masses.trimToSize();
				if (single)
					eqClass.addAllCompounds(names, masses);
				else
					eqClass.addAllPairs(names, masses);
			}
			
			// add class to hash map
			if (!single) {
				if (key == null)
					throw new RuntimeException("No key generated for equivalence class: " + classes.size());
//				key.trimToSize();
				classes.put(key, eqClass);
			} 
			
			single = !single;
		}
		
		// give some time to interrupt this thread
		synchronized(t) {
			t.wait(1);
		}
		
		if (Utilities.DEFAULT_NUMBER_OF_CLASSES < classes.size())
			System.out.print("Initial classes hash capacity ("+Utilities.DEFAULT_NUMBER_OF_CLASSES+") smaller than final map size ("+classes.size()+"). Increase DEFAULT_NUMBER_OF_CLASSES ("+Utilities.DEFAULT_NUMBER_OF_CLASSES+") for better performance.");
		
		
		System.out.print(counter + " annotated compounds in " + classes.size() + " classes");
		reader.close();
		
		return classes; 
	}
	
	/**
	 * Extracts the compound names and masses from the equivalence classes and returns them
	 * in a hash map.
	 * 
	 * @param classes HashMap of basis mass vectors and their equivalence classes.
	 * @return HashMap of compound names and their mass vectors.
	 */
	public static HashMap<String,int[]> getMasses(HashMap<ArrayList<Integer>, EquivalenceClass> classes) {
		int initialCapacity = (int)(1.1*classes.size()/Utilities.RATIO_CLASSES_COMPOUNDS);
		HashMap<String, int[]> massesMap = new HashMap<String, int[]>((int)(initialCapacity/0.75f), 0.75f);
		
		for (Iterator<ArrayList<Integer>> keys = classes.keySet().iterator(); keys.hasNext();) {
			EquivalenceClass eqClass = classes.get(keys.next());
			for (int i=0; i<eqClass.getNames().size(); i++)
				massesMap.put(eqClass.getName(i), eqClass.getMass(i));
		}
		
		if (initialCapacity < massesMap.size())
			System.out.println("Initial masses hash capacity ("+initialCapacity+") smaller than final map size ("+massesMap.size()+"). Decrease RATIO_CLASSES_COMPOUNDS ("+Utilities.RATIO_CLASSES_COMPOUNDS+") to at most "+(int)(1.1*classes.size()/massesMap.size())+" for better performance.");
		
		return massesMap;
	}
	
	/**
	 * Output equivalence classes as table. An equivalence class spans two lines: the odd line contains the single compounds,
	 * the even line contains pairs of compounds of the class (not counting the initial version string line). Each compound
	 * is printed as tab-separated list of the string identifier followed by the mass elements:
		[class 1: single compounds] id[1] C[1] H[1] N[1] O[1] P[1] S[1] ... id[n] C[n] H[n] N[n] O[n] P[n] S[n]
		[class 1: double compounds] id[1-1] C[1-1] ... S[1-1] id[1-2] C[1-2] ... S[1-2] ... id[1-m] C[1-m] ... S[1-m]
		[class 2: single compounds] ...
	 * @param classes
	 * @param version
	 * @param classesFile
	 * @throws IOException
	 */
	public static void writeClasses(HashMap<ArrayList<Integer>, EquivalenceClass> classes, String version, String classesFile, String infoFile) throws IOException {
		int compounds = 0, sum = 0;
		BufferedWriter writer = new BufferedWriter(new FileWriter(classesFile));
		BufferedWriter singleSizesWriter = new BufferedWriter(new FileWriter(classesFile+".single.sizes"));
		BufferedWriter doubleSizesWriter = new BufferedWriter(new FileWriter(classesFile+".double.sizes"));
		BufferedWriter infoWriter = new BufferedWriter(new FileWriter(infoFile, new File(infoFile).exists()));
		
		if (!version.startsWith("# "))
			version = "# "+version;
		writer.write(version+"\n");
		
		Set<ArrayList<Integer>> keys = classes.keySet(); 
		for (Iterator<ArrayList<Integer>> keysIterator = keys.iterator(); keysIterator.hasNext(); ) {			
			int i=0, singleSize = 0, doubleSize = 0;
			EquivalenceClass eqClass = classes.get(keysIterator.next());
			
			for (Iterator<String> namesIterator = eqClass.getNames().iterator(); namesIterator.hasNext(); i++) {
				// print single compound name
				String name = namesIterator.next();
				writer.write(name + "\t");
				// print single compound mass
				int[] mass = eqClass.getMass(i);
				for (int j=0; j<mass.length; j++) {
					writer.write(mass[j] + "\t");
				}
				compounds++;
				singleSize++;
			}
			
			writer.write("\n");
			
			i=0;
			for (Iterator<String> doubleNamesIterator = eqClass.getExplicitDoubleNames().iterator(); doubleNamesIterator.hasNext(); i++) {
				
				// print double compound name
				writer.write(doubleNamesIterator.next() + "\t");
				
				// print double compound mass
				int[] mass = eqClass.getDoubleMass(i);
				for (int j=0; j<mass.length; j++) {
					writer.write(mass[j] + "\t");
				}
				
				if (i % 2 == 0)
					doubleSize++;
			}
			
			singleSizesWriter.write(String.valueOf(singleSize));
			doubleSizesWriter.write(String.valueOf(doubleSize));
			if (keysIterator.hasNext()) {
				singleSizesWriter.write("\t");
				doubleSizesWriter.write("\t");
			}
			writer.write("\n");
			sum += (singleSize + doubleSize);
		}
		
		String message = compounds + " annotated compounds in " + keys.size() + " classes, average of "+((float)sum/(float)keys.size())+" entries ";
		System.out.print(message);
		infoWriter.write(message+"\n");
		
		writer.close();
		singleSizesWriter.write("\n");
		doubleSizesWriter.write("\n");
		singleSizesWriter.close();
		doubleSizesWriter.close();
		infoWriter.close();
	}
	
	/**
	 * Compares two sets of equivalence classes for equality of the contained compounds. First
	 * compares the number of classes, then iterates over the keys of the first set of classes
	 * and checks whether the second set contains the same keys with the same compounds and
	 * identical mass vectors.  
	 * 
	 * @param classes1
	 * @param classes2
	 * @return True, if both sets of classes contain the same compounds, false otherwise. 
	 */
	public static boolean compareClasses(HashMap<ArrayList<Integer>, EquivalenceClass> classes1, HashMap<ArrayList<Integer>, EquivalenceClass> classes2) {
		
		if (classes1.keySet().size() != classes2.keySet().size() || classes1.values().size() != classes2.values().size()) {
			System.out.println("Unequal number of classes.");
			return false;
		}
		
		// compare each class
		for (Iterator<ArrayList<Integer>> keys = classes1.keySet().iterator(); keys.hasNext();) {
			ArrayList<Integer> key = keys.next();
			EquivalenceClass eqClass1 = classes1.get(key);
			EquivalenceClass eqClass2 = classes2.get(key);
			
			if (eqClass1 == null || eqClass2 == null) {
				System.out.println("Equivalence class in classes1, but not in classes2.");
				return false;
			}
			
			if (eqClass1.getNames().size() != eqClass2.getNames().size() || eqClass1.getExplicitDoubleNames().size() != eqClass2.getExplicitDoubleNames().size()) {
				System.out.println("Unequal number of compounds in equivalence classes.");
				return false;
			}
			
			// check if every compound from class1 is contained identically in class2 
			for (int i=0; i<eqClass1.getNames().size(); i++) {
				String compound1 = eqClass1.getNames().get(i);
				if (!eqClass2.getNames().contains(compound1)) {
					System.out.println("Different single compounds in classes1 and classes2.");
					return false;
				}
				for (int j=0; j<Utilities.ELEMENTS.length; j++) {
					int[] mass = eqClass2.getMass(eqClass2.getNames().indexOf(compound1));
					if (eqClass1.getMass(i)[j] != mass[j]) {
						System.out.println("Different single masses in classes1 and classes2.");
						return false;
					}
				}
			}
			
			// check if every pair of compounds from class1 is contained identically in class2
			for (int i=0; i<eqClass1.getExplicitDoubleNames().size(); i++) {
				String compound1 = eqClass1.getExplicitDoubleNames().get(i);
				if (!eqClass2.getExplicitDoubleNames().contains(compound1)) {
					System.out.println("Different double compounds in classes1 and classes2.");
					return false;
				}
				for (int j=0; j<Utilities.ELEMENTS.length; j++) {
					int[] mass = eqClass2.getDoubleMass(eqClass2.getExplicitDoubleNames().indexOf(compound1));
					if (eqClass1.getDoubleMass(i)[j] != mass[j]) {
						System.out.println("Different double masses in classes1 and classes2.");
						return false;
					}
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Performs a consistency check of the equivalence classes. The following tests are included:
	 * 
	 *  - test every class for non-emptiness
	 *  - validate every key for being a vector of length ELEMENTS.length
	 *  - test equivalence of all single and double compounds in each class
	 *  - if exhaustive==true, every class is also tested for equivalence with any other class.
	 * 
	 * @param classes Equivalence classes to check.
	 * @param exhaustive Indicates to also test equivalences between classes.
	 * @throws RuntimeException if an inconsistency as defined above was found.
	 */
	public static void consistency(HashMap<ArrayList<Integer>, EquivalenceClass> classes, int compounds, boolean exhaustive) throws RuntimeException {
		int counter = 0, decimals = 1;
		ArrayList<int[]> referenceMasses = new ArrayList<int[]>();
		
		if (exhaustive)
			System.out.print("Exhaustively checking consistency (O(classes^2))...");
		else 
			System.out.print("Checking consistency...");
		
		Set<ArrayList<Integer>> keys = classes.keySet();
		
		if (keys.size() < 1)
			throw new RuntimeException("No equivalence class.");
			
		for (Iterator<ArrayList<Integer>> keysIterator = keys.iterator(); keysIterator.hasNext(); ) {
			boolean singles = false, doubles = false;
			
			ArrayList<Integer> key = keysIterator.next();
			
			if (key.size() != Utilities.ELEMENTS.length)
				throw new RuntimeException("Invalid key size: "+ key.size());
				
			EquivalenceClass eqClass = classes.get(key);
			
			// create string from the key
			String keyString = "";
			for (Iterator<Integer> it = key.iterator(); it.hasNext(); )
				keyString += " " + it.next();
			
			// take first mass array as reference
			int[] mass = new int[Utilities.ELEMENTS.length], doubleMass = new int[Utilities.ELEMENTS.length];
			String name = "", doubleName = "";
			if (eqClass.getMasses().size() > 0) {
				mass = eqClass.getMass(0);
				name = eqClass.getName(0);
				singles = true;
			}
			if (eqClass.getExplicitDoubleMasses().size() > 0) {
				doubleMass = sumMasses(eqClass.getDoubleMass(0), eqClass.getDoubleMass(1));
				doubleName = eqClass.getDoubleName(0) + " + " + eqClass.getDoubleName(1);
				doubles = true;					
			}

			if (eqClass.getMasses().size() < 1 && eqClass.getExplicitDoubleMasses().size() < 1)
				throw new RuntimeException("Empty equivalence class:" + keyString);
			
			if (exhaustive) {
				// check this class for equivalence with a previous class
				for (Iterator<int[]> ref = referenceMasses.iterator(); ref.hasNext(); ) {
					int[] referenceMass = ref.next();
					if (singles && isEquivalent(referenceMass, mass)) {
						String refString = "";
						for (int i=0; i<Utilities.ELEMENTS.length; i++)
							refString += " " + referenceMass[i];
						throw new RuntimeException("Class " + keyString + " is equivalent to class" + refString);
					}
					if (doubles && isEquivalent(referenceMass, doubleMass)) {
						String refString = "";
						for (int i=0; i<Utilities.ELEMENTS.length; i++)
							refString += " " + referenceMass[i];
						throw new RuntimeException("Class " + keyString + " is equivalent to class" + refString);
					}
				}
			}
			
			// check equivalence of single and double compounds
			if (singles && doubles && !isEquivalent(mass, doubleMass))
				throw new RuntimeException("Invalid compounds '" + name + "' and '" + doubleName + "' in class:" + keyString);
			
			// keep a reference compound from every class for checking class-equivalences
			if (exhaustive) {
				if (singles)
					referenceMasses.add(mass);
				else
					referenceMasses.add(doubleMass);
			}
			
			// check size of doubleArrays for evenness
			if ((eqClass.getExplicitDoubleNames().size() % 2) != 0 || (eqClass.getExplicitDoubleMasses().size() % 2) != 0)
				throw new RuntimeException("Odd number of double compounds in class:" + keyString);
			
			// check single compounds for equivalence
			int i=0;
			for (Iterator<int[]> massIterator = eqClass.getMasses().iterator(); massIterator.hasNext(); i++) {
				if (!isEquivalent(mass, massIterator.next()))
					throw new RuntimeException("Invalid compounds '" + name + "' and '" + eqClass.getName(i) + "' in class:" + keyString);
			}
			
			// check double compounds for equivalence
			i=0;
			for (Iterator<int[]> doubleMassIterator = eqClass.getExplicitDoubleMasses().iterator(); doubleMassIterator.hasNext(); i += 2) {
				int[] sum = sumMasses(doubleMassIterator.next(), doubleMassIterator.next());
				if (!isEquivalent(doubleMass, sum)) {					
					throw new RuntimeException("Invalid compounds '" + doubleName + "' and '" + eqClass.getDoubleName(i) + " + " + eqClass.getDoubleName(i+1) + "' in class:" + keyString);
				}
			}
			
			counter++;
			if (counter % (keys.size()/10) == 0)
				System.out.print((decimals++)*10 + "%");
			if (counter % compounds == 0)
				System.out.print(".");
		}
		System.out.println(" finished.");
	}
	
	/**
	 * Normalize mass[ELEMENTS.length] by its greatest common divisor.
	 * 
	 * Returns an array containing the normalized mass vector (basis vector of the equivalence class).
	 */
	static ArrayList<Integer> normalize(int[] mass) {
		if (mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + mass.length + ".");
		
		ArrayList<Integer> norm = new ArrayList<Integer>(Utilities.ELEMENTS.length);
		int m = mass[0];
		// calculate the greatest common divisor of x[ELEMENTS.length]
		for (int i=0; i<Utilities.ELEMENTS.length-1; i++) {
			int n = mass[i+1];
			int r;
			while (n != 0) {
				r = m % n;
				m = n;
				n = r;
			}
		}

		// normalize mass vector by the gcd
		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			norm.add(mass[i] / m);
		}
		
		return norm;
	}
	
	
	/**
	 * Compares two mass vectors for mass equivalence. Used for checking the consistency of equivalence classes.
	 * @param mass1
	 * @param mass2
	 * @return
	 */
	private static boolean isEquivalent(int[] mass1, int[] mass2) {
		if ((mass1 != null && mass1.length != Utilities.ELEMENTS.length) || (mass2 != null && mass2.length != Utilities.ELEMENTS.length))
			throw new IllegalArgumentException("Mass vector has invalid size: " + ((mass1.length != Utilities.ELEMENTS.length)? mass1.length : mass2.length) + ".");
		
		boolean first = true;
		float fraction = -1;

		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			if (mass1[i] == 0 && mass2[i] == 0)
				continue;
			if (mass1[i] == 0 ^ mass2[i] == 0)
				return false;
			if (first) {
				fraction = (float)mass1[i] / (float)mass2[i];
				first = false;
			} else {
				if (fraction != (float)mass1[i] / (float)mass2[i])
					return false;
			}
		}
		return true;
	}
	
	static int[] sumMasses(int[] mass1, int[] mass2) {
		if ((mass1 != null && mass1.length != Utilities.ELEMENTS.length) || (mass2 != null && mass2.length != Utilities.ELEMENTS.length))
			throw new IllegalArgumentException("Mass vector has invalid size: " + ((mass1.length != Utilities.ELEMENTS.length)? mass1.length : mass2.length) + ".");
		
		int sum[] = new int[Utilities.ELEMENTS.length];
		
		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			sum[i] = mass1[i] + mass2[i];
		}
		
		return sum;
	}

	static int[] multiplyMass(int[] mass, int factor) {
		if (mass != null && mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + mass.length + ".");

		if (factor < 0)
			throw new IllegalArgumentException("Negative factor: " + factor + ".");
		
		int result[] = new int[Utilities.ELEMENTS.length];
		
		for (int i=0; i<Utilities.ELEMENTS.length; i++) {
			result[i] = mass[i] * factor;
		}
		
		return result;	
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
