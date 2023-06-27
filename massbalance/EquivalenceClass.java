/**
 * 
 */
package massbalance;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class representing the mass equivalence classes used for mass-balanced randomization.
 * 
 * @author Georg Basler
 *
 */
public class EquivalenceClass {
	
	// single compound variables
	/*
	 * The ArrayLists for names and masses are synchronous at all times, 
	 * so that any index refers to the same compound in both ArrayLists.
	 */
	private ArrayList<String> names;
	private ArrayList<int[]> masses;
	private ArrayList<String> doubleNames;
	private ArrayList<int[]> doubleMasses;
	
	/**
	 * Creates a new, empty equivalence class.
	 *  
	 */
	public EquivalenceClass() {
		
		names = new ArrayList<String>();
		masses = new ArrayList<int[]>();
		doubleNames = new ArrayList<String>();
		doubleMasses = new ArrayList<int[]>();
		
	}
	
	/**
	 * Creates a new equivalence class containing one compound with the given name and mass.  
	 * 
	 * @param name The name of the compound.
	 * @param mass The ELEMENTS.length mass vector of the compound. 
	 */
	public EquivalenceClass(String name, int[] mass) {
		if (mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + mass.length + ".");
		
		names = new ArrayList<String>();
		masses = new ArrayList<int[]>();
		doubleNames = new ArrayList<String>();
		doubleMasses = new ArrayList<int[]>();
		
		this.addCompound(name, mass);
	}
	
	/**
	 * Creates a new equivalence class containing a pair of compounds with the given names and masses. 
	 * 
	 * @param name1 The name of the compound.
	 * @param mass1 The ELEMENTS.length mass vector of the compound.
	 * @param name2 The name of the compound.
	 * @param mass2 The ELEMENTS.length mass vector of the compound. 
	 */
	public EquivalenceClass(String name1, int[] mass1, String name2, int[] mass2) {
		if (mass1.length != Utilities.ELEMENTS.length || mass2.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + ((mass1.length != Utilities.ELEMENTS.length)? mass1.length : mass2.length) + ".");
		
		names = new ArrayList<String>();
		masses = new ArrayList<int[]>();
		doubleNames = new ArrayList<String>();
		doubleMasses = new ArrayList<int[]>();
		
		this.addPair(name1, mass1, name2, mass2);
	}
	
	/**
	 * Adds a compound name and its mass to the mass equivalence class.
	 * The mass vector must be linearly dependent on the existing mass vectors of the class.
	 * 
	 * @param name
	 * @param mass
	 */
	public void addCompound(String name, int[] mass) {
		if (mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + mass.length + ".");
		
		// check for linear dependence
		if (masses.size() > 0) {
			if (!Equivalence.normalize(masses.get(0)).equals(Equivalence.normalize(mass)))
				throw new IllegalArgumentException("Mass vector of "+name+" is not linearly dependent on "+names.get(0)+".");
		} else if (doubleMasses.size() > 0) {
			if (!Equivalence.normalize(Equivalence.sumMasses(doubleMasses.get(0), doubleMasses.get(1))).equals(Equivalence.normalize(mass)))
				throw new IllegalArgumentException("Mass vector of "+name+" is not linearly dependent on "+doubleNames.get(0)+" + "+doubleNames.get(1)+".");
		}
		
		names.add(name);
		masses.add(mass);
		
		if (Utilities.DEBUG) {
			ArrayList<Integer> norm0 = Equivalence.normalize(masses.get(0));
			ArrayList<Integer> norm1 = Equivalence.normalize(mass);
			for (int i=0; i<Utilities.ELEMENTS.length; i++) {
				if (!norm0.get(i).equals(norm1.get(i))) {
					System.out.println("Non-equivalent compounds: "+ names.get(0) + " != " + name);
					System.out.println("norm0: ");
					for (int j=0; j<Utilities.ELEMENTS.length; j++) {
						System.out.print(norm0.get(j) + " ");
					}
					System.out.println();
					System.out.println("norm1: ");
					for (int j=0; j<Utilities.ELEMENTS.length; j++) {
						System.out.print(norm1.get(j) + " ");
					}
					System.out.println();
//					System.out.println("Key 0: " + Equivalence.getKey(masses.get(0)) + ", Key 1: " + Equivalence.getKey(mass));
				}
			}
		}
	}
	
	/**
	 * Adds a pair of compound names and masses to the mass equivalence class.
	 * The sum of the given mass vectors must be linearly dependent on the existing mass vectors of the class.
	 * 
	 * @param name1
	 * @param mass1
	 * @param name2
	 * @param mass2
	 */
	public void addPair(String name1, int[] mass1, String name2, int[] mass2) {
		if (mass1.length != Utilities.ELEMENTS.length || mass2.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + ((mass1.length != Utilities.ELEMENTS.length)? mass1.length : mass2.length) + ".");
		
		// check for linear dependence
		if (masses.size() > 0) {
			if (!Equivalence.normalize(Equivalence.sumMasses(mass1, mass2)).equals(Equivalence.normalize(masses.get(0))))
				throw new IllegalArgumentException("Mass vector of "+name1+" + "+name2+" is not linearly dependent on "+names.get(0)+".");
		} else if (doubleMasses.size() > 0) {
			if (!Equivalence.normalize(Equivalence.sumMasses(mass1, mass2)).equals(Equivalence.normalize(Equivalence.sumMasses(doubleMasses.get(0), doubleMasses.get(1)))))
				throw new IllegalArgumentException("Mass vector of "+name1+" + "+name2+" is not linearly dependent on "+doubleNames.get(0)+" + "+doubleNames.get(1)+".");
		}
		
		doubleNames.add(name1);
		doubleMasses.add(mass1);
		doubleNames.add(name2);
		doubleMasses.add(mass2);
		
		if (Utilities.DEBUG) {
			ArrayList<Integer> norm0 = Equivalence.normalize(Equivalence.sumMasses(doubleMasses.get(0), doubleMasses.get(1)));
			ArrayList<Integer> norm1 = Equivalence.normalize(Equivalence.sumMasses(mass1, mass2));
			for (int i = 0; i < Utilities.ELEMENTS.length; i++) {
				if (!norm0.get(i).equals(norm1.get(i))) {
					System.out.println("Non-equivalent compounds: " + doubleNames.get(0) + " + " + doubleNames.get(1) + " != " + name1 + " + " + name2);
//					System.out.println("Key 0: " + Equivalence.getKey(norm0) + ", Key 1: " + Equivalence.getKey(norm1));
				}
			}
		}
		
	}
	
	/**
	 * Adds a set of compounds and their corresponding mass vectors
	 * to the mass equivalence class. The sum of the given mass
	 * vectors must be linearly dependent on the existing mass
	 * vectors of the class.
	 * 
	 * @param names
	 * @param masses
	 */
	public void addAllCompounds(ArrayList<String> names, ArrayList<int[]> masses) {
		if (names.size() != masses.size())
			throw new IllegalArgumentException("Number of names (" + names.size() + ") does not match number of masses (" + masses.size() + "). Invalid classes file?");
		
		// need to call the own method to guarantee linear dependence checking 
		for (int i=0; i<names.size(); i++)
			addCompound(names.get(i), masses.get(i));
	}
	
	/**
	 * Adds a set of compound pairs and their corresponding mass vectors
	 * to the mass equivalence class. Consecutive elements are considered
	 * pairs. The sum of the mass vectors of each pair must be linearly
	 * dependent on the existing mass vectors of the class.
	 * @param names
	 * @param masses
	 */
	public void addAllPairs(ArrayList<String> names, ArrayList<int[]> masses) {
		if (names.size() != masses.size())
			throw new IllegalArgumentException("Number of names (" + names.size() + ") does not match number of masses (" + masses.size() + "). Invalid classes file?");
		else if (names.size() % 2 != 0)
			throw new IllegalArgumentException("Odd number of elements: " + names.size() + ". Invalid classes file?");

		// need to call the own method to guarantee linear dependence checking
		for (int i=0; i<names.size()-1; i+=2)
			addPair(names.get(i), masses.get(i), names.get(i+1), masses.get(i+1));
	}
	
	
	/**
	 * Returns the single compound name at the given index.
	 */
	public String getName(int index) {
		return names.get(index);
	}
	
	/**
	 * Returns the single mass vector at the given index.
	 */
	public int[] getMass(int index) {
		return masses.get(index);
	}
	
	/**
	 * Returns the double compound name at the given index.
	 */
	public String getDoubleName(int index) {
		return doubleNames.get(index);	
	}
	
	/**
	 * Returns the double mass vector at the given index.
	 */
	public int[] getDoubleMass(int index) {
		return doubleMasses.get(index);
	}
	
	/**
	 * Returns the set of all single compound names in the class. 
	 */
	public ArrayList<String> getNames() {
		return names;
	}
	
	/**
	 * Returns the set of all single compound masses in the class. 
	 */
	public ArrayList<int[]> getMasses() {
		return masses;
	}

	
	/**
	 * Returns the set of all explicit double compound names in the class. 
	 */
	public ArrayList<String> getExplicitDoubleNames() {
		return doubleNames;
	}
	
	/**
	 * Returns the set of all explicit double compound masses in the class. 
	 */
	public ArrayList<int[]> getExplicitDoubleMasses() {
		return doubleMasses;
	}
	
	/**
	 * Returns a hash map of all pairs of compound names and masses in the class,
	 * including the implicit pairs from single compound equivalences.
	 */
	public HashMap<String[], int[][]> getPairs() {
		HashMap<String[], int[][]> pairs = new HashMap<String[], int[][]>();
		
		// add the explicit pairs
		for (int i=0; i<doubleNames.size()-1; i+=2) {
			if (pairs.put(new String[]{doubleNames.get(i), doubleNames.get(i+1)}, new int[][]{doubleMasses.get(i), doubleMasses.get(i+1)}) != null)
				System.out.println("WARNING: duplicate pair in equivalence class: "+doubleNames.get(i)+" + "+doubleNames.get(i+1));
		}
		// add the implicit pairs of individual compounds
		for (int i=0; i<names.size()-1; i++) {
			for (int j=i+1; j<names.size(); j++) {
				if (pairs.put(new String[]{names.get(i), names.get(j)}, new int[][]{masses.get(i), masses.get(j)}) != null)
					System.out.println("WARNING: duplicate implicit pair in equivalence class: "+names.get(i)+" + "+names.get(j));
			}
		}
		
		return pairs;
	}
	
}
