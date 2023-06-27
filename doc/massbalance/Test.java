
package massbalance;

import java.util.ArrayList;
import java.util.HashMap;


public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String a = "aaaslkdj$$ sdl $$$ alskj4$askl";
		a = a.replaceAll("(\\$)+", "\\$");
		System.out.println(a);
		
		final String str = "APPLEE";
		String replaced = str.replaceAll("(.)\\1", "$1");
		System.out.println(replaced);

		
		// regression tests for Randomize
		
//		System.out.println("*** TESTING bsubtilis NO ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\bsubtilis2\\bsubtilis-compartments.network", "C:\\Temp\\bsubtilis2"});
//		System.out.println("*** TESTING bsubtilis MASSBALANCE ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\bsubtilis2\\bsubtilis-compartments.network", "C:\\Temp\\bsubtilis2", "massbalance", "0-1"});
//		System.out.println("*** TESTING bsubtilis COMPARTMENTS ***");
//		try {
//			Randomize.main(new String[]{"H:\\MassBalance\\data\\bsubtilis2\\bsubtilis-compartments.network", "C:\\Temp\\bsubtilis2", "massbalance", "1-2", "compartments"});
//		} catch (IllegalArgumentException e) {
//			if (e.getMessage().startsWith("Missing compartment"))
//				System.out.println("Passed with "+e.getMessage());
//			else
//				throw e;
//		}
//		System.out.println("*** TESTING bsubtilis ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\bsubtilis2\\bsubtilis-compartments.network", "C:\\Temp\\bsubtilis2", "massbalance", "2-3", "reversible", "nofix", "noniterative", "depth=2", "p=0.3", "substitutability", "strict"});
//		System.out.println("*** TESTING bsubtilis SWITCH ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\bsubtilis2\\bsubtilis-compartments.network", "C:\\Temp\\bsubtilis2", "switch", "0-0", "reversible", "nofix", "noniterative", "depth=3", "p=0.1", "substitutability", "strict"});
//		System.out.println();
		
//		System.out.println("*** TESTING arabidopsis NO ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\arabidopsis\\1752-0509-4-114-s6.xml", "C:\\Temp\\arabidopsis"});
//		System.out.println("*** TESTING arabidopsis MASSBALANCE ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\arabidopsis\\1752-0509-4-114-s6.xml", "C:\\Temp\\arabidopsis", "massbalance", "0-1"});
//		System.out.println("*** TESTING arabidopsis ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\arabidopsis\\1752-0509-4-114-s6.xml", "C:\\Temp\\arabidopsis",  "massbalance", "1-2", "compartments", "reversible", "nofix", "noniterative", "depth=2", "p=0.3", "substitutability", "strict"});
//		System.out.println("*** TESTING arabidopsis SWITCH ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\arabidopsis\\1752-0509-4-114-s6.xml", "C:\\Temp\\arabidopsis", "switch", "0-0", "reversible", "nofix", "noniterative", "depth=3", "p=0.1", "substitutability", "strict"});
//		System.out.println();
//		
//		System.out.println("*** TESTING ecocyc14.6 NO ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\ecocyc14.6\\reactions.dat", "C:\\Temp\\ecocyc14.6"});
//		System.out.println("*** TESTING ecocyc14.6 MASSBALANCE ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\ecocyc14.6\\reactions.dat", "C:\\Temp\\ecocyc14.6", "massbalance", "0-1"});
//		System.out.println("*** TESTING ecocyc14.6 ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\ecocyc14.6\\reactions.dat", "C:\\Temp\\ecocyc14.6",  "massbalance", "1-2", "compartments", "reversible", "nofix", "noniterative", "depth=2", "p=0.3", "substitutability", "strict"});
//		System.out.println("*** TESTING ecocyc14.6 SWITCH ALL ARGS ***");
//		Randomize.main(new String[]{"H:\\MassBalance\\data\\ecocyc14.6\\reactions.dat", "C:\\Temp\\ecocyc14.6", "switch", "0-0", "reversible", "nofix", "noniterative", "depth=3", "p=0.1", "substitutability", "strict"});		
//		System.out.println();
		
		// regression tests for MetabolicGraphs
		
//		System.out.println("*** TESTING arabidopsis properties MASSBALANCE ***");
//		try {
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0"});
//			System.out.println("*** TESTING arabidopsis properties MATRIX ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "matrix-labelled-sorted"});
//			System.out.println("*** TESTING arabidopsis properties degrees ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "degrees"});		
//			System.out.println("*** TESTING arabidopsis properties weights ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "weights"});
//			System.out.println("*** TESTING arabidopsis properties scopes ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "scopes=5"});
//			System.out.println("*** TESTING arabidopsis properties pathLength ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "pathLength"});
//			System.out.println("*** TESTING arabidopsis properties clustering ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "clustering"});
//			System.out.println("*** TESTING arabidopsis properties path ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "path-directed=Ath_C0644,Ath_C0576,Ath_C0004"});
//			System.out.println("*** TESTING arabidopsis properties connected ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "connected=Ath_C0644>Ath_C0004>Ath_C0074"});
//			System.out.println("*** TESTING arabidopsis properties deltaGr ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "deltaGr"});
//			System.out.println("*** TESTING arabidopsis properties deltaGn ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "deltaGn"});
//			System.out.println("*** TESTING arabidopsis properties assortativity ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "assortativity"});
//			System.out.println("*** TESTING arabidopsis properties transitionDegree ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "transitionDegree"});
//			System.out.println("*** TESTING arabidopsis properties matrixSubstitutions ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "matrixSubstitutions"});
//			System.out.println("*** TESTING arabidopsis properties centrality ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "centrality-reversible-double,d=0.5"});
//			System.out.println("*** TESTING arabidopsis properties localCentrality ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "localCentrality"});
//			System.out.println("*** TESTING arabidopsis properties knockoutPath ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "knockoutPath,e=0.1,Ath_R0025_1"});
//			System.out.println("*** TESTING arabidopsis properties write ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "write"});
//			System.out.println("*** TESTING arabidopsis properties write2 ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "write2"});
//			System.out.println("*** TESTING arabidopsis properties FROM-TO ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "4-7"});
//			System.out.println("*** TESTING arabidopsis properties ALL ARGS ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "switch", "2-2", "matrix", "degrees", "weights", "scopes=2", "pathLength", "clustering", "path=Ath_C0644,Ath_C0576,Ath_C0004", "connected=Ath_C0644>Ath_C0004>Ath_C0074", "deltaGr", "deltaGn", "assortativity", "transitionDegree", "matrixSubstitutions", "centrality,d=0.98", "localCentrality", "knockoutPath,e=0.1,Ath_C0004", "write", "write2"});		
//			System.out.println("*** TESTING arabidopsis properties SWITCH ALL ARGS ***");
//			Properties.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","switch", "2-2", "matrix", "scopes=2", "pathLength", "clustering", "path=Ath_C0644,Ath_C0576,Ath_C0004", "connected=Ath_C0644>Ath_C0004>Ath_C0074", "deltaGr", "deltaGn", "assortativity", "transitionDegree", "matrixSubstitutions", "centrality,d=0.98", "localCentrality", "knockoutPath,e=0.1,Ath_C0004", "write", "write2"});
//		} catch (InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		
//		MetabolicGraphs.main(new String[]{"C:\\Temp\\arabidopsis", "arabidopsis","massbalance", "0-0", "connected=Ath_C0644>Ath_C0004>Ath_C0074"});
		
		
//		Integer i = 3;
//		TestIntegerObject(i);
//		System.out.println("TestIntegerObject: "+i);
//		
//		System.out.println("TestIntegers: ");
//		int a = 0, b = 1, c = 2;
//		TestIntegers(a, b, c);
//		System.out.println(a);
//		System.out.println(b);
//		System.out.println(c);
//		
//		System.out.println("TestIntegerArray: ");
//		int[] d = new int[]{0, 1, 2};
//		TestIntegerArray(d);
//		for (int j=0; j<d.length; j++)
//			System.out.println(d[j]);
//		
//		System.out.println("TestArrayList:");
//		ArrayList<Integer> e = new ArrayList<Integer>();
//		e.add(0);
//		e.add(1);
//		e.add(2);
//		TestArrayList(e);
//		for (Integer f : e)
//			System.out.println(f);
		
		// test object parameter by reference
//		HashMap<Vertex, Integer> map = new HashMap<Vertex, Integer>();
//		Vertex v1 = new Vertex("testReaction1", false);
//		map.put(v1, 1);
//		
//		System.out.println("return value: "+TestHashMap(map));
//		
//		Vertex v3 = new Vertex("testReaction3", false);
//		map.put(v3, 3);
//		
//		for (Vertex v : map.keySet()) {
//			System.out.println(v.getName()+", rev: "+(v.reversedReaction()==null)+", "+map.get(v));
//		}
		
//		new Test();
//		
//		Random random = new Random();
//		int day = random.nextInt(365);
//		int year = random.nextInt(47);
//		
//		System.out.println("day = "+(int)(day+1));
//		System.out.println("year = "+(int)(1945)+year);
	}
	
	public Test() {
		
//		Substitution s = invoke();
//		ArrayList<Vertex> substitutes = s.substitutes;
//		ArrayList<int[]> stoichiometry = s.stoichiometry;
//		ArrayList<Integer> edgeIndices = s.edgeIndices;
//		
//		for (Vertex v : substitutes)
//			System.out.print(v.getName()+" [reversed="+(v.reversedReaction()!=null)+"]; ");
//		System.out.println();
//
//		for (int[] ints : stoichiometry) {
//			for (int i=0; i<ints.length; i++)
//				System.out.print(ints[i]+"; ");
//			System.out.print(" - ");
//		}
//		System.out.println();
//		
//		for (Integer integer : edgeIndices)
//			System.out.print(integer.intValue());
//		System.out.println();
//		
//		System.out.println("Finished.");
	}
	
	public static int TestHashMap(HashMap<Vertex, Integer> map) {
		Vertex v = new Vertex("testReaction2", false);
		map.put(v, 2);
		return 0;
	}
	
	public static int TestIntegerObject(Integer i) {
		i += 1;
		return 0;
	}
	
	public static int TestIntegers(int a, int b, int c) {
		a +=1;
		b++;
		c++;
		return 0;
	}
	
	public static int TestIntegerArray(int[] a) {
		a[0] += 1;
		a[1] += 1;
		a[2] += 1;
		return 0;
	}
	
	public static int TestArrayList(ArrayList<Integer> a) {
		a.set(0, a.get(0)+1);
		a.remove(1);
		return 0;
	}

	public Substitution invoke() {
		
		Substitution s = new Substitution();
		
		s.substitutes.clear();
		s.stoichiometry.clear();
		s.edgeIndices.clear();
		
		Vertex test1 = new Vertex("test1", false);
		Vertex test2 = new Vertex("test2", true);
		s.substitutes.add(test1);
		s.substitutes.add(test2);
		
		s.stoichiometry.add(new int[]{1,2,3,4,5});
		int[] test3 = new int[]{6,7,8,9};
		s.stoichiometry.add(test3);
		
		s.edgeIndices.add(new Integer(3));
		Integer test4 = new Integer(4);
		s.edgeIndices.add(test4);
		
		return s;
	}
	
	private class Substitution {
		ArrayList<Vertex> substitutes = new ArrayList<Vertex>(3);
		ArrayList<int[]> stoichiometry = new ArrayList<int[]>(3);
		ArrayList<Integer> edgeIndices = new ArrayList<Integer>(3);
		
	}
}

