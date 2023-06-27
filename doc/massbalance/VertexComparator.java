/**
 * 
 */
package massbalance;

import java.util.Comparator;

/**
 * Class for comparing and ordering vertices. Compounds are compared by their name and compartment,
 * reaction vertices by name and reversibility.
 * 
 * @author Georg Basler
 *
 */
public class VertexComparator implements Comparator<Vertex> {

	/**
	 * Compound vertices are compared their name and compartment,
	 * reaction vertices by their name and reversibility.
	 * 
	 * @return A negative integer, zero, or a positive integer as the
     * first vertex is considered less than, equal to, or greater
     * than the second.
	 */
	public int compare(Vertex v1, Vertex v2) {

		if (v1.getType() != v2.getType()) {
			// compare unequal vertex types by their type
			return (v1.getType() ? -1 : 1);
			
		} else if (v1.getType() == Vertex.COMPOUND) {
			// compare compound vertices by their name and compartment
   			return ((v1.getName()+"$"+v1.getCompartment())).compareTo(
   					(v2.getName()+"$"+v2.getCompartment()));
   			
		} else if (v1.getType() == Vertex.REACTION) {
			// compare reaction vertices by their name and reversibility
			return ((v1.getName()+(v1.isReversed()?"$rev":"")).compareTo(
					(v2.getName()+(v2.isReversed()?"$rev":""))));
		} else {
			throw new RuntimeException("VertexComparator not implemented for vertices of type "+v1.getType()+".");
		}
	}
}
