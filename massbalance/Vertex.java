/**
 * 
 */
package massbalance;

import java.io.Serializable;

/**
 * Basic vertex object representing compound and reaction vertices.
 * 
 * @author Georg Basler
 *
 */
public class Vertex implements Serializable {

	private static final long serialVersionUID = -2969225075285620605L;
	static final boolean COMPOUND = true;
	static final boolean REACTION = false;
	
	private String name;
	private String compartment;
	private Double deltaG;
	private boolean type;
	private int[] mass;
	private boolean isReversed;
	private Vertex reversedReaction;
	
	/**
	 * Creates a new compound vertex with the given name and given compartment.
	 */
	public Vertex(String name, int[] mass, String compartment) {
		if (mass != null && mass.length != Utilities.ELEMENTS.length)
			throw new IllegalArgumentException("Mass vector has invalid size: " + mass.length + ".");
		
		this.name = name;
		this.mass = mass;
		this.compartment = compartment;
		this.type = COMPOUND;
	}
	
	/**
	 * Creates a new reaction with the given name. If reversible==true, a corresponding reversed reaction vertex is created.
	 * 
	 * @param name
	 * @param reversible
	 */
	public Vertex(String name, boolean reversible) {
		this.name = name;
		this.type = REACTION;
		this.isReversed = false;
		if (reversible)
			// create the corresponding reversed reaction
			this.reversedReaction = new Vertex(name, this);
	}
	
	/**
	 * Private constructor to create and link a reversed reaction.
	 * 
	 * @param name
	 * @param reversedReaction
	 */
	private Vertex(String name, Vertex reversedReaction) {
		this.name = name;
		this.type = REACTION;
		this.isReversed = true;
		// link to forward reaction
		this.reversedReaction = reversedReaction;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getName() {
		return this.name;
	}
	
	/**
	 * Returns the deltaGf of a compound.
	 * 
	 * @return deltaGf.
	 */
	public Double getDeltaGf() {
		if (this.type == REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.REACTION.");
		return this.deltaG;
	}
	
	/**
	 * Sets the deltaGf of a compound.
	 * 
	 * @param deltaGf
	 */
	public void setDeltaGf(Double deltaGf) {
		if (this.type == REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.REACTION.");
		this.deltaG = deltaGf;
	}
	
	public boolean getType() {
		return this.type;
	}
	
	public int[] getMass() {
		if (this.type == REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.REACTION.");
		return this.mass;
	}
	
	public String getCompartment() {
		if (this.type == REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.REACTION.");
		return this.compartment;
	}
	
	public boolean isReversed() {
		if (this.type != REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.COMPOUND.");
		
		return this.isReversed;
	}
	
	public Vertex reversedReaction() {
		if (this.type != REACTION)
			throw new IllegalArgumentException("Cannot invoke method on a Vertex of type Vertex.COMPOUND.");
		
		return this.reversedReaction;
	}
	
	/**
	 * Compares the member variables of this vertex to the member variables of the given vertex. In particular,
	 * compares the name, type, the mass for compounds, isReversed and the reversed reaction for reactions.
	 * 
	 * Does not compare the reversed reaction.
	 * 
	 * @param v
	 * @return true, if the vertices are equal, false otherwise.
	 */
	public boolean compare(Vertex v) {
		if (getType() == REACTION && reversedReaction() != null || v.getType() == REACTION && v.reversedReaction() != null) {
			// compare the reversed reaction
			if (!compare_impl(reversedReaction(), v.reversedReaction()))
				return false;
		}
		// compare the vertices
		return compare_impl(this, v);
	}
	
	/**
	 * Implementation used for comparison of vertices and reversed reactions.
	 * 
	 * @param v1
	 * @param v2
	 * @return
	 */
	private boolean compare_impl(Vertex v1, Vertex v2) {
		if (v1.getName() != v2.getName())
			return false;
		if (v1.getType() != v2.getType())
			return false;
		if (v1.getType() == REACTION)
			if (v1.isReversed() != v2.isReversed())
				return false;
		if (v1.getType() == COMPOUND) {
			if (v1.getCompartment() == null ^ v2.getCompartment() == null)
				return false;
			if (v1.getCompartment() != null && !v1.getCompartment().equals(v2.getCompartment()))
				return false;
			if (v1.getMass() != v2.getMass())
				return false;
			if (v1.getMass() != null) {
				for (int i=0; i<v1.getMass().length; i++)
					if (v1.getMass()[i] != v2.getMass()[i])
						return false;
			}
		}
		return true;
	}
	
	@Override
	public Vertex clone() {
		
		if (this.getType() == Vertex.COMPOUND) {
			Vertex compound = new Vertex(this.getName(), this.getMass(), this.getCompartment());
			compound.setDeltaGf(this.getDeltaGf());
			return compound;
		} else {
			if (this.isReversed())
				throw new RuntimeException("Cannot invoke clone on a reversed reaction vertex.");
			return new Vertex(this.getName(), this.reversedReaction() != null);
		}
	}
}
