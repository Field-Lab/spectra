import java.util.LinkedList;

/**
 * Graph object defines all connected components of an undirected
 * real-valued complete graph, in which all edges above a certain threshold are removed
 * The data structure describes all the currently known connected components in a list
 * 
 * In normal use, all singletons are generated in the list, then edges are processed in
 * increasing order, possibly merging certain connected components, until only one tree
 * remains in the list.
 * 
 * @see AdjacencyTree
 * @author Vincent Deo, Stanford University
 */
public class Graph {
	private LinkedList<AdjacencyTree> connComp;
	public boolean debug = false;
	
	/**
	 * Generate new graph with n singleton connected components
	 */
	public Graph(int n) {
		connComp = new LinkedList<AdjacencyTree>();
		for (int i = 0 ; i < n ; i++) {
			connComp.addLast(new AdjacencyTree(i));
		}
	}
	
	/**
	 * Process all real-valued graph edges
	 * For the structure to be correct, edges MUST be in increasing value order
	 * 
	 * @param l left node of edge
	 * @param r right node of edge
	 * @param w value of edge
	 */
	public void addAllEdges(int[] l, int[] r, double[] w) {
		if (debug) {
			System.out.println("Processing " + w.length + " edges.");
			System.out.println("-------------------- Progress ---------------------");
		}
		for (int i = 0 ; i < w.length ; i++) {
			addEdge(l[i],r[i],w[i]);
			if (debug && w.length >= 50 && i%(w.length/50) == 0)
				System.out.print("*");
		}
		if (debug)
			System.out.println();
	}

	/**
	 * Processes one real-valued edge of the graph
	 */
	private void addEdge(int l, int r, double w) {
		AdjacencyTree lTree = findTree(l); // Find connected component of left node
		if (lTree.contains(r)) {// Both nodes already connected
			connComp.addFirst(lTree); // Push back connected component
			return;
		}
		
		AdjacencyTree rTree = findTree(r); // Find connected component of right node
		connComp.addFirst(new AdjacencyTree(lTree, rTree, w)); // Merge under w value and push in graph list
	}
	
	/**
	 * Retrieves and removes the connected component containing a given node
	 */
	private AdjacencyTree findTree(int node) {
		for (int i = 0; i < connComp.size(); i++) {
			if (connComp.get(i).contains(node)) {
				AdjacencyTree t = connComp.remove(i);
				return t;
			}
		}
		return null;
	}
	
	/**
	 * Generates sequence of connected components of the graph 
	 */
	public String toString() {
		String s = new String();
		for (AdjacencyTree t: connComp) {
			s = s + t.toString() + "\n";
		}
		return s;
	}
	
	/**
	 * Return the final merge herarchy once all edges are processed and only a single connected component remains
	 */
	public AdjacencyTree getFinalTree() {
		if (connComp.size() > 1)
			return null;
		else
			return connComp.getFirst();
	}
	
	/**
	 * Removes all singleton connected components from the Graph
	 * This is useful when singletons remain after processing all non-zero, non-nan edges,
	 * which are singular ellipses to discard.
	 */
	public void removeSingletons() {
		if (connComp.isEmpty())
			return;
		AdjacencyTree t = connComp.pollFirst();
		removeSingletons();
		if (t.innerSize() > 0)
			connComp.addFirst(t);
	}
}
