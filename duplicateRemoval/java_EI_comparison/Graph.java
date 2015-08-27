import java.util.LinkedList;

public class Graph {
	private LinkedList<AdjacencyTree> connComp;
	public boolean debug = false;
	
	public Graph(int n) {
		connComp = new LinkedList<AdjacencyTree>();
		for (int i = 0 ; i < n ; i++) {
			connComp.addLast(new AdjacencyTree(i));
		}
	}
	
	// Pass in list of edges, WITH W increasingly sorted already (ask matlab...)
	public void addAllEdges(int[] l, int[] r, double[] w) {
		if (debug) {
			System.out.println("Processing " + w.length + " edges.");
			System.out.println("-------------------- Progress ---------------------");
		}
		for (int i = 0 ; i < w.length ; i++) {
			addEdge(l[i],r[i],w[i]);
			if (debug && i%(w.length/50) == 0)
				System.out.print("*");
		}
		if (debug)
			System.out.println();
	}

	private void addEdge(int l, int r, double w) {
		AdjacencyTree lTree = findTree(l);
		if (lTree.contains(r)) {// Both nodes already connected 
			connComp.addFirst(lTree);
			return;
		}
		
		AdjacencyTree rTree = findTree(r);
		connComp.addFirst(new AdjacencyTree(lTree, rTree, w));
	}
	
	private AdjacencyTree findTree(int node) {
		for (int i = 0; i < connComp.size(); i++) {
			if (connComp.get(i).contains(node)) {
				AdjacencyTree t = connComp.remove(i);
				return t;
			}
		}
		return null;
	}
	
	public String toString() {
		String s = new String();
		for (AdjacencyTree t: connComp) {
			s = s + t.toString() + "\n";
		}
		return s;
	}
	
	public AdjacencyTree getFinalTree() {
		if (connComp.size() > 1)
			return null;
		else
			return connComp.getFirst();
	}
}
