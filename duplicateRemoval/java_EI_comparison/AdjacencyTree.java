import java.util.LinkedList;

public class AdjacencyTree {
	AdjacencyTree left;
	AdjacencyTree right;

	int[] sortedSubNodes;

	double thresholdValue;

	// Build Leaf
	public AdjacencyTree(int node) {
		left = null;
		right = null;
		sortedSubNodes = new int[] { node };
		thresholdValue = 0;
	}

	public AdjacencyTree(AdjacencyTree a, AdjacencyTree b, double dist) {
		// a and b MUST be index-exclusive
		// smaller indexes go
		left = a;
		right = b;
		sortedSubNodes = merge(a.sortedSubNodes, b.sortedSubNodes);
		thresholdValue = dist;
	}

	private int[] merge(int[] a, int[] b) {
		int[] out = new int[a.length + b.length];
		int aIn = 0;
		int bIn = 0;
		for (int k = 0; k < out.length; k++) {
			if (aIn < a.length && bIn < b.length) {
				if (a[aIn] < b[bIn])
					out[k] = a[aIn++];
				else
					out[k] = b[bIn++];
			} else if (aIn < a.length) {
				out[k] = a[aIn++];
			} else {
				out[k] = b[bIn++];
			}
		}
		return out;
	}

	public boolean contains(int node) {
		return recDychoSearch(node, 0, sortedSubNodes.length - 1);
	}

	private boolean recDychoSearch(int node, int start, int end) {
		if (start == end)
			return (sortedSubNodes[start] == node);

		int m = (start + end) / 2;
		return (recDychoSearch(node, start, m) || recDychoSearch(node, m + 1,
				end));
	}

	public String toString() {
		String s = new String();
		for (int n = 0; n < sortedSubNodes.length; n++) {
			s = s + " " + sortedSubNodes[n];
		}
		return s;
	}

	public int depth() {
		if (left == null)
			if (right == null)
				return 1;
			else
				return 1 + right.depth();
		else if (right == null)
			return 1 + left.depth();
		else
			return 1 + Math.max(left.depth(), right.depth());
	}

	public int innerSize() {
		if (left == null)
			if (right == null)
				return 0;
			else
				return 1 + right.innerSize();
		else if (right == null)
			return 1 + left.innerSize();
		else
			return 1 + left.innerSize() + right.innerSize();
	}
	
	public double[] thresholdList() {
		LinkedList<Double> thrList = new LinkedList<Double>();
		thresholdListRec(thrList);
		int n = thrList.size();
		double[] thrArr = new double[n];
		for (int k = 0 ; k < n ; k++) {
			thrArr[k] = thrList.poll();
		}
		return thrArr;
	}
	
	private void thresholdListRec(LinkedList<Double> baseList) {
		if (left == null)
			if (right == null) // Leaf case
				return;
			else { // Single son
				baseList.addLast(thresholdValue);
				right.thresholdListRec(baseList);
			}
		else if (right == null) { // Single son
			baseList.addLast(thresholdValue);
			left.thresholdListRec(baseList);
		} else { // Double son
			baseList.addLast(thresholdValue);
			left.thresholdListRec(baseList);
			right.thresholdListRec(baseList);
		}
	}
	
	public int[][] partitioning(double threshold) {
		LinkedList<int[]> connComp = new LinkedList<int[]>();
		partitioningRec(threshold, connComp);
		
		int[][] connCompArr = new int[connComp.size()][];
		for ( int n = 0 ; n < connCompArr.length ; n++ ) {
			connCompArr[n] = connComp.poll().clone(); // Clone not necessary as done by matlab. Better clone for passing to further java.
		}
		return connCompArr;
		
	}
	
	private void partitioningRec(double threshold, LinkedList<int[]> connComp) {
		// Leaf case
		if (left == null && right == null) {
			connComp.addLast(sortedSubNodes);
			return;
		}
		// Subthreshold case
		if (threshold > thresholdValue) {
			connComp.addLast(sortedSubNodes);
			return;
		}
		// Above Threshold case, recursive descent.
		// Should be 2 sons exactly
		if (left != null)
			left.partitioningRec(threshold, connComp);
		if (right != null)
			right.partitioningRec(threshold, connComp);
	}
	
	public double[][] enumRectangles() {
		LinkedList<Rectangle> rectList = new LinkedList<Rectangle>();
		
		rectList.add(new Rectangle(thresholdValue, thresholdValue + 1, 0, sortedSubNodes.length));
		
		enumRectanglesRec(rectList,0);
		
		int n = rectList.size();
		
		double[][] rectArr = new double[n][4];
		for (int k = 0 ; k < n ; k++) {
			rectArr[k] = rectList.poll().toArr();
		}
		return rectArr;
	}
	
	private void enumRectanglesRec(LinkedList<Rectangle> baseList,int offset) {
		if (left == null)
			if (right == null) // Leaf case
				return;
			else { // Single son
				baseList.addLast(new Rectangle(right.thresholdValue, thresholdValue, offset, offset + right.sortedSubNodes.length));
				right.enumRectanglesRec(baseList, offset + right.sortedSubNodes.length);
				System.out.println("Err: single son should not happen.");
			}
		else if (right == null) { // Single son
			baseList.addLast(new Rectangle(left.thresholdValue, thresholdValue, offset, offset + left.sortedSubNodes.length));
			left.enumRectanglesRec(baseList, offset + left.sortedSubNodes.length);
			System.out.println("Err: single son should not happen.");
		} else { // Double son
			baseList.addLast(new Rectangle(left.thresholdValue, thresholdValue, offset, offset + left.sortedSubNodes.length));
			baseList.addLast(new Rectangle(right.thresholdValue, thresholdValue, offset + left.sortedSubNodes.length, offset + sortedSubNodes.length));
			left.enumRectanglesRec(baseList, offset);
			right.enumRectanglesRec(baseList, offset + left.sortedSubNodes.length);
		}
	}
	
	class Rectangle {
		double lowerThresh;
		double upperThresh;
		int startNum;
		int endNum;
		
		Rectangle(double low, double high, int s, int e) {
			lowerThresh = low;
			upperThresh = high;
			startNum = s;
			endNum = e;
		}
		
		public double[] toArr() {
			return new double[] {lowerThresh, upperThresh, (double) startNum, (double) endNum};
		}
	}
}
