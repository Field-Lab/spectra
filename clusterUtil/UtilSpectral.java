/**
 * Class of math/algo utility method used at some point in spectral clustering
 * or clustering post-analysis
 * 
 * @author Vincent Deo, Stanford University
 */
public class UtilSpectral {

	/**
	 * Determines the index of the n-th smallest value in each row of a 2D
	 * matrix. Could be done in matlab, but blazingly faster in Java except for
	 * n = 1 Used when determining the spectral clustering local sigma scale,
	 * given by the distance from a point to its n-th closest neighbor
	 * 
	 * @param n
	 *            number of smallest value to find in each row
	 * @param dist
	 *            matrix of values
	 * @return array of n-th smallest value in each row of dist
	 */
	public static double[] nthSmallest(int n, double[][] dist) {
		// output alloc
		double[] sigmas = new double[dist.length];

		/*
		 * Done by starting an insertion sort at the start of the array, then
		 * stop This is enough for n small. But will be badly quadratic for n
		 * larger, with n of the order of row length OPTIMIZATION n small: store
		 * sorted buffer of n lowest values, go through and update OPTIMIZATION
		 * n large: quicksort everything first...
		 */
		for (int row = 0; row < dist.length; row++) {
			int l = dist[row].length;
			for (int k = 0; k < n; k++) {
				swap(dist[row], k, findMin(dist[row], k, l - 1));
			}
			sigmas[row] = dist[row][n - 1];
		}
		return sigmas;
	}

	/**
	 * Determines the index of the n-th smallest value in an array. Done by
	 * casting into single row 2d matrix and calling surcharged function aboved
	 * 
	 * @param n
	 *            number of smallest value to find in each row
	 * @param dist
	 *            array of values
	 * @return size 1 array containing n-th smallest value of dist
	 */
	public static double[] nthSmallest(int n, double[] dist) {
		double[][] pass = new double[1][];
		pass[0] = dist;
		return nthSmallest(n, pass);
	}

	/*
	 * Array index swapping function
	 */
	private static void swap(double[] arr, int index, int index2) {
		double temp = arr[index];
		arr[index] = arr[index2];
		arr[index2] = temp;
	}

	/*
	 * Find index of minimum of subarray
	 */
	private static int findMin(double[] arr, int start, int end) {
		double min = Double.MAX_VALUE;
		int index = -1;

		for (int i = start; i <= end; i++) {
			if (arr[i] <= min) {
				index = i;
				min = arr[i];
			}
		}
		return index;
	}

	/*
	 * Shortcut utility function to compute affinity/laplacian matrices here
	 * rather than in Matlab RAM management is more convenient by reference
	 * Timing is roughly equivalent with gaussian affinity Timing is much faster
	 * with power decrease affinities - clustering quality untested
	 * 
	 * @param nth n-th neighbor ot consider for local scaling
	 * 
	 * @param thrVal threshold value at which to cut affinity metric to 0
	 * 
	 * @param refSpikes PCs of reference spikes, used for computing the affinity
	 * Laplacian
	 * 
	 * @param refSigmas computed local scales for the reference spikes
	 * 
	 * @param allSpikes PCs of a bigger subset of spikes on electrode, possibly
	 * all
	 * 
	 * @param eigenVectors precomputed lowest eigenvalue eigenvectors of
	 * affinity laplacian
	 */
	public static double[][] computeSpectralProjections(int nth, double thrVal, double[][] refSpikes,
			double[] refSigmas, double[][] allSpikes, double[][] eigenVectors) {
		int outDims = eigenVectors[0].length; // Laplacian eigenspace dimension
		int nRef = refSpikes.length; // number of reference spikes
		int nTot = allSpikes.length; // number of spikes total
		int PCDims = refSpikes[0].length; // dimensions in PC space

		// nRef x 1 affinity computation buffer
		double[] dstMatrix = new double[nRef]; 

		// nTot x outDims output
		double[][] outputs = new double[nTot][outDims];


		for (int i = 0; i < nTot; i++) { // Loop on all spikes
			for (int j = 0; j < nRef; j++) { // Loop on ref spikes
				dstMatrix[j] = 0;
				// Fill dstMatrix[j] with quadratic euclidian distance between ref spike j and (all) spike i
				for (int d = 0; d < PCDims; d++) {
					dstMatrix[j] += (allSpikes[i][d] - refSpikes[j][d]) * (allSpikes[i][d] - refSpikes[j][d]);
				} // d
			} // j

			// Inverse sqrt of local sigma for spike i
			double sigma = 1 / Math.sqrt(nthSmallest(nth, dstMatrix)[0]);

			double expVal;
			
			// Loop on ref spikes
			for (int j = 0; j < nRef; j++) {
				// Compute normalized rows and columns of extended all-to-ref laplacian
				expVal = 1 / (1 + dstMatrix[j] * sigma * refSigmas[j]); // This is quad law, not exp law
				
				if (expVal < thrVal) // threshold
					expVal = 0;
				// expVal = -dstMatrix[j] * sigma * refSigmas[j];
				
				// Compute matrix multiplication with eigenVectors step by step
				for (int out = 0; out < outDims; out++) {
					outputs[i][out] += eigenVectors[j][out] * expVal;
				}
			}
		} // i
		
		return outputs;
	}
}