public class UtilSpectral {
	
	//
	public static double[] nthSmallest(int n, double[][] dist) {
		
		double[] sigmas = new double[dist.length];
		
		for (int row = 0 ; row < dist.length ; row++) {
			int l = dist[row].length;
			for (int k = 0 ; k < n ; k++) {
				swap(dist[row],k,findMin(dist[row],k,l-1));
			}
			sigmas[row] = dist[row][n-1];
		}
		return sigmas;
	}
	
	public static double[] nthSmallest(int n, double[] dist) {
		double[][] pass = new double[1][];
		pass[0] = dist;
		return nthSmallest(n,pass);
	}
	
	private static void swap(double[] arr, int index, int index2) {
		double temp = arr[index];
		arr[index] = arr[index2];
		arr[index2] = temp;
	}
	
	private static int findMin(double[] arr, int start, int end) {
		double min = Double.MAX_VALUE;
		int index = -1;
		
		for (int i = start; i <= end ; i++) {
			if (arr[i] <= min) {
				index = i;
				min = arr[i];
			}
		}
		return index;
	}
	
	public static double[][] computeSpectralProjections(int nth, double thrVal, double[][] refSpikes, double[] refSigmas,
			double[][] allSpikes, double[][] eigenVectors) {
		int outDims = eigenVectors[0].length;
		int nRef = refSpikes.length;
		int nTot = allSpikes.length;
		int PCDims = refSpikes[0].length;
		
		double[] dstMatrix = new double[nRef];
		
		double[][] outputs = new double[nTot][outDims];
		
		for (int i = 0; i < nTot ; i++) {
			for (int j = 0 ; j < nRef ; j++) {
				dstMatrix[j] = 0;
				for (int d = 0 ; d < PCDims ; d++) {
					dstMatrix[j] += (allSpikes[i][d] - refSpikes[j][d]) * (allSpikes[i][d] - refSpikes[j][d]);
				}
			}
			
			double sigma = 1/Math.sqrt(nthSmallest(nth,dstMatrix)[0]);
			
			double expVal;
			for (int j = 0 ; j < nRef ; j++) {
				expVal = 1/(1 + dstMatrix[j] * sigma * refSigmas[j]);
				if (expVal < thrVal)
					expVal = 0;
				//expVal = -dstMatrix[j] * sigma * refSigmas[j];
				for (int out = 0 ; out < outDims ; out++) {
					outputs[i][out] += eigenVectors[j][out] * expVal;
				}
			}
		}
		return outputs;
	}
}