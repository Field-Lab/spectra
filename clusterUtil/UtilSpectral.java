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
	
	static void swap(double[] arr, int index, int index2) {
		double temp = arr[index];
		arr[index] = arr[index2];
		arr[index2] = temp;
	}
	
	static int findMin(double[] arr, int start, int end) {
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
}