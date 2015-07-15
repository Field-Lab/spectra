public class neuronCompLoop {

	public static int[][] computeIntersects (int[] neurArray1, int[] neurArray2, int[] length1, int[] length2) {
		int m = length1.length;
		int n = length2.length;
		
		int[][] result = new int[m][n];
		
		int offset1 = 0;
		int offset2 = 0;
		
		top:
		for (int i = 0 ; i < m ; i++) {
			for (int j = 0 ; j < n ; j++) {
				int k = 0;
				int l = 0;
				while (k < length1[i] && l < length2[j]) {
					try {
					if (neurArray1[offset1 + k] == neurArray2[offset2 + l]) {
						k++;
						l++;
						result[i][j]++;
						continue;
					}
					if (neurArray1[offset1 + k] < neurArray2[offset2 + l]) {
						k++;
					} else {
						l++;
					}
				} catch (Exception e) {
					System.out.println(i + " " + j + " " + k + " " + l + " " + offset1 + " " + offset2);
					break top;
				}
				}
				offset2 += length2[j];
			}
			offset1 += length1[i];
			offset2 = 0;
		}
		return result;
	}
}