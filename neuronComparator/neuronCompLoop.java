import java.util.Arrays;

/**
 * Algo utility class for comparing two different clusterings,
 * ie two different sets of neurons for the same data
 * 
 * @author Vincent Deo, Stanford University
 */
public class neuronCompLoop {
    
    /**
     * Optimizes the permutation and subsets yielding the maximum trace of the
     * intersection score matrix.
     */
    public static double maxMetric(double[][] intersect) {
        double maxTemp = 0;
        boolean rows = intersect.length <= intersect[0].length;
        int lowLength = Math.min(intersect.length, intersect[0].length);
        int highLength = Math.max(intersect.length, intersect[0].length);
        
        int[] subset = new int[lowLength];
        int[] perm = new int[lowLength];
        for(int i = 0 ; i < lowLength ; i++)
            subset[i] = i;
        do {
            for(int i = 0 ; i < lowLength ; i++)
                perm[i] = i;
            do {
                double sum = 0;
                for (int n = 0 ; n < lowLength ; n++) {
                    if (rows)
                        sum += intersect[n][subset[perm[n]]];
                    else
                        sum += intersect[subset[perm[n]]][n];
                    if (sum > maxTemp)
                    maxTemp = sum;
                }
            } while (nextPerm(perm));
        } while (nextSubset(subset,highLength));
        
        return maxTemp;
    }
    
    /**
     * Compute the cardinality of all intersections between neurons
     * 
     * @param neurArray1 Concatenated spike times of all neurons in clustering 1
     * @param neurArray2 Concatenated spike times of all neurons in clustering 2
     * @param length1 array of length nNeurons1, length1[i] = nSpikes(neuron_i) in clustering 1
     * @param length2 array of length nNeurons2, length2[i] = nSpikes(neuron_i) in clustering 2
     */
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
    
    /**
     *  Next permutation function
     * Operates on starting array 1 ... n
     * And iterates permutations by side effects returning true
     * until n ... 1 is reached, then returns false
     */
    public static boolean nextPerm(int[] perm) {
        int i = 0;
        int j = 0;
        //From right to left, find the first one that is not in ascending order.
        for (i = perm.length - 2; i >= 0; i--) {
            if (perm[i] >= perm[i + 1])
                continue;
            //From right to left, find the first one that is larger than num[i]
            for (j = perm.length - 1; j > i; j--) {
                if (perm[j] > perm[i])
                    break;
            }
            break;
        }
        //If we find i, swap the number on position i and j
        if (i >= 0) {
            int tmp = perm[i];
            perm[i] = perm[j];
            perm[j] = tmp;
        } else
            return false;
        //Reverse the numbers which are on the right of i
        int start = i + 1;
        int end = perm.length - 1;
        while (start < end) {
            int tmp = perm[start];
            perm[start] = perm[end];
            perm[end] = tmp;
            start++;
            end--;
        }
        return true;
    }
    
    /**
     *  Next subset function
     * Operates on starting array 1 ... k
     * And iterates subsets of size k by side effects returning true
     * until n-k+1 ... n is reached, then returns false
     */
    public static boolean nextSubset(int[] subset, int totSize) {
        // From right: find first gap.
        int val = totSize-1;
        int index = subset.length-1;
        while (index >= 0 && subset[index] == val) {
            index--;
            val--;
        }
        if (index == -1)
            return false;
        subset[index]++;
        for(int i = index + 1 ; i < subset.length ; i++) {
            subset[i] = subset[index] + i - index;
        }
        return true;
    }
}