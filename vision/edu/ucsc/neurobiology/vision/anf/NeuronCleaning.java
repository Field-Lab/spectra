package edu.ucsc.neurobiology.vision.anf;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.analysis.AutocorrelationCalculator;
import edu.ucsc.neurobiology.vision.analysis.CrossCorrelationCalculator;
import edu.ucsc.neurobiology.vision.calculations.AbstractCalculation;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.VisionHeader;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionParams;

/**
 * Calculation to clean the raw neurons found in a data set. First, all neurons
 * which do not spike much are removed. Second, all contaminated neurons are
 * removed. Last, duplicate neurons are remoced.
 * 
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, UCSC
 */
public class NeuronCleaning extends AbstractCalculation {

	private String neuronFileName;
	private NeuronFile neuronFile;
	private int minSpikesCount;
	private double maxContamination;
	private int coincidenceTime; // in samples. All correlated spikes within
									// this range are considered contamination.
									// A range of phase shifts are tested, and
									// the worst is kept.
	private double maxCorrelation;
	static final int deltaT = 100; // samples for cross correlation plot. Must
									// be larger than coincidenceTime

	public void startCalculation() throws Exception {
		Vision app = Vision.getInstance();
		app.sendMessage("Neuron Cleaning...");

		long t1 = System.currentTimeMillis();

		neuronFile = new NeuronFile(neuronFileName);
		cleanNeurons(neuronFile, minSpikesCount, maxContamination, coincidenceTime);
		neuronFile.undeleteNeurons();
		neuronFile.close();

		long t2 = System.currentTimeMillis();
		app.sendMessage("Done in: " + (t2 - t1) / 1000. + " s.");

		Vision.getInstance().getCalculationManager().calculationDone();
	}

	public void cleanNeurons(NeuronFile neuronFile, int minSpikesCount, double maxContamination,
			final int coincidenceTime) throws IOException {

		neuronFile.undeleteNeurons();

		IntegerList neuronsList = new IntegerList();
		int[] neurons = neuronFile.getIDList();
		for (int i = 0; i < neurons.length; i++) {
			neuronsList.add(neurons[i]);
		}
		IntegerList[] neuronsLists = removeLowSpikeAndContam(neuronFile, neuronsList, minSpikesCount, maxContamination);

		for (int i = 0; i < neuronsLists[0].size(); i++) {
			neuronFile.deleteNeuron(neuronsLists[0].get(i));
		}

		neuronsLists = removeCorrelated(neuronFile, neuronsLists[1], coincidenceTime, maxCorrelation);

		for (int i = 0; i < neuronsLists[0].size(); i++) {
			neuronFile.deleteNeuron(neuronsLists[0].get(i));
		}

		// generated the final neuron file
		File f = new File(StringUtil.removeExtension(neuronFileName) + ".neurons");
		if (f.exists() && f.isFile()) {
			f.delete();
		}

		// create the new neuron file
		VisionHeader h = neuronFile.getHeader();
		h.minNeuronSpikes = minSpikesCount;
		h.acfT1 = VisionParams.ACFT1;
		h.acfT2 = VisionParams.ACFT2;
		h.maxContamination = maxContamination;
		h.coincidenceTime = coincidenceTime;
		h.maxCorrelation = maxCorrelation;
		h.removeDuplicates = 1; // shlens
		NeuronFile finalNeuronFile = new NeuronFile(f.getAbsolutePath(), h, neuronFile.getHeaderCapacity(),
				neuronFile.getTTLTimes());

		int[] idList = neuronFile.getIDList();
		for (int i = 0; i < idList.length; i++) {
			finalNeuronFile.addNeuron(neuronFile.getElectrode(idList[i]), idList[i],
					neuronFile.getSpikeTimes(idList[i]), neuronFile.getSpikeCount(idList[i]));
		}
		finalNeuronFile.close();

	}

	/**
	 * @param neuronFile
	 * @param idList
	 * @param minSpikesCount
	 * @param maxContamination
	 * @return
	 * @throws IOException
	 * 
	 *             Do first pass cleaning of neurons to reduce memory
	 *             constraints.
	 */
	public static IntegerList[] removeLowSpikeAndContam(NeuronFile neuronFile, IntegerList idList, int minSpikesCount,
			double maxContamination) throws IOException {

		Vision app = Vision.getInstance();
		IntegerList toRemoveList = new IntegerList();
		IntegerList toKeepList = new IntegerList();

		int rawNeuronCount = idList.size();
		int afterSpikeCount = rawNeuronCount;
		int afterContam = rawNeuronCount;
		for (int i = 0; i < idList.size(); i++) {
			app.sendMessage("Checking for spike count and contamination: " + i + "(" + idList.get(i) + ") ");

			// remove slow neurons
			if (neuronFile.getSpikeCount(idList.get(i)) < minSpikesCount) {
				toRemoveList.add(idList.get(i));
				afterSpikeCount--;
				afterContam--;
				continue;
			} else if (AutocorrelationCalculator.getContamination(idList.get(i), neuronFile) >= maxContamination) {
				toRemoveList.add(idList.get(i));
				afterContam--;
				continue;
			} else {
				toKeepList.add(idList.get(i));
			}

		}
		System.out.println("Raw Neurons: " + rawNeuronCount);
		System.out.println("After spike count cut: " + afterSpikeCount);
		System.out.println("After contamination cut: " + afterContam);

		return new IntegerList[] { toRemoveList, toKeepList };
	}

	public static IntegerList[] removeCorrelated(NeuronFile neuronFile, IntegerList idList, int coincidenceTime,
			double maxCorrelation) throws IOException {
		Vision app = Vision.getInstance();
		IntegerList toKeep = new IntegerList();
		IntegerList toRemove = new IntegerList();

		// Now load all remaining neurons.
		HashMap<Integer, int[]> times = new HashMap<Integer, int[]>();
		for (int i = 0; i < idList.size(); i++) {
			times.put(idList.get(i), neuronFile.getSpikeTimes(idList.get(i)));
		}

		for (int i = 0; i < idList.size(); i++) {
			toKeep.add(idList.get(i));
		}

		// Remove duplicated neurons
		if (coincidenceTime != 0) {

			DoubleHistogram pHist = new DoubleHistogram("", 0, 1.1, 0.02);
			DoubleHistogram ccH = null;

			for (int i = 0; i < toKeep.size(); i++) {

				app.sendMessage("Checking for duplication: " + i + "(" + toKeep.get(i) + ") ");
				int[] t1 = times.get(toKeep.get(i));

				for (int j = i + 1; j < toKeep.size(); j++) {

					int[] t2 = times.get(toKeep.get(j));

					ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(t1, t2, 1, deltaT, ccH);

					// sums all bins within coincidence time * 2, to get forward
					// and
					// backward in time
					// keep the maximum
					double[] cc = ccH.toArray();
					double maxV = Double.NEGATIVE_INFINITY;
					for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
						double v = 0;
						for (int m = 0; m < 2 * coincidenceTime; m++) {
							v += cc[k + m];
						}
						if (v > maxV) {
							maxV = v;
						}
					}
					double p = maxV / Math.min(t1.length, t2.length);

					pHist.fill(p, 1);
					if (p >= maxCorrelation) {
						if (t1.length < t2.length) {
							toRemove.add(toKeep.get(i));
							toKeep.remove(i);
							i--;
							break;
						} else {
							toRemove.add(toKeep.get(j));
							toKeep.remove(j);
							j--;
						}
					}
				}
			}
		} // PlotUtilities.showData("", pHist, new HistogramStyle());

		System.out.println("Unique Neurons: " + toKeep.size());
		return new IntegerList[] { toRemove, toKeep };
	}

	public void setParameters(HashMap parameters) {
		neuronFileName = (String) parameters.get("Neuron_File");
		minSpikesCount = Integer.parseInt((String) parameters.get("Minimun Number of Spikes"));
		maxContamination = Double.parseDouble((String) parameters.get("Maximum Contamination"));
		coincidenceTime = Integer.parseInt((String) parameters.get("Coincidence Time"));
		maxCorrelation = Double.parseDouble((String) parameters.get("Maximum Correlation"));
	}

	// ----------------------------------------------------------------------------------------
	// ------------------ Custom additions - Matlab Compatibility and debug
	// -------------------
	// ----------------------------------------------------------------------------------------

	/*
	 * Custom function to analyze correlation value repartitions and feed to
	 * Matlab for analysis
	 * 
	 * @return matrix of spike train correlation coefficients for all neurons
	 * 
	 * @author Vincent Deo - Stanford University - 09/29/2015
	 */
	public static double[][] correlationMatrix(NeuronFile neuronFile, int coincidenceTime) throws IOException {
		int[] ids = neuronFile.getIDList();
		return correlationMatrix(neuronFile, ids, coincidenceTime);
	}

	/*
	 * Custom function to analyze correlation value repartitions and feed to
	 * Matlab for analysis
	 * 
	 * @return matrix of spike train correlation coefficients for selected
	 * neurons
	 * 
	 * @author Vincent Deo - Stanford University - 09/29/2015
	 */
	public static double[][] correlationMatrix(NeuronFile neuronFile, int[] ids, int coincidenceTime)
			throws IOException {
		Vision app = Vision.getInstance();

		double[][] corrMatrix = new double[ids.length][ids.length];

		HashMap<Integer, int[]> times = new HashMap<Integer, int[]>();
		for (int i = 0; i < ids.length; i++) {
			times.put(ids[i], neuronFile.getSpikeTimes(ids[i]));
		}

		// Remove duplicated neurons
		if (coincidenceTime != 0) {

			DoubleHistogram pHist = new DoubleHistogram("", 0, 1.1, 0.02);
			DoubleHistogram ccH = null;

			for (int i = 0; i < ids.length; i++) {
				int[] t1 = times.get(ids[i]);

				for (int j = i; j < ids.length; j++) {

					int[] t2 = times.get(ids[j]);

					ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(t1, t2, 1, deltaT, ccH);

					// sums all bins within coincidence time * 2, to get forward
					// and backward in time
					// keep the maximum
					double[] cc = ccH.toArray();
					double maxV = Double.NEGATIVE_INFINITY;
					for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
						double v = 0;
						for (int m = 0; m < 2 * coincidenceTime; m++) {
							v += cc[k + m];
						}
						if (v > maxV) {
							maxV = v;
						}
					}

					double p = maxV / Math.min(t1.length, t2.length);

					pHist.fill(p, 1);

					corrMatrix[i][j] = p;
					corrMatrix[j][i] = p;
				}
			}
		}
		return corrMatrix;
	}

	public static double[] getContam(NeuronFile neuronFile, int[] ids) throws IOException {
		double[] contamVal = new double[ids.length];
		for (int i = 0; i < ids.length; i++) {
			contamVal[i] = getContam(neuronFile, ids[i]);
		}
		return contamVal;
	}

	private static double getContam(NeuronFile neuronFile, int id) throws IOException {
		return AutocorrelationCalculator.getContamination(id, neuronFile);
	}

	public static double getContam(int[] spikes, int nSamples) {
		return AutocorrelationCalculator.getContamination(spikes, nSamples);
	}

	public static double[] getCorrVal(NeuronFile neuronFile, int[] id1, int[] id2, int coincidenceTime)
			throws IOException {
		double[] corrVal = new double[id1.length];
		for (int i = 0; i < id1.length; i++) {
			corrVal[i] = getCorrVal(neuronFile, id1[i], id2[i], coincidenceTime);
		}
		return corrVal;
	}

	private static double getCorrVal(NeuronFile neuronFile, int id1, int id2, int coincidenceTime) throws IOException {
		DoubleHistogram ccH = null;

		int[] t1 = neuronFile.getSpikeTimes(id1);
		int[] t2 = neuronFile.getSpikeTimes(id2);
		ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(t1, t2, 1, deltaT, ccH);
		double[] cc = ccH.toArray();
		double maxV = Double.NEGATIVE_INFINITY;
		for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
			double v = 0;
			for (int m = 0; m < 2 * coincidenceTime; m++) {
				v += cc[k + m];
			}
			if (v > maxV) {
				maxV = v;
			}
		}
		return maxV / Math.min(t1.length, t2.length);
	}

	public static double getCorrVal(int[] sp1, int[] sp2, int coincidenceTime) {
		DoubleHistogram ccH = null;
		ccH = CrossCorrelationCalculator.getCrossCorrelationHistogram(sp1, sp2, 1, deltaT, ccH);
		double[] cc = ccH.toArray();
		double maxV = Double.NEGATIVE_INFINITY;
		for (int k = 0; k < cc.length - 2 * coincidenceTime; k++) {
			double v = 0;
			for (int m = 0; m < 2 * coincidenceTime; m++) {
				v += cc[k + m];
			}
			if (v > maxV) {
				maxV = v;
			}
		}
		return maxV / Math.min(sp1.length, sp2.length);
	}

}
