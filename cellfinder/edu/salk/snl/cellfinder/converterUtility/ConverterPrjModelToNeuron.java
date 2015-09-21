package edu.salk.snl.cellfinder.converterUtility;

import java.io.File;
import java.util.HashMap;

import edu.salk.snl.cellfinder.cellfinder.CellFinderData;
import edu.salk.snl.cellfinder.cluster.ClusterData;
import edu.salk.snl.cellfinder.cluster.ClusterModel;
import edu.salk.snl.cellfinder.data.NeuronsData;
import edu.salk.snl.cellfinder.data.NeuronsModel;
import edu.salk.snl.cellfinder.sort.SpikeSorter;

/**
 * Custom class that allows shell-style conversion of .prj files to .neurons
 * files without recomputing the clustering, ie according to an existing neuron
 * file
 * 
 * The resulting neuroning might only consist in only the cores of the clusters
 * given that Vision's initial clustering forces a full partition of spikes
 * whereas CellFinder's reattribution requires high belonging probability.
 * 
 * @author Vincent Deo, Stanford University
 */
public class ConverterPrjModelToNeuron {

	/**
	 * Shell starter
	 * 
	 * @param args
	 *            {<prj-file>, <model-file>, <neur-raw-file>} OR {<root-path>}
	 *            in which case <root-path> must be a folder containing all 33
	 *            files.
	 */
	public static void main(String[] args) {
		if (args.length != 3 || args.length != 1) {
			System.out.println("usage: <prj-file> <model-file> <neur-raw-file>");
			System.exit(-1);
		}
		String prjFileName, modelFileName, neuronFileToWrite;
		if (args.length == 3) { // {<prj-file>, <model-file>, <neur-raw-file>
								// case
			// Assign names
			prjFileName = args[0];
			modelFileName = args[1];
			neuronFileToWrite = args[2];
		} else { // {<root-path>} case
			String argPath = args[0];
			// Check for ending / or \
			if (args[0].charAt(args[0].length() - 1) != File.separatorChar)
				argPath = args[0] + File.separatorChar;
			// Extract folder name
			File name = new File(argPath);
			String dataName = name.getName();
			// Assign file names
			prjFileName = argPath + dataName + ".prj";
			modelFileName = argPath + dataName + ".model";
			neuronFileToWrite = argPath + dataName + ".neurons-raw";
		}

		// Create cell finder data structure
		CellFinderData cellFinderData = new CellFinderData();

		// load up files - expect path / format / other errors to come up from
		// CellFinder here
		cellFinderData.loadData(prjFileName);
		cellFinderData.loadModel(modelFileName);

		// load up data
		NeuronsModel neuronsModel = cellFinderData.getNeuronsModel();
		NeuronsData neuronsData = cellFinderData.getNeuronsData();

		// Calculate the statistics of each electrode (e.g. contamination,
		// number
		// of spikes) and use this for filtering bad spikes.
		neuronsModel.updateModels(neuronsData);

		for (int electrode = 1; electrode < cellFinderData.getNumberElectrodes(); electrode++) {
			// determine if bad electrode
			if (neuronsModel.isBadElectrode(electrode)) {
				System.out.println("Unable to export sorted spikes. Bad electrode.");
				continue;
			}

			// load up electrode spikes
			ClusterData clusterData = neuronsData.loadElectrode(electrode);

			// select clustering model
			ClusterModel clusterModel = neuronsModel.loadModel(electrode);
			SpikeSorter spikeSorter = neuronsModel.loadSorter(electrode);

			// load clustering model
			clusterData.loadModel(clusterModel, spikeSorter);

			// Spike sort based on model
			clusterData.updateClusters();
		}

		// Do necessary CellFinder stuff
		HashMap<String, String> temp = new HashMap<String, String>();
		temp.put("minSpikeRate", "0.0");
		temp.put("maxContamination", "0.0");

		// Write and exit
		cellFinderData.saveNeurons(neuronFileToWrite, temp);
		System.out.println("Done saving neurons-raw file.");
		System.exit(0);
	}

}
