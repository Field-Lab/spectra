package edu.salk.snl.cellfinder.converterUtility;

import java.util.HashMap;

import edu.salk.snl.cellfinder.cellfinder.CellFinderData;
import edu.salk.snl.cellfinder.cluster.ClusterData;
import edu.salk.snl.cellfinder.cluster.ClusterModel;
import edu.salk.snl.cellfinder.data.NeuronsData;
import edu.salk.snl.cellfinder.data.NeuronsModel;
import edu.salk.snl.cellfinder.sort.SpikeSorter;

public class ConverterPrjModelToNeuron {

	public static void main(String[] args) {
		if (args.length != 3) {
			System.out.println("usage: <prj-file> <model-file> <neur-raw-file>");
			System.exit(-1);
		}

		String prjFileName = args[0];
		String modelFileName = args[1];
		String neuronFileToWrite = args[2];

		// create cell finder data structure
		CellFinderData cellFinderData = new CellFinderData();

		// load up files
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

		HashMap<String, String> temp = new HashMap<String, String>();
		temp.put("minSpikeRate", "0.0");
		temp.put("maxContamination", "0.0");

		cellFinderData.saveNeurons(neuronFileToWrite, temp);
		System.out.println("Done saving neurons-raw file.");
		System.exit(0);
	}

}
