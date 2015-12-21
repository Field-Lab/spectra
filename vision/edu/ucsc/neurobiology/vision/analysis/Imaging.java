package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.text.*;
import java.util.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Imaging
    extends AbstractCalculation implements SampleListener {

    private String neuronFileName, rawDataSource;
    private int nlPoints, nrPoints, nPoints;
    private int nElectrodes;
    private long nSpikes;
    private long sampleIndex;
    private LinkedList<short[]> sampleList;
    private DataContainer spikeData;
    private NeuronFile.SpikeStream1 spikeIterator;
    private NeuronFile.ExtendedSpike currentSpike;
    private float[] mean;
    private double alpha;
    private long spikeNumber = 1; // DO NOT CHANGE, the algorithm needs a starting value of 1
    private int spikesToUse;
    private MultipleCompressedSampleInputStream sampleInputStream;
    private RawDataHeader header;
    private NeuronFile neuronFile;
    private double[] step, remainder;
    private Vision app;
    private long t1;
    private JPanel diagnosticPanel;
    private java.util.Timer timer;
    private JLabel samplesReadLabel;
    private int numThreads;
    private int[] edgeIDs;
    private EICalculator[] eiCalcs;
    
    private PhysiologicalImagingFile eiFile;


    public Imaging() {
        if (Vision.isGUIBased()) {
            diagnosticPanel = new JPanel(new TableLayout(0, 0));
        }
    }


    public void startCalculation() throws IOException {
        app = Vision.getInstance();
        app.sendMessage("Preparing...");



        // get the raw data file
        sampleInputStream = new MultipleCompressedSampleInputStream(rawDataSource);
        header = sampleInputStream.getHeader();
        nElectrodes = header.getNumberOfElectrodes();

        mean = new float[nElectrodes];

        neuronFile = new NeuronFile(neuronFileName);
        spikeIterator = neuronFile.getSpikeStream();
        nSpikes = spikeIterator.getNumberOfSpikes();

        try {
            String eiFileName = StringUtil.removeExtension(neuronFileName) + VisionParams.EI_FILE_EXTENSION;
            eiFile = new PhysiologicalImagingFile(eiFileName, nlPoints, nrPoints, header.getArrayID());
        } catch (IOException e) {
            e.printStackTrace();
        }

        // set up the step for selecting spikes
        int[] idList = neuronFile.getIDList();
        int maxID = MathUtil.max(idList) + 1;  // Why not just get the last one; it's already sorted?
        step = new double[maxID];
        remainder = new double[maxID];
        for (int i = 0; i < idList.length; i++) {
            step[idList[i]] = neuronFile.getSpikeCount(idList[i]) / spikesToUse;
        }

        // load the first spike
        currentSpike = spikeIterator.next();
        spikeNumber++;

        spikeData = new DataContainer();
        spikeData.spike = new short[nPoints + 2][nElectrodes];

        // Flexible threading, could rewrite to use abstracted ParallelUtils for edgeIDs
        int i1 = 0;
        int i2;
        int istep = idList.length / numThreads;
        edgeIDs = new int[numThreads];
        eiCalcs = new EICalculator[numThreads];
        for (int i = 0; i < numThreads; i++) {
            i2 = istep * (i+1);
            if (i == (numThreads - 1)) i2 = idList.length;
            
            int[] l = new int[i2 - i1];
            System.arraycopy(idList, i1, l, 0, l.length);
            edgeIDs[i] = l[l.length-1];
            eiCalcs[i] = new EICalculator(nlPoints, nrPoints, l, nElectrodes);
            
            i1 = i2;
        }

        // create the sample buffer
        sampleList = new LinkedList<short[]> ();
        for (int i = 0; i < nPoints; i++) {
            sampleList.add(new short[nElectrodes]);
        }

        app.startProgressBar();
        timer = new java.util.Timer(true);
        timer.scheduleAtFixedRate(timerTask, 0, 100);

        // start
        t1 = System.currentTimeMillis();
        sampleInputStream.addSampleListener(this);
        sampleInputStream.start();
        for (EICalculator e : eiCalcs) e.start();

        app.sendMessage("Calculating EI...");
    }


    TimerTask timerTask = new TimerTask() {
        public void run() {
            int samplesRead = sampleInputStream.getSamplesRead();
            app.setProgress( (int) (100.0 * samplesRead / header.getNumberOfSamples()));

            // show GUI info if appropriate
            if (Vision.isGUIBased()) {
                NumberFormat f = DecimalFormat.getInstance();
                samplesReadLabel.setText(f.format(sampleInputStream.getSamplesRead()));
            }

            // finish the calculation
            if (sampleInputStream.isFinished()) {
                for (EICalculator e : eiCalcs) e.finish();
                while (anyEICalcsAlive()) {
                    try {
                        Thread.sleep(100);
                    } catch (InterruptedException e) {}
                }

                // save the covariance matrices
                try {
//                    String eiFileName = StringUtil.removeExtension(neuronFileName) +
//                                        ".ei";
//                    PhysiologicalImagingFile eiFile = new PhysiologicalImagingFile(
//                        eiFileName, nlPoints, nrPoints,
//                        sampleInputStream.getHeader().getArrayID());

                    for (EICalculator e : eiCalcs) {
                        for (Integer id : e.calc.keySet()) {
                            WaveformCalculator calc = e.calc.get(id);
                            calc.finish();
                            eiFile.appendImage(id, calc.getCount(), calc.getImage());
                        }
                    }
                    eiFile.close();

                } catch (IOException e) {
                    e.printStackTrace();
                }

                long t2 = System.currentTimeMillis();
                app.sendMessage("Done in: " + (t2 - t1) / 1000. + " s.");

                timer.cancel();
                app.endProgressBar();
                Vision.getInstance().getCalculationManager().calculationDone();
            }
        }
    };

    
    private boolean anyEICalcsAlive() {
        for (EICalculator e : eiCalcs) {
            if (e.isAlive()) return true;
        }
        return false;
    }

    
    public void processSample(short[] sample) {
        if (spikeNumber >= nSpikes) { // make sure we don't get a infinite loop
            return;
        }

        // filter the sample and add it into the sample queue
        short[] s = sampleList.removeFirst();
        System.arraycopy(sample, 0, s, 0, nElectrodes);
        for (int i = 0; i < nElectrodes; i++) {
            mean[i] += alpha * (s[i] - mean[i]);
            s[i] -= mean[i];
        }
        sampleList.addLast(s);

        // process possible spikes
        while (sampleIndex - currentSpike.time == nrPoints) {
            for (int i = 0; i < nPoints; i++) {
                s = sampleList.get(i);
                for (int e = 0; e < nElectrodes; e++) {
                    spikeData.spike[i][e] = s[e];
                }
            }

            spikeData.neuronID = currentSpike.neuronID;
            spikeData.time = currentSpike.time;

            for (int i = 0; i < edgeIDs.length; i++) {
                if (spikeData.neuronID <= edgeIDs[i]) {
                    eiCalcs[i].add(spikeData);
                    break;
                }
            }

            nextSpike();
            if (spikeNumber >= nSpikes) { // make sure we don't get a infinite loop
                break;
            }
        }

        sampleIndex++;
    }


    // select the next spike in a uniform fashion
    private void nextSpike() {
        while (spikeIterator.hasNext()) {
            currentSpike = spikeIterator.next();
            spikeNumber++;
            int id = currentSpike.neuronID;

            remainder[id]++;
            if (remainder[id] > step[id]) {
                remainder[id] -= step[id];
                return;
            }
        }
    }


    public void finishSampleProcessing() {}


    static class DataContainer {
        int neuronID;
        int time;
        short[][] spike;
    }
    

    public void setParameters(HashMap<String, String> parameters) {
        String datasetFolder = (String) parameters.get("Dataset Folder");
        rawDataSource = (String) parameters.get("Raw Data File");
        neuronFileName = datasetFolder + File.separator + new File(rawDataSource).getName() +
                         ".neurons";
        nlPoints = Integer.parseInt( (String) parameters.get("Left Samples"));
        nrPoints = Integer.parseInt( (String) parameters.get("Right Samples"));
        this.nPoints = nlPoints + nrPoints + 1;

        double timeConstant = Double.parseDouble( (String) parameters.get(
            "Mean Time Constant"));
        alpha = (float) (1.0 / (timeConstant * 20000.0));

        spikesToUse = Integer.parseInt( (String) parameters.get("Spikes To Average"));

        if (parameters.get("nThreads") == null) {
            numThreads = 2;
        } else {
            numThreads = Integer.parseInt( (String) parameters.get("nThreads"));
        }
        
        if (Vision.isGUIBased()) {
            // create the diagnostic panel
            diagnosticPanel.add(new JLabel("Samples read:"));
            samplesReadLabel = new JLabel("0", JLabel.RIGHT);
            diagnosticPanel.add(samplesReadLabel);
        }
    }


    public JComponent getDiagnosticPanel() {
        return diagnosticPanel;
    }

}
