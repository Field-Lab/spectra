package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.LinkedList;

import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionParams;

/**
 * A class that allows access to the stored spike information (and TTL
 * information). There is no assumption made about the file size. The file is
 * not entirely read into memory, only its header is read. The file can be read
 * in two ways: 1) reading all spike times on a given electrode or 2) iterating
 * over all the spikes.
 * 
 * @see #getSpikeTimes(int electrode)
 * @see #iterator()
 * @author Charles A. Loomis, University of California, Santa Cruz<br>
 *         Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, UCSC
 */
public class SpikeFile {
	private static final int MAGIC_NUMBER = 0xABCDEF;
	private OpenMode accessMode;

	private final int spikeSize = 2 + 4;
	private int headerSize;
	private int nElectrodes;

	private VisionHeader header;

	private long ttlLocation;
	private int[] spikesPerElectrode;
	private int nSpikes;
	private IntegerList ttlArray;
	private String fileNameNoExtension;
	private int[][] spikeTimes;
	private RandomAccessFile f;

	/**
	 * Creates s SpikeFile which allows access to the spike information. There
	 * is no assumption made about the file size. The file is not entirely read
	 * into memory but only it's header.
	 */
	public SpikeFile(String fileName, int arrayID, float meanTimeConstant, float threshold, int nSamples,
			int samplingFrequency) throws IOException {

		accessMode = OpenMode.WRITE;

		ElectrodeMap m = ElectrodeMapFactory.getElectrodeMap(arrayID);
		this.nElectrodes = m.getNumberOfElectrodes();

		header = new VisionHeader();
		header.magic = MAGIC_NUMBER;
		header.headerVersion = 1; // -- shlens HACK
		header.version = 1;
		header.meanTimeConstant = meanTimeConstant;
		header.threshold = threshold;
		header.arrayID = arrayID;
		header.nSamples = nSamples;
		header.samplingFrequency = samplingFrequency;
		header.visionVersion = VisionParams.VERSION;

		// System.out.println("header size " + header.size()); //-- shlens

		this.ttlLocation = -1;
		this.spikesPerElectrode = new int[nElectrodes];
		this.nSpikes = 0;
		this.ttlArray = new IntegerList();

		File file = new File(fileName);
		this.fileNameNoExtension = StringUtil.removeExtension(file.getAbsolutePath());

		f = new RandomAccessFile(file, "rw");
		writeHeader();
		f.seek(headerSize);
	}

	/**
	 * Creates a SpikeFile which allows access to the spike information. There
	 * is no assumption made about the file size. The file is not entirely read
	 * into memory but only it's header.
	 */
	public SpikeFile(String fileName) throws IOException {
		accessMode = OpenMode.READ;
		File file = new File(fileName);
		this.fileNameNoExtension = StringUtil.removeExtension(file.getAbsolutePath());
		if (file.length() == 0)
			throw new IOException("Sorry, the file " + file.getName() + " is empty.");

		f = new RandomAccessFile(file, "r");
		readHeader();
		nElectrodes = ElectrodeMapFactory.getElectrodeMap(header.arrayID).getNumberOfElectrodes();
	}

	private void readHeader() throws IOException {
		f.seek(0);
		header = new VisionHeader();
		IOUtil.readPublicFields(header, f);

		if (header.magic != MAGIC_NUMBER)
			throw new IOException("The spike file does not start with the magic number");

		ttlLocation = f.readLong();
		// System.err.println("TTL location " + ttlLocation);

		// read the spike counts
		int nElectrodes = f.readInt();
		spikesPerElectrode = new int[nElectrodes];
		nSpikes = 0;
		for (int e = 0; e < nElectrodes; e++) {
			spikesPerElectrode[e] = f.readInt();
			if (e != 0)
				nSpikes += spikesPerElectrode[e];
		}

		headerSize = (int) f.getFilePointer();
		// System.err.println("Header Size: " + headerSize);

		// now read the TTLs
		ttlArray = new IntegerList();
		f.seek(ttlLocation);
		for (int i = 0; i < spikesPerElectrode[0]; i++) {
			ttlArray.add(f.readInt());
		}
	}

	private void writeHeader() throws IOException {
		f.seek(0);
		IOUtil.writePublicFields(header, f);

		spikesPerElectrode[0] = ttlArray.size();

		ttlLocation = f.length();
		f.writeLong(ttlLocation);

		// write the spike counts
		f.writeInt(nElectrodes);
		for (int e = 0; e < nElectrodes; e++) {
			f.writeInt(spikesPerElectrode[e]);
		}

		headerSize = (int) f.getFilePointer();

		// now write the TTLs
		if (ttlLocation != -1) {
			f.seek(ttlLocation);
			for (int i = 0; i < ttlArray.size(); i++) {
				f.writeInt(ttlArray.get(i));
			}
		}
	}

	public void addSpike(short electrode, int time) throws IOException {
		if (accessMode != OpenMode.WRITE)
			throw new Error("Cannot add spikes in READ mode");

		if (electrode == 0) {
			ttlArray.add(time);
		} else {
			f.writeShort(electrode);
			f.writeInt(time);
			// f.writeFloat(spike.amplitude);
			// f.writeFloat(spike.width);

			spikesPerElectrode[electrode]++;
		}
	}

	public int getArrayID() {
		return header.arrayID;
	}

	public int getNumberOfSamples() {
		return header.nSamples;
	}

	public int[] getTTLTimes() {
		return ttlArray.toArray();
	}

	public int getFirstTTL() {
		return ttlArray.get(0);
	}

	public int getLastTTL() {
		return ttlArray.get(ttlArray.size() - 1);
	}

	public int getSpikesCount() {
		return nSpikes;
	}

	public int getSpikesCount(int electrode) {
		return spikesPerElectrode[electrode];
	}

	/**
	 * Return the fileName as the String for this object.
	 */
	public String getFileNameNoExtension() {
		return fileNameNoExtension;
	}

	public double getThreshold() {
		return header.threshold;
	}

	public double getMeanTimeConstant() {
		return header.meanTimeConstant;
	}

	public float getSamplingFrequency() {
		return header.samplingFrequency;
	}

	public int[] getSpikeTimes(int electrode) throws IOException {
		int[] times = new int[getSpikesCount(electrode)];
		SpikeIterator iter = iterator();
		int spikeIndex = 0;

		while (iter.hasNext()) {
			Spike spike = iter.next();
			if (spike.electrode != electrode) {
				continue;
			}
			times[spikeIndex] = spike.time;
			spikeIndex++;
		}
		return times;
	}

	public int[][] getSpikeTimes() throws IOException {
		int[][] times = new int[nElectrodes][];
		for (int i = 0; i < times.length; i++) {
			times[i] = new int[getSpikesCount(i)];
		}
		int[] index = new int[times.length];

		SpikeIterator iter = iterator();
		while (iter.hasNext()) {
			Spike s = iter.next();
			times[s.electrode][index[s.electrode]] = s.time;
			index[s.electrode]++;
		}

		return times;
	}

	// ----------------------- Vincent Deo - Stanford University - 04/30/2015
	public int[] getSpikesTimesUntil(int startTime, int stopTime, int electrode) throws IOException {
		LinkedList<Integer> times = new LinkedList<Integer>();
		SpikeIterator iter = iterator();
		int spikeIndex = 0;

		while (iter.hasNext()) {
			Spike spike = iter.next();
			if (spike.electrode != electrode || spike.time < startTime)
				continue;
			if (spike.time >= stopTime)
				break;
			times.addLast((Integer) spike.time);
			spikeIndex++;
		}
		int[] output = new int[spikeIndex];
		for (int i = 0; i < spikeIndex; i++) {
			output[i] = times.removeFirst().intValue();
		}
		return output;
	}

	public int[][] getSpikesTimesUntil(int startTime, int stopTime) throws IOException {
		@SuppressWarnings("unchecked")
		LinkedList<Integer>[] times = (LinkedList<Integer>[]) new LinkedList<?>[nElectrodes];
		for (int i = 0; i < nElectrodes; i++) {
			times[i] = new LinkedList<Integer>();
		}
		SpikeIterator iter = iterator();
		int[] spikeIndex = new int[nElectrodes];

		while (iter.hasNext()) {
			Spike spike = iter.next();
			if (spike.time < startTime)
				continue;
			if (spike.time >= stopTime)
				break;
			times[spike.electrode].addLast((Integer) spike.time);
			spikeIndex[spike.electrode]++;
		}
		int[][] output = new int[nElectrodes][];
		for (int el = 0; el < nElectrodes; el++) {
			output[el] = new int[spikeIndex[el]];
			for (int i = 0; i < spikeIndex[el]; i++) {
				output[el][i] = times[el].removeFirst().intValue();
			}
		}
		return output;
	}

	public void pushMatlabSpikes(int[][] spikes) throws IOException {
		for (int i = 0; i < spikes.length; i++) {
			addSpike((short) spikes[i][2], spikes[i][1]);
		}
	}

	// -----------------------------------------------------------------------

	public SpikeIterator iterator(int startSample) throws IOException {
		SpikeIterator i = new SpikeIteratorImp();

		while (i.hasNext()) {
			Spike s = i.next();
			if (s.time >= startSample) {
				break;
			}
		}

		return i;
	}

	public SpikeIterator iterator() throws IOException {
		if (accessMode != OpenMode.READ) {
			throw new Error("Cannot get iterator in WRITE mode");
		}
		return new SpikeIteratorImp();
	}

	/**
	 * This iterator returns all of the spikes. Be careful, the spike returned
	 * from this iterator are reused. You should not keep a reference to a
	 * returned spike after the subsequent call to the next() method.
	 */
	private class SpikeIteratorImp implements SpikeIterator {

		private byte[] buffer = new byte[4096 * spikeSize];
		private int index, bytesRead;
		private Spike spike = new Spike(0, 0, 0);
		private int spikeIndex = 0;

		public SpikeIteratorImp() throws IOException {
			f.seek(headerSize);
		}

		public boolean hasNext() {
			return spikeIndex < nSpikes;
		}

		public Spike next() {
			if (spikeIndex >= nSpikes) { // no more spikes !
				return null;
			}

			if (index == bytesRead) {
				try {
					bytesRead = f.read(buffer);
				} catch (IOException ex) {
					ex.printStackTrace();
					return null;
				}

				// System.err.println("read " + bytesRead + ", spikeIndex " +
				// spikeIndex);
				if (bytesRead == -1 || bytesRead == 0) {
					return null;
				}
				index = 0;
			}

			int ch1, ch2, ch3, ch4;

			// read the electrode
			ch1 = buffer[index++] & 0xff;
			ch2 = buffer[index++] & 0xff;
			short el = (short) ((ch1 << 8) + (ch2 << 0));

			// read the time
			ch1 = buffer[index++] & 0xff;
			ch2 = buffer[index++] & 0xff;
			ch3 = buffer[index++] & 0xff;
			ch4 = buffer[index++] & 0xff;
			int t = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);

			spike.setValues(t, el, 0);
			spikeIndex++;
			return spike;
		}

		public Spike current() {
			return spike;
		}
	}

	public int getNumberOfElectrodes() {
		return nElectrodes;
	}

	public void close() throws IOException {
		if (accessMode == OpenMode.WRITE) {
			writeHeader();
		}
		f.close();
	}

	public VisionHeader getHeader() {
		return (VisionHeader) header.clone();
	}
}

enum OpenMode {
	READ, WRITE
};