package edu.ucsc.neurobiology.vision.io;

import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;

import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap;
import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory;

/**
 * This class is an interface to the projections file created by
 * PCANFProjections. It contains the spike data projected on the PCA
 * eigenvectors obtained from all spikes on a given electrode. For an example of
 * using this class look into PCANFClustering.
 * 
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ProjectionsFile {
	private RandomAccessFile projectionsFile;

	private final static int MAGIC = 0xBBBBBB;
	public final int HEADER_SIZE;
	public final static int ENTRY_SIZE = 8 + 4;
	public final static int writeBufferSize = 5000;
	public int maxSpikesPerElectrode;

	private VisionHeader header;
	private int nElectrodes;
	private byte[] readBuffer;
	private byte[] writeBuffer;
	private long[] startLocation;
	private int[] writtenSpikes;
	private int[] ttlTimes;

	public ProjectionsFile(String prjFileName) throws IOException {
		projectionsFile = new RandomAccessFile(prjFileName, "r");

		header = new VisionHeader();
		IOUtil.readPublicFields(header, projectionsFile);

		// check the header
		if (header.magic != MAGIC) {
			throw new IOException("This .prj file does not contain the magic number");
		}

		// read the TTLs
		int nTTL = projectionsFile.readInt();
		ttlTimes = new int[nTTL];
		for (int i = 0; i < nTTL; i++) {
			ttlTimes[i] = projectionsFile.readInt();
		}
		HEADER_SIZE = (int) projectionsFile.getFilePointer();

		// FIX1
		// nElectrodes =
		// ElectrodeMapFactory.getNumberOfElectrodes(header.arrayID);
		ElectrodeMap m = ElectrodeMapFactory.getElectrodeMap(header.arrayID);
		nElectrodes = m.getNumberOfElectrodes();

		int maxSpikesPerElectrode = Integer.MIN_VALUE;
		for (int electrode = 0; electrode < nElectrodes; electrode++) {
			projectionsFile.seek(HEADER_SIZE + electrode * ENTRY_SIZE);
			long seekLocation = projectionsFile.readLong(); // do not remove
			int nSpikes = projectionsFile.readInt();
			if (nSpikes > maxSpikesPerElectrode) {
				maxSpikesPerElectrode = nSpikes;
			}
		}
		this.maxSpikesPerElectrode = maxSpikesPerElectrode;

		readBuffer = new byte[maxSpikesPerElectrode * (header.nDimensions + 1) * 4];
	}

	public ProjectionsFile(String prjFileName, VisionHeader header, int[] spikeCount, int[] ttlTimes)
			throws IOException {
		// check the header
		header.magic = MAGIC;
		header.version = 1;

		this.header = header;
		this.projectionsFile = new RandomAccessFile(prjFileName, "rw");
		IOUtil.writePublicFields(header, projectionsFile);

		// write the TTLs
		projectionsFile.writeInt(ttlTimes.length);
		for (int i = 0; i < ttlTimes.length; i++)
			projectionsFile.writeInt(ttlTimes[i]);

		HEADER_SIZE = (int) projectionsFile.getFilePointer();

		ElectrodeMap m = ElectrodeMapFactory.getElectrodeMap(header.arrayID);
		nElectrodes = m.getNumberOfElectrodes();

		startLocation = new long[nElectrodes];
		long index = HEADER_SIZE + ENTRY_SIZE * nElectrodes;
		for (int i = 0; i < nElectrodes; i++) {
			if (spikeCount[i] != 0) {
				startLocation[i] = index;
				index += spikeCount[i] * (header.nDimensions + 1) * 4;
			} else {
				startLocation[i] = -1;
			}
		}

		writtenSpikes = new int[nElectrodes];
		this.maxSpikesPerElectrode = -1;
		writeBuffer = new byte[writeBufferSize * (header.nDimensions + 1) * 4];
	}

	public int[] getTTLTimes() {
		return ttlTimes.clone();
	}

	// -----------------------------------------------------------------
	// Vincent Deo - Stanford University - 05/19/2015
	// for Matlab compatibility
	public class ProjectionsToMatlab {
		public float[][] data;
		public int[] times;

		ProjectionsToMatlab(int dims, int max) {
			data = new float[dims][max];
			times = new int[max];
		}
	}

	public ProjectionsToMatlab readProjections(int electrode) throws IOException {
		ProjectionsToMatlab output = new ProjectionsToMatlab(header.nDimensions, maxSpikesPerElectrode);
		readProjections(electrode, output.data, output.times);
		return output;
	}

	// -----------------------------------------------------------------

	synchronized public int readProjections(int electrode, float[][] data, int[] times) throws IOException {
		projectionsFile.seek(HEADER_SIZE + electrode * ENTRY_SIZE);
		long seekLocation = projectionsFile.readLong();
		if (seekLocation == -1)
			return 0; // This is the no spikes case. see how seekLocation is
						// set.
		int nSpikes = projectionsFile.readInt();

		projectionsFile.seek(seekLocation);
		projectionsFile.readFully(readBuffer, 0, nSpikes * (header.nDimensions + 1) * 4);

		int byteIndex = 0, ch1, ch2, ch3, ch4;
		for (int i = 0; i < nSpikes; i++) {
			// reconstruct the spike time
			ch1 = readBuffer[byteIndex++] & 0xff;
			ch2 = readBuffer[byteIndex++] & 0xff;
			ch3 = readBuffer[byteIndex++] & 0xff;
			ch4 = readBuffer[byteIndex++] & 0xff;
			if ((ch1 | ch2 | ch3 | ch4) < 0) {
				throw new EOFException();
			}
			times[i] = (ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0);

			// reconstruct the projections
			for (int d = 0; d < header.nDimensions; d++) {
				ch1 = readBuffer[byteIndex++] & 0xff;
				ch2 = readBuffer[byteIndex++] & 0xff;
				ch3 = readBuffer[byteIndex++] & 0xff;
				ch4 = readBuffer[byteIndex++] & 0xff;
				if ((ch1 | ch2 | ch3 | ch4) < 0) {
					throw new EOFException();
				}
				data[d][i] = Float.intBitsToFloat((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
			}
		}

		return nSpikes;
	}

	synchronized public void saveData(final int electrode, int nSpikesToSave, int[][] spikeTimes,
			float[][][] spikeProjections) throws IOException {

		// prepare the data to be written
		int byteIndex = 0, intBits;
		for (int i = 0; i < nSpikesToSave; i++) {
			// convert to bytes the spike time
			intBits = spikeTimes[electrode][i];
			writeBuffer[byteIndex++] = (byte) ((intBits >>> 24) & 0xFF);
			writeBuffer[byteIndex++] = (byte) ((intBits >>> 16) & 0xFF);
			writeBuffer[byteIndex++] = (byte) ((intBits >>> 8) & 0xFF);
			writeBuffer[byteIndex++] = (byte) ((intBits >>> 0) & 0xFF);

			// convert to bytes the projections
			for (int d = 0; d < header.nDimensions; d++) {
				intBits = Float.floatToIntBits(spikeProjections[electrode][i][d]);
				writeBuffer[byteIndex++] = (byte) ((intBits >>> 24) & 0xFF);
				writeBuffer[byteIndex++] = (byte) ((intBits >>> 16) & 0xFF);
				writeBuffer[byteIndex++] = (byte) ((intBits >>> 8) & 0xFF);
				writeBuffer[byteIndex++] = (byte) ((intBits >>> 0) & 0xFF);
			}
		}

		// write the data
		projectionsFile.seek(startLocation[electrode] + writtenSpikes[electrode] * (header.nDimensions + 1) * 4);
		projectionsFile.write(writeBuffer, 0, byteIndex);
		writtenSpikes[electrode] += nSpikesToSave;
	}

	public void close() throws IOException {
		// write the header
		if (writtenSpikes != null) {
			projectionsFile.seek(HEADER_SIZE);
			for (int i = 0; i < nElectrodes; i++) {
				projectionsFile.writeLong(startLocation[i]);
				projectionsFile.writeInt(writtenSpikes[i]);
			}
		}
		projectionsFile.close();
	}

	public VisionHeader getHeader() {
		return (VisionHeader) header.clone();
	}

}