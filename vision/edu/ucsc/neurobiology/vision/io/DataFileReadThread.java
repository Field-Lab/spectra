package edu.ucsc.neurobiology.vision.io;

import java.io.IOException;
import java.util.LinkedList;

/**
 * A data input thread class Written only for a Matlab backline purpose to
 * accelerate data access and prefiltering, which are done in this thread and
 * not in the Matlab main thread.
 * 
 * Requires Vision's RawDataFile class.
 * 
 * @author Vincent Deo - Stanford University - 07/08/2015
 */
public class DataFileReadThread extends Thread {

	// Input management
	private final int nElectrodes;
	private RawDataFile rawDataFile;

	// Buffer management
	private short[][] dataBuffer = null;
	private float[][] filteredBuffer = null;

	// Filtering
	private float[] filterState;
	private final float alpha;

	// Matlab command pile
	// IO thread pops commands and executes them
	// Matlab thread adds commands at the tail
	// Only getter commands need to wait until execution
	private LinkedList<Command> commandList = new LinkedList<Command>();

	/**
	 * Constructor associates basic parameters only
	 * 
	 * @param rawDataFile
	 *            wrapper of data source
	 * @param nEl
	 *            number of electrodes (also available with
	 *            "rawDataFile.getHeader().getNumberOfElectrodes()" )
	 * @param alpha
	 *            DC removal filter parameters
	 */
	public DataFileReadThread(RawDataFile rawDataFile, int nEl, float alpha) {
		super();
		this.rawDataFile = rawDataFile;
		nElectrodes = nEl;
		filterState = new float[nEl];
		this.alpha = alpha;

		this.setPriority(MAX_PRIORITY);
	}

	/**
	 * Safe-typing of possible commands
	 */
	private enum CommandType {
		LOAD, FILTER, RETURN, RETURNFILTERED, SETFILTER, QUIT, CLEAR
	}

	/**
	 * Wrapper of Commands passed by Matlab
	 */
	public class Command {
		public final CommandType commandType;

		public int startSample;
		public int bufferLength;

		public float[] filterStateToSet;

		public Command(CommandType type) {
			commandType = type;
		}
	}

	/**
	 * Thread run method Waits for Matlab's command to be available in the
	 * command structure Then executes and wait until a "quit" command is passed
	 */
	@Override
	public void run() {

		loop: while (true) {
			Command c;
			synchronized (commandList) {
				if (commandList.isEmpty())
					try {
						commandList.wait();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				c = commandList.pollFirst();
			} // synchronized (commandList)

			switch (c.commandType) {
			case LOAD:
				if (dataBuffer == null || dataBuffer.length != c.bufferLength
						|| dataBuffer[1].length != nElectrodes)
					dataBuffer = new short[c.bufferLength][nElectrodes];
				try {
					rawDataFile.getData(c.startSample, dataBuffer);
				} catch (IOException e) {
					e.printStackTrace();
				}
				break;
			case FILTER:
				if (dataBuffer != null) {
					if (filteredBuffer == null
							|| filteredBuffer.length != dataBuffer.length
							|| filteredBuffer[1].length != nElectrodes)
						filteredBuffer = new float[dataBuffer.length][nElectrodes];
					filter();
				}
				break;
			case RETURN:
				// System.out.println("IOT: trying to catch lock on c.");
				synchronized (c) {
					// System.out.println("IOT: lock on c catched.");
					c.notifyAll();
				}
				break;
			case RETURNFILTERED:
				// System.out.println("IOT: trying to catch lock on c.");
				synchronized (c) {
					// System.out.println("IOT: lock on c catched.");
					c.notifyAll();
				}
				break;
			case SETFILTER:
				filterState = c.filterStateToSet;
				break;
			case CLEAR:
				dataBuffer = null;
				filteredBuffer = null;
				break;
			case QUIT:
				break loop;
			default:
				throw new Error("Unrecognized command at switch");
			}
		} // while(true)
			// System.out.println("IOT: End of thread.");
	}

	// -------------------------- COMMAND SETTERS --------------------------
	// ------------------------- MATLAB CALL ONLY -------------------------
	public void loadNextBuffer(int startSample, int bufferLength,
			boolean filterTag) {
		Command c = new Command(CommandType.LOAD);
		c.startSample = startSample;
		c.bufferLength = bufferLength;

		Command f = null;
		if (filterTag)
			f = new Command(CommandType.FILTER);

		synchronized (commandList) {
			commandList.addLast(c);
			if (f != null)
				commandList.addLast(f);
			commandList.notifyAll();
		}
	}

	public void filterBuffer() {
		Command c = new Command(CommandType.FILTER);

		synchronized (commandList) {
			commandList.addLast(c);
			commandList.notifyAll();
		}
	}

	// Getters == wait on Command, IOThread will notify
	// Non returning commands - let it go in the command pile.
	public short[][] getShortBuffer() throws InterruptedException {
		Command c = new Command(CommandType.RETURN);
		if (this.isAlive()) {
			synchronized (c) {
				synchronized (commandList) {
					commandList.addLast(c);
					commandList.notifyAll();
				}
				// System.out.println("MT: Waiting on command");
				// long s = System.currentTimeMillis();
				c.wait();
				// long f = System.currentTimeMillis();
				// System.out.println("Waited for buffer " + (f - s) + " ms.");
			}
			return dataBuffer;
		} else
			throw new Error("MT: Buffer request on dead IO Thread");
	}

	public float[][] getFilteredBuffer() throws InterruptedException {
		Command c = new Command(CommandType.RETURNFILTERED);
		if (this.isAlive()) {
			synchronized (c) {
				synchronized (commandList) {
					commandList.addLast(c);
					commandList.notifyAll();
				}
				// System.out.println("MT: Waiting on command");
				// long s = System.currentTimeMillis();
				c.wait();
				// long f = System.currentTimeMillis();
				// System.out.println("Waited for buffer " + (f - s) + " ms.");
			}
			return filteredBuffer;
		} else
			throw new Error("MT: Buffer request on dead IO Thread");

	}

	public void setFilterState(float[] state) {
		Command c = new Command(CommandType.SETFILTER);
		c.filterStateToSet = state.clone();
		synchronized (commandList) {
			commandList.addLast(c);
			commandList.notifyAll();
		}
	}

	public void setQuit() {
		Command c = new Command(CommandType.QUIT);
		synchronized (commandList) {
			commandList.addLast(c);
			commandList.notifyAll();
		}
	}

	public void clearBuffers() {
		Command c = new Command(CommandType.CLEAR);
		synchronized (commandList) {
			commandList.addLast(c);
			commandList.notifyAll();
		}
	}

	// -------------------------- INTERNAL METHODS --------------------------
	// ------------------------- IOTHREAD CALL ONLY -------------------------

	private void filter() {
		for (int el = 0; el < nElectrodes; el++) {
			for (int s = 0; s < dataBuffer.length; s++) {
				filterState[el] += alpha
						* ((float) dataBuffer[s][el] - filterState[el]);
				filteredBuffer[s][el] = (float) dataBuffer[s][el]
						- filterState[el];
			}
		}
	}
	// ----------------------------------------------------------------------
}