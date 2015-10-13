package edu.ucsc.neurobiology.vision.io;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.analysis.MakeParametersFile;
import edu.ucsc.neurobiology.vision.analysis.PCACalculator;
import edu.ucsc.neurobiology.vision.io.tags.DoubleArrayTag;
import edu.ucsc.neurobiology.vision.io.tags.DoubleTag;
import edu.ucsc.neurobiology.vision.io.tags.NeuroInputStream;
import edu.ucsc.neurobiology.vision.io.tags.NeuroOutputStream;
import edu.ucsc.neurobiology.vision.io.tags.ObjectEmbeddingTag;
import edu.ucsc.neurobiology.vision.io.tags.StringTag;
import edu.ucsc.neurobiology.vision.math.ArrayMeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.Num;
import edu.ucsc.neurobiology.vision.math.expressions.Expression;
import edu.ucsc.neurobiology.vision.neuronviewer.NeuronViewer;
import edu.ucsc.neurobiology.vision.util.IntegerList;

/**
 * This class is used to create and read parameter files. Each neuron is
 * associated with an ID (same as in the neurons file) and with a list of
 * parameters. The parameters can be strings, doubles and double[]. Parameter
 * files are created by MakeParametersFile and used primarily by the
 * NeuronViewer.
 * 
 * @see MakeParametersFile
 * @see NeuronViewer
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, The Salk Institute
 */
public class ParametersFile {
	private ArrayList<ObjectEmbeddingTag[]> table;
	private RandomAccessFile file;
	private NeuroInputStream nInput;
	private NeuroOutputStream nOutput;
	private int nColumns;
	private int maxRows;
	private String[][] columnNamesAndTypes;
	private LinkedHashMap<String, Integer> columnNamesMap;
	private boolean needsFlush;

	/**
	 * Create a new file.
	 * 
	 * @param fileName
	 * @param columnNamesAndTypes
	 * @param maxRows
	 */
	public ParametersFile(String fileName, String[][] columnNamesAndTypes, int maxRows) throws IOException {
		needsFlush = true;
		File f = new File(fileName);
		if (f.exists()) {
			f.delete();
		}

		this.columnNamesAndTypes = columnNamesAndTypes;
		this.nColumns = columnNamesAndTypes.length;
		this.maxRows = maxRows;
		this.columnNamesMap = new LinkedHashMap<String, Integer>();
		for (int i = 0; i < nColumns; i++) {
			columnNamesMap.put(columnNamesAndTypes[i][0], new Integer(i));
		}

		openFile(fileName);
	}

	public ParametersFile(String fileName) throws IOException {
		this(new File(fileName));
	}

	/**
	 * Opens an existing params file.
	 * 
	 * @param f
	 *            File
	 * @throws IOException
	 */
	public ParametersFile(File f) throws IOException {
		needsFlush = false;
		openFile(f.getAbsolutePath());

		this.nColumns = file.readInt();
		int nRows = file.readInt();
		this.maxRows = file.readInt();
		columnNamesAndTypes = new String[nColumns][2];
		columnNamesMap = new LinkedHashMap<String, Integer>();
		for (int i = 0; i < nColumns; i++) {
			columnNamesAndTypes[i][0] = nInput.readString();
			columnNamesAndTypes[i][1] = nInput.readString();
			columnNamesMap.put(columnNamesAndTypes[i][0], new Integer(i));

			// System.err.println(columnNamesAndTypes[i][0] + " - " +
			// columnNamesAndTypes[i][1]);
		}

		// Reading as a block here significantly improves performance
		// in OS 10.5, Java 1.5, when reading over NFS. OS 10.5 with
		// NFS appears to not have sufficient buffering.
		int[][] seekLocations = new int[nRows][nColumns];
		byte[] buffer = new byte[4 * nRows * nColumns];
		file.readFully(buffer);

		DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
		for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nColumns; j++) {
				seekLocations[i][j] = dis.readInt();
			}
		}

		file.seek(seekLocations[0][0]);
		buffer = new byte[(int) (file.length() - seekLocations[0][0])];
		file.readFully(buffer);
		ByteArrayInputStream bb = new ByteArrayInputStream(buffer);
		NeuroInputStream input = new NeuroInputStream(bb);

		for (int i = 0; i < nRows; i++) {
			ObjectEmbeddingTag[] row = new ObjectEmbeddingTag[nColumns];
			for (int j = 0; j < nColumns; j++) {
				row[j] = (ObjectEmbeddingTag) input.readTag();
			}
			table.add(row);
		}
	}

	public boolean hasParameter(String paramName) {
		for (int i = 0; i < nColumns; i++) {
			if (columnNamesAndTypes[i][0].equals(paramName)) {
				return true;
			}
		}
		return false;
	}

	public double[] getArrayAverage(String expression, String className) {
		HashMap<Integer, double[]> xMap = null;
		try {
			xMap = evaluate(expression, "classID==\"" + className + "\"");
		} catch (CannotEvaluateException ex) {
			return null;
		}

		int n = 0;
		double[] average = null;
		for (Integer id : xMap.keySet()) {
			double[] x = xMap.get(id);
			if (x != null) {
				if (average == null) {
					average = x;
				} else {
					MathUtil.add(average, x);
				}
				n++;
			}
		}
		MathUtil.divide(average, n);

		return average;
	}

	public double[][] getArrayAverage1(String expression, String className) {
		HashMap<Integer, double[]> xMap = null;
		try {
			xMap = evaluate(expression, "classID==\"" + className + "\"");
		} catch (CannotEvaluateException ex) {
			return null;
		}

		ArrayMeanVarianceCalculator mvc = null;
		for (Integer id : xMap.keySet()) {
			double[] x = xMap.get(id);
			if (x != null) {
				if (mvc == null) {
					mvc = new ArrayMeanVarianceCalculator(x.length);
				}

				mvc.add(x);
			}
		}

		return new double[][] { mvc.getMean(), mvc.getMeanVariance() };
	}

	public double getAverage(String expression, String className) {
		return getAverage1(expression, className).x;
	}

	public Num getAverage1(String expression, String... className) {
		String cut = "";
		for (int i = 0; i < className.length; i++) {
			cut += "classID==\"" + className[i] + "\"";
			if (i != className.length - 1) {
				cut += "|";
			}
		}

		HashMap<Integer, Double> xMap = null;
		try {
			xMap = evaluate(expression, cut);
		} catch (CannotEvaluateException ex) {
			return new Num(Double.NaN, Double.NaN);
		}

		MeanVarianceCalculator mvc = new MeanVarianceCalculator();
		for (Integer id : xMap.keySet()) {
			Double x = xMap.get(id);
			if (x != null) {
				mvc.add(x);
			}
		}

		return new Num(mvc.getMean(), mvc.getMeanVariance());
	}

	synchronized public String[][] getColumnNamesAndTypes() {
		return columnNamesAndTypes;
	}

	synchronized private void openFile(String fileName) throws IOException {
		table = new ArrayList<ObjectEmbeddingTag[]>();
		file = new RandomAccessFile(fileName, "rw");
		nInput = new NeuroInputStream(new RandomAccessInputStream(file));
		nOutput = new NeuroOutputStream(new RandomAccessOutputStream(file));
	}

	synchronized public void addRow(Object[] row) {
		if (table.size() == maxRows) {
			throw new IllegalStateException("The maximum number of rows (" + maxRows + ") is exceded");
		}
		if (row.length != nColumns) {
			throw new IllegalArgumentException("Wrong row length. Received " + row.length + " expected " + nColumns);
		}

		ObjectEmbeddingTag[] tagRow = new ObjectEmbeddingTag[row.length];
		try {
			for (int i = 0; i < nColumns; i++) {
				String cName = "edu.ucsc.neurobiology.vision.io.tags." + columnNamesAndTypes[i][1] + "Tag";
				tagRow[i] = (ObjectEmbeddingTag) Class.forName(cName).newInstance();
				tagRow[i].setEmbeddedObject(row[i]);
			}
		} catch (ClassNotFoundException e) {
			Vision.reportFatalException("", e);
		} catch (InstantiationException e1) {
			Vision.reportFatalException("", e1);
		} catch (IllegalAccessException e2) {
			Vision.reportFatalException("", e2);
		}

		table.add(tagRow);
	}

	synchronized public Object[] getRow(int id) {
		int rowID = getRowForID(id);
		if (rowID == -1) {
			throw new IllegalArgumentException("No such neuron ID " + id);
		}

		ObjectEmbeddingTag[] row = (ObjectEmbeddingTag[]) table.get(rowID);
		Object[] obj = new Object[row.length];
		for (int i = 0; i < row.length; i++) {
			obj[i] = row[i].getEmbeddedObject();
		}
		return obj;
	}

	synchronized private int getRowForID(int id) {
		for (int i = 0; i < table.size(); i++) {
			DoubleTag t = (DoubleTag) table.get(i)[0];
			if (t.intValue() == id) {
				return i;
			}
		}

		return -1;
	}

	synchronized public int[] getIDList() {
		int[] idList = new int[table.size()];
		for (int i = 0; i < table.size(); i++) {
			DoubleTag idTag = (DoubleTag) table.get(i)[0];
			idList[i] = (int) idTag.doubleValue();
		}
		Arrays.sort(idList);

		return idList;
	}

	synchronized public final ObjectEmbeddingTag getCellByRow(int row, int column) {
		return ((ObjectEmbeddingTag[]) table.get(row))[column];
	}

	synchronized public final ObjectEmbeddingTag getCellByID(int id, int column) {
		return ((ObjectEmbeddingTag[]) table.get(getRowForID(id)))[column];
	}

	synchronized public final ObjectEmbeddingTag getCell(int id, String columnName) {
		if (!columnNamesMap.containsKey(columnName)) {
			return null;
		} else {
			int columnIndex = ((Integer) columnNamesMap.get(columnName)).intValue();
			return ((ObjectEmbeddingTag[]) table.get(getRowForID(id)))[columnIndex];
		}
	}

	synchronized public final void setCell(int id, String columnName, Object value) {
		int columnIndex = ((Integer) columnNamesMap.get(columnName)).intValue();
		ObjectEmbeddingTag tag = ((ObjectEmbeddingTag[]) table.get(getRowForID(id)))[columnIndex];
		tag.setEmbeddedObject(value);

		needsFlush = true;
	}

	synchronized public final double getDoubleCell(int id, String columnName) {
		ObjectEmbeddingTag t = getCell(id, columnName);
		return ((DoubleTag) t).doubleValue();
	}

	synchronized public final String getStringCell(int id, String columnName) {
		ObjectEmbeddingTag t = getCell(id, columnName);
		return (t != null) ? (String) ((StringTag) t).getEmbeddedObject() : null;
	}

	synchronized public final double[] getArrayCell(int id, String columnName) {
		ObjectEmbeddingTag t = getCell(id, columnName);
		return ((DoubleArrayTag) t).getArray();
	}

	synchronized public final double[] getDoubleArrayCell(int id, int column) {
		ObjectEmbeddingTag t = getCellByID(id, column);
		return ((DoubleArrayTag) t).getArray();
	}

	/**
	 * 
	 * 
	 * @param expressionString
	 * @param cutString
	 * @return
	 * @throws CannotEvaluateException
	 * 
	 *             Function assumes that pca has already been removed with
	 *             evaluateEx
	 */
	synchronized public HashMap evaluate(String expressionString, String cutString) throws CannotEvaluateException {

		HashMap<Integer, Object> columnMap = new HashMap<Integer, Object>();
		HashMap<String, Object> variables = new HashMap<String, Object>();
		Expression expression = new Expression(expressionString);
		if (cutString == null || cutString.trim().length() == 0) {
			cutString = "true";
		}
		Expression cutExpression = new Expression(cutString);

		for (int row = 0; row < table.size(); row++) {
			int id = ((DoubleTag) getCellByRow(row, 0)).intValue();

			variables.clear();
			for (int col = 0; col < nColumns; col++) {
				Object obj = getCellByRow(row, col).getEmbeddedObject();
				variables.put(columnNamesAndTypes[col][0], obj);
			}

			Object cut = null;
			try {
				cut = cutExpression.evaluate(variables);
			} catch (Exception ex) {
				throw new CannotEvaluateException("Cannot Evaluate Cut: " + ex.getMessage());
			}

			if (!(cut instanceof Boolean)) {
				throw new CannotEvaluateException("The Cut should be Boolean");
			}

			if (((Boolean) cut).booleanValue() == true) {
				Object result = null;
				try {
					result = expression.evaluate(variables);
				} catch (Exception ex1) {
					throw new CannotEvaluateException("Cannot Evaluate Expression: " + ex1.getMessage());
				}
				columnMap.put(new Integer(id), result);
			}
		}

		return columnMap;
	}

	synchronized public Object evaluate(String expressionString, int id) throws CannotEvaluateException {

		HashMap<String, Object> variables = new HashMap<String, Object>();
		Expression expression = new Expression(expressionString);
		int row = getRowForID(id);

		variables.clear();
		for (int col = 0; col < nColumns; col++) {
			Object obj = getCellByRow(row, col).getEmbeddedObject();
			variables.put(columnNamesAndTypes[col][0], obj);
		}

		try {
			return expression.evaluate(variables);
		} catch (Exception ex1) {
			throw new CannotEvaluateException("Cannot Evaluate Expression: " + ex1.getMessage());
		}
	}

	public int[] getNeuronsInClass(String className) {
		int[] id = getIDList();
		IntegerList idList = new IntegerList();
		for (int i = 0; i < id.length; i++) {
			if (getCellByID(id[i], 1).getEmbeddedObject().equals(className)) {
				idList.add(id[i]);
			}
		}
		return idList.toArray();
	}

	synchronized public LinkedHashMap<Integer, ? extends Object> getClassIDs() {
		LinkedHashMap<Integer, Object> columnMap = new LinkedHashMap<Integer, Object>();

		int columnIndex = ((Integer) columnNamesMap.get("classID")).intValue();
		// System.out.println(columnIndex);
		// if (columnIndex == -1) {
		// throw new Error("columnIndex is -1");
		// }

		for (int row = 0; row < table.size(); row++) {
			int id = ((DoubleTag) getCellByRow(row, 0)).intValue();
			Object obj = getCellByID(id, columnIndex).getEmbeddedObject();
			columnMap.put(new Integer(id), obj);
		}

		return columnMap;
	}

	public HashMap<Integer, Double> evaluateEx(String complexExpression, String cutExpression)
			throws CannotEvaluateException {

		if (complexExpression.startsWith("pca")) {
			if (complexExpression.charAt(3) != '(') {
				throw new CannotEvaluateException("Missing ( at character 3");
			}
			if (complexExpression.charAt(complexExpression.length() - 1) != ')') {
				throw new CannotEvaluateException("Missing ) at the end of expression");
			}
			int comma = complexExpression.lastIndexOf(',');
			if (comma == -1) {
				throw new CannotEvaluateException("pca() takes 2 parameters");
			}
			String expr = complexExpression.substring(4, comma);
			String eigenString = complexExpression.substring(comma + 1, complexExpression.length() - 1).trim();

			// System.err.println(expr);
			// System.err.println(eigenString);

			int eigen;
			try {
				eigen = Integer.parseInt(eigenString);
			} catch (NumberFormatException ex) {
				throw new CannotEvaluateException("the second parameter should be an integer: " + eigenString);
			}

			HashMap<Integer, double[]> vMap = evaluate(expr, cutExpression);
			if (vMap.isEmpty() || containsNulls(vMap)) {
				return null;
			} else {
				PCACalculator pca = new PCACalculator(vMap);
				return pca.getPCAComponent(eigen);
			}
		} else {
			return evaluate(complexExpression, cutExpression);
		}
	}

	synchronized public HashMap evaluate(String expressionString, int[] idList) throws CannotEvaluateException {
		if (expressionString.startsWith("pca")) {
			if (expressionString.charAt(3) != '(') {
				throw new CannotEvaluateException("Missing ( at character 3");
			}
			if (expressionString.charAt(expressionString.length() - 1) != ')') {
				throw new CannotEvaluateException("Missing ) at the end of expression");
			}
			int comma = expressionString.lastIndexOf(',');
			if (comma == -1) {
				throw new CannotEvaluateException("pca() takes 2 parameters");
			}
			String expr = expressionString.substring(4, comma);
			String eigenString = expressionString.substring(comma + 1, expressionString.length() - 1).trim();

			// System.err.println(expr);
			// System.err.println(eigenString);

			int eigen;
			try {
				eigen = Integer.parseInt(eigenString);
			} catch (NumberFormatException ex) {
				throw new CannotEvaluateException("the second parameter should be an integer: " + eigenString);
			}

			HashMap<Integer, double[]> vMap = evaluate(expr, idList);
			if (vMap.isEmpty() || containsNulls(vMap)) {
				return null;
			} else {
				PCACalculator pca = new PCACalculator(vMap);
				return pca.getPCAComponent(eigen);
			}
		} else {
			return evaluateSubfunc(expressionString, idList);
		}
	}

	/**
	 * 
	 * 
	 * @param expressionString
	 * @param idList
	 * @return
	 * @throws CannotEvaluateException
	 * 
	 *             Function assumes that pca has already been removed by
	 *             evaluate().
	 */
	synchronized public HashMap<Integer, Object> evaluateSubfunc(String expressionString, int[] idList)
			throws CannotEvaluateException {

		HashMap<Integer, Object> columnMap = new HashMap<Integer, Object>();
		HashMap<String, Object> variables = new HashMap<String, Object>();
		Expression expression = new Expression(expressionString);

		int maxID = MathUtil.max(getIDList());
		int[] id2rowMap = new int[maxID + 1];
		final int nRows = table.size();
		for (int row = 0; row < nRows; row++) {
			DoubleTag idTag = (DoubleTag) table.get(row)[0];
			id2rowMap[idTag.intValue()] = row;
		}

		for (int i = 0; i < idList.length; i++) {
			int id = idList[i];

			variables.clear();
			int row = id2rowMap[id];
			for (int column = 0; column < nColumns; column++) {
				Object obj = table.get(row)[column].getEmbeddedObject();
				variables.put(columnNamesAndTypes[column][0], obj);
			}

			Object result = null;
			try {
				result = expression.evaluate(variables);
			} catch (Exception ex) {
				throw new CannotEvaluateException("Cannot Evaluate Expression: " + ex.getMessage());
			}
			columnMap.put(new Integer(id), result);
		}

		return columnMap;
	}

	synchronized public void flush() throws IOException {
		file.seek(0);

		int nRows = table.size();

		// write the header
		file.writeInt(nColumns);
		file.writeInt(nRows);
		file.writeInt(maxRows);

		// write the column names
		for (int i = 0; i < nColumns; i++) {
			nOutput.writeString(columnNamesAndTypes[i][0]);
			nOutput.writeString(columnNamesAndTypes[i][1]);
		}
		final int headerSize = (int) file.getFilePointer();

		// write the cells
		file.seek(headerSize + nColumns * maxRows * 4);
		int[][] seekLocations = new int[nRows][nColumns];
		for (int i = 0; i < nRows; i++) {
			ObjectEmbeddingTag[] row = (ObjectEmbeddingTag[]) table.get(i);
			for (int j = 0; j < nColumns; j++) {
				seekLocations[i][j] = (int) file.getFilePointer();
				nOutput.writeTag(row[j]);
			}
		}

		// write the table of links
		file.seek(headerSize);
		for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nColumns; j++) {
				file.writeInt(seekLocations[i][j]);
			}
		}

		needsFlush = false;
	}

	synchronized public void printTable() {
		for (int i = 0; i < table.size(); i++) {
			for (int j = 0; j < nColumns; j++) {
				System.out.print(getCellByRow(i, j));
			}
			System.out.println();
		}
	}

	/**
	 * This should be designed so that it is safe to run it more than once.
	 * file.close() implements Closeable, which is specified to have no effect
	 * if run more than once.
	 * 
	 * @param save
	 * @throws IOException
	 */
	synchronized public void close(boolean save) throws IOException {
		if (save && needsFlush())
			flush();
		file.close();
	}

	synchronized public boolean needsFlush() {
		return needsFlush;
	}

	public static void profileEvaluate() throws IOException {
		final ParametersFile p = new ParametersFile("c:\\data003\\data003.params");
		final int[] idList = p.getIDList();

		Thread t = new Thread() {
			public void run() {
				double t1 = System.currentTimeMillis();
				for (int i = 0; i < 100; i++) {
					try {
						p.evaluate("norm(RedTimeCourse#GreenTimeCourse#BlueTimeCourse)", idList);
					} catch (CannotEvaluateException ex) {
						ex.printStackTrace();
					}
				}
				double t2 = System.currentTimeMillis();
				System.out.println((t2 - t1) / 1000.0);
			}
		};
		t.start();
	}

	public static boolean containsNulls(HashMap _m) {
		if (_m == null) {
			return true;
		} else {
			HashMap<Object, Object> m = _m;
			boolean allNulls = true;
			for (Object key : m.values()) {
				if (key != null) {
					allNulls = false;
				}
			}
			return allNulls;
		}
	}

	@Override
	/**
	 * This is a backup safeguard, mostly here because people seemed to forget to close STAFiles and 
	 * EI files appropriately from Matlab.  Do not rely on this to close files; you should always run close
	 * explicitly.
	 */
	protected void finalize() throws Throwable {
		try {
			close(false);
		} catch (Exception e) {
		} finally {
			super.finalize();
		}
	}

	// ------------------------------------------------------------------------------------------------------

	public double[][] gatherGaussians(int[] neurons) {

		double[][] gaussians = new double[neurons.length][];

		for (int i = 0; i < neurons.length; i++) {
			gaussians[i] = new double[] { getDoubleCell(neurons[i], "x0"), getDoubleCell(neurons[i], "y0"),
					getDoubleCell(neurons[i], "SigmaX"), getDoubleCell(neurons[i], "SigmaY"),
					getDoubleCell(neurons[i], "Theta") };
		}
		return gaussians;

	}
}
