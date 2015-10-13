package edu.ucsc.neurobiology.vision.io;

import java.io.*;
import java.util.*;


/**
 * Implements the 512 style header. Provides methods for reading, writing and accessing
 * header fields.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataHeader512 extends RawDataHeader {
    public static final boolean DEBUG = false;
    
    public static final int DEFAULT_TIME_BASE = 1904;
    
    public static final int HEADER_LENGTH_TAG = 0;
    public static final int TIME_TAG = 1;
    public static final int COMMENT_TAG = 2;
    public static final int FORMAT_TAG = 3;
    public static final int ARRAY_ID_TAG = 4;
    public static final int FREQUENCY_TAG = 5;
    public static final int TRIGGER_TAG = 6; //deprecated on 7/3/2015, mgrivich
    public static final int DATASET_IDENTIFIER_TAG = 7;
    public static final int TRIGGER_TAG_V2 = 8; //new on 7/3/2015
    public static final int DATA_TAG = 499;
    public static final int FILE_TYPE = 0x512;
    
    private int headerLength = -1;
    private int timeBase = -1;
    private long secondsTime = -1;
    private String comment = null;
    private String datasetIdentifier = "";
    private int format = -1;
    private int nElectrodes = -1;
    private int arrayID = -1;
    private int frequency = -1;
    private int nSamples = -1;
    
    
    public RawDataHeader512(long secondsTime, int nElectrodes, int frequency, int nSamples, int arrayID,
    		int format, String datasetIdentifier, String comment) throws IOException {
    	this(DEFAULT_TIME_BASE, secondsTime, nElectrodes, frequency, nSamples, arrayID, format, datasetIdentifier, comment);
    }
    
    public RawDataHeader512(int timeBase, long secondsTime, int nElectrodes, int frequency, int nSamples, 
    		int arrayID, int format, String datasetIdentifier, String comment) throws IOException {
    	this.timeBase = timeBase;
        this.secondsTime = secondsTime;
        this.comment = comment;
        this.format = format;
        this.nElectrodes = nElectrodes;
        this.arrayID = arrayID;
        this.frequency = frequency;
        this.nSamples = nSamples;
        this.datasetIdentifier = datasetIdentifier;

        this.headerLength = getBinaryRepresentation().length;
    }
    
    public RawDataHeader512(InputStream input) throws IOException {
        DataInputStream stream = new DataInputStream(input);

        boolean tagsLeft = true;
        while (tagsLeft) {
            final int tag = stream.readInt();
            final int length = stream.readInt();
            byte[] b;
            switch (tag) {
                case HEADER_LENGTH_TAG:
                    this.headerLength = stream.readInt();
                    if (DEBUG) System.out.println("headerLength: " + headerLength);
                    break;

                case TIME_TAG:
                    this.timeBase = stream.readInt();
                    this.secondsTime = stream.readLong();
                    if (DEBUG) {
                        System.out.println("timeBase: " + timeBase);
                        System.out.println("secondsTime: " + secondsTime);
                        System.out.println("printTime: " + printTime());
                    }
                    break;

                case COMMENT_TAG:
                    b = new byte[length];
                    for (int i = 0; i < length; i++) b[i] = stream.readByte();
                    this.comment = new String(b);
                    if (DEBUG) System.out.println("comment: " + comment);
                    break;

                case FORMAT_TAG:
                    this.format = stream.readInt();
                    if (DEBUG) System.out.println("format: " + format);
                    break;

                case ARRAY_ID_TAG:
                    this.nElectrodes = stream.readInt();
                    this.arrayID = stream.readInt();
                    if (DEBUG) System.out.println("nElectrodes: " + nElectrodes);
                    break;

                case FREQUENCY_TAG:
                    this.frequency = stream.readInt();
                    if (DEBUG) System.out.println("frequency: " + frequency);
                    break;

                case TRIGGER_TAG:
                    stream.readInt();
                    stream.readInt();
                    if (DEBUG) System.out.println("trigger: ");
                    break;
                    
                case TRIGGER_TAG_V2:
                	stream.readInt();
                	stream.readInt();
                	stream.readInt();
                	stream.readInt();
                	if (DEBUG) System.out.println("trigger_v2: ");
                	break;

                case DATASET_IDENTIFIER_TAG:
                    b = new byte[length];
                    for (int i = 0; i < length; i++) b[i] = stream.readByte();
                    this.datasetIdentifier = new String(b);
                    if (DEBUG) System.out.println("datasetIdentifier: " + datasetIdentifier);
                    break;

                case DATA_TAG:
                    this.nSamples = stream.readInt();
                    if (DEBUG) System.out.println("nSamples: " + nSamples);
                    tagsLeft = false; // DATA_TAG signals the end of the header
                    break;

                default: 
                    b = new byte[length];
                    for (int i = 0; i < length; i++) b[i] = stream.readByte();
                    this.comment = new String(b);
                    System.out.println("Warning: Unknown Tag:" + tag + " Length:" + length);
                    break;
            } // switch(tag)
        } // while(readingHeader)
        
        if ((headerLength == -1) || (timeBase == -1) || (secondsTime == -1) ||
            (comment == null) || (format == -1) || (nElectrodes == -1) ||
            (arrayID == -1) || (frequency == -1) || (nSamples == -1)) {
            throw new IOException("there are missing tags");
        }

    }

    
    public void setFormat(int format) {
        this.format = format;
    }

    public void setNumberOfSamples(int nSamples) {
        this.nSamples = nSamples;
    }

    public void setDatasetIdentifier(String datasetIdentifier) {
        this.datasetIdentifier = datasetIdentifier;
    }

    public int    getFormat()             { return format;            }
    public int    getNumberOfElectrodes() { return nElectrodes;       }
    public int    getArrayID()            { return arrayID;           }
    public int    getSamplingFrequency()  { return frequency;         }
    public int    getNumberOfSamples()    { return nSamples;          }
    public String getComment()            { return comment;           }
    public String getDatasetIdentifier()  { return datasetIdentifier; }
    public int    getTimeBase()           { return timeBase;          }
    public long   getTime()               { return secondsTime;       }
    public int getHeaderSize() { return headerLength; }

    
    public int getSampleSize() {
        return 2 + (nElectrodes - 1) * 3 / 2;
    }


    public byte[] getBinaryRepresentation() throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        DataOutputStream data = new DataOutputStream(bos);
        int tag;
        try {
            // read the header length tag
            data.writeInt(HEADER_LENGTH_TAG);
            data.writeInt(4);
            data.writeInt(headerLength);

            // read the time tag
            data.writeInt(TIME_TAG);
            data.writeInt(12);
            data.writeInt(timeBase);
            data.writeLong(secondsTime);

            // read the user comment tag
            data.writeInt(COMMENT_TAG);
            data.writeInt(comment.length());
            byte[] b = comment.getBytes();
            for (int i = 0; i < b.length; i++) {
                data.writeByte(b[i]);
            }

            // reading the format tag
            data.writeInt(FORMAT_TAG);
            data.writeInt(4);
            data.writeInt(format);

            // reading electrode ID tag
            data.writeInt(ARRAY_ID_TAG);
            data.writeInt(8);
            data.writeInt(nElectrodes);
            data.writeInt(arrayID);

            // reading the frequency tag
            data.writeInt(FREQUENCY_TAG);
            data.writeInt(4);
            data.writeInt(frequency);

            // reading the trigger params
            data.writeInt(TRIGGER_TAG);
            data.writeInt(8);
            data.writeInt(0);
            data.writeInt(0);

            // reading the dataset id tag
            data.writeInt(DATASET_IDENTIFIER_TAG);
            data.writeInt(datasetIdentifier.length());
            b = datasetIdentifier.getBytes();
            for (int i = 0; i < b.length; i++)
                data.writeByte(b[i]);

            // writing the data tag
            data.writeInt(DATA_TAG);
            data.writeInt(4);
            data.writeInt(nSamples);

            data.close();
        } catch (IOException e) {
            throw new IOException("Could not generate binary header.");
        }

        return bos.toByteArray();
    }

    
    /**
     * Contract of this method is that it processes the timeBase and secondsTime and 
     * returns an object that represents the proper time in a Java-friendly format
     * that will be printed correctly with toString.
     */
    public Date printTime() {
    	Calendar c = Calendar.getInstance(TimeZone.getTimeZone("GMT"));
    	c.clear();
    	c.set(Calendar.YEAR, timeBase);
    	long baseMillis = c.getTimeInMillis();
    	return new Date(baseMillis + secondsTime*1000);
    }
    
    
    @Override
    public String toString() {
        String s = "";        
        s += "\n timeBase = " + timeBase;
        s += "\n secondsTime = " + secondsTime;
        s += "\n printTime = " + printTime();
        
        s += "\n dataset ID = " + datasetIdentifier;
        s += "\n format = " + format;
        s += "\n nSamples = " + nSamples;
        s += "\n nElectrodes = " + nElectrodes;
        s += "\n arrayID = " + extractArrayID(arrayID) +
            " (" + extractElectrodeMapPart(arrayID) + "/" +
            extractElectrodeMapPartsCount(arrayID) + ")";
        s += "\n scanRate = " + frequency;
        s += "\n headerSize = " + headerLength;
        s += "\n comment = " + comment;

//        System.out.println(getExperimentIdentifier());
//        System.out.println(getDatasetName());

        return s;
    }

}