package sc.fiji.pQCT;

import java.io.*;
import javax.activation.*; //UnsupportedDataTypeException

public class CTHeaderReader{
	File fileIn;
	public double percent;
	public String site;
	CTHeaderReader(String fileName){
		fileIn = new File(fileName);
	}
	
	//Read the file to memory
	public void read() throws Exception{
		long fileLength = fileIn.length();
		byte[] fileData;
		try {
			BufferedInputStream inputStream = new BufferedInputStream(
					new FileInputStream(fileIn));
			DataInputStream dataInputStream = new DataInputStream(inputStream);
			// Allocate memory for reading the file into memory
			fileData = new byte[(int) fileLength];
			// Read the data to memory
			dataInputStream.read(fileData, 0, (int) fileLength);
			// Close the file after reading
			dataInputStream.close();
			readHeader(fileData);
		} catch (Exception e) {
			throw new UnsupportedDataTypeException("Could not read input file or file did not exist.");
		}
	}
	
	//Extract the header
	private void readHeader(byte[] fileData){
		int offset = 334;
		percent = Double
				.longBitsToDouble((long) (((long) (fileData[offset + 7] & 0xFF)) << 56
						| ((long) (fileData[offset + 6] & 0xFF)) << 48
						| ((long) (fileData[offset + 5] & 0xFF)) << 40
						| ((long) (fileData[offset + 4] & 0xFF)) << 32
						| ((long) (fileData[offset + 3] & 0xFF)) << 24
						| ((long) (fileData[offset + 2] & 0xFF)) << 16
						| ((long) (fileData[offset + 1] & 0xFF)) << 8 | ((long) (fileData[offset + 0] & 0xFF)) << 0));
		offset = 1316;
		site = new String(fileData, offset + 1, (int) fileData[offset]);
	}
}