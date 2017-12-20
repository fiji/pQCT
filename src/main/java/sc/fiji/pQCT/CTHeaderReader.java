
package sc.fiji.pQCT;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import javax.activation.UnsupportedDataTypeException;

class CTHeaderReader {

	private static final int PERCENT_OFFSET = 334;
	private static final int SITE_LENGTH_OFFSET = 1316;
	private final File fileIn;
	double percent;
	String site;

	CTHeaderReader(final String fileName) {
		fileIn = new File(fileName);
	}

	private void readHeader(final DataInputStream input) throws IOException {
		input.skipBytes(PERCENT_OFFSET);
		percent = input.readDouble();
		input.skipBytes(SITE_LENGTH_OFFSET - PERCENT_OFFSET);
		final int siteLength = input.readByte();
		final byte[] siteBytes = new byte[siteLength];
		input.read(siteBytes, 0, siteLength);
		site = new String(siteBytes);
	}

	// Read the file to memory
	void read() throws Exception {
		try (BufferedInputStream inputStream = new BufferedInputStream(
			new FileInputStream(fileIn)))
		{
			final DataInputStream dataInputStream = new DataInputStream(inputStream);
			readHeader(dataInputStream);
		}
		catch (final Exception e) {
			throw new UnsupportedDataTypeException(
				"Could not read input file or file did not exist.");
		}
	}
}
