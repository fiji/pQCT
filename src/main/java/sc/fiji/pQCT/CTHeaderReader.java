/*
BSD 2-Clause License

Copyright (c) 2018, Timo Rantalainen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
