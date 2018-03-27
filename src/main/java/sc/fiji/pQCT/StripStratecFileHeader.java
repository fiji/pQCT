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
import java.io.FileOutputStream;
import java.util.Arrays;

import javax.activation.UnsupportedDataTypeException;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class StripStratecFileHeader implements PlugIn {

	@Override
	public void run(final String arg) {
		final GenericDialog dialog = new GenericDialog("Strip Stratec Header");
		// MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment
		dialog.addCheckbox("Strip_MeasInfo", false);
		dialog.addCheckbox("Strip_PatBirth", false);
		dialog.addCheckbox("Strip_PatMenoAge", false);
		dialog.addCheckbox("Strip_PatName", true);
		dialog.addCheckbox("Strip_PatTitleAndComment", false);
		dialog.addStringField("Stratec_file_to_strip", Prefs.getDefaultDirectory() +
			"I0020001.m01", 60);
		dialog.addStringField("File_save_name", Prefs.getDefaultDirectory() +
			"I0020001.m01", 60);
		dialog.showDialog();
		if (!dialog.wasOKed()) {
			return;
		}

		final boolean[] toStrip = new boolean[5];
		for (int i = 0; i < toStrip.length; ++i) {
			toStrip[i] = dialog.getNextBoolean();
		}
		final String fileIn = dialog.getNextString();
		final String fileOut = dialog.getNextString();
		if (fileIn == null || fileOut == null) {
			IJ.log("Give both input and output file as parameters");
			return;
		}
		final File test = new File(fileIn);
		if (!test.exists()) {
			IJ.error("Input file didn't exist");
			return;
		}

		try {
			stripFile(fileIn, fileOut, toStrip);
		}
		catch (final Exception err) {
			IJ.error("Stratec file header stripping failed", err.getMessage());
		}
	}

	private static void stripFile(final String fileNameIn,
		final String fileNameOut, final boolean[] toStrip) throws Exception
	{
		final File fileIn = new File(fileNameIn);
		final long fileLength = fileIn.length();
		final byte[] fileData;
		try (BufferedInputStream inputStream = new BufferedInputStream(
			new FileInputStream(fileIn)))
		{
			final DataInputStream dataInputStream = new DataInputStream(inputStream);
			fileData = new byte[(int) fileLength];
			dataInputStream.read(fileData, 0, (int) fileLength);
		}
		catch (final Exception e) {
			throw new UnsupportedDataTypeException("Could not read input file.");
		}
		stripHeader(fileData, toStrip); // Strip header
		writeFile(fileNameOut, fileData);
	}

	// Writing dummy header containing sufficient details for the
	// PqctAnalysis imageJ plugin, might not suffice for Geanie or
	// Stratec software
	private static void stripHeader(final byte[] data,
		final boolean[] toStrip)
	{
		// MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment
		final int[] fromIndices = { 662, 1091, 1095, 1099, 1141 };
		final int[] stripLengths = { 324, 4, 4, 41, 124 };
		for (int s = 0; s < fromIndices.length; ++s) {
			if (toStrip[s]) {
				final int toIndex = fromIndices[s] + stripLengths[s];
				Arrays.fill(data, fromIndices[s], toIndex, (byte) 0);
			}
		}
	}

	private static void writeFile(final String fileName, final byte[] fileData) {
		try (FileOutputStream writer = new FileOutputStream(fileName)) {
			writer.write(fileData);
		}
		catch (final Exception err) {
			IJ.error("Saving failed");
		}
	}
}
