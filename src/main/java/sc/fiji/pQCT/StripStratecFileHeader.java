/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

	ImageJ Stratec file header stripper. Files stripped with this plug-in still contain a header,
	just not data identifying the patient without having some additional knowledge.
    Copyright (C) 2012 Timo Rantalainen
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
		/*MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment*/
		dialog.addCheckbox("Strip_MeasInfo", false);
		dialog.addCheckbox("Strip_PatBirth", false);
		dialog.addCheckbox("Strip_PatMenoAge", false);
		dialog.addCheckbox("Strip_PatName", true);
		dialog.addCheckbox("Strip_PatTitleAndComment", false);
		dialog.addStringField("Stratec_file_to_strip", Prefs.getDefaultDirectory() +
			"I0020001.M01", 60);
		dialog.addStringField("File_save_name", Prefs.getDefaultDirectory() +
			"I0020001.M01", 60);
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

	// Writing dummy header containing sufficient details for
	// Distribution_Analysis imageJ plugin, might not suffice for Geanie or
	// Stractec software
	private static void stripHeader(final byte[] data,
		final boolean[] toStrip)
	{
		/*MeasInfo,PatBirth,PatMenoAge,PatName,PatTitle&Comment*/
		final int[] offsetsToStrip = { 662, 1091, 1095, 1099, 1141 };
		final int[] stripLengths = { 324, 4, 4, 41, 124 };
		for (int s = 0; s < offsetsToStrip.length; ++s) {
			if (toStrip[s]) {
				Arrays.fill(data, offsetsToStrip[s], stripLengths[s], (byte) 0);
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
