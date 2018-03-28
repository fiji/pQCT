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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import javax.activation.UnsupportedDataTypeException;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

// This file format is supported in SCIFIO already,
// but we'll keep this plugin around for people who don't have it enabled.
// Also SCIFIO is as of yet experimental code.
// TODO Remove when ImageJ/Fiji comes with SCIFIO enabled by default
public class ReadStratecFile extends ImagePlus implements PlugIn {

	private static final int HEADER_LENGTH = 1609;
	private String PatName;
	private long PatNo;
	private int PatMeasNo;
	private long PatBirth;
	private long MeasDate;
	private double VoxelSize;
	private int PicX0;
	private int PicY0;
	private int PicMatrixX;
	private int PicMatrixY;
	private String MeasInfo;
	private String Device;
	private String PatID;
	private double ObjLen;
	private String fileName;
	private String properties;

	@Override
	public void run(final String arg) {
		final String path;
		if (!arg.isEmpty()) {
			// Called by HandleExtraFileTypes
			final File file = new File(arg);
			path = file.getParent() + "/";
			fileName = file.getName();
		}
		else {
			// Select file manually
			final OpenDialog od = new OpenDialog("Select stratec image (I*.M*)", arg);
			path = od.getDirectory();
			fileName = od.getFileName();
		}
		if (fileName == null) return;
		readFile(path);
		fileInfo();
		if (arg.isEmpty() && getHeight() > 0) {
			show();
		}
	}

	private void fileInfo() {
		FileInfo fi = getFileInfo();
		if (fi == null) {
			fi = new FileInfo();
		}
		fi.pixelWidth = VoxelSize;
		fi.pixelHeight = VoxelSize;
		fi.width = PicMatrixX;
		fi.height = PicMatrixY;
		fi.valueUnit = "mm";
		fi.fileName = fileName;
		fi.info = properties;
		fi.fileFormat = FileInfo.RAW;
		fi.compression = FileInfo.COMPRESSION_NONE;
		fi.fileType = FileInfo.GRAY16_SIGNED;
		setFileInfo(fi);
	}

	private static String getNByteString(final ByteBuffer buffer, final int pos) {
		buffer.position(pos);
		final byte n = buffer.get();
		final byte[] bytes = new byte[n];
		buffer.get(bytes);
		return new String(bytes);
	}

	private void readFile(final String path) {
		final File file = new File(path + fileName);
		final int bytes = (int) file.length();
		if (bytes < HEADER_LENGTH) {
			IJ.error("Reading the Stratec file failed: File length < 1609 bytes.");
			return;
		}
		try (final DataInputStream dataInputStream = new DataInputStream(
			new BufferedInputStream(new FileInputStream(file))))
		{
			// Allocate memory for reading the file into memory
			final byte[] data = new byte[bytes];
			// Read the data to memory
			dataInputStream.read(data, 0, bytes);
			final ByteBuffer buffer = ByteBuffer.wrap(data).order(
				ByteOrder.LITTLE_ENDIAN);
			readHeader(buffer);
			readImage(buffer, path);
		}
		catch (final IOException e) {
			IJ.error("Reading the Stratec file failed: " + e.getMessage());
		}
	}

	private void readHeader(final ByteBuffer buffer)
		throws UnsupportedDataTypeException
	{
		Device = getNByteString(buffer, 1050);
		if (!Device.toLowerCase().contains(".typ")) {
			throw new UnsupportedDataTypeException("Device string not found.");
		}
		VoxelSize = buffer.getDouble(12);
		ObjLen = buffer.getDouble(318);
		MeasInfo = getNByteString(buffer, 662);
		MeasDate = buffer.getInt(986);
		PatMeasNo = buffer.getShort(1085);
		PatNo = buffer.getInt(1087);
		PatBirth = buffer.getInt(1091);
		PatName = getNByteString(buffer, 1099);
		PatID = getNByteString(buffer, 1282);
		PicX0 = buffer.getShort(1525);
		PicY0 = buffer.getShort(1527);
		PicMatrixX = buffer.getShort(1529);
		PicMatrixY = buffer.getShort(1531);
	}

	private void readImage(final ByteBuffer buffer, final String path) {
		final ImagePlus tempImage = NewImage.createShortImage(fileName + " " +
			Double.toString(VoxelSize), PicMatrixX, PicMatrixY, 1,
			NewImage.FILL_BLACK);
		setImage(tempImage.getImage());
		setProcessor(fileName, tempImage.getProcessor());
		setProperties(path);
		final short[] pixels = (short[]) getProcessor().getPixels();
		final int size = PicMatrixX * PicMatrixY;
		int min = Short.MAX_VALUE;
		int max = Short.MIN_VALUE;
		buffer.position(HEADER_LENGTH);
		for (int i = 0; i < size; i++) {
			final int pixel = readSignedShort(buffer);
			final int unsignedShort = pixel & 0xFFFF;
			min = Math.min(min, unsignedShort);
			max = Math.max(max, unsignedShort);
			pixels[i] = (short) pixel;
		}
		setDisplayRange(min, max);
		final Calibration cal = getCalibration();
		final double[] coefficients = { -32.768, 0.001 };
		cal.setFunction(Calibration.STRAIGHT_LINE, coefficients, "1/cm");
		cal.setUnit("mm");
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = VoxelSize;
	}

	private static int readSignedShort(final ByteBuffer buffer) {
		final int bitMask = 0x8000;
		final short pixel = buffer.getShort();
		return (pixel >= 0 ? -bitMask : bitMask - 1) + pixel;
	}

	private void setProperties(final String directory) {
		final String[] propertyNames = { "File Name", "File Path", "Pixel Spacing",
			"ObjLen", "MeasInfo", "Acquisition Date", "Device", "PatMeasNo", "PatNo",
			"Patient's Birth Date", "Patient's Name", "Patient ID", "PicX0", "PicY0",
			"Width", "Height", "Stratec File" };
		final String[] propertyValues = { fileName, directory, Double.toString(
			VoxelSize), Double.toString(ObjLen), MeasInfo, Long.toString(MeasDate),
			Device, Integer.toString(PatMeasNo), Long.toString(PatNo), Long.toString(
				PatBirth), PatName, PatID, Integer.toString(PicX0), Integer.toString(
					PicY0), Integer.toString(PicMatrixX), Integer.toString(PicMatrixY),
			"1" };
		final StringBuilder builder = new StringBuilder();
		for (int i = 0; i < propertyNames.length; ++i) {
			builder.append(propertyNames[i]).append(": ").append(propertyValues[i])
				.append("\n");
		}
		properties = builder.toString();
		setProperty("Info", properties);
	}
}
