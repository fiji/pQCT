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

import java.util.StringTokenizer;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.ImageInfo;
import ij.plugin.PlugIn;
import ij.text.TextPanel;

public class ExportHeader implements PlugIn {

	private String imageInfo;
	private String site;
	private double percent;

	@Override
	public void run(final String arg) {
		final ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) return;
		if (imp.getType() != ImagePlus.GRAY16) {
			IJ.error("Distribution analysis expects 16-bit greyscale data");
			return;
		}
		imageInfo = new ImageInfo().getImageInfo(imp);
		final String imageName;
		if (getInfoProperty(imageInfo, "File Name") != null) {
			imageName = getInfoProperty(imageInfo, "File Name");
		}
		else {
			if (imp.getImageStackSize() == 1) {
				imageName = imp.getTitle();
			}
			else {
				imageName = imageInfo.substring(0, imageInfo.indexOf("\n"));
			}
			imageInfo += "File Name:" + imageName + "\n";
		}

		// Try to get CT file additional header info
		String filePath = "";
		if (getInfoProperty(imageInfo, "File Path") != null) {
			filePath = getInfoProperty(imageInfo, "File Path");
		}
		final String nameSub = imageName != null ? imageName.substring(1, imageName
			.length()) : "";
		final String cName = "C" + nameSub;
		try {
			final CTHeaderReader ct = new CTHeaderReader(filePath + cName);
			ct.read();
			site = ct.site;
			percent = ct.percent;
		}
		catch (final Exception err) {
			IJ.log("Failed to read CT header");
			site = "Unknown";
			percent = Double.NaN;
		}

		TextPanel textPanel = IJ.getTextPanel();
		if (textPanel == null) {
			textPanel = new TextPanel();
		}
		if (textPanel.getLineCount() == 0) {
			writeHeader(textPanel);
		}

		String results = "";
		results = printResults(results, imp);
		textPanel.appendLine(results);
		textPanel.updateDisplay();
	}

	private static String getInfoProperty(final String properties,
		final CharSequence propertyToGet)
	{
		final StringTokenizer st = new StringTokenizer(properties, "\n");
		String currentToken;
		while (st.hasMoreTokens()) {
			currentToken = st.nextToken();
			if (currentToken.contains(propertyToGet)) {
				final String[] returnValue = currentToken.split(":", 2);
				return returnValue[1].trim();
			}
		}
		return null;
	}

	private String printResults(String results, final ImagePlus imp) {
		if (imp == null) {
			return results;
		}
		final String[] propertyNames = { "File Name", "File Path", "Patient's Name",
			"Patient ID", "Patient's Birth Date", "Acquisition Date", "Pixel Spacing",
			"ObjLen" };
		final String fileNameProperty = getInfoProperty(imageInfo, "File Name");
		if (fileNameProperty != null) {
			results += fileNameProperty + "\t";
		}
		else if (imp.getImageStackSize() == 1) {
			results += getInfoProperty(imageInfo, "Title") + "\t";
		}
		else {
			results += imageInfo.substring(0, imageInfo.indexOf("\n")) + "\t";
		}
		final StringBuilder resultsBuilder = new StringBuilder(results);
		for (int i = 1; i < propertyNames.length; ++i) {
			resultsBuilder.append(getInfoProperty(imageInfo, propertyNames[i]))
				.append("\t");
		}
		resultsBuilder.append(site).append("\t");
		resultsBuilder.append(percent).append("\t");
		return resultsBuilder.toString();
	}

	private static void writeHeader(final TextPanel textPanel) {
		final String[] propertyNames = { "File Name", "File Path", "Patient's Name",
			"Patient ID", "Patient's Birth Date", "Acquisition Date", "Pixel Spacing",
			"Object Length", "Measurement Mask", "Percent Length" };
		final StringBuilder headings = new StringBuilder();
		for (final String propertyName : propertyNames) {
			headings.append(propertyName).append("\t");
		}
		textPanel.setColumnHeadings(headings.toString());
	}
}
