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

	ImageJ density distribution analysis plugin
    Copyright (C) 2011 Timo Rantalainen
*/

package sc.fiji.pQCT;

import java.util.StringTokenizer;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.plugin.filter.Info;
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
		// TODO Fix deprecated call
		imageInfo = new Info().getImageInfo(imp, imp.getChannelProcessor());
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
			IJ.log("CT header " + err.toString());
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
