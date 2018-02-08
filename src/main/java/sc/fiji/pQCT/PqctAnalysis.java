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

import java.awt.Color;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.stream.DoubleStream;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.filter.Info;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.text.TextPanel;
import sc.fiji.pQCT.analysis.ConcentricRingAnalysis;
import sc.fiji.pQCT.analysis.CorticalAnalysis;
import sc.fiji.pQCT.analysis.DetermineAlpha;
import sc.fiji.pQCT.analysis.DistributionAnalysis;
import sc.fiji.pQCT.analysis.MassDistribution;
import sc.fiji.pQCT.analysis.SoftTissueAnalysis;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.io.ScaledImageData;
import sc.fiji.pQCT.selectroi.RoiSelector;
import sc.fiji.pQCT.selectroi.SelectROI;
import sc.fiji.pQCT.selectroi.SelectSoftROI;
import sc.fiji.pQCT.utils.ResultsWriter;

public class PqctAnalysis implements PlugIn {

	@Override
	public void run(final String arg) {
		IJ.log("Launching plugin");
		final ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null) return;
		if (imp.getType() != ImagePlus.GRAY16) {
			IJ.error("Distribution analysis expects 16-bit greyscale data");
			return;
		}
		// Set sector widths and division numbers
		// Distribution analysis sectorWidth, Distribution analysis sectors,
		// Concentric distribution analysis sectorWidth, Concentric distribution
		// analysis sectors
		final int[] sectorsAndDivisions = { 10, 3, 10, 10 };

		// TODO replace calls to deprecated code
		String imageInfo = new Info().getImageInfo(imp, imp.getChannelProcessor());
		// Check image calibration
		final Calibration cal = imp.getCalibration();
		double[] calibrationCoefficients = { 0, 1 };
		if (getInfoProperty(imageInfo, "Stratec File") == null) {
			if (cal != null && cal.getCoefficients() != null) {
				calibrationCoefficients = cal.getCoefficients();
			}
		}
		else {
			calibrationCoefficients = new double[2];
			/*Read calibration from TYP file database*/
			final String typFileName = getInfoProperty(imageInfo, "Device");
			try {
				final ClassLoader loader = getClass().getClassLoader();
				final InputStream ir = loader.getResourceAsStream("typ/" + typFileName);
				final byte[] typFileData = new byte[ir.available()];
				ir.read(typFileData);
				ir.close();
				final String typFiledDataString = new String(typFileData, "ISO-8859-1");
				// break the typFileDataString into lines
				final StringTokenizer st = new StringTokenizer(typFiledDataString,
					"\n");
				final List<String> typFileLines = new Vector<>();
				while (st.hasMoreTokens()) {
					typFileLines.add(st.nextToken());
				}
				// Search for XSlope and XInter
				final String[] searchFor = { "XInter", "XSlope" };
				for (int i = 0; i < searchFor.length; ++i) {
					int index = 0;
					String temp = typFileLines.get(index);
					while (!temp.contains(searchFor[i]) && index < typFileLines.size()) {
						++index;
						temp = typFileLines.get(index);
					}
					if (temp.contains(searchFor[i])) { // Found line
						final StringTokenizer st2 = new StringTokenizer(temp, "=");
						final List<String> typFileLineTokens = new Vector<>();
						while (st2.hasMoreTokens()) {
							typFileLineTokens.add(st2.nextToken().trim());
						}
						calibrationCoefficients[i] = Double.valueOf(typFileLineTokens.get(
							1));
					}
					else {
						calibrationCoefficients[i] = i * 1000.0;
					}
				}
				calibrationCoefficients[1] /= 1000.0; // 1.495
			}
			catch (final NullPointerException npe) {
				IJ.error(".TYP file not found");
			}
			catch (final IOException e) {
				IJ.error(".TYP file could not be read");
			}
		}
		double resolution = cal.pixelWidth;
		if (getInfoProperty(imageInfo, "Pixel Spacing") != null) {
			String temp = getInfoProperty(imageInfo, "Pixel Spacing");
			if (temp.contains("\\")) {
				temp = temp.substring(0, temp.indexOf("\\"));
			}
			resolution = Double.valueOf(temp);
		}
		// Get parameters for scaling the image and for thresholding
		final GenericDialog dialog = new GenericDialog("Analysis parameters");
		final String[] topLabels = { "Flip_horizontal", "Flip_vertical",
			"No_filtering", "Measurement_tube" };
		final boolean[] defaultTopValues = new boolean[4];
		dialog.addCheckboxGroup(1, 4, topLabels, defaultTopValues);
		dialog.addNumericField("Air_threshold", -40, 4, 8, null);
		dialog.addNumericField("Fat threshold", 40, 4, 8, null);
		dialog.addNumericField("Muscle_threshold", 40, 4, 8, null);
		dialog.addNumericField("Edge_divisions", 45, 4, 8, null);
		dialog.addNumericField("Marrow_threshold", 80, 4, 8, null);
		dialog.addNumericField("Soft_tissue_threshold", 200.0, 4, 8, null);
		dialog.addNumericField("Rotation_threshold", 200.0, 4, 8, null);
		dialog.addNumericField("Area threshold", 550.0, 4, 8, null); // 550.0
		dialog.addNumericField("bMD threshold", 690.0, 4, 8, null); // 690.0
		dialog.addNumericField("Scaling_coefficient (slope)",
			calibrationCoefficients[1], 4, 8, null);
		dialog.addNumericField("Scaling_constant (intercept)",
			calibrationCoefficients[0], 4, 8, null);

		// Get ROI selection
		final String[] choiceLabels = { "Bigger", "Smaller", "Left", "Right", "Top",
			"Bottom", "Central", "Peripheral", "SecondLargest", "TwoLargestLeft",
			"TwoLargestRight", "FirstFromLeft", "SecondFromLeft", "ThirdFromLeft",
			"FourthFromLeft", "FifthFromLeft", "FirstFromTop", "SecondFromTop",
			"ThirdFromTop", "FourthFromTop", "FifthFromTop" };
		dialog.addChoice("Roi_selection", choiceLabels, choiceLabels[0]);
		dialog.addChoice("Soft_Tissue_Roi_selection", choiceLabels,
			choiceLabels[0]);
		final String[] rotationLabels = { "According_to_Imax/Imin",
			"Furthest_point", "All_Bones_Imax/Imin", "Not_selected_to_right",
			"Selected_to_right" };
		dialog.addChoice("Rotation_selection", rotationLabels, rotationLabels[0]); // "According_to_Imax/Imin"

		final String[] middleLabels = { "Analyse_cortical_results",
			"Analyse_mass_distribution", "Analyse_concentric_density_distribution",
			"Analyse_density_distribution", "Analyse_soft_tissues",
			"Prevent_peeling_PVE_pixels", "Allow_cleaving", "Suppress_result_image",
			"Limit_ROI_search_to_manually_selected",
			"Set_distribution_results_rotation_manually" };
		final boolean[] middleDefaults = new boolean[10];
		middleDefaults[4] = true;
		dialog.addCheckboxGroup(4, 3, middleLabels, middleDefaults);

		dialog.addNumericField("Manual_rotation_[+-_180_deg]", 0.0, 4, 8, null);

		final String[] bottomLabels = new String[8];
		final boolean[] bottomDefaults = new boolean[8];
		bottomLabels[0] = "Guess_flip";
		bottomLabels[1] = "Guess_right";
		bottomLabels[2] = "Guess_larger";
		bottomLabels[3] = "Stacked_bones";
		bottomLabels[4] = "Guess_stacked";
		bottomLabels[5] = "Invert_flip_guess";
		bottomLabels[6] = "Flip_distribution_results";
		bottomLabels[7] = "Save_visual_result_image_on_disk";
		dialog.addCheckboxGroup(2, 5, bottomLabels, bottomDefaults);

		dialog.addStringField("Image_save_path", Prefs.getDefaultDirectory(), 40);
		// TODO Change help URL
		dialog.addHelp("http://bonej.org/densitydistribution");
		dialog.showDialog();
		if (!dialog.wasOKed()) {
			return;
		}
		for (int i = 0; i < defaultTopValues.length; ++i) {
			defaultTopValues[i] = dialog.getNextBoolean();
		}
		final double[] thresholdsAndScaling = new double[11];
		for (int i = 0; i < thresholdsAndScaling.length; ++i) {
			thresholdsAndScaling[i] = dialog.getNextNumber();
		}
		final String[] alignmentStrings = new String[3];
		for (int i = 0; i < alignmentStrings.length; ++i) {
			alignmentStrings[i] = dialog.getNextChoice();
		}
		for (int i = 0; i < middleDefaults.length; ++i) {
			middleDefaults[i] = dialog.getNextBoolean();
		}
		final double manualAlpha = dialog.getNextNumber() * Math.PI / 180.0;
		for (int i = 0; i < bottomDefaults.length; ++i) {
			bottomDefaults[i] = dialog.getNextBoolean();
		}
		final String imageSavePath = dialog.getNextString();
		final ScaledImageData scaledImageData;

		String imageName = getInfoProperty(imageInfo, "File Name");
		if (imageName == null) {
			if (imp.getImageStackSize() == 1) {
				imageName = imp.getTitle();
			}
			else {
				imageName = imageInfo.substring(0, imageInfo.indexOf("\n"));
			}
			imageInfo += "File Name:" + imageName + "\n";
		}

		final short[] tempPointer = (short[]) imp.getProcessor().getPixels();
		final int[] signedShort = new int[tempPointer.length];
		final float[] floatPointer = (float[]) imp.getProcessor().toFloat(1, null)
			.getPixels();
		if (imp.getOriginalFileInfo().fileType == ij.io.FileInfo.GRAY16_SIGNED ||
			cal.isSigned16Bit())
		{
			for (int i = 0; i < tempPointer.length; ++i) {
				signedShort[i] = (int) (floatPointer[i] - Math.pow(2.0, 15.0));
			}
		}
		else {
			/*
			Apply the original calibration of the image prior to applying the calibration got from the user
			-> enables using ImageJ for figuring out the calibration without too much fuss.
			*/
			try {
				double[] origCalCoeffs = imp.getOriginalFileInfo().coefficients;
				if (origCalCoeffs == null) {
					origCalCoeffs = cal.getCoefficients();
				}
				for (int i = 0; i < tempPointer.length; ++i) {
					signedShort[i] = (int) (floatPointer[i] * origCalCoeffs[1] +
						origCalCoeffs[0]);
				}
			}
			catch (final Exception err) {
				for (int i = 0; i < tempPointer.length; ++i) {
					signedShort[i] = tempPointer[i];
				}
			}
		}
		final ImageAndAnalysisDetails details = new ImageAndAnalysisDetails(
			defaultTopValues, thresholdsAndScaling, alignmentStrings, choiceLabels,
			rotationLabels, middleDefaults, manualAlpha, bottomDefaults,
			sectorsAndDivisions);
		// Scale and 3x3 median filter the data
		scaledImageData = new ScaledImageData(signedShort, imp.getWidth(), imp
			.getHeight(), resolution, details.scalingFactor, details.constant,
			details.flipHorizontal, details.flipVertical, details.noFiltering);
		RoiSelector roi = null;
		RoiSelector softRoi = null;

		try {
			if (details.cOn || details.mOn || details.conOn || details.dOn) {
				roi = new SelectROI(scaledImageData, details, imp,
					details.boneThreshold, true);
			}
			if (details.stOn) {
				// An ROI appears on the image every now and then, haven't figured out
				// why -> remove any unwanted rois prior to soft-tissue analysis
				imp.setRoi(null, false);
				softRoi = new SelectSoftROI(scaledImageData, details, imp);
				if (roi == null) {
					roi = softRoi;
				}
			}
		}
		catch (final ExecutionException err) {
			IJ.log("Caught sieve error " + err.toString());
			return;
		}

		if (roi == null) {
			IJ.log("No analysis was selected.");
			return;
		}
		boolean alphaOn = false;
		DetermineAlpha determineAlpha = null;
		if (details.cOn || details.mOn || details.conOn || details.dOn) {
			determineAlpha = new DetermineAlpha((SelectROI) roi, details);
			alphaOn = true;
		}

		details.flipDistribution = roi.details.flipDistribution;
		details.stacked = roi.details.stacked;

		TextPanel textPanel = IJ.getTextPanel();
		if (textPanel == null) {
			textPanel = new TextPanel();
		}
		final ResultsWriter resultsWriter = new ResultsWriter(imageInfo, alphaOn);

		if (textPanel.getLineCount() == 0) {
			resultsWriter.writeHeader(textPanel, details);
		}

		String results = "";
		results = resultsWriter.printResults(results, details, imp);
		if (determineAlpha != null) {
			results = printAlpha(results, determineAlpha);
		}

		ImagePlus resultImage = null;
		boolean makeImage = true;
		if (details.suppressImages && !details.saveImageOnDisk && roi != null) {
			makeImage = false;
		}
		else {
			resultImage = getRGBResultImage(roi.scaledImage, roi.width, roi.height,
				imageSavePath);
			resultImage.setTitle(imp.getTitle() + "-result");
		}
		if (details.stOn) {
			final SoftTissueAnalysis softTissueAnalysis = new SoftTissueAnalysis(
				(SelectSoftROI) softRoi);
			results = printSoftTissueResults(results, softTissueAnalysis);
			if (makeImage && resultImage != null) {
				resultImage = tintSoftTissue(resultImage, softRoi.softSieve);
			}
		}
		if (details.cOn) {
			final CorticalAnalysis cortAnalysis = new CorticalAnalysis(
				(SelectROI) roi);
			results = printCorticalResults(results, cortAnalysis);
			if (makeImage && resultImage != null) {
				resultImage = tintBoneStratec(resultImage, roi.sieve, roi.scaledImage,
					roi.details.marrowThreshold, cortAnalysis.cortexSieve);
			}

		}
		if (details.mOn) {
			final MassDistribution massDistribution = new MassDistribution(
				(SelectROI) roi, details, determineAlpha);
			results = printMassDistributionResults(results, massDistribution,
				details);
		}
		if (details.conOn) {
			final ConcentricRingAnalysis concentricRingAnalysis =
				new ConcentricRingAnalysis((SelectROI) roi, details, determineAlpha);
			results = printConcentricRingResults(results, concentricRingAnalysis,
				details);
			if (!details.dOn && makeImage && resultImage != null) {
				resultImage = drawPeriRadii(resultImage,
					concentricRingAnalysis.boneCenter, determineAlpha.pindColor,
					concentricRingAnalysis.rU, concentricRingAnalysis.theta);
				resultImage = drawMarrowCenter(resultImage, determineAlpha.alpha /
					Math.PI * 180.0, concentricRingAnalysis.boneCenter);
			}
		}

		if (details.dOn) {
			final DistributionAnalysis distributionAnalysis =
				new DistributionAnalysis((SelectROI) roi, details, determineAlpha);
			results = printDistributionResults(results, distributionAnalysis,
				details);
			if (makeImage && resultImage != null) {
				resultImage = drawRadii(resultImage, distributionAnalysis.marrowCenter,
					determineAlpha.pindColor, distributionAnalysis.r,
					distributionAnalysis.r2, distributionAnalysis.theta);
				resultImage = drawMarrowCenter(resultImage, determineAlpha.alpha /
					Math.PI * 180.0, distributionAnalysis.marrowCenter);
			}
		}

		if ((details.dOn || details.conOn) && makeImage && resultImage != null) {
			resultImage = drawRotated(resultImage, determineAlpha.alpha / Math.PI *
				180.0);
		}

		if (!details.suppressImages && resultImage != null) {
			resultImage = drawScale(resultImage, roi.pixelSpacing);
			resultImage.show();
		}
		if (details.saveImageOnDisk && resultImage != null) {
			if (details.suppressImages) {
				resultImage = drawScale(resultImage, roi.pixelSpacing);
			}
			final FileSaver fSaver = new FileSaver(resultImage);
			fSaver.saveAsPng(imageSavePath + imageName + ".png");
		}
		textPanel.appendLine(results);
		textPanel.updateDisplay();
	}

	public static String getInfoProperty(final String properties,
		final CharSequence propertyToGet)
	{
		final StringTokenizer st = new StringTokenizer(properties, "\n");
		String currentToken = null;
		while (st.hasMoreTokens()) {
			currentToken = st.nextToken();
			if (currentToken.contains(propertyToGet)) {
				break;
			}
		}
		if (currentToken == null) {
			return null;
		}

		final StringTokenizer st2 = new StringTokenizer(currentToken, ":");
		String token2 = null;
		while (st2.hasMoreTokens()) {
			token2 = st2.nextToken();
		}
		return token2 != null ? token2.trim() : null;
	}

	private static ImagePlus drawMarrowCenter(final ImagePlus tempImage,
		final double aplha, final double[] marrowCenter)
	{
		for (int i = 0; i < 10; i++) {// 45;i++) {//
			int x = ((int) (marrowCenter[0] + i));
			int y = ((int) (marrowCenter[1]));
			final ImageProcessor processor = tempImage.getProcessor();
			processor.setColor(new Color(0, 255, 255));
			processor.drawPixel(x, y);
			x = (int) (marrowCenter[0]);
			y = (int) (marrowCenter[1] + i);
			processor.setColor(new Color(255, 0, 255));
			processor.drawPixel(x, y);
			// Plot rotated axes...
			final double cos = i * Math.cos(-aplha / 180 * Math.PI);
			final double sin = i * Math.sin(-aplha / 180 * Math.PI);
			x = ((int) (marrowCenter[0] + cos));
			y = ((int) (marrowCenter[1] + sin));
			processor.setColor(new Color(0, 255, 0));
			processor.drawPixel(x, y);
			x = ((int) (marrowCenter[0] - sin));
			y = ((int) (marrowCenter[1] + cos));
			processor.setColor(new Color(0, 0, 255));
			processor.drawPixel(x, y);
		}
		return tempImage;
	}

	/*Concentric rings distribution result image*/
	private static ImagePlus drawPeriRadii(final ImagePlus tempImage,
		final double[] marrowCenter, final Vector<Integer> pindColor,
		final double[] r, final double[] theta)
	{
		// Draw unrotated radii
		for (int i = 0; i < theta.length; i++) {
			final int x = ((int) (marrowCenter[0] + r[i] * Math.cos(theta[i])));
			final int y = ((int) (marrowCenter[1] + r[i] * Math.sin(theta[i])));
			final double colorScale = pindColor.get(i) / 359.0;
			tempImage.getProcessor().setColor(new Color(0, (int) (255.0 * colorScale),
				(int) (255.0 * (1.0 - colorScale))));
			tempImage.getProcessor().drawPixel(x, y);
		}
		return tempImage;
	}

	private static ImagePlus drawRadii(final ImagePlus tempImage,
		final double[] marrowCenter, final Vector<Integer> pindColor,
		final double[] r, final double[] r2, final double[] theta)
	{
		// Draw unrotated radii
		for (int i = 0; i < 360; i++) {
			int x = ((int) (marrowCenter[0] + r[i] * Math.cos(theta[i])));
			int y = ((int) (marrowCenter[1] + r[i] * Math.sin(theta[i])));
			final double colorScale = pindColor.get(i) / 359.0;
			tempImage.getProcessor().setColor(new Color((int) (255.0 * colorScale), 0,
				(int) (255.0 * (1.0 - colorScale))));
			tempImage.getProcessor().drawPixel(x, y);
			x = ((int) (marrowCenter[0] + r2[i] * Math.cos(theta[i])));
			y = ((int) (marrowCenter[1] + r2[i] * Math.sin(theta[i])));
			tempImage.getProcessor().setColor(new Color(0, (int) (255.0 * colorScale),
				(int) (255.0 * (1.0 - colorScale))));
			tempImage.getProcessor().drawPixel(x, y);
		}
		return tempImage;
	}

	private static ImagePlus drawRotated(final ImagePlus tempImage,
		final double alpha)
	{
		tempImage.getProcessor().setBackgroundValue(0.0);
		tempImage.getProcessor().setInterpolationMethod(ImageProcessor.BICUBIC);
		final int width = tempImage.getWidth();
		final int height = tempImage.getHeight();
		final int hypot = (int) Math.sqrt(width * width + height * height);
		final ImageProcessor tIP;
		final int nW = (int) Math.abs(Math.ceil(Math.sin((alpha - 45.0) / 180.0 *
			Math.PI) * hypot));
		final int nH = (int) Math.abs(Math.ceil(Math.cos((alpha - 45.0) / 180.0 *
			Math.PI) * hypot));
		int nSize = Math.max(nW, nH);
		int offs = nSize - width;
		if (offs % 2 != 0) {
			offs++;
		}
		else {
			nSize = nSize + 1;
		}
		offs = offs / 2;
		tIP = expandImage(tempImage.getProcessor(), nSize, nSize, offs, offs);
		tempImage.setProcessor(null, tIP);
		tempImage.getProcessor().rotate(alpha);
		return tempImage;
	}

	private static ImagePlus drawScale(final ImagePlus tempImage,
		final double pixelSpacing)
	{
		final Calibration cal = new Calibration();
		cal.setUnit("mm");
		cal.pixelWidth = cal.pixelHeight = pixelSpacing;
		tempImage.setCalibration(cal);
		tempImage.getProcessor().setColor(new Color(255, 0, 0));
		// System.out.println("drawLine");
		tempImage.getProcessor().drawLine(5, 5, (int) (5.0 + 10.0 / pixelSpacing),
			5);
		tempImage.getProcessor().drawString("1 cm", 5, 20);
		return tempImage;
	}

	// TODO Can be called from there?
	/*Function taken from ij.plugin.CanvasResizer*/
	private static ImageProcessor expandImage(final ImageProcessor ipOld,
		final int wNew, final int hNew, final int xOff, final int yOff)
	{
		final ImageProcessor ipNew = ipOld.createProcessor(wNew, hNew);
		ipNew.setColor(new Color(0.0f, 0.0f, 0.0f));
		ipNew.setBackgroundValue(0.0);
		ipNew.fill();
		ipNew.insert(ipOld, xOff, yOff);
		return ipNew;
	}

	/*Get image into which we'll start adding stuff*/
	private static ImagePlus getRGBResultImage(final double[] values,
		final int width, final int height, final String path)
	{
		final ImagePlus tempImage = new ImagePlus(path + "Visual results");
		tempImage.setProcessor(new FloatProcessor(width, height, values));
		new ImageConverter(tempImage).convertToRGB();
		return tempImage;
	}

	private static String printAlpha(String results,
		final DetermineAlpha determineAlpha)
	{
		results += Double.toString(determineAlpha.alpha * 180 / Math.PI) + "\t";
		results += Double.toString(determineAlpha.rotationCorrection) + "\t";
		results += Double.toString(determineAlpha.distanceBetweenBones) + "\t";
		return results;
	}

	private static String printConcentricRingResults(final String results,
		final ConcentricRingAnalysis ringAnalysis,
		final ImageAndAnalysisDetails details)
	{
		final int limit = 360 / details.concentricSector;
		final StringBuilder resultsBuilder = new StringBuilder(results);
		for (int i = 0; i < limit; ++i) {
			resultsBuilder.append(ringAnalysis.pericorticalRadii[i]).append("\t");
		}
		for (int j = 0; j < details.concentricDivisions; ++j) {
			for (int i = 0; i < limit; ++i) {
				resultsBuilder.append(ringAnalysis.BMDs.get(j)[i]).append("\t");
			}
		}
		return resultsBuilder.toString();
	}

	private static String printCorticalResults(final String results,
		final CorticalAnalysis cortAnalysis)
	{
		final StringBuilder builder = new StringBuilder(results);
		DoubleStream.of(cortAnalysis.maMassD, cortAnalysis.stratecMaMassD,
			cortAnalysis.marrowDensity, cortAnalysis.marrowArea, cortAnalysis.bMD,
			cortAnalysis.area, cortAnalysis.CoD, cortAnalysis.CoA, cortAnalysis.sSI,
			cortAnalysis.sSIMax, cortAnalysis.sSIMin, cortAnalysis.iPo,
			cortAnalysis.iMax, cortAnalysis.iMin, cortAnalysis.dwIPo,
			cortAnalysis.dwIMax, cortAnalysis.dwIMin, cortAnalysis.ToD,
			cortAnalysis.ToA, cortAnalysis.medullaryArea, cortAnalysis.bSId).mapToObj(
				Double::toString).forEach(s -> builder.append(s).append("\t"));
		return builder.toString();
	}

	private static String printDistributionResults(final String results,
		final DistributionAnalysis distributionAnalysis,
		final ImageAndAnalysisDetails details)
	{
		final StringBuilder resultsBuilder = new StringBuilder(results);
		resultsBuilder.append(distributionAnalysis.peeledBMD).append("\t");
		// Radial distribution
		for (int i = 0; i < details.divisions; ++i) {
			resultsBuilder.append(distributionAnalysis.radialDistribution[i]).append(
				"\t");
		}
		final int iterations = 360 / details.sectorWidth;
		// Polar distribution
		for (int i = 0; i < iterations; ++i) {
			resultsBuilder.append(distributionAnalysis.polarDistribution[i]).append(
				"\t");
		}

		for (int pp = 0; pp < iterations; ++pp) {
			resultsBuilder.append(distributionAnalysis.endocorticalRadii[pp]).append(
				"\t");
		}
		for (int pp = 0; pp < iterations; ++pp) {
			resultsBuilder.append(distributionAnalysis.pericorticalRadii[pp]).append(
				"\t");
		}
		// Cortex bMD values
		for (int pp = 0; pp < iterations; ++pp) {
			resultsBuilder.append(distributionAnalysis.endoCorticalBMDs[pp]).append(
				"\t");
		}
		for (int pp = 0; pp < iterations; ++pp) {
			resultsBuilder.append(distributionAnalysis.midCorticalBMDs[pp]).append(
				"\t");
		}
		for (int pp = 0; pp < iterations; ++pp) {
			resultsBuilder.append(distributionAnalysis.periCorticalBMDs[pp]).append(
				"\t");
		}
		return resultsBuilder.toString();
	}

	private static String printMassDistributionResults(final String results,
		final MassDistribution massDistribution,
		final ImageAndAnalysisDetails details)
	{
		final StringBuilder resultsBuilder = new StringBuilder(results);
		for (int pp = 0; pp < (360 / details.sectorWidth); pp++) {
			resultsBuilder.append(massDistribution.bMCs[pp]).append("\t");
		}
		return resultsBuilder.toString();
	}

	private static String printSoftTissueResults(String results,
		final SoftTissueAnalysis softTissueAnalysis)
	{
		results += softTissueAnalysis.totalMuD + "\t";
		results += softTissueAnalysis.totalMuA + "\t";
		results += softTissueAnalysis.muD + "\t";
		results += softTissueAnalysis.muA + "\t";
		results += softTissueAnalysis.intraMuFatD + "\t";
		results += softTissueAnalysis.intraMuFatA + "\t";
		results += softTissueAnalysis.fatD + "\t";
		results += softTissueAnalysis.fatA + "\t";
		results += softTissueAnalysis.subCutFatDMedian + "\t";
		results += softTissueAnalysis.subCutFatD + "\t";
		results += softTissueAnalysis.subCutFatA + "\t";

		results += softTissueAnalysis.meD + "\t";
		results += softTissueAnalysis.meA + "\t";
		results += softTissueAnalysis.boneD + "\t";
		results += softTissueAnalysis.boneA + "\t";
		results += softTissueAnalysis.peeledD + "\t";
		results += softTissueAnalysis.peeledA + "\t";

		results += softTissueAnalysis.limbD + "\t";
		results += softTissueAnalysis.limbA + "\t";
		results += softTissueAnalysis.fatPercentage + "\t";
		return results;
	}

	/*Add bone sieve Stratec*/
	private static ImagePlus tintBoneStratec(final ImagePlus tempImage,
		final byte[] sieve, final double[] scaledImage,
		final double marrowThreshold, final byte[] stratecSieve)
	{
		for (int y = 0; y < tempImage.getHeight(); ++y) {
			for (int x = 0; x < tempImage.getWidth(); ++x) {
				final int value = tempImage.getProcessor().getPixel(x, y);
				final int[] rgb = new int[3];
				for (int i = 0; i < 3; ++i) {
					rgb[i] = (value >> (i * 8)) & 0XFF;
				}
				final int index = x + y * tempImage.getWidth();
				if (sieve[index] == 1) {
					// Tint bone area with purple
					tempImage.getProcessor().setColor(new Color(rgb[2], 0, rgb[0]));
					if (scaledImage[index] <= marrowThreshold) {
						// Tint marrow area with green
						if (rgb[0] < 255 - 50) {
							rgb[0] += 50;
						}
						tempImage.getProcessor().setColor(new Color(0, 0, rgb[0]));
					}
					tempImage.getProcessor().drawPixel(x, y);
				}
				if (stratecSieve[index] == 1) {
					// Tint stratec bone area with cyan
					tempImage.getProcessor().setColor(new Color(0, rgb[0], rgb[0]));
					tempImage.getProcessor().drawPixel(x, y);
				}
			}
		}
		return tempImage;
	}

	private static ImagePlus tintSoftTissue(final ImagePlus tempImage,
		final byte[] sieve)
	{
		for (int y = 0; y < tempImage.getHeight(); ++y) {
			for (int x = 0; x < tempImage.getWidth(); ++x) {
				final int value = tempImage.getProcessor().getPixel(x, y);
				final int[] rgb = new int[3];
				for (int i = 0; i < 3; ++i) {
					rgb[i] = (value >> (i * 8)) & 0XFF;
				}
				final byte pixel = sieve[x + y * tempImage.getWidth()];
				final Color sieveColor;
				switch (pixel) {
					case 2:
						sieveColor = new Color(rgb[2], rgb[1], 0);
						break;
					case 3:
						sieveColor = new Color(rgb[2], 0, 0);
						break;
					case 4:
						sieveColor = new Color(0, rgb[1], 0);
						break;
					case 5:
						sieveColor = new Color(rgb[2], 0, rgb[0]);
						break;
					default:
						continue;
				}
				tempImage.getProcessor().setColor(sieveColor);
				tempImage.getProcessor().drawPixel(x, y);
			}
		}
		return tempImage;
	}
}
