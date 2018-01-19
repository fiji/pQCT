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

package sc.fiji.pQCT.utils;

import ij.ImagePlus;
import ij.text.TextPanel;
import sc.fiji.pQCT.DistributionAnalysisPlugIn;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;

//TODO Refactor into Distribution_Analysis (make a subpackage where the massive class is chopped up?)
public class ResultsWriter {

	private final String imageInfo;
	private final boolean alphaOn;

	public ResultsWriter(final String imageInfo, final boolean alphaOn) {
		this.imageInfo = imageInfo;
		this.alphaOn = alphaOn;
	}

	public String printResults(final String results,
		final ImageAndAnalysisDetails details, final ImagePlus imp)
	{
		final String[] propertyNames = { "File Name", "Patient's Name",
			"Patient ID", "Patient's Birth Date", "Acquisition Date", "Pixel Spacing",
			"ObjLen" };
		final String[] parameters = { Double.toString(details.airThreshold), Double
			.toString(details.fatThreshold), Double.toString(details.muscleThreshold),
			Double.toString(details.marrowThreshold), Double.toString(
				details.softThreshold), Double.toString(details.rotationThreshold),
			Double.toString(details.areaThreshold), Double.toString(
				details.bMDThreshold), Double.toString(details.scalingFactor), Double
					.toString(details.constant) };
		final StringBuilder resultsBuilder = new StringBuilder(results);
		if (imp != null) {
			if (DistributionAnalysisPlugIn.getInfoProperty(imageInfo,
				"File Name") != null)
			{
				resultsBuilder.append(DistributionAnalysisPlugIn.getInfoProperty(
					imageInfo, "File Path"));
				resultsBuilder.append(DistributionAnalysisPlugIn.getInfoProperty(
					imageInfo, "File Name")).append("\t");
			}
			else {
				if (imp.getImageStackSize() == 1) {
					resultsBuilder.append(DistributionAnalysisPlugIn.getInfoProperty(
						imageInfo, "Title")).append("\t");
				}
				else {
					resultsBuilder.append(imageInfo.substring(0, imageInfo.indexOf("\n")))
						.append("\t");
				}
			}
			for (int i = 1; i < propertyNames.length; ++i) {
				resultsBuilder.append(DistributionAnalysisPlugIn.getInfoProperty(
					imageInfo, propertyNames[i])).append("\t");
			}
		}

		for (final String parameter : parameters) {
			resultsBuilder.append(parameter).append("\t");
		}

		resultsBuilder.append(details.manualRotation).append("\t");
		resultsBuilder.append(details.flipDistribution).append("\t");
		resultsBuilder.append(details.guessFlip).append("\t");
		resultsBuilder.append(details.guessLarger).append("\t");
		resultsBuilder.append(details.stacked).append("\t");
		resultsBuilder.append(details.invertGuess).append("\t");
		resultsBuilder.append(details.allowCleaving).append("\t");
		resultsBuilder.append(details.preventPeeling).append("\t");
		resultsBuilder.append(details.roiChoice).append("\t");
		resultsBuilder.append(details.rotationChoice).append("\t");
		resultsBuilder.append(details.flipHorizontal).append("\t");
		resultsBuilder.append(details.flipVertical).append("\t");
		return resultsBuilder.toString();
	}

	public void writeHeader(final TextPanel textPanel,
		final ImageAndAnalysisDetails details)
	{
		final StringBuilder headings = new StringBuilder(String.join("\t",
			"File Name", "Patient's Name", "Patient ID", "Patient's Birth Date",
			"Acquisition Date", "Pixel Spacing", "Object Length", "Air Threshold",
			"Fat Threshold", "Muscle Threshold", "Marrow Threshold", "Soft Threshold",
			"Rotation Threshold", "Area Threshold", "bMD Threshold",
			"Scaling Coefficient", "Scaling Constant", "Manual Rotation",
			"Flip Distribution", "Guess right", "Guess larger", "Stacked bones",
			"Invert guess", "Allow Cleaving", "Prevent PVE peeling", "Roi choice",
			"Rotation choice", "Flip Horizontal", "Flip Vertical"));
		if (alphaOn) {
			headings.append("\t").append(String.join("\t", "Alpha [deg]",
				"Rotation correction [deg]", "Distance between bones[mm]"));
		}

		if (details.stOn) {
			headings.append("\t").append(String.join("\t", "muD [mg/cm3]",
				"muA [cm2]", "LeanMuD [mg/cm3]", "LeanMuA [cm2]", "IntraFatD [mg/cm3]",
				"IntraFatA [cm2]", "fatD [mg/cm3]", "fatA [cm2]",
				"subCutFatDMedian [mg/cm3]", "subCutFatD [mg/cm3]", "subCutFatA [cm2]",
				"MedD [mg/cm3]", "MedA [cm2]", "boneD [mg/cm3]", "boneA [cm2]",
				"peeledD [mg/cm3]", "peeledA [cm2]", "limbD [mg/cm3]", "limbA [cm2]",
				"Density weighted fat percentage [%]"));
		}

		if (details.cOn) {
			headings.append("\t").append(String.join("\t", "maMassD [g/cm3]",
				"stratecMaMassD [g/cm3]", "marrowDensity [mg/cm3]", "marrowArea [mm2]",
				"CoD [mg/cm3]", "CoA [mm2]", "Stratec CoD [mg/cm3]",
				"Stratec CoA [mm2]", "sSI [mm3]", "SSImax [mm3]", "SSImin [mm3]",
				"iPo [mm4]", "Imax [mm4]", "Imin [mm4]", "dwIPo [mg/cm]",
				"dwImax [mg/cm]", "dwImin [mg/cm]", "ToD [mg/cm3]", "ToA[mm2]",
				"medullaryArea [mm2]", "bSId[g/cm4]"));
		}
		if (details.mOn) {
			for (int i = 0; i < (360 / details.sectorWidth); ++i) {
				headings.append("\t").append(i * details.sectorWidth).append("?- ")
					.append((i + 1) * details.sectorWidth).append("?mineral mass [mg]");
			}
		}

		if (details.conOn) {
			for (int i = 0; i < (360 / details.concentricSector); ++i) {
				headings.append("\t").append(i * details.concentricSector).append("?- ")
					.append((i + 1) * details.concentricSector).append(
						"?concentric analysis pericortical radius [mm]");
			}
			for (int j = 0; j < details.concentricDivisions; ++j) {
				for (int i = 0; i < (360 / details.concentricSector); ++i) {
					headings.append("\t").append("Division ").append(j + 1).append(
						" sector ").append(i * details.concentricSector).append("?- ")
						.append((i + 1) * details.concentricSector).append(
							"?vBMD [mg/cm3]");
				}
			}
		}

		if (details.dOn) {
			headings.append("\t").append("Peeled mean vBMD [mg/cm3]");
			// Radial distribution
			for (int i = 0; i < details.divisions; ++i) {
				headings.append("\t").append("Radial division ").append(i).append(
					" vBMD [mg/cm3]");
			}

			final int iterations = (360 / details.sectorWidth);
			for (int i = 0; i < iterations; ++i) {
				final String rowStart = (i * details.sectorWidth) + "?- " + ((i + 1) *
					details.sectorWidth);
				headings.append("\t").append("Polar sector ").append(i).append(
					" vBMD [mg/cm3]");
				headings.append("\t").append(rowStart).append(
					" ?endocortical radius [mm]");
				headings.append("\t").append(rowStart).append(
					" ?pericortical radius [mm]");
				headings.append("\t").append(rowStart).append(
					" ?endocortical vBMD [mg/cm3]");
				headings.append("\t").append(rowStart).append(
					" ?midcortical vBMD [mg/cm3]");
				headings.append("\t").append(rowStart).append(
					" ?pericortical vBMD [mg/cm3]");
			}
		}
		textPanel.setColumnHeadings(headings.toString());
	}
}
