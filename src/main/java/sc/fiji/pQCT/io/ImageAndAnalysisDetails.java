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

package sc.fiji.pQCT.io;

public class ImageAndAnalysisDetails {

	public final boolean flipHorizontal;
	public final boolean flipVertical;
	public final boolean noFiltering;
	public final boolean sleeveOn;
	public final double scalingFactor;
	public final double constant;

	public final double airThreshold; // Fat lower threshold
	public final double fatThreshold; // Fat higher threshold
	public final double muscleThreshold; // Muscle lower threshold
	// Used with livewire to include intermuscular fat
	public final double edgeDivisions;
	public final double marrowThreshold; // Marrow higher threshold
	public final double softThreshold; // Soft tissues higher threshold
	// For cortical area analyses (CoA, sSI, I) + peeling distal pixels
	public final double areaThreshold;
	public final double rotationThreshold;
	public final double bMDThreshold; // For cortical bMD analyses
	// Thresholding bone from the rest and cortical area analyses (CoA, sSI, I)
	public final double boneThreshold;

	public final boolean cOn; // Basic analyses
	public final boolean mOn; // Mass distribution
	public final boolean conOn; // Concentric rings analysis
	public final boolean dOn; // Distribution analysis
	public final boolean stOn; // Soft tissue analysis

	public final int sectorWidth;
	public final int divisions;
	public final int concentricSector;
	public final int concentricDivisions;
	public final String roiChoice;
	public final String roiChoiceSt;
	public final String rotationChoice;
	public final String[] choiceLabels;
	public final String[] rotationLabels;
	public final boolean preventPeeling;
	public final boolean allowCleaving;
	public final boolean suppressImages;
	public final boolean manualRoi;
	public final boolean manualRotation;
	public final double manualAlpha;
	public final boolean guessFlip;
	public final boolean guessLarger;
	public final boolean guessStacked;
	public final boolean invertGuess;
	public final boolean saveImageOnDisk;
	public boolean flipDistribution;
	public boolean stacked;

	// ImageJ plugin constructor
	public ImageAndAnalysisDetails(final boolean[] defaultTopValues,
		final double[] thresholdsAndScaling, final String[] alignmentStrings,
		final String[] choiceLabels, final String[] rotationLabels,
		final boolean[] middleDefaults, final double manualAlpha,
		final boolean[] bottomDefaults, final int[] sectorsAndDivisions)
	{
		flipHorizontal = defaultTopValues[0];
		flipVertical = defaultTopValues[1];
		noFiltering = defaultTopValues[2];
		sleeveOn = defaultTopValues[3];

		airThreshold = thresholdsAndScaling[0];
		fatThreshold = thresholdsAndScaling[1];
		muscleThreshold = thresholdsAndScaling[2];
		edgeDivisions = thresholdsAndScaling[3];
		marrowThreshold = thresholdsAndScaling[4];
		softThreshold = thresholdsAndScaling[5];
		rotationThreshold = thresholdsAndScaling[6];
		areaThreshold = thresholdsAndScaling[7];
		bMDThreshold = thresholdsAndScaling[8];
		scalingFactor = thresholdsAndScaling[9];
		constant = thresholdsAndScaling[10];
		boneThreshold = areaThreshold;

		roiChoice = alignmentStrings[0];
		roiChoiceSt = alignmentStrings[1];
		rotationChoice = alignmentStrings[2];
		this.choiceLabels = choiceLabels;
		this.rotationLabels = rotationLabels;

		cOn = middleDefaults[0];
		mOn = middleDefaults[1];
		conOn = middleDefaults[2];
		dOn = middleDefaults[3];
		stOn = middleDefaults[4];
		preventPeeling = middleDefaults[5];
		allowCleaving = middleDefaults[6];
		suppressImages = middleDefaults[7];
		manualRoi = middleDefaults[8];
		manualRotation = middleDefaults[9];

		this.manualAlpha = manualAlpha;

		guessFlip = bottomDefaults[0];
		// Index 1 omitted on purpose
		guessLarger = bottomDefaults[2];
		stacked = bottomDefaults[3];
		guessStacked = bottomDefaults[4];
		invertGuess = bottomDefaults[5];
		flipDistribution = bottomDefaults[6];
		saveImageOnDisk = bottomDefaults[7];

		sectorWidth = sectorsAndDivisions[0];
		divisions = sectorsAndDivisions[1];
		concentricSector = sectorsAndDivisions[2];
		concentricDivisions = sectorsAndDivisions[3];
	}
}
