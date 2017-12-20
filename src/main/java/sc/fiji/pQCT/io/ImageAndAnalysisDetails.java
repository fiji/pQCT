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

    Copyright (C) 2011 Timo Rantalainen
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
