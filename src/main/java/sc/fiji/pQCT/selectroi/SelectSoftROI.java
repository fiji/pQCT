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

package sc.fiji.pQCT.selectroi;

import java.awt.Polygon;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.stream.IntStream;

import ij.ImagePlus;
import ij.gui.Roi;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.io.ScaledImageData;

//Debugging

public class SelectSoftROI extends RoiSelector {
	
	private final int RADIAL_DIVISIONS = 720;
	public byte[] eroded;

	// ImageJ constructor
	public SelectSoftROI(final ScaledImageData dataIn,
		final ImageAndAnalysisDetails detailsIn, final ImagePlus imp)
		throws ExecutionException
	{
		super(dataIn, detailsIn, imp);
		// Soft tissue analysis
		softSieve = null;
		if (details.stOn) {
			// Get rid of measurement tube used at the UKK institute
			final byte[] sleeve;
			if (details.sleeveOn) {
				sleeve = removeSleeve(softScaledImage, 25.0);
				final int size = width * height;
				IntStream.range(0, size).filter(i -> sleeve[i] == 1).forEach(
					i -> softScaledImage[i] = minimum);
			}

			// Ignore data outside manually selected ROI, if manualRoi has been
			// selected
			final Roi ijROI = imp.getRoi();
			if (ijROI != null && details.manualRoi) {
				// Set pixels outside the manually selected ROI to zero
				final double[] tempScaledImage = clone(softScaledImage);
				// Check whether pixel is within ROI, mark with bone threshold
				for (int j = 0; j < height; j++) {
					for (int i = 0; i < width; i++) {
						if (!ijROI.contains(i, j)) {
							softScaledImage[i + j * width] = minimum;
						}
					}
				}
				// Check whether a polygon can be acquired and include polygon points
				// too
				final Polygon polygon = ijROI.getPolygon();
				if (polygon != null) {
					for (int j = 0; j < polygon.npoints; j++) {
						final int index = polygon.xpoints[j] + polygon.ypoints[j] * width;
						softScaledImage[index] = tempScaledImage[index];
					}
				}
			}

			final Vector<Object> masks = getSieve(softScaledImage, airThreshold,
				details.roiChoiceSt, details.guessStacked, details.stacked, false,
				true);
			softSieve = (byte[]) masks.get(0);

			// Erode three layers of pixels from the fat sieve to get rid of higher
			// density layer (i.e. skin) on top of fat to enable finding muscle border
			byte[] muscleSieve = clone(softSieve);
			final double[] muscleImage = clone(softScaledImage);

			// Remove skin by eroding three layers of pixels
			for (int i = 0; i < 3; ++i) {
				muscleSieve = erode(muscleSieve);
			}
			// The three layers of skin removed
			final byte[] subCutaneousFat = clone(muscleSieve);

			// Remove everything other than the selected limb from the image
			for (int i = 0; i < muscleSieve.length; ++i) {
				if (muscleSieve[i] < 1) {
					muscleImage[i] = minimum;
				}
			}
			// Look for muscle outline
			final Vector<Object> muscleMasks = getSieve(muscleImage,
				details.muscleThreshold, "Bigger", details.guessStacked,
				details.stacked, false, false);
			final List<DetectedEdge> muscleEdges = (Vector<DetectedEdge>) muscleMasks
				.get(2);
			muscleEdges.sort(Collections.reverseOrder());
			int tempMuscleArea = 0;
			muscleSieve = new byte[softSieve.length];
			int areaToAdd = 0;
			// Include areas that contribute more than 1.0% on top of what is already
			// included
			while (areaToAdd < muscleEdges.size() && tempMuscleArea *
				0.01 < muscleEdges.get(areaToAdd).area)
			{
				final byte[] tempMuscleSieve = fillSieve(muscleEdges.get(areaToAdd).iit,
					muscleEdges.get(areaToAdd).jiit, width, height, muscleImage,
					details.muscleThreshold);
				for (int i = 0; i < tempMuscleSieve.length; ++i) {
					if (tempMuscleSieve[i] > 0) {
						muscleSieve[i] = tempMuscleSieve[i];
					}
				}
				tempMuscleArea += muscleEdges.get(areaToAdd).area;
				areaToAdd++;
			}

			// Dilate the sieve to include all muscle pixels
			byte[] tempMuscleSieve = clone(muscleSieve);
			tempMuscleSieve = dilate(tempMuscleSieve, (byte) 1, (byte) 0, (byte) 2);
			
			eroded = new byte[softSieve.length];
			for (int i = 0; i < tempMuscleSieve.length; ++i) {
				if (tempMuscleSieve[i] == 1) {
					subCutaneousFat[i] = 0;
				}
			}

			// create temp boneResult to wipe out bone and marrow
			final Vector<Object> masks2 = getSieve(softScaledImage, softThreshold, details.roiChoiceSt,
					details.guessStacked, details.stacked, false, false);
			final byte[] boneResult = (byte[]) masks2.get(1);
			

			for (int i = 0; i < softSieve.length; ++i) {
				if (softSieve[i] == 1 && softScaledImage[i] >= airThreshold &&
					softScaledImage[i] < fatThreshold)
				{
					// Fat
					softSieve[i] = 2;
				}
				if (muscleSieve[i] == 1 && boneResult[i] == 0) {
					if (softScaledImage[i] >= muscleThreshold &&
						softScaledImage[i] < softThreshold)
					{
						// Muscle
						softSieve[i] = 3;
					}
					if (softScaledImage[i] >= airThreshold &&
						softScaledImage[i] < muscleThreshold)
					{
						// Intra/Intermuscular fat
						softSieve[i] = 4;
					}
				}
				if (subCutaneousFat[i] == 1) {
					// Subcut fat
					softSieve[i] = 5;
				}
				if (boneResult[i] == 1) {
					if (softScaledImage[i] >= fatThreshold) {
						// Bone & marrow
						softSieve[i] = 6;
					}
					else {
						// Marrow fat
						softSieve[i] = 7;
					}
				}
				if (softSieve[i] > 0 && subCutaneousFat[i] == 0 &&
					tempMuscleSieve[i] == 0)
				{
					// Skin eroded pixels
					eroded[i] = 1;
				}
			}
		}
	}
}
