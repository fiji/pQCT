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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ij.ImagePlus;
import ij.gui.Roi;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.io.ScaledImageData;

public class SelectSoftROI extends RoiSelector {

	private static final int RADIAL_DIVISIONS = 720;
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
			/*Get rid of measurement tube used at the UKK institute*/
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
				final double[] tempScaledImage = softScaledImage.clone();
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

			byte[] muscleSieve = softSieve.clone();
			final double[] muscleImage = softScaledImage.clone();
			final byte[] subCutaneousFat = softSieve.clone();
			// Remove skin by eroding three layers of pixels
			for (int i = 0; i < 3; ++i) {
				erode(muscleSieve);
			}
			// Remove everything other than the selected limb from the image
			for (int i = 0; i < muscleSieve.length; ++i) {
				if (muscleSieve[i] < 1) {
					muscleImage[i] = minimum;
				}
			}
			/*Look for muscle outline*/
			final Vector<Object> muscleMasks = getSieve(muscleImage,
				details.muscleThreshold, "Bigger", details.guessStacked,
				details.stacked, false, false);
			final List<DetectedEdge> muscleEdges = (Vector<DetectedEdge>) muscleMasks
				.get(2);
			muscleEdges.sort(Collections.reverseOrder());
			int tempMuscleArea = 0;
			muscleSieve = new byte[softSieve.length];
			int areaToAdd = 0;
			/*Include areas that contribute more than 1.0% on top of what is already included*/
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

			// Re-segment soft-tissues using livewire based on the muscleSieve:
			// 1) bring rays back from image edges to centre of soft-tissue mask 1 deg
			// apart
			// 2) use livewire on the 360 edge pixels
			// 3) rotate livewire init pixels around a few times to get the
			// segmentation to go through subcut/intramuscular fat
			final Vector<Object> masks2 = getSieve(softScaledImage, softThreshold,
				details.roiChoiceSt, details.guessStacked, details.stacked, false,
				false);
			final byte[] boneResult = (byte[]) masks2.get(1);
			final Collection<Double> ii = new ArrayList<>();
			final Collection<Double> jj = new ArrayList<>();
			final double[] softCentre = new double[2];

			for (int i = 0; i < width; ++i) {
				for (int j = 0; j < height; ++j) {
					if (boneResult[i + j * width] > 0) {
						ii.add((double) i);
						jj.add((double) j);
						softCentre[0] += i;
						softCentre[1] += j;
					}
				}
			}

			softCentre[0] /= ii.size();
			softCentre[1] /= jj.size();

			final double maxR = Math.sqrt(Math.pow(Math.max(Math.abs(softCentre[0] -
				width), softCentre[0]), 2.0) + Math.pow(Math.max(Math.abs(
					softCentre[1] - height), softCentre[1]), 2.0));
			final double[] rs = new double[RADIAL_DIVISIONS];
			double r;
			double t;
			final double[] theta = new double[RADIAL_DIVISIONS];
			final int[][] edgeCoords = new int[RADIAL_DIVISIONS][2];
			// Get the extremes of muscle area with 1 deg increments in polar
			// coordinates
			for (int i = 0; i < RADIAL_DIVISIONS; ++i) {
				r = maxR;
				t = i / (RADIAL_DIVISIONS * 2.0 * Math.PI);
				while (Math.round(r * Math.cos(t) + softCentre[0]) < 0 || Math.round(r *
					Math.sin(t) + softCentre[1]) < 0 || Math.round(r * Math.cos(t) +
						softCentre[0]) >= width || Math.round(r * Math.sin(t) +
							softCentre[1]) >= height)
				{
					r -= 0.1;
				}
				while (Math.round(r * Math.cos(t) + softCentre[0]) >= 0 && Math.round(
					r * Math.sin(t) + softCentre[1]) >= 0 && Math.round(r * Math.cos(t) +
						softCentre[0]) < width && Math.round(r * Math.sin(t) +
							softCentre[1]) < height && muscleSieve[(int) (Math.round(r * Math
								.cos(t) + softCentre[0]) + Math.round(r * Math.sin(t) +
									softCentre[1]) * width)] < 1)
				{
					r -= 0.1;
				}
				rs[i] = r;
				theta[i] = t;
				edgeCoords[i][0] = (int) (Math.round(rs[i] * Math.cos(theta[i]) +
					softCentre[0]));
				edgeCoords[i][1] = (int) (Math.round(rs[i] * Math.sin(theta[i]) +
					softCentre[1]));
			}

			// Arraylists for edge, and livewire seed coordinates
			ArrayList<Integer> edgeii = new ArrayList<>();
			ArrayList<Integer> edgejj = new ArrayList<>();

			// Create list of seed coordinates
			for (int i = 0; i < RADIAL_DIVISIONS; ++i) {
				edgeii.add(edgeCoords[i][0]);
				edgejj.add(edgeCoords[i][1]);
			}
			final Vector<Object> tempVO = getLassoEdge(edgeii, edgejj, softCentre,
				muscleSieve);
			edgeii = (ArrayList<Integer>) tempVO.get(0);
			edgejj = (ArrayList<Integer>) tempVO.get(1);

			// Fill in muscle mask with inter-muscular fat
			muscleSieve = getByteMask(width, height, edgeii, edgejj);
			// Dilate the sieve to include all muscle pixels
			muscleSieve = dilateMuscleMask(muscleSieve, softScaledImage, width,
				height, muscleThreshold);
			/*Re-segmenting done*/

			// Wipe muscle area +3 layer of pixels away from subcut.
			final byte[] tempMuscleSieve = Arrays.copyOf(muscleSieve,
				muscleSieve.length);
			eroded = new byte[softSieve.length];
			for (int i = 0; i < tempMuscleSieve.length; ++i) {
				if (tempMuscleSieve[i] == 1) {
					subCutaneousFat[i] = 0;
				}
			}

			for (int i = 0; i < softSieve.length; ++i) {
				if (softSieve[i] == 1 && softScaledImage[i] >= airThreshold &&
					softScaledImage[i] < fatThreshold)
				{
					softSieve[i] = 2; // Fat
				}
				if (muscleSieve[i] == 1 && boneResult[i] == 0) {
					if (softScaledImage[i] >= muscleThreshold &&
						softScaledImage[i] < softThreshold)
					{
						softSieve[i] = 3; // Muscle
					}
					if (softScaledImage[i] >= airThreshold &&
						softScaledImage[i] < muscleThreshold)
					{
						softSieve[i] = 4; // Intra/Intermuscular fat
					}
				}
				if (subCutaneousFat[i] == 1) {
					softSieve[i] = 5; // Subcut fat
				}
				if (boneResult[i] == 1) {
					if (softScaledImage[i] >= fatThreshold) {
						softSieve[i] = 6; // Bone & marrow
					}
					else {
						softSieve[i] = 7; // Marrow fat
					}
				}
				if (softSieve[i] > 0 && subCutaneousFat[i] == 0 &&
					tempMuscleSieve[i] == 0)
				{
					eroded[i] = 1; // Skin eroded pixels
				}
			}
		}
	}

	private static byte[] dilateMuscleMask(final byte[] mask,
		final double[] softScaledImage, final int width, final int height,
		final double threshold)
	{
		final ArrayList<Integer> initialI = new ArrayList<>();
		final ArrayList<Integer> initialJ = new ArrayList<>();
		int i;
		int j;
		for (i = 0; i < width; ++i) {
			for (j = 0; j < height; ++j) {
				if (mask[i + j * width] == 1) {
					initialI.add(i);
					initialJ.add(j);
				}
			}
		}

		while (!initialI.isEmpty() && initialI.get(initialI.size() - 1) > 0 &&
			initialI.get(initialI.size() - 1) < width - 1 && initialJ.get(initialJ
				.size() - 1) > 0 && initialJ.get(initialJ.size() - 1) < height - 1)
		{
			i = initialI.get(initialI.size() - 1);
			j = initialJ.get(initialJ.size() - 1);
			initialI.remove(initialI.size() - 1);
			initialJ.remove(initialJ.size() - 1);

			if (mask[i + j * width] == 0 && softScaledImage[i + j *
				width] >= threshold)
			{
				mask[i + j * width] = 1;
			}

			if (mask[i - 1 + j * width] == 0 && softScaledImage[i + j *
				width] >= threshold)
			{
				initialI.add(i - 1);
				initialJ.add(j);
			}

			if (mask[i + 1 + j * width] == 0 && softScaledImage[i + j *
				width] >= threshold)
			{
				initialI.add(i + 1);
				initialJ.add(j);
			}

			if (mask[i + (j - 1) * width] == 0 && softScaledImage[i + j *
				width] >= threshold)
			{
				initialI.add(i);
				initialJ.add(j - 1);
			}

			if (mask[i + (j + 1) * width] == 0 && softScaledImage[i + j *
				width] >= threshold)
			{
				initialI.add(i);
				initialJ.add(j + 1);
			}

		}
		return mask;
	}

	private static byte[] fillBorder(final byte[] mask, final int width) {
		final ArrayList<Integer> initialI = new ArrayList<>();
		final ArrayList<Integer> initialJ = new ArrayList<>();
		initialI.add(0);
		initialJ.add(0);
		int i;
		int j;
		while (!initialI.isEmpty()) {
			final int lastI = initialI.size() - 1;
			final int lastJ = initialJ.size() - 1;
			i = initialI.get(lastI);
			j = initialJ.get(lastJ);
			initialI.remove(lastI);
			initialJ.remove(lastJ);
			final int offset = j * width;
			if (mask[i + offset] == 0) {
				mask[i + offset] = 1;
			}
			// TODO Loop neighbourhood with "inBounds"
			if (i - 1 >= 0 && mask[i - 1 + offset] == 0) {
				initialI.add(i - 1);
				initialJ.add(j);
			}

			if (i + 1 < width && mask[i + 1 + offset] == 0) {
				initialI.add(i + 1);
				initialJ.add(j);
			}

			if (j - 1 >= 0 && mask[i + (j - 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j - 1);
			}

			if (j + 1 < width && mask[i + (j + 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j + 1);
			}

		}
		return mask;
	}

	// TODO Refactor so that fillBorder and -Mask are one method
	private static byte[] fillMask(int i, int j, final byte[] mask,
		final int width, final int height)
	{
		final ArrayList<Integer> initialI = new ArrayList<>();
		final ArrayList<Integer> initialJ = new ArrayList<>();
		initialI.add(i);
		initialJ.add(j);
		while (!initialI.isEmpty() && initialI.get(initialI.size() - 1) > 0 &&
			initialI.get(initialI.size() - 1) < width - 1 && initialJ.get(initialJ
				.size() - 1) > 0 && initialJ.get(initialJ.size() - 1) < height - 1)
		{
			i = initialI.get(initialI.size() - 1);
			j = initialJ.get(initialJ.size() - 1);
			initialI.remove(initialI.size() - 1);
			initialJ.remove(initialJ.size() - 1);

			if (mask[i + j * width] == 0) {
				mask[i + j * width] = 1;
			}

			if (mask[i - 1 + j * width] == 0) {
				initialI.add(i - 1);
				initialJ.add(j);
			}

			if (mask[i + 1 + j * width] == 0) {
				initialI.add(i + 1);
				initialJ.add(j);
			}

			if (mask[i + (j - 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j - 1);
			}

			if (mask[i + (j + 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j + 1);
			}

		}
		return mask;
	}

	private static int[] findMaskFillInit(final byte[] mask, final int width,
		final int height, final ArrayList<Integer> edgeii,
		final ArrayList<Integer> edgejj)
	{
		byte[] tempMask = mask.clone();
		tempMask = fillBorder(tempMask, width);
		final int[] returnCoordinates = new int[2];
		final int[] steer = new int[2];
		for (int j = 0; j < edgeii.size() - 1; ++j) {
			returnCoordinates[0] = edgeii.get(j);
			returnCoordinates[1] = edgejj.get(j);
			double direction = Math.atan2(edgejj.get(j + 1) - returnCoordinates[1],
				edgeii.get(j + 1) - returnCoordinates[0]);
			direction += Math.PI / 4.0;
			for (int i = 0; i < 3; ++i) {
				steer[0] = (int) Math.round(Math.cos(direction));
				steer[1] = (int) Math.round(Math.sin(direction));
				/*Handle OOB*/
				while ((returnCoordinates[0] + steer[0]) < 0 || (returnCoordinates[0] +
					steer[0]) >= width || (returnCoordinates[1] + steer[1]) < 0 ||
					(returnCoordinates[1] + steer[1]) >= height)
				{
					direction += Math.PI / 4.0;
					steer[0] = (int) Math.round(Math.cos(direction));
					steer[1] = (int) Math.round(Math.sin(direction));
				}
				if (tempMask[returnCoordinates[0] + steer[0] + (returnCoordinates[1] +
					steer[1]) * width] == 0)
				{
					returnCoordinates[0] += steer[0];
					returnCoordinates[1] += steer[1];
					return returnCoordinates;
				}
				direction += Math.PI / 4.0;
			}
		}
		return null;
	}

	private static byte[] getByteMask(final int width, final int height,
		final ArrayList<Integer> edgeii, final ArrayList<Integer> edgejj)
	{
		final byte[] mask = new byte[width * height];
		for (int i = 0; i < edgeii.size(); ++i) {
			mask[edgeii.get(i) + edgejj.get(i) * width] = 1;
		}
		final int[] fillInitCoords = findMaskFillInit(mask, width, height, edgeii,
			edgejj);
		if (fillInitCoords == null) {
			return mask;
		}
		return fillMask(fillInitCoords[0], fillInitCoords[1], mask, width, height);
	}

	private Vector<Object> getLassoEdge(final ArrayList<Integer> edgeii,
		final ArrayList<Integer> edgejj, final double[] softCentre,
		final byte[] image)
	{

		final int sectorToConsider = (int) (details.edgeDivisions / 360.0 *
			(RADIAL_DIVISIONS));

		// Pop edge into DetectedRadialEdgeTheta, sort by incrementing radius
		// Start the algorithm from the most distant point from the centre of area
		final List<DetectedRadialEdgeTheta> radialEdge = new Vector<>();
		double ii;
		double jj;
		for (int i = 0; i < edgeii.size(); ++i) {
			ii = edgeii.get(i) - softCentre[0];
			jj = edgejj.get(i) - softCentre[1];
			radialEdge.add(new DetectedRadialEdgeTheta(Math.sqrt(Math.pow(ii, 2.0) +
				Math.pow(jj, 2.0)), i));
		}
		// Get maximal radius, and it's location
		Collections.sort(radialEdge);
		final int maxIndex = radialEdge.get(radialEdge.size() - 1).index;
		final int[] indices = new int[RADIAL_DIVISIONS];
		int count = 0;
		for (int i = maxIndex; i < RADIAL_DIVISIONS; ++i) {
			indices[count] = i;
			++count;
		}
		for (int i = 0; i < maxIndex; ++i) {
			indices[count] = i;
			++count;
		}
		int currentI = 0;
		final ArrayList<Integer> boundaryIndices = new ArrayList<>();
		final ArrayList<Integer> currentIndices = new ArrayList<>();
		// Loop through the boundary pixels. Look for the furthest point which can
		// be seen from the current point by checking if a line between the points
		// has any roi points in it
		boundaryIndices.add(0);
		while (currentI < (RADIAL_DIVISIONS - 1)) {
			// Go through all point pairs from the furthest, select the first, which
			// can be seen. Check the next 179 points
			currentIndices.clear();
			for (int i = currentI; i < Math.min(currentI + (sectorToConsider - 1),
				RADIAL_DIVISIONS); ++i)
			{
				currentIndices.add(i);
			}

			int sightIndex = 1;
			for (int i = currentIndices.size() - 1; i > 0; --i) {
				final List<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii
					.get(indices[currentI]), edgejj.get(indices[currentI])),
					new Coordinate(edgeii.get(indices[currentIndices.get(i)]), edgejj.get(
						indices[currentIndices.get(i)])));
				if (isPathFree(image, tempCoordinates)) {
					sightIndex = i;
					break;
				}
			}
			currentI = currentI + sightIndex;
			if (currentI <= (RADIAL_DIVISIONS - 1)) {
				boundaryIndices.add(currentI);
			}
		}

		// Create the boundary
		final Collection<Integer> returnii = new ArrayList<>();
		final Collection<Integer> returnjj = new ArrayList<>();
		for (int i = 1; i < boundaryIndices.size(); ++i) {
			final List<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii
				.get(indices[boundaryIndices.get(i - 1)]), edgejj.get(
					indices[boundaryIndices.get(i - 1)])), new Coordinate(edgeii.get(
						indices[boundaryIndices.get(i)]), edgejj.get(indices[boundaryIndices
							.get(i)])));
			for (int j = 1; j < tempCoordinates.size(); ++j) {
				returnii.add((int) tempCoordinates.get(j).ii);
				returnjj.add((int) tempCoordinates.get(j).jj);
			}
		}

		// Add the final missing bit
		final List<Coordinate> tempCoordinates = getLine(new Coordinate(edgeii.get(
			indices[boundaryIndices.get(boundaryIndices.size() - 1)]), edgejj.get(
				indices[boundaryIndices.get(boundaryIndices.size() - 1)])),
			new Coordinate(edgeii.get(indices[boundaryIndices.get(0)]), edgejj.get(
				indices[boundaryIndices.get(0)])));
		for (int j = 1; j < tempCoordinates.size(); ++j) {
			returnii.add((int) tempCoordinates.get(j).ii);
			returnjj.add((int) tempCoordinates.get(j).jj);
		}
		final Vector<Object> returnV = new Vector<>();
		returnV.add(returnii);
		returnV.add(returnjj);
		return returnV;
	}

	/*
	    Get the line coordinates from origin to target
	    Digital Differential Analyzer (DDA) algorithm for line
	    http://www.tutorialspoint.com/computer_graphics/line_generation_algorithm.htm
	*/
	private static List<Coordinate> getLine(final Coordinate origin,
		final Coordinate target)
	{
		final Coordinate difference = target.subtract(origin);
		final double steps = difference.maxVal();
		final int points = (int) (steps + 1);
		final double iiIncrement = difference.ii / steps;
		final double jjIncrement = difference.jj / steps;
		return IntStream.range(0, points).mapToObj(i -> {
			final double ii = Math.round(origin.ii + i * iiIncrement);
			final double jj = Math.round(origin.jj + i * jjIncrement);
			return new Coordinate(ii, jj);
		}).collect(Collectors.toList());
	}

	// Helper function to check whether a path is free or blocked
	private boolean isPathFree(final byte[] image,
		final List<Coordinate> pathCoordinates)
	{
		final int size = pathCoordinates.size();
		if (size < 2) {
			return true;
		}
		final List<Coordinate> coordinates = pathCoordinates.subList(1, size - 1);
		return coordinates.stream().mapToInt(c -> (int) (c.ii + c.jj * width))
			.filter(i -> image[i] == 1).count() <= 2;
	}

}
