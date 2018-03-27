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

import static java.util.stream.Collectors.toList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

import ij.ImagePlus;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.io.ScaledImageData;

public abstract class RoiSelector {

	public final ImageAndAnalysisDetails details;
	public final double[] scaledImage;
	public final double[] softScaledImage;
	public final double minimum;
	public final int height;
	public final int width;
	// For cortical area analyses (CoA, sSI, I), and peeling distal pixels
	public final double areaThreshold;
	// For cortical bMD analyses
	public final double BMDthreshold;
	public final ScaledImageData scaledImageData;
	public final ImagePlus imp;
	final double airThreshold;
	final double fatThreshold;
	final double muscleThreshold;
	// Thresholding soft tissues + marrow from bone
	final double softThreshold;
	public double[] cortexROI;
	public Vector<Integer> boneMarrowRoiI;
	public Vector<Integer> boneMarrowRoiJ;
	// For bMD analyses
	public Vector<Integer> cortexRoiI;
	// For bMD analyses
	public Vector<Integer> cortexRoiJ;
	// For area analyses
	public Vector<Integer> cortexAreaRoiI;
	// For area analyses
	public Vector<Integer> cortexAreaRoiJ;
	public Vector<Integer> area;
	public Vector<Integer> length;
	public int selection;
	public double pixelSpacing;
	public byte[] sieve;
	// Mask for soft tissues
	public byte[] softSieve;
	// Will contain filled bones
	byte[] result;

	RoiSelector(final ScaledImageData dataIn,
		final ImageAndAnalysisDetails detailsIn, final ImagePlus imp)
	{
		scaledImageData = dataIn;
		this.imp = imp;
		details = detailsIn;
		scaledImage = dataIn.scaledImage.clone();
		softScaledImage = dataIn.softScaledImage.clone();
		pixelSpacing = dataIn.pixelSpacing;
		width = dataIn.width;
		height = dataIn.height;
		airThreshold = details.airThreshold;
		fatThreshold = details.fatThreshold;
		muscleThreshold = details.muscleThreshold;
		// For cortical area analyses (CoA, sSI, I) + peeling distal pixes
		areaThreshold = details.areaThreshold;
		// For cortical bMD analyses
		BMDthreshold = details.bMDThreshold;
		// Thresholding soft tissues + marrow from bone
		softThreshold = details.softThreshold;
		minimum = dataIn.minimum;
	}

	public byte[] fillSieve(final Vector<Integer> roiI,
		final Vector<Integer> roiJ, final int width, final int height,
		final double[] scaledImage, final double threshold)
	{

		final int[][] fourconnectedNHood = { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0,
			1 } };

		// Fill the area enclosed by the traced edge contained in roiI,roiJ
		// beginning needs to be within the traced edge
		byte[] sieveTemp = new byte[width * height];
		int x;
		int y;
		for (int z = 0; z < roiI.size(); ++z) {
			sieveTemp[roiI.get(z) + roiJ.get(z) * width] = 1;
		}

		// Determine the flood fill init
		int[] tempCoordinates;
		while (true) {
			tempCoordinates = findFillInit(sieveTemp, roiI, roiJ, scaledImage,
				threshold);
			if (tempCoordinates == null) {
				return sieveTemp;
			}
			x = tempCoordinates[0];
			y = tempCoordinates[1];

			final Vector<Integer> initialX = new Vector<>();
			final Vector<Integer> initialY = new Vector<>();
			initialX.add(x);
			initialY.add(y);
			sieveTemp[x + y * width] = 1;
			final byte[] sieveTemp2 = sieveTemp.clone();
			boolean noLeak = true;
			while (!initialX.isEmpty()) {

				x = initialX.remove(initialX.size() - 1);
				y = initialY.remove(initialY.size() - 1);

				final int index = x + y * width;
				if (sieveTemp2[index] == 0) {
					sieveTemp2[index] = 1;
				}
				if (x < 1 || x >= width - 1 || y < 1 || y >= height - 1) {
					noLeak = false;
					break;
				}
				// Check 4-connected neighbours
				for (final int[] aFourconnectedNHood : fourconnectedNHood) {
					if (sieveTemp2[x + aFourconnectedNHood[0] + (y +
						aFourconnectedNHood[1]) * width] == 0)
					{
						initialX.add(x + aFourconnectedNHood[0]);
						initialY.add(y + aFourconnectedNHood[1]);
					}

				}

			}
			if (noLeak) {
				sieveTemp = sieveTemp2.clone();
			}
		}
	}

	// DetectedEdge
	public static int[] twoLargestBonesDetectedEdges(
		final List<DetectedEdge> edges)
	{
		// Identify the two longest circumferences
		final List<DetectedEdge> temp3 = new ArrayList<>(edges);
		Collections.sort(temp3);
		int counter = 0;
		final int[] twoLongest = new int[2];
		while (edges.get(counter).area != temp3.get(temp3.size() - 1).area) {
			++counter;
		}
		twoLongest[0] = counter;
		counter = 0;
		if (temp3.size() > 1) {
			while (edges.get(counter).area != temp3.get(temp3.size() - 2).area ||
				counter == twoLongest[0])
			{
				++counter;
			}
		}
		twoLongest[1] = counter;
		return twoLongest;
	}

	// DetectedEdge
	private double[] calcDistancesFromCentreOfLimb(final List<DetectedEdge> edges,
		final double[] tempScaledImage, final double fatThreshold)
	{
		final List<double[]> bones = new ArrayList<>(edges.size());
		for (int i = 0; i < edges.size(); ++i) {
			bones.add(new double[3]);
			for (int j = 0; j < 3; ++j) {
				bones.get(i)[j] = 0;
			}
		}
		// Find the centre of area of the limb
		final int maxIndice = selectRoiBiggestBoneDetectedEdges(edges);
		final byte[] limbSieve = new byte[tempScaledImage.length];
		limbSieve[edges.get(maxIndice).iit.get(0) + edges.get(maxIndice).jiit.get(
			0) * width] = 1;
		// Dilate muscleSieve, into neighbouring fat pixels
		int tempDil = 1;
		while (tempDil > 0) {
			tempDil = dilateLimb(limbSieve, (byte) 1, (byte) 0, (byte) 4,
				fatThreshold, tempScaledImage);
		}
		double limbCenterX = 0.0;
		double limbCenterY = 0.0;
		double limbPoints = 0.0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i) {
				if (limbSieve[i + j * width] == (byte) 1) {
					limbCenterX += i;
					limbCenterY += j;
					limbPoints += 1;
				}
			}
		}
		limbCenterX /= limbPoints;
		limbCenterY /= limbPoints;
		// Find the centres of circumference of the bones
		final double[] distanceFromCentreOfLimb = new double[edges.size()];
		for (int i = 0; i < edges.size(); ++i) {
			for (int j = 0; j < edges.get(i).iit.size(); j++) {
				bones.get(i)[0] += edges.get(i).iit.get(j);
				bones.get(i)[1] += edges.get(i).jiit.get(j);
				bones.get(i)[2] += 1;
			}
			bones.get(i)[0] /= bones.get(i)[2];
			bones.get(i)[1] /= bones.get(i)[2];
			// Square root omitted, as it does not affect the order...
			distanceFromCentreOfLimb[i] = Math.pow(limbCenterX - bones.get(i)[0],
				2.0) + Math.pow(limbCenterY - bones.get(i)[1], 2.0);
		}
		return distanceFromCentreOfLimb;
	}

	// Remove the extra part from vectors and replace with a straight line
	private Vector<Vector<Integer>> cleave(final byte[] result,
		final Vector<Integer> fatRoiI, final Vector<Integer> fatRoiJ,
		final int[] cleavingIndices)
	{
		final int initI = fatRoiI.get(cleavingIndices[0]);
		final int initJ = fatRoiJ.get(cleavingIndices[0]);
		final int targetI = fatRoiI.get(cleavingIndices[1]);
		final int targetJ = fatRoiJ.get(cleavingIndices[1]);
		// remove cleaved elements
		int replacementI = fatRoiI.get(cleavingIndices[0]);
		int replacementJ = fatRoiJ.get(cleavingIndices[0]);
		// the elements to be cleaved
		final List<Integer> cleavedI = new Vector<>(fatRoiI.subList(
			cleavingIndices[0] + 1, cleavingIndices[1] + 1));
		final List<Integer> cleavedJ = new Vector<>(fatRoiJ.subList(
			cleavingIndices[0] + 1, cleavingIndices[1] + 1));
		for (int i = cleavingIndices[0]; i < cleavingIndices[1]; ++i) {
			// Remove the elements to be cleaved
			fatRoiI.removeElementAt(cleavingIndices[0]);
			fatRoiJ.removeElementAt(cleavingIndices[0]);
		}
		// Insert replacement line
		final double replacementLength = cleavingIndices[1] - cleavingIndices[0];
		final double repILength = targetI - initI;
		final double repJLength = targetJ - initJ;
		double relativeLength;
		final Vector<Integer> insertionI = new Vector<>();
		final Vector<Integer> insertionJ = new Vector<>();
		insertionI.add(replacementI);
		insertionJ.add(replacementJ);
		for (int k = cleavingIndices[0]; k < cleavingIndices[1]; ++k) {
			relativeLength = k - cleavingIndices[0];
			replacementI = ((int) (repILength * (relativeLength /
				replacementLength))) + initI;
			replacementJ = ((int) (repJLength * (relativeLength /
				replacementLength))) + initJ;
			if (replacementI != insertionI.lastElement() || replacementJ != insertionJ
				.lastElement())
			{
				insertionI.add(replacementI);
				insertionJ.add(replacementJ);
				result[replacementI + replacementJ * width] = 1;
			}
		}
		fatRoiI.addAll(cleavingIndices[0], insertionI);
		fatRoiJ.addAll(cleavingIndices[0], insertionJ);
		Collections.reverse(insertionI);
		Collections.reverse(insertionJ);
		cleavedI.addAll(0, insertionI);
		cleavedJ.addAll(0, insertionJ);
		final Vector<Vector<Integer>> returnVectorPair = new Vector<>();
		returnVectorPair.add(new Vector<>());
		returnVectorPair.add(new Vector<>());
		returnVectorPair.get(0).addAll(cleavedI);
		returnVectorPair.get(1).addAll(cleavedJ);
		return returnVectorPair;
	}

	/*Cleaving is made by looking at the ratios of
	distances between two points along the edge and the shortest distance
	between the points. If the maximum of the  ratio is big enough, the
	highest ratio points will be connected with a straigth
	line and the edge with higher indices will be removed. E.g.
	for a circle, the maximum ratio is (pi/2)/d ~= 1.57 and for square
	it is 2/sqrt(2) = sqrt(2) ~= 1.41.*/
	private Vector<Vector<Vector<Integer>>> cleaveEdge(final byte[] result,
		final Vector<Integer> fatRoiI, final Vector<Integer> fatRoiJ,
		final double minRatio, final double minLength)
	{
		double distanceAlongTheEdge;
		double distance;
		double ratio;
		final double minEdge = fatRoiI.size() / minLength;
		final int[] cleavingIndices = new int[2];
		final Vector<Vector<Vector<Integer>>> returnVectorVectorPointer =
			new Vector<>();
		while (true) {
			double highestRatio = minRatio - 0.1;
			/*Go through all point pairs*/
			for (int i = 0; i < fatRoiI.size() - 11; ++i) {
				for (int j = i + 10; j < fatRoiI.size(); ++j) {
					distance = Math.sqrt(Math.pow(fatRoiI.get(j) - fatRoiI.get(i), 2.0) +
						Math.pow(fatRoiJ.get(j) - fatRoiJ.get(i), 2.0));
					distanceAlongTheEdge = Math.min((j - i), fatRoiI.size() - j + i);
					ratio = distanceAlongTheEdge / distance;
					if (ratio > highestRatio && distanceAlongTheEdge > minEdge) {
						highestRatio = ratio;
						cleavingIndices[0] = i;
						cleavingIndices[1] = j;
					}

				}
			}
			/*If ratio is high enough, cleave at the highest ratio point pair*/
			if (highestRatio < minRatio) {
				break;
			}
			returnVectorVectorPointer.add(cleave(result, fatRoiI, fatRoiJ,
				cleavingIndices));
		}
		/*Insert the last retained part to first index.*/
		final Vector<Vector<Integer>> returnVectorPair = new Vector<>();
		returnVectorPair.add(new Vector<>());
		returnVectorPair.add(new Vector<>());
		returnVectorPair.get(0).addAll(fatRoiI);
		returnVectorPair.get(1).addAll(fatRoiJ);
		if (returnVectorVectorPointer.size() < 1) {
			returnVectorVectorPointer.add(returnVectorPair);
		}
		else {
			returnVectorVectorPointer.insertElementAt(returnVectorPair, 0);
		}
		return returnVectorVectorPointer;
	}

	private int dilateLimb(final byte[] data, final byte dilateVal,
		final byte min, final byte temp, final double threshold,
		final double[] scaledImage)
	{
		// Dilate algorithm
		// Best dilate by one solution taken from
		// http://ostermiller.org/dilate_and_erode.html
		int dilated = 0;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				if (data[i * width + j] == dilateVal) {
					if (i > 0 && data[(i - 1) * width + j] == min && scaledImage[(i - 1) *
						width + j] >= threshold)
					{
						data[(i - 1) * width + j] = temp;
					}
					if (j > 0 && data[(i) * width + j - 1] == min && scaledImage[(i) *
						width + j - 1] >= threshold)
					{
						data[(i) * width + j - 1] = temp;
					}
					if (i + 1 < height && data[(i + 1) * width + j] == min &&
						scaledImage[(i + 1) * width + j] >= threshold)
					{
						data[(i + 1) * width + j] = temp;
					}
					if (j + 1 < width && data[(i) * width + j + 1] == min &&
						scaledImage[(i) * width + j + 1] >= threshold)
					{
						data[(i) * width + j + 1] = temp;
					}
				}
			}
		}
		for (int i = 0; i < width * height; i++) {
			if (data[i] == temp) {
				// Set to proper value here...
				data[i] = dilateVal;
				++dilated;
			}
		}
		return dilated;
	}

	// DetectedEdge version
	private Vector<Object> fillResultEdge(byte[] result,
		final Vector<Integer> iit, final Vector<Integer> jiit,
		final double[] scaledImage, final double threshold)
	{
		if (iit.isEmpty()) {
			return null;
		}
		Vector<Object> results = null;
		int pixelsFilled = 0;
		// Set initial fill pixel to the first pixel above threshold not on the
		// border
		boolean possible = true;
		final byte[] tempResult = result.clone();
		int[] tempCoordinates = findFillInit(tempResult, iit, jiit, scaledImage,
			threshold);
		while (possible && tempCoordinates != null) {
			final Vector<Object> returned = resultFill(tempCoordinates[0],
				tempCoordinates[1], tempResult);
			possible = (Boolean) returned.get(0);
			pixelsFilled += (Integer) returned.get(1);
			tempCoordinates = findFillInit(tempResult, iit, jiit, scaledImage,
				threshold);
		}
		if (possible) {
			results = new Vector<>();
			result = tempResult;
			results.add(result);
			results.add(iit);
			results.add(jiit);
			results.add(pixelsFilled);
		}
		return results;
	}

	// DetectEdge
	private Vector<Object> findEdge(final double[] scaledImage,
		final double threshold, final boolean allowCleaving)
	{
		int i = 0;
		int j = 0;
		int tempI;
		int tempJ;
		byte[] result = new byte[scaledImage.length];
		final Collection<DetectedEdge> edges = new Vector<>();
		while ((i < (width - 1)) && (j < (height - 1))) {
			while (j < height - 1 && i < width && scaledImage[i + j *
				width] < threshold)
			{
				i++;
				if (result[i + j * width] == 1) {
					while (j < height - 1 && result[i + j * width] > 0) {
						i++;
						if (i == width && j < height - 2) {
							i = 0;
							j++;
						}

					}
				}
				if (i == width) {
					j++;
					if (j >= height - 1) break;
					i = 0;
				}
			}
			tempI = i;
			tempJ = j;

			if (i >= width - 1 && j >= height - 1) {
				break; /*Go to end...*/
			}
			result[i + j * width] = 1;

			// Tracing algorithm DetectedEdge
			final Vector<Object> returned = traceEdge(scaledImage, result, threshold,
				i, j);
			result = (byte[]) returned.get(0);
			final Vector<Integer> newIit = (Vector<Integer>) returned.get(1);
			final Vector<Integer> newJiit = (Vector<Integer>) returned.get(2);
			// Tracing algorithm done...

			if (allowCleaving) {
				final Vector<Vector<Vector<Integer>>> returnedVectors = cleaveEdge(
					result, newIit, newJiit, 3.0, 6.0);
				for (final Vector<Vector<Integer>> returnedVector : returnedVectors) {
					// Fill edge within result..
					final Vector<Integer> iit = new Vector<>();
					final Vector<Integer> jiit = new Vector<>();
					for (int ii = 0; ii < returnedVector.get(0).size(); ++ii) {
						iit.add(returnedVector.get(0).get(ii));
						jiit.add(returnedVector.get(1).get(ii));
					}
					final Vector<Object> results = fillResultEdge(result, iit, jiit,
						scaledImage, threshold);
					if (results != null) {
						result = (byte[]) results.get(0);
						edges.add(new DetectedEdge((Vector<Integer>) results.get(1),
							(Vector<Integer>) results.get(2), (Integer) results.get(3)));
					}
				}
			}
			else {
				// Fill edge within result..
				final Vector<Integer> iit = new Vector<>();
				final Vector<Integer> jiit = new Vector<>();
				for (int ii = 0; ii < newIit.size(); ++ii) {
					iit.add(newIit.get(ii));
					jiit.add(newJiit.get(ii));
				}
				final Vector<Object> results = fillResultEdge(result, iit, jiit,
					scaledImage, threshold);
				if (results != null) {
					result = (byte[]) results.get(0);
					edges.add(new DetectedEdge((Vector<Integer>) results.get(1),
						(Vector<Integer>) results.get(2), (Integer) results.get(3)));
				}
			}
			// Find next empty spot
			i = tempI;
			j = tempJ;
			while (j < height && scaledImage[i + j * width] >= threshold) {
				i++;
				if (i == width) {
					i = 0;
					j++;
				}
			}
		}

		final Vector<Object> returnVector = new Vector<>();
		returnVector.add(result);
		returnVector.add(edges);
		return returnVector;
	}

	// DetectedEdge. Find fill init by steering clockwise from next to previous
	private int[] findFillInit(final byte[] result, final Vector<Integer> iit,
		final Vector<Integer> jiit, final double[] scaledImage,
		final double threshold)
	{
		final int[] returnCoordinates = new int[2];
		final int[] steer = new int[2];
		for (int j = 0; j < iit.size() - 1; ++j) {
			returnCoordinates[0] = iit.get(j);
			returnCoordinates[1] = jiit.get(j);
			double direction = Math.atan2(jiit.get(j + 1) - returnCoordinates[1], iit
				.get(j + 1) - returnCoordinates[0]);
			for (int i = 0; i < 8; ++i) {
				direction += Math.PI / 4.0;
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

				if (result[returnCoordinates[0] + steer[0] + (returnCoordinates[1] +
					steer[1]) * width] == 0 && scaledImage[returnCoordinates[0] +
						steer[0] + (returnCoordinates[1] + steer[1]) * width] >= threshold)
				{
					returnCoordinates[0] += steer[0];
					returnCoordinates[1] += steer[1];
					return returnCoordinates;
				}
				if (result[returnCoordinates[0] + steer[0] + (returnCoordinates[1] +
					steer[1]) * width] == 1)
				{
					break;
				}
			}
		}
		return null;
	}

	// DetectedEdge
	private static boolean guessFlipLarger(final List<DetectedEdge> edges,
		final boolean stacked)
	{
		final List<Integer> temp = new Vector<>();
		final List<Integer> temp2 = new Vector<>();
		for (final DetectedEdge edge : edges) {
			temp.add(edge.area);
			temp2.add(edge.area);
		}
		Collections.sort(temp);

		final int[] counter = { 0, 0 };

		while (!temp2.get(counter[0]).equals(temp.get(temp.size() - 1))) {
			++counter[0];
		}
		boolean returnValue = false;
		if (temp.size() > 1) {
			while (!temp2.get(counter[1]).equals(temp.get(temp.size() - 2))) {
				++counter[1];
			}
			if (stacked) {
				returnValue = edges.get(counter[0]).jiit.get(0) >= edges.get(
					counter[1]).jiit.get(0);
			}
			else {
				returnValue = edges.get(counter[0]).iit.get(0) >= edges.get(
					counter[1]).iit.get(0);
			}
		}
		return returnValue;
	}

	// DetectedEdge Only two biggest bone will be considered..
	private static boolean guessFlipSelection(final List<DetectedEdge> edges,
		final int selection, final boolean stacked)
	{

		final int[] considered = twoLargestBonesDetectedEdges(edges);
		if (selection != considered[0] && selection != considered[1]) {
			// selection is not the biggest or the second biggest bone -> can't make a
			// guess, return false
			return false;
		}

		final int[] possibleCoords;
		final int selectionCoord;
		final DetectedEdge first = edges.get(considered[0]);
		final DetectedEdge second = edges.get(considered[1]);
		if (stacked) {
			selectionCoord = edges.get(selection).jiit.get(0);
			possibleCoords = new int[] { first.jiit.get(0), second.jiit.get(0) };
		}
		else {
			selectionCoord = edges.get(selection).iit.get(0);
			possibleCoords = new int[] { first.iit.get(0), second.iit.get(0) };
		}
		return (selection == considered[0] && selectionCoord > possibleCoords[1]) ||
			(selection == considered[1] && selectionCoord > possibleCoords[0]);
	}

	private Vector<Object> resultFill(int i, int j, final byte[] tempResult) {
		final Vector<Integer> initialI = new Vector<>();
		final Vector<Integer> initialJ = new Vector<>();
		initialI.add(i);
		initialJ.add(j);
		int pixelsFilled = 0;
		while (!initialI.isEmpty() && initialI.lastElement() > 0 && initialI
			.lastElement() < width - 1 && initialJ.lastElement() > 0 && initialJ
				.lastElement() < height - 1)
		{
			i = initialI.lastElement();
			j = initialJ.lastElement();
			initialI.remove(initialI.size() - 1);
			initialJ.remove(initialJ.size() - 1);

			if (tempResult[i + j * width] == 0) {
				tempResult[i + j * width] = 1;
				++pixelsFilled;
			}

			if (tempResult[i - 1 + j * width] == 0) {
				initialI.add(i - 1);
				initialJ.add(j);
			}

			if (tempResult[i + 1 + j * width] == 0) {
				initialI.add(i + 1);
				initialJ.add(j);
			}

			if (tempResult[i + (j - 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j - 1);
			}

			if (tempResult[i + (j + 1) * width] == 0) {
				initialI.add(i);
				initialJ.add(j + 1);
			}

		}
		final Vector<Object> returnValue = new Vector<>();
		returnValue.add(initialI.isEmpty() && initialJ.isEmpty());
		returnValue.add(pixelsFilled);
		return returnValue;
	}

	// DetectedEdge
	private static int selectRoiBiggestBoneDetectedEdges(
		final List<DetectedEdge> edges)
	{
		int maxArea = 0;
		int maxPos = -1;
		for (int i = 0; i < edges.size(); i++) {
			final int area = edges.get(i).area;
			if (area > maxArea) {
				maxArea = area;
				maxPos = i;
			}
		}
		return maxPos;
	}

	// DetectedEdge
	private static int selectRoiBottomBone(final Collection<DetectedEdge> edges) {
		return selectRoiFirstNthFromTop(edges, edges.size() - 1);
	}

	// DetectedEdge
	private int selectRoiCentralBone(final List<DetectedEdge> edges,
		final double[] tempScaledImage, final double fatThreshold)
	{
		final double[] distanceFromCentreOfLimb = calcDistancesFromCentreOfLimb(
			edges, tempScaledImage, fatThreshold);
		final double[] temp = Arrays.copyOf(distanceFromCentreOfLimb,
			distanceFromCentreOfLimb.length);
		Arrays.sort(temp);
		int counter = 0;
		while (Double.compare(distanceFromCentreOfLimb[counter], temp[0]) != 0) {
			++counter;
		}
		return counter;
	}

	// DetectedEdge , indexing from 0
	private static int selectRoiFirstNthFromLeft(
		final Collection<DetectedEdge> edges, final int nth)
	{
		final List<Integer> indices = edges.stream().map(e -> e.iit.get(0)).collect(
			toList());
		final Integer nthLeast = indices.stream().sorted().skip(nth).findFirst()
			.orElse(-1);
		return indices.indexOf(nthLeast);
	}

	// DetectedEdge, indexing from 0
	private static int selectRoiFirstNthFromTop(
		final Collection<DetectedEdge> edges, final int nth)
	{
		final List<Integer> indices = edges.stream().map(e -> e.jiit.get(0))
			.collect(toList());
		final Integer nthLeast = indices.stream().sorted().skip(nth).findFirst()
			.orElse(-1);
		return indices.indexOf(nthLeast);
	}

	// DetectedEdge
	private static int selectRoiLeftmostBone(
		final Collection<DetectedEdge> edges)
	{
		return selectRoiFirstNthFromLeft(edges, 0);
	}

	// DetectedEdge
	private int selectRoiPeripheralBone(final List<DetectedEdge> edges,
		final double[] tempScaledImage, final double fatThreshold)
	{
		final double[] distanceFromCentreOfLimb = calcDistancesFromCentreOfLimb(
			edges, tempScaledImage, fatThreshold);
		final double[] temp = Arrays.copyOf(distanceFromCentreOfLimb,
			distanceFromCentreOfLimb.length);
		Arrays.sort(temp);
		int counter = 0;
		while (distanceFromCentreOfLimb[counter] != temp[temp.length - 1]) {
			++counter;
		}
		return counter;
	}

	// DetectedEdge
	private static int selectRoiRightmostBone(
		final Collection<DetectedEdge> edges)
	{
		return selectRoiFirstNthFromLeft(edges, edges.size() - 1);
	}

	// DetectedEdge
	private static int selectRoiSecondLargestBoneDetectedEdges(
		final List<DetectedEdge> edges)
	{
		if (edges.size() < 2) {
			return -1;
		}
		final List<DetectedEdge> sortedEdges = edges.stream().sorted().collect(
			toList());
		return edges.indexOf(sortedEdges.get(edges.size() - 2));
	}

	// DetectedEdge
	private static int selectRoiSmallestBoneDetectedEdges(
		final List<DetectedEdge> edges)
	{
		final DetectedEdge leastAreaEdge = edges.stream().sorted().findFirst()
			.orElse(null);
		return edges.indexOf(leastAreaEdge);
	}

	// DetectedEdge
	private static int selectRoiTopBone(final Collection<DetectedEdge> edges) {
		return selectRoiFirstNthFromTop(edges, 0);
	}

	// DetectedEdge
	private static int selectRoiTwoLargestLeft(final List<DetectedEdge> edges) {
		if (edges.size() < 2) {
			// In case only one ROI has been found..
			return -1;
		}
		final int[] twoBones = twoLargestBonesRetainOrderDetectedEdges(edges);
		final Collection<DetectedEdge> tempEdges = new Vector<>();
		for (final int twoBone : twoBones) {
			tempEdges.add(edges.get(twoBone));
		}
		final int tempSelection = selectRoiLeftmostBone(tempEdges);
		return twoBones[tempSelection];
	}

	// DetectedEdge
	private static int selectRoiTwoLargestRight(final List<DetectedEdge> edges) {
		if (edges.size() < 2) {
			return -1;
		} // In case only one ROI has been found..
		final int[] twoBones = twoLargestBonesRetainOrderDetectedEdges(edges);
		final Collection<DetectedEdge> tempEdges = new Vector<>();
		for (final int twoBone : twoBones) {
			tempEdges.add(edges.get(twoBone));
		}
		final int tempSelection = selectRoiRightmostBone(tempEdges);
		return twoBones[tempSelection];
	}

	/*	Edge Tracing DetectedEdge
	trace edge by advancing according to the previous direction
	if above threshold, turn to negative direction
	if below threshold, turn to positive direction
	Idea taken from http://www.math.ucla.edu/~bertozzi/RTG/zhong07/report_zhong.pdf
	The paper traced continent edges on map/satellite image
	*/
	private Vector<Object> traceEdge(final double[] scaledImage,
		final byte[] result, final double threshold, int i, int j)
	{
		final Collection<Integer> iit = new Vector<>();
		final Collection<Integer> jiit = new Vector<>();
		iit.add(i);
		jiit.add(j);
		// begin by advancing right. Positive angles rotate the direction clockwise.
		double direction = 0;
		double previousDirection;
		final int initI;
		final int initJ;
		initI = i;
		initJ = j;
		while (true) {
			int counter = 0;
			previousDirection = direction;
			// Handle going out of bounds by considering out of bounds to be less than
			// threshold
			if ((i + ((int) Math.round(Math.cos(direction)))) >= 0 && (i + ((int) Math
				.round(Math.cos(direction))) < width) && (j + ((int) Math.round(Math
					.sin(direction))) >= 0) && (j + ((int) Math.round(Math.sin(
						direction))) < height) && scaledImage[i + ((int) Math.round(Math
							.cos(direction))) + (j + ((int) Math.round(Math.sin(
								direction)))) * width] > threshold)
			{
				// Rotate counter clockwise
				while (counter < 8 && i + ((int) Math.round(Math.cos(direction -
					Math.PI / 4.0))) >= 0 && i + ((int) Math.round(Math.cos(direction -
						Math.PI / 4.0))) < width && j + ((int) Math.round(Math.sin(
							direction - Math.PI / 4.0))) >= 0 && j + ((int) Math.round(Math
								.sin(direction - Math.PI / 4.0))) < height && scaledImage[i +
									((int) Math.round(Math.cos(direction - Math.PI / 4.0))) + (j +
										((int) Math.round(Math.sin(direction - Math.PI / 4.0)))) *
										width] > threshold)
				{
					direction -= Math.PI / 4.0;
					++counter;
					if (Math.abs(direction - previousDirection) >= 180) {
						break;
					}
				}
			}
			else {
				// Rotate clockwise
				while (counter < 8 && (i + ((int) Math.round(Math.cos(
					direction))) < 0 || i + ((int) Math.round(Math.cos(
						direction))) >= width || j + ((int) Math.round(Math.sin(
							direction))) < 0 || j + ((int) Math.round(Math.sin(
								direction))) >= height || scaledImage[i + ((int) Math.round(Math
									.cos(direction))) + (j + ((int) Math.round(Math.sin(
										direction)))) * width] < threshold))
				{
					direction += Math.PI / 4.0;
					++counter;
					if (Math.abs(direction - previousDirection) >= 180) {
						break;
					}
				}

			}
			i += (int) Math.round(Math.cos(direction));
			j += (int) Math.round(Math.sin(direction));
			if ((i == initI && j == initJ) || counter > 7 || scaledImage[i + j *
				width] < threshold || result[i + j * width] == 1 || result[i + j *
					width] > 3)
			{
				for (int ii = 0; ii < result.length; ++ii) {
					if (result[ii] > 1) {
						result[ii] = 1;
					}
				}
				final Vector<Object> returnVector = new Vector<>();
				returnVector.add(result);
				returnVector.add(iit);
				returnVector.add(jiit);
				return returnVector;
			}
			else {
				if (result[i + j * width] == 0) {
					result[i + j * width] = 2;
				}
				else if (result[i + j * width] != 1) {
					result[i + j * width]++;
				}
				iit.add(i);
				jiit.add(j);

			}
			// Keep steering counter clockwise not to miss single pixel structs...
			direction -= Math.PI / 2.0;
		}
	}

	// DetectedEdge
	private static int[] twoLargestBonesRetainOrderDetectedEdges(
		final List<DetectedEdge> edges)
	{
		// Identify the two longest circumferences
		final List<DetectedEdge> temp3 = new Vector<>(edges);
		Collections.sort(temp3);
		int counter = 0;
		final int[] twoLongest = new int[2];
		while (edges.get(counter).area != temp3.get(temp3.size() - 1).area) {
			++counter;
		}
		twoLongest[0] = counter;
		counter = 0;
		if (temp3.size() > 1) {
			while (edges.get(counter).area != temp3.get(temp3.size() - 2).area ||
				counter == twoLongest[0])
			{
				++counter;
			}
		}
		twoLongest[1] = counter;
		Arrays.sort(twoLongest);
		return twoLongest;
	}

	byte[] dilate(final byte[] data, final byte dilateVal, final byte min,
		final byte temp)
	{
		// Dilate algorithm
		// Best dilate by one solution taken from
		// http://ostermiller.org/dilate_and_erode.html
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				if (data[i * width + j] == dilateVal) {
					if (i > 0 && data[(i - 1) * width + j] == min) {
						data[(i - 1) * width + j] = temp;
					}
					if (j > 0 && data[(i) * width + j - 1] == min) {
						data[i * width + j - 1] = temp;
					}
					if (i + 1 < height && data[(i + 1) * width + j] == min) {
						data[(i + 1) * width + j] = temp;
					}
					if (j + 1 < width && data[i * width + j + 1] == min) {
						data[i * width + j + 1] = temp;
					}
				}
			}
		}
		for (int i = 0; i < data.length; i++) {
			if (data[i] == temp) {
				// Set to proper value here...
				data[i] = dilateVal;
			}
		}
		return data;
	}

	/*A function to get rid of the measurement tube used at UKK-institute
	with Stratex XCT3000 device. Needed for soft tissue analysis*/
	byte[] removeSleeve(final double[] scaledImage,
		final double sleeveThreshold)
	{
		int i = 10;
		int j = 10;
		while ((j < height - 12 && i < width - 11 && scaledImage[i + j *
			width] < sleeveThreshold) || scaledImage[i + j * width] == 0)
		{
			i++;
			if (i == width - 11) {
				j++;
				if (j >= height - 12) break;
				i = 10;
			}
		}
		// Sleeve found
		byte[] sleeve = new byte[width * height];
		final Vector<Integer> initialI = new Vector<>();
		final Vector<Integer> initialJ = new Vector<>();
		initialI.add(i);
		initialJ.add(j);
		while (!initialI.isEmpty() && initialI.lastElement() > 0 && initialI
			.lastElement() < width - 1 && initialJ.lastElement() > 0 && initialJ
				.lastElement() < height - 1)
		{
			i = initialI.lastElement();
			j = initialJ.lastElement();
			initialI.remove(initialI.size() - 1);
			initialJ.remove(initialJ.size() - 1);
			if (scaledImage[i + j * width] > sleeveThreshold && sleeve[i + j *
				width] == 0)
			{
				sleeve[i + j * width] = 1;
			}
			if (scaledImage[i - 1 + j * width] > sleeveThreshold && sleeve[i - 1 + j *
				width] == 0)
			{
				initialI.add(i - 1);
				initialJ.add(j);
			}
			if (scaledImage[i + 1 + j * width] > sleeveThreshold && sleeve[i + 1 + j *
				width] == 0)
			{
				initialI.add(i + 1);
				initialJ.add(j);
			}
			if (scaledImage[i + (j - 1) * width] > sleeveThreshold && sleeve[i + (j -
				1) * width] == 0)
			{
				initialI.add(i);
				initialJ.add(j - 1);
			}
			if (scaledImage[i + (j + 1) * width] > sleeveThreshold && sleeve[i + (j +
				1) * width] == 0)
			{
				initialI.add(i);
				initialJ.add(j + 1);
			}

		}
		sleeve = dilate(sleeve, (byte) 1, (byte) 0, (byte) 2);
		sleeve = dilate(sleeve, (byte) 1, (byte) 0, (byte) 2);
		return sleeve;
	}

	byte[] erode(final byte[] data) {
		// Erode algorithm
		// Modified from the best dilate by one solution taken from
		// http://ostermiller.org/dilate_and_erode.html
		for (int i = 1; i < height - 1; i++) {
			for (int j = 1; j < width - 1; j++) {
				if (data[i * width + j] > 0) {
					if (data[(i - 1) * width + j] == 0 | data[(i) * width + j - 1] == 0 |
						data[(i + 1) * width + j] == 0 | data[(i) * width + j + 1] == 0)
					{
						data[i * width + j] = -1;
					} // Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i = 0; i < data.length; i++) {
			data[i] = (byte) Math.max(0, data[i]);
		}
		return data;
	}

	// DetectedEdges
	Vector<Object> getSieve(final double[] tempScaledImage,
		final double boneThreshold, final String roiChoice,
		final boolean guessStacked, final boolean stacked, final boolean guessFlip,
		final boolean allowCleaving) throws ExecutionException
	{
		// Trace bone edges
		final Vector<?> results = findEdge(tempScaledImage, boneThreshold,
			allowCleaving);
		result = (byte[]) results.get(0);
		@SuppressWarnings("unchecked")
		final List<DetectedEdge> edges = (Vector<DetectedEdge>) results.get(1);
		if (edges.size() < 1) {
			final double minValue = Arrays.stream(tempScaledImage).min().orElse(0);
			final double maxValue = Arrays.stream(tempScaledImage).max().orElse(0);
			throw new ExecutionException(
				"Couldn't find a bone. The range of intensities in the file is " +
					minValue + " to " + maxValue +
					". Set the Thresholds between these values", new Throwable());
		}

		// Select correct bone outline
		int selection = 0;
		final int choiceIndex = Arrays.asList(details.choiceLabels).indexOf(
			roiChoice);
		switch (choiceIndex) {
			case 0:
				selection = selectRoiBiggestBoneDetectedEdges(edges);
				break;
			case 1:
				selection = selectRoiSmallestBoneDetectedEdges(edges);
				break;
			case 2:
				selection = selectRoiLeftmostBone(edges);
				break;
			case 3:
				selection = selectRoiRightmostBone(edges);
				break;
			case 4:
				selection = selectRoiTopBone(edges);
				break;
			case 5:
				selection = selectRoiBottomBone(edges);
				break;
			case 6:
				selection = selectRoiCentralBone(edges, tempScaledImage,
					details.fatThreshold);
				break;
			case 7:
				selection = selectRoiPeripheralBone(edges, tempScaledImage,
					details.fatThreshold);
				break;
			case 8:
				selection = selectRoiSecondLargestBoneDetectedEdges(edges);
				break;
			case 9:
				selection = selectRoiTwoLargestLeft(edges);
				break;
			case 10:
				selection = selectRoiTwoLargestRight(edges);
				break;
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
				selection = selectRoiFirstNthFromLeft(edges, choiceIndex - 11);
				break;
			case 16:
			case 17:
			case 18:
			case 19:
			case 20:
				selection = selectRoiFirstNthFromTop(edges, choiceIndex - 16);
				break;
		}

		if (guessStacked) {
			final int[] guessingStack = twoLargestBonesDetectedEdges(edges);
			final DetectedEdge edge = edges.get(guessingStack[0]);
			final DetectedEdge edge2 = edges.get(guessingStack[1]);
			final double stackedThreshold = 1.1 * Math.abs(edge.iit.get(0) - edge2.iit
				.get(1));
			details.stacked = Math.abs(edge.jiit.get(0) - edge2.jiit.get(
				1)) > stackedThreshold;
		}

		// Try to guess whether to flip the distribution
		if (guessFlip) {
			if (details.guessLarger) {
				details.flipDistribution = guessFlipLarger(edges, stacked);
			}
			else {
				details.flipDistribution = guessFlipSelection(edges, selection,
					stacked);
			}
			// Flip flip, if roiChoice is smaller or second largest
			if (details.invertGuess) {
				details.flipDistribution = !details.flipDistribution;
			}
		}

		final byte[] tempSieve = fillSieve(edges.get(selection).iit, edges.get(
			selection).jiit, width, height, tempScaledImage, boneThreshold);
		final Vector<Object> returnVector = new Vector<>();
		returnVector.add(tempSieve);
		returnVector.add(result);
		returnVector.add(edges);
		returnVector.add(selection);
		return returnVector;
	}
}
