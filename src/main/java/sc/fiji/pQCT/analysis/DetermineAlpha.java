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

package sc.fiji.pQCT.analysis;

import static java.util.stream.Collectors.toCollection;
import static java.util.stream.IntStream.range;

import java.util.Arrays;
import java.util.Collections;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.selectroi.DetectedEdge;
import sc.fiji.pQCT.selectroi.RoiSelector;
import sc.fiji.pQCT.selectroi.SelectROI;

public class DetermineAlpha {

	private final ImageAndAnalysisDetails details;
	public Vector<Integer> pindColor;
	public double alpha;
	public double rotationCorrection;
	public double distanceBetweenBones;
	Vector<Integer> pind;

	public DetermineAlpha(final SelectROI roi,
		final ImageAndAnalysisDetails details)
	{
		this.details = details;
		// Calculate CSMIs and rotation angle to align maximal and minimal bending
		// axes with X and Y axes
		rotationCorrection = details.sectorWidth / 2.0;
		// Rotation according to Imax/Imin for bone of interest or according to all
		// bones
		final String choice = details.rotationChoice;
		final String[] labels = details.rotationLabels;
		if (choice.equals(labels[0])) {
			final double[] csmiValues = csmi(roi.sieve, roi.width, roi.height);
			determineMomentAlpha(csmiValues);
		}
		if (choice.equals(labels[2])) {
			final int voxels = roi.width * roi.height;
			final byte[] tempCsmiSieve = new byte[voxels];
			final double[] image = roi.scaledImage;
			range(0, voxels).filter(i -> image[i] >= details.rotationThreshold)
				.forEach(i -> tempCsmiSieve[i] = 1);
			final double[] csmiValues = csmi(tempCsmiSieve, roi.width, roi.height);
			determineMomentAlpha(csmiValues);
		}
		/*Rotation according to the furthest point*/
		if (choice.equals(labels[1])) {
			/*Calculate alpha from periosteal radii*/
			final double[] marrowCenter = new double[2];
			final IntStream iValues;
			final IntStream jValues;
			if (roi.boneMarrowRoiI.isEmpty()) {
				iValues = roi.cortexAreaRoiI.stream().mapToInt(i -> i);
				jValues = roi.cortexAreaRoiJ.stream().mapToInt(j -> j);
			}
			else {
				iValues = roi.boneMarrowRoiI.stream().mapToInt(i -> i);
				jValues = roi.boneMarrowRoiJ.stream().mapToInt(j -> j);
			}
			marrowCenter[0] = iValues.average().orElse(0.0);
			marrowCenter[1] = jValues.average().orElse(0.0);
			final DetectedEdge edge = roi.edges.get(roi.selection);
			final double[] radii = new double[edge.length];
			for (int i = 0; i < edge.length; ++i) {
				final double x = edge.iit.get(i) - marrowCenter[0];
				final double y = edge.jiit.get(i) - marrowCenter[1];
				radii[i] = Math.sqrt(x * x + y * y);
			}
			final double[] sumRadii = new double[radii.length];
			for (int i = 5; i < radii.length - 6; ++i) {
				for (int j = -5; j < 6; ++j) {
					sumRadii[i] += radii[i + j];
				}
			}
			final double greatestR = Arrays.stream(sumRadii).max().orElse(0);
			int index = 0;
			while (Double.compare(sumRadii[index], greatestR) != 0)
			{
				++index;
			}
			final double x = edge.iit.get(index) - marrowCenter[0];
			final double y = edge.jiit.get(index) - marrowCenter[1];
			alpha = Math.PI - Math.atan2(y, x);
		}

		/*Rotate unselected bone to right*/
		if (choice.equals(labels[3]) || choice.equals(labels[4])) {
			/*Create temp roi for rotating using rotationThreshold..*/
			final SelectROI tempRoi;
			try {
				tempRoi = new SelectROI(roi.scaledImageData, roi.details, roi.imp,
					details.rotationThreshold, false);
			}
			catch (final ExecutionException e) {
				e.printStackTrace();
				return;
			}
			/*Find the second biggest bone (could be bigger than the selected roi...*/
			final int[] twoBones = RoiSelector.twoLargestBonesDetectedEdges(
				tempRoi.edges);
			final int otherBoneSelection;
			if (tempRoi.selection == twoBones[0]) {
				otherBoneSelection = twoBones[1];
			}
			else {
				otherBoneSelection = twoBones[0];
			}
			/*Fill a sieve with a second bone and acquire coordinates...*/
			final Vector<Integer> sRoiI = tempRoi.edges.get(otherBoneSelection).iit;
			final Vector<Integer> sRoiJ = tempRoi.edges.get(otherBoneSelection).jiit;
			final byte[] secondBoneSieve = tempRoi.fillSieve(sRoiI, sRoiJ,
				tempRoi.width, tempRoi.height, tempRoi.scaledImage,
				details.rotationThreshold);
			final double[] selectedBoneCenter = calculateCenter(tempRoi.sieve,
				tempRoi.width, tempRoi.height); /*Calculate selected bone centre*/
			final double[] otherBoneCenter = calculateCenter(secondBoneSieve,
				tempRoi.width, tempRoi.height); /*Calculate other bone centre*/
			final double x;
			final double y;
			if (choice.equals(labels[3])) {
				// Rotate unselected bone to right
				// Use the selected bone as origin for rotation
				x = otherBoneCenter[0] - selectedBoneCenter[0];
				y = otherBoneCenter[1] - selectedBoneCenter[1];
			}
			else {
				// Rotate selected bone to right
				// Use the selected bone as origin for rotation
				x = selectedBoneCenter[0] - otherBoneCenter[0];
				y = selectedBoneCenter[1] - otherBoneCenter[1];
			}
			alpha = -Math.atan2(y, x);
			distanceBetweenBones = Math.sqrt(Math.pow(x, 2.0) + Math.pow(y, 2.0)) *
				roi.pixelSpacing;
		}

		/*Manual rotation*/
		if (details.manualRotation) {
			alpha = details.manualAlpha;
		}

		/*Flip distribution*/
		if (details.flipDistribution) {
			rotationCorrection = -rotationCorrection;
		}

		final int rotationIndex = (int) (alpha / Math.PI * 180.0 +
			rotationCorrection);

		// Calculate CSMIs and rotation angle to align maximal and minimal bending
		// axes with X and Y axes
		pind = rotateIndex(rotationIndex);

		pindColor = details.flipDistribution ? rotateIndex(rotationIndex)
			: rotateIndex(-rotationIndex);
	}

	private static double[] calculateCenter(final byte[] sieve, final int width,
		final int height)
	{
		final double[] originBone = new double[3];
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i) {
				if (sieve[i + j * width] > 0) {
					originBone[0] += i;
					originBone[1] += j;
					originBone[2] += 1;
				}
			}
		}
		originBone[0] /= originBone[2];
		originBone[1] /= originBone[2];
		return originBone;
	}

	private static double[] csmi(final byte[] sieve, final int width,
		final int height)
	{
		final double[] cortexCenter = new double[2];
		double points = 0;
		final Vector<Integer> bmcI = new Vector<>();
		final Vector<Integer> bmcJ = new Vector<>();
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (sieve[i + j * width] > 0) {
					cortexCenter[0] += i;
					cortexCenter[1] += j;
					bmcI.add(i);
					bmcJ.add(j);
					++points;
				}
			}
		}
		cortexCenter[0] /= points;
		cortexCenter[1] /= points;

		final double[] returnValues = new double[3];
		// Calculating cross-sectional moment of inertia in the original image
		// orientation
		for (int i = 0; i < bmcI.size(); i++) {
			final double x = bmcI.get(i) - cortexCenter[0];
			final double y = bmcJ.get(i) - cortexCenter[1];
			returnValues[0] += x * x;
			returnValues[1] += y * y;
			returnValues[2] += x * y;
		}
		return returnValues;
	}

	private void determineMomentAlpha(final double[] csmiValues) {
		final double xmax = csmiValues[0];
		final double ymax = csmiValues[1];
		// Calculate rotation required to align rotation axes
		if (Double.compare(xmax, ymax) == 0) {
			alpha = 0;
			return;
		}
		final double moment = csmiValues[2];
		final double vali1;
		final double vali2;
		alpha = Math.atan(2.0 * moment / (ymax - xmax)) / 2.0;
		// Calculate the maximal and minimial cross-sectional moments of inertia
		vali1 = (ymax + xmax) / 2 + (ymax - xmax) / 2 * Math.cos(2 * (-alpha)) -
			moment * Math.sin(2 * (-alpha));
		vali2 = (ymax + xmax) / 2 - (ymax - xmax) / 2 * Math.cos(2 * (-alpha)) +
			moment * Math.sin(2 * (-alpha));
		// The according to Imax/Imin alpha may align rotation axis
		// corresponding
		// to maximal CSMI with either horizontal
		// or vertical axis, whichever rotation is smaller...
		// Always rotate towards horizontal axis... maximal bending axis will be
		// aligned with horizontal axis
		// Note that e.g. tibial mid-shaft rotation is completely different if
		// only tibia or if both tibia and fibula
		// are consireder!!!
		if (vali1 > vali2) {
			if (alpha < 0) {
				alpha = alpha + Math.PI / 2.0;
			}
			else {
				alpha = alpha - Math.PI / 2.0;
			}
		}
	}

	private Vector<Integer> rotateIndex(final int rotationAngle) {
		final int initialIndex = 360 - Math.abs(rotationAngle);
		final Vector<Integer> rotateIndexVector = Stream.iterate(initialIndex,
			i -> i++ % 360).limit(360).collect(toCollection(Vector::new));
		/*Flip rotateIndexVector, for e.g. comparing left to right*/
		if (details.flipDistribution) {
			Collections.reverse(rotateIndexVector);
		}
		return rotateIndexVector;
	}
}
