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

import java.util.concurrent.ExecutionException;
import java.util.stream.IntStream;

import sc.fiji.pQCT.selectroi.SelectROI;

public class CorticalAnalysis {

	// Medullary area = ToA-CoA
	public final double medullaryArea;
	public final double maMassD;
	public final double stratecMaMassD;
	public final double sSIMax;
	public final double sSIMin;
	public final double bSId;
	public final double iPo;
	public final double iMax;
	public final double iMin;
	public final double dwIPo;
	public final double dwIMax;
	public final double dwIMin;
	public double bMD;
	public double area;
	public double marrowArea;
	public double marrowDensity;
	public double ToA;
	public double ToD;
	public double sSI;
	// Stratec/Geanie compatible CoD and CoA
	public byte[] cortexSieve;
	public double CoD;
	public double CoA;

	public CorticalAnalysis(final SelectROI roi) {
		ToA = 0;
		ToD = 0;
		marrowArea = 0;
		marrowDensity = 0;
		final int iterations = roi.width * roi.height;
		for (int i = 0; i < iterations; i++) {
			if (roi.sieve[i] <= 0) {
				continue;
			}
			final double value = roi.scaledImage[i];
			ToA++;
			ToD += value;
			// Marrow analysis
			if (value < roi.details.marrowThreshold) {
				marrowArea++;
				marrowDensity += value;
			}
		}
		ToD /= ToA;
		final double spacingSq = roi.pixelSpacing * roi.pixelSpacing;
		ToA *= spacingSq;
		marrowDensity /= marrowArea;
		marrowArea *= spacingSq;

		// Mass density is calculated by converting the bMD to Hounsfield Units, and
		// scaling the HUs to comparable HUs between machinesHUs are then scaled to
		// mass density as described in Schneider et al. Phys. Med. Biol. 45 (2000)
		// 459 -- 478.
		double H; // Machine comparable Hounsfield Unit
		double mu; // Med bMD as attenuation coefficient
		double muH2O; // Water as attenuation coefficient
		muH2O = (0.0 - roi.details.constant) / roi.details.scalingFactor;
		mu = (marrowDensity - roi.details.constant) / roi.details.scalingFactor;
		// Equation 6 in Schneider et al. 2000 *1000 omitted
		H = mu / muH2O - 1.0;
		// Equation 21 in Schneider et al. 2000 *10^-3 omitted
		maMassD = 1.018 + 0.893 * H;

		// Stratec pQCT is calibrated so that fat is 0 vBMD and water is 50 vBMD,
		// Sievanen J Bone Miner Res. 1998 May;13(5):871-82.
		muH2O = (50.0 - roi.details.constant) / roi.details.scalingFactor;
		mu = (marrowDensity - roi.details.constant) / roi.details.scalingFactor;
		H = mu / muH2O - 1.0; // Equation 6 in Schneider et al. 2000 *1000 omitted
		// Equation 21 in Schneider et al. 2000 *10^-3 omitted
		stratecMaMassD = 1.018 + 0.893 * H;
		// To make it look nicer, we'll use a unit of g^2/cm^4
		bSId = ToD * ToD * ToA / 100000000.0;
		bMD = 0;
		area = 0;
		for (int j = 0; j < roi.cortexRoiI.size(); j++) {
			bMD += roi.cortexROI[roi.cortexRoiI.get(j) + roi.cortexRoiJ.get(j) *
				roi.width];
		}
		bMD /= roi.cortexRoiI.size();
		final double[] cortexCenter = new double[2];
		cortexCenter[0] = roi.cortexAreaRoiI.stream().mapToInt(i -> i).average()
			.orElse(-1.0);
		cortexCenter[1] = roi.cortexAreaRoiJ.stream().mapToInt(j -> j).average()
			.orElse(-1.0);
		// Calculate cortical area from 550 threshold...
		final int cortexROIs = roi.cortexAreaRoiI.size();
		area = cortexROIs * roi.pixelSpacing * roi.pixelSpacing;
		medullaryArea = ToA - area;

		final double cortexX = cortexCenter[0];
		final double cortexY = cortexCenter[1];
		// y for cortical pixels. used for BSI calculations, i.e. density weighted
		// section modulus
		final double maxRadiusY = IntStream.range(0, cortexROIs).mapToDouble(i -> {
			final double x = roi.cortexAreaRoiI.get(i) - cortexX;
			final double y = roi.cortexAreaRoiJ.get(i) - cortexY;
			return Math.sqrt(x * x + y * y);
		}).max().orElse(0.0);

		// Calculate CSMIs and rotation angle to align maximal and minimal bending
		// axes with X and Y axes
		double ssimo = 0;
		double ssixmax = 0;
		double ssiymax = 0;
		double moment = 0;
		double xmax = 0;
		double ymax = 0;
		double dwmo = 0;
		double dwxmax = 0;
		double dwymax = 0;
		sSI = 0;
		// Calculating cross-sectional moment of inertia in the original image
		// orientation
		for (int i = 0; i < cortexROIs; i++) {
			final int roiX = roi.cortexAreaRoiI.get(i);
			final int roiY = roi.cortexAreaRoiJ.get(i);
			final double x = roiX - cortexX;
			final double y = roiY - cortexY;
			final double scale = roi.scaledImage[roiX + roiY * roi.width];
			final double pX = x * roi.pixelSpacing;
			final double pY = y * roi.pixelSpacing;
			xmax = xmax + pX * pX * spacingSq;
			ymax = ymax + pY * pY * spacingSq;
			moment = moment + pX * pY * spacingSq;
			final double dwSpacing = (roi.pixelSpacing / 10) * (roi.pixelSpacing /
				10);
			dwxmax = dwxmax + (pX / 10) * (pX / 10) * dwSpacing * scale;
			dwymax = dwymax + (pY / 10) * (pY / 10) * dwSpacing * scale;
			dwmo = dwmo + (pX / 10) * (pY / 10) * dwSpacing * scale;
			final double ssiScale = spacingSq * (scale / 1200);
			final double ssiMaxR = maxRadiusY * roi.pixelSpacing;
			ssixmax += pX * pX * ssiScale / ssiMaxR;
			ssiymax += pY * pY * ssiScale / ssiMaxR;
			ssimo += x * roi.pixelSpacing * y * roi.pixelSpacing * ssiScale / ssiMaxR;
			sSI = sSI + (pX * pX + pY * pY) * ssiScale / ssiMaxR;
		}

		iPo = xmax + ymax;
		dwIPo = dwxmax + dwymax;
		// Ipolar caclulated
		// Calculation of Imax and Imin
		// Calculate rotation required to align rotation axes
		final double alpha = Math.atan(2 * moment / (ymax - xmax)) / 2;
		// Calculate the maximal and minimial cross-sectional moments of inertia
		final double[] moments = momentOfInertiaMinMax(xmax, ymax, alpha, moment);
		iMin = moments[0];
		iMax = moments[1];
		// Calculate the maximal and minimal density weighted cross-sectional
		// moments of inertia
		final double[] dwIMoments = momentOfInertiaMinMax(dwxmax, dwymax, alpha,
			dwmo);
		dwIMin = dwIMoments[0];
		dwIMax = dwIMoments[1];
		final double[] ssi = momentOfInertiaMinMax(ssixmax, ssiymax, alpha, ssimo);
		sSIMin = ssi[0];
		sSIMax = ssi[1];

		// Calculate Stratec/Geanie compatible CoA and CoD, i.e. define a ROI larger
		// than the bone and calculate
		// CoD and CoA from the ROI independent of whether the cortex is continuous.
		final SelectROI tempRoi;
		try {
			tempRoi = new SelectROI(roi.scaledImageData, roi.details, roi.imp,
				roi.details.rotationThreshold, false);
		}
		catch (final ExecutionException e) {
			e.printStackTrace();
			return;
		}
		CoD = 0;
		CoA = 0;
		int CoDcounter = 0;
		cortexSieve = new byte[roi.scaledImage.length];
		for (int j = 0; j < roi.scaledImage.length; ++j) {
			if (tempRoi.sieve[j] > 0 && roi.scaledImage[j] >= roi.BMDthreshold) {
				CoD += roi.scaledImage[j];
				++CoDcounter;
				cortexSieve[j] = 1;
			}
			if (tempRoi.sieve[j] > 0 && roi.scaledImage[j] >= roi.areaThreshold) {
				CoA += 1.0;
			}
		}
		CoD /= CoDcounter;
		CoA *= spacingSq;
	}

	private static double[] momentOfInertiaMinMax(final double xMax,
		final double yMax, final double alpha, final double moment)
	{
		final double angle = 2.0 * (-alpha);
		final double maxAvg = (yMax + xMax) / 2.0;
		final double minAvg = (yMax - xMax) / 2.0;
		final double max = maxAvg + minAvg * Math.cos(angle) - moment * Math.sin(
			angle);
		final double min = maxAvg - minAvg * Math.cos(angle) + moment * Math.sin(
			angle);
		if (min > max) {
			return new double[] { max, min };
		}
		return new double[] { min, max };
	}
}
