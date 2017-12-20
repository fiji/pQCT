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

import java.util.Vector;

import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.selectroi.SelectROI;

public class MassDistribution {

	public final int height;
	public final int width;
	public final double[] bMCs;
	private final double sectorWidth;
	private final double[] boneCenter;
	private final Vector<Integer> pind;
	private final SelectROI roi;
	private final double[] bMC = new double[360];

	public MassDistribution(final SelectROI roi,
		final ImageAndAnalysisDetails details, final DetermineAlpha determineAlpha)
	{
		pind = determineAlpha.pind;
		this.roi = roi;
		sectorWidth = details.sectorWidth;
		height = roi.height;
		width = roi.width;
		boneCenter = new double[2];
		double points = 0;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				if (roi.sieve[i + j * width] > 0) {
					boneCenter[0] += i;
					boneCenter[1] += j;
					++points;
				}
			}
		}
		boneCenter[0] /= points;
		boneCenter[1] /= points;
		final int sectors = (int) (360 / sectorWidth);
		bMCs = new double[sectors];
		calculateDistribution();
		rotateResults();
	}

	private void calculateDistribution() {
		// Calculate radii in polar coordinate system originating from bone marrow
		// center of mass
		final double rIncrement = 0.1;
		// Finding endocortical and pericortical
		// borders uMath.sing polar coordinates
		for (int et = 0; et < 360; et++) {
			final double theta = Math.PI / 180.0 * et;
			bMC[et] = 0;
			double R = 0;
			while (roi.sieve[(int) (boneCenter[0] + R * Math.cos(theta)) +
				((int) ((boneCenter[1] + R * Math.sin(theta))) * width)] > 0 ||
				roi.sieve[(int) (boneCenter[0] + (R + 0.5) * Math.cos(theta)) +
					((int) ((boneCenter[1] + (R + 0.5) * Math.sin(theta))) *
						width)] > 0 || roi.sieve[(int) (boneCenter[0] + (R + 1) * Math.cos(
							theta)) + ((int) ((boneCenter[1] + (R + 1) * Math.sin(theta))) *
								width)] > 0 || roi.sieve[(int) (boneCenter[0] + (R + 2) * Math
									.cos(theta)) + ((int) ((boneCenter[1] + (R + 2) * Math.sin(
										theta))) * width)] > 0 || roi.sieve[(int) (boneCenter[0] +
											(R + 3) * Math.cos(theta)) + ((int) ((boneCenter[1] + (R +
												3) * Math.sin(theta))) * width)] > 0 ||
				roi.sieve[(int) (boneCenter[0] + (R + 4) * Math.cos(theta)) +
					((int) ((boneCenter[1] + (R + 4) * Math.sin(theta))) * width)] > 0 ||
				roi.sieve[(int) (boneCenter[0] + (R + 6) * Math.cos(theta)) +
					((int) ((boneCenter[1] + (R + 6) * Math.sin(theta))) * width)] > 0)
			{
				// Calculate bMC rho*dV, dV=dA*slice_thickness
				// dA=pi*((r(et)*resolution)^2-((r(et)-0.1)*resolution)^2),
				// slice_thickness = 1 mm
				// (could be set to actual slice thickness, but makes no
				// difference for comparisons -> 1 mm is used bMD divided by
				// 1000, because unit is mg/cm3 and area is mm2
				final double tempBMD = roi.scaledImage[(int) (boneCenter[0] + R * Math
					.cos(theta)) + ((int) ((boneCenter[1] + R * Math.sin(theta))) *
						width)];
				bMC[et] += tempBMD / 1000.0 * Math.PI / 360.0 * ((R *
					roi.pixelSpacing) * (R * roi.pixelSpacing) - ((R - rIncrement) *
						roi.pixelSpacing) * ((R - rIncrement) * roi.pixelSpacing));
				R += rIncrement;
			}
		}
	}

	private void rotateResults() {
		// Calculate the division and sector values of vBMD
		for (int pp = 0; pp < bMCs.length; pp++) {
			for (int dd = 0; dd < sectorWidth; dd++) {
				bMCs[pp] += bMC[pind.get((int) (pp * sectorWidth + dd))];
			}
		}
	}
}
