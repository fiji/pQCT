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

import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.selectroi.SelectROI;

public class ConcentricRingAnalysis {

	private final double pixelSpacing;
	public final int height;
	public final int width;
	public final double[] boneCenter;
	// Variables for radii calculations
	public final double[] theta = new double[360];
	public final double[] rU = new double[360];
	public final double[] pericorticalRadii;
	public final Vector<double[]> BMDs = new Vector<>();
	private final double sectorWidth;
	private final double divisions;
	// Variables for moment calculations
	private final Vector<Integer> pind;
	private final SelectROI roi;
	private final List<double[]> bMDJ = new Vector<>();

	public ConcentricRingAnalysis(final SelectROI roi,
		final ImageAndAnalysisDetails details, final DetermineAlpha determineAlpha)
	{
		pind = determineAlpha.pind;
		this.roi = roi;
		sectorWidth = details.concentricSector;
		divisions = details.concentricDivisions;
		height = roi.height;
		width = roi.width;
		pixelSpacing = roi.pixelSpacing;
		boneCenter = new double[2];
		int points = 0;
		for (int j = 0; j < height; j++) {
			final int offset = j * width;
			for (int i = 0; i < width; i++) {
				if (roi.sieve[offset + i] > 0) {
					boneCenter[0] += i;
					boneCenter[1] += j;
					++points;
				}
			}
		}
		boneCenter[0] = boneCenter[0] / points;
		boneCenter[1] = boneCenter[1] / points;
		calculateRadii();
		rotateResults();
		final int size = (int) (360.0 / sectorWidth);
		pericorticalRadii = new double[size];
	}

	private void calculateRadii() {
		// Calculate radii in polar coordinate system originating from bone marrow
		// center of mass
		for (int i = 0; i < divisions; ++i) {
			bMDJ.add(new double[360]);
		}
		final double rIncrement = 0.1;
		// Finding endocortical and pericortical
		// borders uMath.sing polar coordinates
		for (int et = 0; et < 360; ++et) {
			final Vector<Double> BMD_temp = new Vector<>();
			theta[et] = Math.PI / 180.0 * et;
			final double cosTheta = Math.cos(theta[et]);
			final double sinTheta = Math.sin(theta[et]);
			final double cX = boneCenter[0];
			final double cY = boneCenter[1];
			double r = 0;
			while (roi.sieve[(int) (cX + r * cosTheta) + ((int) ((cY + r *
				sinTheta)) * width)] > 0 || roi.sieve[(int) (cX + (r + 0.5) *
					cosTheta) + ((int) ((cY + (r + 0.5) * sinTheta)) * width)] > 0 ||
				roi.sieve[(int) (cX + (r + 1) * cosTheta) + ((int) ((cY + (r + 1) *
					sinTheta)) * width)] > 0 || roi.sieve[(int) (cX + (r + 2) *
						cosTheta) + ((int) ((cY + (r + 2) * sinTheta)) * width)] > 0 ||
				roi.sieve[(int) (cX + (r + 3) * cosTheta) + ((int) ((cY + (r + 3) *
					sinTheta)) * width)] > 0 || roi.sieve[(int) (cX + (r + 4) *
						cosTheta) + ((int) ((cY + (r + 4) * sinTheta)) * width)] > 0 ||
				roi.sieve[(int) (cX + (r + 6) * cosTheta) + ((int) ((cY + (r + 6) *
					sinTheta)) * width)] > 0)
			{
				// Calculate BMC rho*dV, dV=dA*slice_thickness
				// dA=pi*((r(et)*resolution)^2-((r(et)-0.1)*resolution)^2),
				// slice_thickness = 1 mm
				// (could be set to actual slice thickness, but makes no
				// difference for comparisons -> 1 mm is used bMD divided by
				// 1000, because unit is mg/cm3 and area is mm2
				BMD_temp.add(roi.scaledImage[(int) (cX + r * cosTheta) + ((int) ((cY +
					r * sinTheta)) * width)]);
				r += rIncrement;
			}
			rU[et] = r;
			final int analysisThickness = BMD_temp.size();
			// Dividing the cortex to three divisions -> save the mean vBMD for each
			// division
			if (analysisThickness < divisions) {
				return;
			}
			// cortex
			for (int div = 0; div < divisions; ++div) {
				double mo = 0;
				for (int ka = (int) (analysisThickness * div /
					divisions); ka < (int) (analysisThickness * (div + 1.0) /
						divisions); ka++)
				{
					bMDJ.get(div)[et] += BMD_temp.get(ka);
					mo++;
				}
				bMDJ.get(div)[et] /= mo;
			}
		}
	}

	private void rotateResults() {
		final double[] pRad = Arrays.stream(rU).map(r -> r * pixelSpacing)
			.toArray();
		final int size = (int) (360.0 / sectorWidth);
		for (int div = 0; div < divisions; ++div) {
			BMDs.add(new double[size]);
		}
		// Calculate the endocortical and pericortical radii along with the
		// corresponding radii after peeling one layer of pixels.
		// Calculate the division and sector values of vBMD
		for (int pp = 0; pp < size; pp++) {
			for (int dd = 0; dd < (int) sectorWidth; dd++) {
				final int index = (int) (pp * sectorWidth + dd);
				pericorticalRadii[pp] += pRad[pind.get(index)] / sectorWidth;
				for (int div = 0; div < divisions; ++div) {
					BMDs.get(div)[pp] += bMDJ.get(div)[pind.get(index)] / sectorWidth;
				}
			}
		}
	}
}
