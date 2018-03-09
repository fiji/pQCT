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

package sc.fiji.pQCT.analysis;

import static java.util.Arrays.stream;

import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.selectroi.SelectROI;

//Debugging
//import ij.IJ;

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
		final int size = (int) (360.0 / sectorWidth);
		pericorticalRadii = new double[size];
		calculateRadii();
		rotateResults();
		
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
		//Calculate pericortical radii
		//Double[] pRad = DistributionAnalysis.clone(rU);
		//IJ.log(String.format("rU length %d",rU.length));
		//IJ.log(String.format("pRad length %d",pRad.length));
		double[] pRad =(double []) stream(rU).map(r -> {return r *= pixelSpacing;}).toArray();
		//IJ.log(String.format("pRad length %d",pRad.length));
		final int size = (int) (360.0 / sectorWidth);
		//IJ.log(String.format("size %d",size));
		for (int div = 0; div < divisions; ++div) {
			BMDs.add(new double[size]);
		}
		// Calculate the endocortical and pericortical radii along with the
		// corresponding radii after peeling one layer of pixels.
		// Calculate the division and sector values of vBMD
		for (int pp = 0; pp < size; pp++) {
			for (int dd = 0; dd < (int) sectorWidth; dd++) {
				int index = pind.get((int) (pp * sectorWidth + dd));
				//IJ.log(String.format("pp %d dd %d index %d length prRad %d length pRad %d",pp,dd,index,
				//			pericorticalRadii.length, pRad.length));
				pericorticalRadii[pp] += pRad[index] / sectorWidth;
				for (int div = 0; div < divisions; ++div) {
					BMDs.get(div)[pp] += bMDJ.get(div)[index] / sectorWidth;
				}
			}
		}
		
	}
}
