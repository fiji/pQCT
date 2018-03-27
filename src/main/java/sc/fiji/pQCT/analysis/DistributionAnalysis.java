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
import static java.util.stream.IntStream.range;

import java.util.List;
import java.util.Vector;

import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.selectroi.SelectROI;

public class DistributionAnalysis {

	public final int height;
	public final int width;
	public final double[] marrowCenter = new double[2];
	public final double[] theta = new double[360];
	public final double[] r = new double[360];
	public final double[] r2 = new double[360];
	public final double[] endocorticalRadii;
	public final double[] pericorticalRadii;
	public final double[] endoCorticalBMDs;
	public final double[] midCorticalBMDs;
	public final double[] periCorticalBMDs;
	public final double[] radialDistribution;
	public final double[] polarDistribution;
	public final double peeledBMD;
	private final double pixelSpacing;
	private final double[] originalROI;
	private double[] peeledROI;
	private final double sectorWidth;
	private final double divisions;
	private final double minimum;
	private final double threshold;
	private final double[] rS = new double[360];
	private final double[] rU = new double[360];
	private final List<double[]> bMDJ = new Vector<>();
	private final Vector<Integer> pInd;
	private double maxRadius;

	public DistributionAnalysis(final SelectROI roi,
		final ImageAndAnalysisDetails details, final DetermineAlpha determineAlpha)
	{
		pInd = determineAlpha.pind;
		sectorWidth = details.sectorWidth;
		final int size = (int) (360.0 / sectorWidth);
		endocorticalRadii = new double[size];
		pericorticalRadii = new double[size];
		endoCorticalBMDs = new double[size];
		midCorticalBMDs = new double[size];
		periCorticalBMDs = new double[size];
		polarDistribution = new double[size];
		divisions = details.divisions;
		radialDistribution = new double[(int) divisions];
		final boolean preventPeeling = details.preventPeeling;
		threshold = details.bMDThreshold;
		minimum = roi.minimum;
		final Vector<Integer> marrowI = roi.boneMarrowRoiI;
		final Vector<Integer> marrowJ = roi.boneMarrowRoiJ;
		height = roi.height;
		width = roi.width;
		pixelSpacing = roi.pixelSpacing;
		originalROI = roi.cortexROI.clone();
		peeledROI = roi.cortexROI.clone();
		
		//Test peeledROI min and max values
		final int peeledSize = width * height;
		peeledROI = erode(peeledROI,width,height,minimum);
		
		for (int i = 0; i < marrowI.size(); i++) {
			marrowCenter[0] += (double) marrowI.get(i);
			marrowCenter[1] += (double) marrowJ.get(i);
		}
		marrowCenter[0] /= marrowI.size();
		marrowCenter[1] /= marrowJ.size();
		
		peeledBMD = range(0, peeledSize).filter(i -> peeledROI[i] >= threshold).mapToDouble(ii -> peeledROI[ii]).average().orElse(0.0);

		//Try old implementation here
		final Vector<Integer> cortexI = new Vector<>();
		final Vector<Integer> cortexJ = new Vector<>();
		double maxRadiusY = 0;
		for (int j = 0; j< height;j++){
			for (int i = 0; i<width;i++){
				if (peeledROI[i+j*width] >= threshold){
					if (Math.sqrt((i-marrowCenter[0])*(i-marrowCenter[0])+(j-marrowCenter[1])*(j-marrowCenter[1])) > maxRadiusY){
						maxRadiusY = Math.sqrt((i-marrowCenter[0])*(i-marrowCenter[0])+(j-marrowCenter[1])*(j-marrowCenter[1]));
					}
				}
				if (originalROI[i+j*width] >= threshold){
					cortexI.add(i);
					cortexJ.add(j);
				}
			}
		}
		final double[] cortexCenter = new double[2];
		for (int i = 0; i< cortexI.size();i++){
			cortexCenter[0]+=(double)cortexI.get(i);
			cortexCenter[1]+=(double)cortexJ.get(i);

		}
		cortexCenter[0] /= cortexI.size();
		cortexCenter[1] /= cortexJ.size();

		// y for cortical pixels. used for BSI calculations, i.e. density weighted section modulus
		maxRadiusY = 0;
		for (int i = 0; i< cortexI.size();i++){
			if (Math.sqrt((cortexI.get(i) -cortexCenter[0])*(cortexI.get(i) -cortexCenter[0])
				+(cortexJ.get(i) -cortexCenter[1])*(cortexJ.get(i) -cortexCenter[1])) > maxRadiusY){
				maxRadiusY = Math.sqrt((cortexI.get(i) -cortexCenter[0])*(cortexI.get(i) -cortexCenter[0])
				+(cortexJ.get(i) -cortexCenter[1])*(cortexJ.get(i) -cortexCenter[1]));
			}
		}

		maxRadius = range(0, peeledSize).filter(i -> originalROI[i] >= threshold)
			.mapToDouble(index -> {
				int i = index % width;
				int j = (index-i) / width;
				double x = i - marrowCenter[0];
				double y = j - marrowCenter[1];
				return Math.sqrt(x * x + y * y);
			}).max().orElse(0.0);

		//Needs to be rounded to 0.1
		maxRadius = Math.round(maxRadius*10.0)/10.0;

		calculateRadii(preventPeeling);
		rotateResults();
	}

	// TODO Add a boolean parameter preventPeeling, and combine method with
	// calculateRadiiNoPeeling
	private void calculateRadii(final boolean preventPeeling) {
		// Calculate radii in polar coordinate system originating from bone marrow
		// center of mass
		for (int i = 0; i < divisions; ++i) {
			bMDJ.add(new double[360]);
		}
		// Finding endocortical and pericortical borders uMath.sing polar
		// coordinates
		final double x = marrowCenter[0];
		final double y = marrowCenter[1];
		
		for (int et = 0; et < 360; ++et) {
			final Vector<Double> BMD_temp = new Vector<>();
			theta[et] = Math.PI / 180.0 * et;
			
			if (et > 0) {
				r[et] = Math.round((rS[et - 1] / 2.0)* 10.0)/ 10.0;
			}
			
			// Anatomical endosteal border
			final double sinTheta = Math.sin(theta[et]);
			final double cosTheta = Math.cos(theta[et]);
			
			r[et] = expandRadius(originalROI, threshold, r[et], x, y, cosTheta,sinTheta);
			rS[et] = r[et];
			if (preventPeeling){
				r2[et] = r[et];
			}else{
				r[et] = expandRadius(peeledROI, 1.0, r[et], x, y, cosTheta, sinTheta);
				r2[et] = r[et];
				r[et] = r[et] + 0.1;
			}
			

			// Return from rMax to identify periosteal border
			double rTemp = maxRadius;
			final double[] roiToObserve = preventPeeling ? originalROI : peeledROI;

			while (	rTemp > r2[et]){
				final int index = (int) (x+rTemp*cosTheta)+ (((int) (y+rTemp*sinTheta))*width);
				if (roiToObserve[index]>0){
					// The loop went until no longer on bone
					rTemp+= 0.1;
					break;
				}
				rTemp -= 0.1;				
			}
			
			// Identify anatomical periosteal border
			if (preventPeeling){
				rU[et] = rTemp;
			}else{
				rU[et] = expandRadiusMulti(originalROI, threshold, rTemp, x, y, cosTheta,sinTheta);
			}

			// Get BMD through the cortex by repeating the incrementing
			while (r[et]<rTemp){
				r[et] = r[et] + 0.1;
				final int index = (int) (x+r[et]*cosTheta)+ (((int) (y+r[et]*sinTheta))*width);
				if (roiToObserve[index] > 0){
					BMD_temp.add(originalROI[index]);
				}
			}
			
			// Get the BMDs here
			
			// Dividing the cortex to three divisions -> save the mean vBMD for each
			// division
			final double analysisThickness = BMD_temp.size();
			if (analysisThickness < divisions) {
				break;
			}
			for (int div = 0; div < divisions; ++div) {
				int mo = 0;
				for (int ka = (int) (analysisThickness * div / divisions); ka <(int) (analysisThickness *(div +1.0)/ divisions); ka++){
					bMDJ.get(div)[et] += BMD_temp.get(ka);
					mo++;
				}
				bMDJ.get(div)[et] /= mo;
				
			}
		}
	}

	// TODO Refactor into a static utility method for all classes instead of
	// repeating code
	public static double[] erode(final double[] data, final int width, final int height, final double bgVal) {
		// Erode algorithm
		// Modified from the best dilate by one solution taken from
		// http://ostermiller.org/dilate_and_erode.html
		for (int i = 1; i < height-1; i++) {
			for (int j = 1; j < width-1; j++) {
				final int index = i * width + j;
				if (data[index] > bgVal) {
					if (data[(i - 1) * width + j] == bgVal | 
						data[(i) * width + j - 1] == bgVal |
						data[(i + 1) * width + j] == bgVal |
						data[(i) * width +j + 1] == bgVal)
					{
						// Erode the pixel if any of the neighborhood pixels is background
						data[index] = bgVal - 1;
					}
				}
			}
		}
		for (int i = 0; i < width * height; i++) {
			if (data[i] < bgVal) {
				data[i] = bgVal;
			}
		}
		return data;
	}

	// TODO Replace while(roi[x + r * cos(theta) ...]) loops with similar methods
	private double expandRadius(final double[] roi, final double threshold,
		final double radius, final double x, final double y, final double cos,
		final double sin)
	{
		double expandedR = radius;
		final double maxR = maxRadius;
		while (true) {
			final int index = (int) (x+expandedR*cos)+ (((int) (y+expandedR*sin))*width);
			if (roi[index] >= threshold | expandedR >= maxR) {
				break;
			}
			expandedR += 0.1;
		}
		return expandedR;
	}

	// TODO Replace while(roi[x + r * cos(theta) ...] || roi[x + r + 2 *
	// cos(theta) ...]) loops with similar methods
	private double expandRadiusMulti(final double[] roi, final double threshold,
		final double radius, final double x, final double y, final double cos,
		final double sin)
	{
		double expandedR = radius;
		final double maxR = maxRadius / pixelSpacing;
		while (true) {
			final double radii[] = { expandedR, expandedR + 0.5, expandedR + 1.0,
				expandedR + 2.0, expandedR + 3.0, expandedR + 4.0, expandedR + 6.0 };
			final int[] indices = stream(radii).mapToInt(r -> (int) ((x + r * cos) +
				((int) (y + r * sin)) * width)).toArray();
			if (stream(indices).noneMatch(i -> roi[i] > threshold) |
				expandedR >= maxR)
			{
				break;
			}
			expandedR += 0.1;
		}
		return expandedR;
	}

	// TODO Refactor into a static utility method for all classes instead of
	// repeating code
	private void rotateResults() {
		// Calculate the endocortical and pericortical radii along with the
		// corresponding radii after peeling one layer of pixels
		final double[] pRad = stream(rU).map(r -> r *= pixelSpacing).toArray();
		final double[] eRad = stream(rS).map(r -> r *= pixelSpacing).toArray();
		final int size = (int) (360.0 / sectorWidth);
		final double[][] corticalDensity = new double[(int) divisions][size];
		// Calculate the division and sector values of vBMD
		for (int pp = 0; pp < size; ++pp) {
			for (int dd = 0; dd < (int) sectorWidth; ++dd) {
				
				
				
				final int index = pInd.get((int) (pp * sectorWidth + dd));
				
					
				
				
				
				endocorticalRadii[pp] += eRad[index] / sectorWidth;
				pericorticalRadii[pp] += pRad[index] / sectorWidth;
				
				
				// Cortex
				endoCorticalBMDs[pp] += bMDJ.get(0)[index] / sectorWidth;
				midCorticalBMDs[pp] += bMDJ.get(1)[index] / sectorWidth;
				periCorticalBMDs[pp] += bMDJ.get(2)[index] / sectorWidth;
			}
			corticalDensity[0][pp] = endoCorticalBMDs[pp];
			corticalDensity[1][pp] = midCorticalBMDs[pp];
			corticalDensity[2][pp] = periCorticalBMDs[pp];
		}

		// Radial distribution
		for (int i = 0; i < divisions; ++i) {
			for (int j = 0; j < size; ++j) {
				radialDistribution[i] += corticalDensity[i][j];
			}
			radialDistribution[i] /= size;
		}
		
		// Polar distribution
		for (int j = 0; j < size; ++j) {
			for (int i = 0; i < divisions; ++i) {
				polarDistribution[j] += corticalDensity[i][j];
			}
			polarDistribution[j] /= divisions;
		}
	}
}
