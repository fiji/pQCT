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
	private final double[] peeledROI;
	private final double sectorWidth;
	private final double divisions;
	private final double minimum;
	private final double threshold;
	private final double[] rS = new double[360];
	private final double[] rU = new double[360];
	private final List<double[]> bMDJ = new Vector<>();
	private final Vector<Integer> pInd;
	private final double maxRadius;

	public DistributionAnalysis(final SelectROI roi,
		final ImageAndAnalysisDetails details, final DetermineAlpha determineAlpha)
	{
		pInd = determineAlpha.pind;
		sectorWidth = details.sectorWidth;
		final int size = (int) (360 / sectorWidth);
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
		erode(peeledROI);
		for (int i = 0; i < marrowI.size(); i++) {
			marrowCenter[0] += marrowI.get(i);
			marrowCenter[1] += marrowJ.get(i);
		}
		marrowCenter[0] /= marrowI.size();
		marrowCenter[1] /= marrowJ.size();
		final int peeledSize = width * height;
		peeledBMD = range(0, peeledSize).filter(i -> peeledROI[i] >= threshold)
			.average().orElse(0.0);
		maxRadius = range(0, peeledSize).filter(i -> peeledROI[i] >= threshold)
			.mapToDouble(index -> {
				final int j = index / width;
				final int i = index - width * j;
				final double x = i - marrowCenter[0];
				final double y = j - marrowCenter[1];
				return Math.sqrt(x * x + y * y);
			}).max().orElse(0.0);
		if (preventPeeling) {
			calculateRadiiNoPeeling();
		}
		else {
			calculateRadii();
		}
		rotateResults();
	}

	// TODO Add a boolean parameter preventPeeling, and combine method with
	// calculateRadiiNoPeeling
	private void calculateRadii() {
		// Calculate radii in polar coordinate system originating from bone marrow
		// center of mass
		for (int i = 0; i < divisions; ++i) {
			bMDJ.add(new double[360]);
		}
		// Finding endocortical and pericortical borders uMath.sing polar
		// coordinates
		for (int et = 0; et < 360; ++et) {
			final Vector<Double> BMD_temp = new Vector<>();
			theta[et] = Math.PI / 180.0 * et;
			if (et > 0) {
				r[et] = r[et - 1] / 2.0;
			}

			// Anatomical endosteal border
			final double sinTheta = Math.sin(theta[et]);
			final double cosTheta = Math.cos(theta[et]);
			final double x = marrowCenter[0];
			final double y = marrowCenter[1];
			r[et] = expandRadius(originalROI, threshold, r[et], x, y, cosTheta,
				sinTheta);
			rS[et] = r[et];
			r[et] = expandRadius(peeledROI, 1.0, r[et], x, y, cosTheta, sinTheta);
			r2[et] = r[et];
			r[et] = r[et] + 0.1;

			// Peeled periosteal border
			r[et] = expandRadiusMulti(peeledROI, 0.0, r[et], x, y, cosTheta,
				sinTheta);
			final long bMDIterations = (long) ((r[et] - r2[et]) / 1.0);
			for (int i = 1; i < bMDIterations; i++) {
				final double radius = r2[et] + i * 0.1;
				final int index = (int) ((x + radius * cosTheta) + (y + radius *
					sinTheta) * width);
				if (peeledROI[index] > 0) {
					BMD_temp.add(originalROI[index]);
				}
			}

			// Anatomical periosteal border
			rU[et] = r[et];
			rU[et] = expandRadiusMulti(originalROI, threshold, rU[et], x, y, cosTheta,
				sinTheta);

			// Dividing the cortex to three divisions -> save the mean vBMD for each
			// division
			final double analysisThickness = BMD_temp.size();
			if (analysisThickness < divisions) {
				break;
			}
			for (int div = 0; div < divisions; ++div) {
				int mo = 0;
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

	private void calculateRadiiNoPeeling() {
		// Calculate radii in polar coordinate system originating from bone marrow
		// center of mass
		for (int i = 0; i < divisions; ++i) {
			bMDJ.add(new double[360]);
		}
		// Finding endocortical and pericortical
		// borders uMath.sing polar coordinates
		for (int et = 0; et < 360; ++et) {
			theta[et] = Math.PI / 180.0 * et;
			final Vector<Double> BMD_temp = new Vector<>();
			if (et > 0) {
				r[et] = r[et - 1] / 2.0;
			}
			final double x = marrowCenter[0];
			final double y = marrowCenter[1];
			final double cos = Math.cos(theta[et]);
			final double sin = Math.sin(theta[et]);
			// Anatomical endosteal border
			r[et] = expandRadius(originalROI, threshold, r[et], x, y, cos, sin);
			r2[et] = r[et];
			rS[et] = r[et];
			// Anatomical periosteal border
			r[et] = expandRadiusMulti(originalROI, 0.0, r[et], x, y, cos, sin);
			final long bMDIterations = (long) ((r[et] - r2[et]) / 1.0);
			for (int i = 1; i < bMDIterations; i++) {
				final double radius = r2[et] + i * 0.1;
				final int index = (int) ((x + radius * cos) + (y + radius * sin) *
					width);
				if (originalROI[index] > 0) {
					BMD_temp.add(originalROI[index]);
				}
			}
			rU[et] = r[et];
			final double analysisThickness = BMD_temp.size();
			// Dividing the cortex to three divisions -> save the mean vBMD for each
			// division
			if (analysisThickness < divisions) {
				break;
			}
			for (int div = 0; div < divisions; ++div) {
				int mo = 0;
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

	// TODO Refactor into a static utility method for all classes instead of
	// repeating code
	private void erode(final double[] data) {
		// Erode algorithm
		// Modified from the best dilate by one solution taken from
		// http://ostermiller.org/dilate_and_erode.html
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				final int index = i * width + j;
				if (data[index] > minimum) {
					if (i > 0 && data[(i - 1) * width + j] == minimum || j > 0 &&
						data[(i) * width + j - 1] == minimum || i + 1 < height && data[(i +
							1) * width + j] == minimum || j + 1 < width && data[(i) * width +
								j + 1] == minimum)
					{
						data[index] = minimum - 1;
					} // Erode the pixel if any of the neighborhood pixels is background
				}
			}
		}
		for (int i = 0; i < width * height; i++) {
			if (data[i] < minimum) {
				data[i] = minimum;
			}
		}
	}

	// TODO Replace while(roi[x + r * cos(theta) ...]) loops with similar methods
	private double expandRadius(final double[] roi, final double threshold,
		final double radius, final double x, final double y, final double cos,
		final double sin)
	{
		double expandedR = radius;
		final double maxR = maxRadius / pixelSpacing;
		while (true) {
			final int index = (int) (x + expandedR * cos + (y + expandedR * sin) *
				width);
			if (roi[index] >= threshold || expandedR >= maxR) {
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
				(y + r * sin) * width)).toArray();
			if (stream(indices).noneMatch(i -> roi[i] > threshold) ||
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
		final double[] pRad = rU.clone();
		stream(pRad).forEach(r -> r *= pixelSpacing);
		final double[] eRad = rS.clone();
		stream(eRad).forEach(r -> r *= pixelSpacing);
		final double[] pPRad = r.clone();
		stream(pPRad).forEach(r -> r *= pixelSpacing);
		final double[] pERad = r2.clone();
		stream(pERad).forEach(r -> r *= pixelSpacing);
		final int size = (int) (360 / sectorWidth);
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
