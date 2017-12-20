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

import java.util.DoubleSummaryStatistics;
import java.util.function.Function;
import java.util.stream.IntStream;

import sc.fiji.pQCT.selectroi.SelectSoftROI;

public class SoftTissueAnalysis {

	public double muA;
	public double intraMuFatA;
	public double totalMuA;
	public double fatA;
	public double subCutFatA;
	public double limbA;
	public double muD;
	public double intraMuFatD;
	public double totalMuD;
	public double fatD;
	public double subCutFatD;
	public double subCutFatDMedian;
	public double limbD;
	public double fatPercentage;
	public double meA;
	public double meD;
	public double boneA;
	public double boneD;
	public double peeledA;
	public double peeledD;

	public SoftTissueAnalysis(final SelectSoftROI roi) {
		final int roiSize = roi.width * roi.height;
		final Function<Integer, DoubleSummaryStatistics> typeStatistics =
			type -> IntStream.range(0, roiSize).filter(i -> roi.softSieve[i] == type)
				.mapToDouble(i -> roi.softScaledImage[i]).summaryStatistics();
		final double[][] stats = IntStream.range(1, 8).mapToObj(
			typeStatistics::apply).map(summary -> new double[] { summary.getCount(),
				summary.getAverage() }).toArray(double[][]::new);
		final double areaScale = roi.pixelSpacing * roi.pixelSpacing / 100.0;
		limbA = stats[0][0] * areaScale;
		limbD = stats[0][1];
		final double weightedLimbArea = stats[0][0] + 1000 * limbA;
		final double totalFatSum = stats[1][0] + stats[3][0] + stats[4][0];
		final double totalFatCount = stats[1][1] + stats[3][1] + stats[4][1];
		fatA = totalFatCount * areaScale;
		fatD = totalFatSum / totalFatCount;
		final double weightedFatArea = totalFatSum + totalFatCount * 1000;
		muA = stats[2][0] * areaScale;
		muD = stats[2][1];
		intraMuFatA = stats[3][0];
		intraMuFatD = stats[3][1];
		final double totalMuscleCount = stats[2][0] + stats[3][0];
		final double totalMuscleSum = stats[2][1] + stats[3][1];
		totalMuA = totalMuscleCount * areaScale;
		totalMuD = totalMuscleSum / totalMuscleCount;
		subCutFatA = stats[4][0] * areaScale;
		subCutFatD = stats[4][1];
		boneA = stats[5][0] * areaScale;
		boneD = stats[5][1];
		meA = stats[6][0] * areaScale;
		meD = stats[6][1];
		final DoubleSummaryStatistics peeledStats = IntStream.range(0, roiSize)
			.filter(i -> roi.eroded[i] == 1).mapToDouble(i -> roi.softScaledImage[i])
			.summaryStatistics();
		peeledD = peeledStats.getAverage();
		peeledA = peeledStats.getCount() * areaScale;
		final double[] subCFatPixels = IntStream.range(0, roiSize).filter(
			i -> roi.softSieve[i] == 5).mapToDouble(i -> roi.softScaledImage[i])
			.sorted().toArray();
		subCutFatDMedian = median(subCFatPixels);
		fatPercentage = (weightedFatArea / weightedLimbArea) * 100.0;
	}

	// TODO Use a pre-existing method
	private static double median(final double[] a) {
		if (a.length < 0) {
			return 0.0;
		}
		if (a.length == 1) {
			return a[0];
		}
		final int middle = a.length / 2;
		if (a.length % 2 == 0) {
			return (a[middle - 1] + a[middle]) / 2;
		}
		else {
			return a[middle];
		}
	}
}
