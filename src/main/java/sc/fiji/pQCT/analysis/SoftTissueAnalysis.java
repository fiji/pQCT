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

import java.util.DoubleSummaryStatistics;
import java.util.function.Function;
import java.util.stream.IntStream;

import sc.fiji.pQCT.selectroi.SelectSoftROI;

//Debugging
import ij.IJ;		//log

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
				
		for (int i = 0;i<stats.length;++i){
			IJ.log(String.format("Stats %d one %.2f two %.2f",i,stats[i][0],stats[i][1]));
		}
		
		//Revert these back to the old ones...
				
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
