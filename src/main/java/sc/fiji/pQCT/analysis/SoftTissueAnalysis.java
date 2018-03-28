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

		double weightedFatArea = 0;
		double weightedLimbArea = 0;
		for (int i = 0; i < roi.width * roi.height; i++) {
			if (roi.softSieve[i] > 0) {
				// Bone & Marrow not excluded!!
				limbA += 1;
				limbD += roi.softScaledImage[i];
				weightedLimbArea += roi.softScaledImage[i] + 1000.0;
			}
			if (roi.softSieve[i] == 2 || roi.softSieve[i] == 4 ||
				roi.softSieve[i] == 5)
			{ // Fat
				fatA += 1;
				fatD += roi.softScaledImage[i];
				weightedFatArea += roi.softScaledImage[i] + 1000.0;
			}
			if (roi.softSieve[i] == 3) {
				// Muscle no IntraFat
				muA += 1;
				muD += roi.softScaledImage[i];
				totalMuA += 1;
				totalMuD += roi.softScaledImage[i];
			}
			if (roi.softSieve[i] == 4) {
				// IntraFat
				intraMuFatA += 1;
				intraMuFatD += roi.softScaledImage[i];
				totalMuA += 1;
				totalMuD += roi.softScaledImage[i];
			}
			if (roi.softSieve[i] == 5) {
				// subCutFat
				subCutFatA += 1;
				subCutFatD += roi.softScaledImage[i];
			}
			if (roi.softSieve[i] == 6) {
				// Bone area
				boneA += 1;
				boneD += roi.softScaledImage[i];
			}
			if (roi.softSieve[i] == 7) {
				// MedFat
				meA += 1;
				meD += roi.softScaledImage[i];
			}
			if (roi.eroded[i] == 1) {
				// PeeledA
				peeledA += 1;
				peeledD += roi.softScaledImage[i];
			}

		}

		final double areaScale = roi.pixelSpacing * roi.pixelSpacing / 100.0;
		limbD /= limbA;
		limbA *= areaScale;
		fatD /= fatA;
		fatA *= areaScale;

		meD /= meA;
		meA *= areaScale;
		boneD /= boneA;
		boneA *= areaScale;
		peeledD /= peeledA;
		peeledA *= areaScale;

		// Added SubCutFatDMedian 2016/01/08
		final double[] subCFatPixels = new double[(int) subCutFatA];
		int cnt = 0;
		for (int i = 0; i < roi.width * roi.height; i++) {
			if (roi.softSieve[i] == 5) {
				// subCutFat
				subCFatPixels[cnt] += roi.softScaledImage[i];
				++cnt;
			}
		}
		subCutFatDMedian = median(subCFatPixels);
		subCutFatD /= subCutFatA;
		subCutFatA *= areaScale;
		muD /= muA;
		muA *= areaScale;
		totalMuD /= totalMuA;
		totalMuA *= areaScale;
		intraMuFatD /= intraMuFatA;
		intraMuFatA *= areaScale;
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
