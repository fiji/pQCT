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

package sc.fiji.pQCT.selectroi;

import java.awt.Polygon;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import sc.fiji.pQCT.io.ImageAndAnalysisDetails;
import sc.fiji.pQCT.io.ScaledImageData;

public class SelectROI extends RoiSelector {

	public final Vector<DetectedEdge> edges;

	public SelectROI(final ScaledImageData dataIn,
		final ImageAndAnalysisDetails detailsIn, final ImagePlus imp,
		final double boneThreshold, final boolean setRoi) throws ExecutionException
	{
		super(dataIn, detailsIn, imp);
		/*Select ROI and set everything else than the roi to minimum*/
		cortexROI = new double[width * height];
		cortexRoiI = new Vector<>();
		cortexRoiJ = new Vector<>();
		cortexAreaRoiI = new Vector<>();
		cortexAreaRoiJ = new Vector<>();
		boneMarrowRoiI = new Vector<>();
		boneMarrowRoiJ = new Vector<>();
		Roi ijROI = imp.getRoi();
		final double[] tempScaledImage = scaledImage.clone();
		// scaledImage.clone();
		if (ijROI != null && details.manualRoi) {
			// Set pixels outside the manually selected ROI to zero
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					if (!ijROI.contains(i, j)) {
						// Check whether pixel is within ROI, mark with bone threshold
						tempScaledImage[i + j * width] = minimum;
					}
				}
			}
			final Polygon polygon = ijROI.getPolygon();
			if (polygon != null) {
				// Check whether a polygon can be acquired and include its points
				for (int j = 0; j < polygon.npoints; j++) {
					final int index = polygon.xpoints[j] + polygon.ypoints[j] * width;
					tempScaledImage[index] = scaledImage[index];
				}
			}
		}
		final Vector<Object> boneMasks = getSieve(tempScaledImage, boneThreshold,
			details.roiChoice, details.guessStacked, details.stacked,
			details.guessFlip, details.allowCleaving);
		sieve = (byte[]) boneMasks.get(0);
		result = (byte[]) boneMasks.get(1);
		final Vector<DetectedEdge> boneEdges = (Vector<DetectedEdge>) boneMasks.get(
			2);
		selection = (Integer) boneMasks.get(3);
		/*Add the roi to the image*/
		if (setRoi) {
			final int[] xcoordinates = new int[boneEdges.get(selection).iit.size()];
			final int[] ycoordinates = new int[boneEdges.get(selection).iit.size()];
			for (int i = 0; i < boneEdges.get(selection).iit.size(); ++i) {
				xcoordinates[i] = boneEdges.get(selection).iit.get(i);
				ycoordinates[i] = boneEdges.get(selection).jiit.get(i);
			}
			/*Flip the original image prior to adding the ROI, if scaled image is flipped*/
			if ((details.flipHorizontal || details.flipVertical) && imp.getRoi() != null){
				// Remove existing ROIs in order to flip the whole image...
				IJ.run(imp, "Select None", "");
			}
			if (details.flipHorizontal) {
				imp.getProcessor().flipVertical();
				imp.updateAndDraw();
			}
			if (details.flipVertical) {
				imp.getProcessor().flipHorizontal();
				imp.updateAndDraw();
			}
			ijROI = new PolygonRoi(xcoordinates, ycoordinates, boneEdges.get(
				selection).iit.size(), Roi.POLYGON);
			imp.setRoi(ijROI);
		}
		
		//Visualise scaledImage
		/*
		ImagePlus tempImage = new ImagePlus("peeledROI");
		tempImage.setProcessor(new FloatProcessor(width,height,scaledImage));
		tempImage.show();
		
		ImagePlus tempImage2 = new ImagePlus("sieve");
		tempImage2.setProcessor(new ByteProcessor(width,height,sieve));
		tempImage2.setDisplayRange(0,1);
		tempImage2.show();
		*/

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++) {
				int index = i + j * width;
				if (scaledImage[index]<areaThreshold & sieve[index] > 0){
					boneMarrowRoiI.add(i);
					boneMarrowRoiJ.add(j);
				}
				if (scaledImage[index]>=areaThreshold & sieve[index] > 0){
					cortexAreaRoiI.add(i);
					cortexAreaRoiJ.add(j);
				}
				if (scaledImage[index]>=BMDthreshold & sieve[index] > 0){
					cortexROI[index] = scaledImage[index];				
					cortexRoiI.add(i);
					cortexRoiJ.add(j);
				} else {
					cortexROI[index] = minimum;
				}
			}
		}
		edges = boneEdges;
	}
}
