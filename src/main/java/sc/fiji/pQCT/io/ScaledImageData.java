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

package sc.fiji.pQCT.io;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class ScaledImageData {

	public final double[] scaledImage;
	public final double[] softScaledImage;
	public final double minimum;
	public final int width;
	public final int height;
	public final double pixelSpacing;

	// Constructor
	public ScaledImageData(final int[] data, final int widthIn,
		final int heightIn, final double voxelSize, final double scalingFactor,
		final double constant, final boolean flipHorizontal,
		final boolean flipVertical, final boolean noFiltering)
	{
		height = heightIn;
		width = widthIn;
		pixelSpacing = voxelSize;
		final int filterSize = 3;
		final int size = width * height;
		final double[] unFiltered = IntStream.range(0, size).mapToDouble(
			i -> data[i]).map(x -> x * scalingFactor).map(d -> d + constant)
			.toArray();
		minimum = Arrays.stream(unFiltered).min().orElse(Double.POSITIVE_INFINITY);
		softScaledImage = medianFilter(unFiltered, width, height, 7); // Median
		if (noFiltering) {
			scaledImage = unFiltered;
		}
		else {
			scaledImage = medianFilter(unFiltered, width, height, filterSize); // Median
		}
		if (flipHorizontal) {// Flip the image around the horizontal axis...
			flipHorizontally();
		}
		if (flipVertical) {// Flip the image around the horizontal axis...
			flipVertically();
		}
	}

	private void flipHorizontally() {
		final long midW = (long) (width / 2.0);
		for (int j = 0; j < height; ++j) {
			final int offset = j * height;
			for (int i = 0; i < midW; ++i) {
				final int sourceIndex = offset + i;
				final int targetIndex = offset - 1 - i;
				scaledImage[targetIndex] = scaledImage[sourceIndex];
				softScaledImage[targetIndex] = softScaledImage[sourceIndex];
			}
		}
	}

	private void flipVertically() {
		final long midH = (long) (height / 2.0);
		for (int j = 0; j < midH; ++j) {
			for (int i = 0; i < width; ++i) {
				final int sourceIndex = j * width + i;
				final int targetIndex = (height - j - 1) * width + i;
				scaledImage[targetIndex] = scaledImage[sourceIndex];
				softScaledImage[targetIndex] = softScaledImage[sourceIndex];
			}
		}
	}

	private double[] medianFilter(final double[] data, final int width,
		final int height, final int filterSize)
	{
		/*Fill filtered with min value to get the frame from messing up with edge detection*/
		final double[] filtered = DoubleStream.generate(() -> minimum).limit(width *
			height).toArray();
		final double[] toMedian = new double[filterSize * filterSize];
		final int noGo = (int) Math.floor(filterSize / 2.0);
		final int median = (int) Math.floor(filterSize * filterSize / 2.0);
		for (int row = noGo; row < height - noGo; row++) {
			for (int col = noGo; col < width - noGo; col++) {
				int newPixel = 0;
				for (int rowOffset = -noGo; rowOffset <= noGo; rowOffset++) {
					for (int colOffset = -noGo; colOffset <= noGo; colOffset++) {
						final int rowTotal = row + rowOffset;
						final int colTotal = col + colOffset;
						toMedian[newPixel] = data[rowTotal * width + colTotal];
						newPixel++;
					}
				}
				Arrays.sort(toMedian);
				filtered[row * width + col] = toMedian[median];
			}
		}
		return filtered;
	}
}
