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

package sc.fiji.pQCT.selectroi.liveWireEngine;

import static java.util.Arrays.stream;

import java.util.Arrays;
import java.util.PriorityQueue;

/**
 * Modified by Timo Rantalainen 2012 - 2014 from IvusSnakes
 * (http://ivussnakes.sourceforge.net/) ImageJ plugin A Class to calculate
 * LiveWire paths.
 * <p>
 * Changed the implementation back to the one suggested in Barret &amp;
 * Mortensen 1997. Interactive live-wire boundary extraction. Medical Image
 * Analysis (1996/7) volume 1, number 4, pp 331-341.
 * </p>
 */
public class LiveWireCosts implements Runnable {

	private static final int[][] LAPLACIAN_NEIGHBOUR_INDICES = { { -1, -1 }, { -1,
		0 }, { -1, 1 }, { 0, 1 }, { 1, 1 }, { 1, 0 }, { 1, -1 }, { 0, -1 } };
	private final double[][] imagePixels; // stores Pixels from original image
	private final PriorityQueue<PixelNode> pixelCosts;
	private final double[][] gradientRows; // stores image gradient modulus
	private final double[][] gradientColumns; // stores image gradient modulus
	// it is oriented: X = LEFT TO RIGHT
	// Y = UP TO DOWN
	private final double[][] gradientr; // stores image gradient RESULTANT modulus
	private final double[][] laplacian;
	private final int[][][] whereFrom; // stores where from path started
	// stores whether the nodes were marked or not
	private final boolean[][] visited;
	private final int rows;
	private final int columns;
	private final double gw;// Gradient Magnitude Weight
	private final double dw;// Gradient Direction Weight
	private final double zw;// Binary Laplacian Weight
	private int sr; // seed x and seed y, weight zero for this point
	private int sc;
	private int tr; // thread x and y passed as parameters
	private int tc;
	private Thread myThread;
	private boolean myThreadRuns;// flag for thread state

	/**
	 * Constructor
	 *
	 * @param imagePixels 2D gray scale image in
	 */
	// initializes Dijkstra with the image
	public LiveWireCosts(final double[][] imagePixels) {

		// initializes weights for edge cost taken from Barret 1997
		// these are default values
		gw = 0.43;
		zw = 0.43;
		dw = 0.13;
		// initializes all other matrices
		rows = imagePixels.length;
		columns = imagePixels[0].length;
		this.imagePixels = imagePixels;
		pixelCosts = new PriorityQueue<>();
		whereFrom = new int[rows][columns][2];
		visited = new boolean[rows][columns];
		gradientRows = new double[rows][columns];
		gradientColumns = new double[rows][columns];
		gradientr = new double[rows][columns];
		initGradient();
		laplacian = new double[rows][columns];
		initLaplacian();
	}

	/**
	 * Returns the path from seed point to point r,c
	 *
	 * @param r x-coordinate of the target
	 * @param c y-coordinate of the target
	 * @return m x 2 array of the m-length path with x-, and y-coordinates
	 */
	public int[][] returnPath(final int r, final int c) {
		// returns the path given mouse position

		final int[][] pathCoordinates = new int[rows * columns][];

		if (!visited[r][c]) {
			// TODO refactor so that method blocks (waits) while !visited
			// attempt to get path before creating it
			// this might occur because of the thread
			return null;
		}
		int length = 0;
		int myr = r;
		int myc = c;
		pathCoordinates[length] = new int[] { r, c };
		do { // while we haven't found the seed
			++length;
			myr = whereFrom[myr][myc][0];
			myc = whereFrom[myr][myc][1];
			pathCoordinates[length][0] = myr;
			pathCoordinates[length][1] = myc;
		}
		while (!(myr == sr && myc == sc));

		// path is from last point to first
		// we need to invert it
		final int[][] pathToReturn = new int[length + 1][2];
		for (int i = 0; i <= length; i++) {
			pathToReturn[i][0] = pathCoordinates[length - i][0];
			pathToReturn[i][1] = pathCoordinates[length - i][1];
		}
		return pathToReturn;
	}

	/*set the seed point to start Dijkstra
	@param r x-coordinate of the seed
	@param c y-coordinate of the seed
	*/
	public void setSeed(final int r, final int c) {
		myThreadRuns = false;
		if (myThread != null) {
			try {
				myThread.join();
			}
			catch (final InterruptedException ignored) {}
		}
		tr = r;
		tc = c;
		myThreadRuns = true;
		myThread = new Thread(this);
		myThread.start();
	}

	/** Implement the Runnable interface */
	public void run() {
		// runs set point in parallel
		final int r = tr;
		final int c = tc;
		int[] nextIndex;
		sr = r;
		sc = c;

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				visited[i][j] = false;
			}
		}

		visited[r][c] = true;
		int[] coordinates = { r, c };
		whereFrom[r][c] = coordinates;

		updateCosts(r, c, 0);

		while ((pixelCosts.peek() != null) && (myThreadRuns)) {
			nextIndex = pixelCosts.peek().getIndex();
			whereFrom[nextIndex[0]][nextIndex[1]] = pixelCosts.peek().getWhereFrom();
			updateCosts(nextIndex[0], nextIndex[1], pixelCosts.peek().getDistance());

			// removes pixels that are already visited and went to the queue
			while (true) {
				if (pixelCosts.peek() == null) break;
				coordinates = pixelCosts.peek().getIndex();
				if (!visited[coordinates[0]][coordinates[1]]) break;
				pixelCosts.poll();
			}
		}
		pixelCosts.clear();
	}

	/**
	 * Check neighbours for Laplacian zero-crossing
	 * 
	 * @param tempLap mask of zero-crossings
	 * @param neighbourhood the neighbourhood to check
	 * @param centre the coordinates of the centre pixel to check the
	 *          neighbourhood for
	 * @return tempLap the zero-crossings mask
	 */
	private double[][] checkNeighbours(final double[][] tempLap,
		final int[][] neighbourhood, final int[] centre)
	{
		for (final int[] coordinates : neighbourhood) {
			final double atCoordinates = laplacian[coordinates[0]][coordinates[1]];
			final double atCenter = laplacian[centre[0]][centre[1]];
			if (Math.signum(atCenter) != Math.signum(atCoordinates) && Math.abs(
				atCenter) < Math.abs(atCoordinates))
			{
				/*zero-crossing detected, change to 0 to disable laplacian*/
				tempLap[centre[0]][centre[1]] = 0;
			}
		}
		return tempLap;
	}

	// returns the edge cost of going from sx,sy to dx,dy
	private double edgeCost(final int sr, final int sc, final int dr,
		final int dc)
	{
		// fg is the Gradient Magnitude
		double edgeCostSum = gw * gradientr[dr][dc] + zw * laplacian[dr][dc];
		edgeCostSum += edgeDirectionCost(sr, sc, dr, dc) * dw;
		return edgeCostSum;
	}

	/**
	 * Calculates edge direction cost
	 * 
	 * @param sr source pixel x-coordinate
	 * @param sc source pixel y-coordinate
	 * @param dr destination pixel x-coordinate
	 * @param dc destination pixel y-coordinate
	 * @return edgeDirectionCostValue edge direction cost
	 */
	private double edgeDirectionCost(final int sr, final int sc, final int dr,
		final int dc)
	{
		final Vector2d Dp = (new Vector2d(gradientRows[sr][sc],
			-gradientColumns[sr][sc])).getUnit();
		final Vector2d Dq = (new Vector2d(gradientRows[dr][dc],
			-gradientColumns[dr][dc])).getUnit();
		final Vector2d p = new Vector2d(sr, sc);
		final Vector2d q = new Vector2d(dr, dc);
		final Vector2d L;
		if (Dp.dotProduct(q.sub(p)) >= 0) {
			L = q.sub(p).getUnit();
		}
		else {
			L = p.sub(q).getUnit();
		}
		/*Barret 1996/1997 eqs 3 & 4*/
		return 2.0 / (3.0 * Math.PI) * (Math.acos(Dp.dotProduct(L)) + Math.acos(L
			.dotProduct(Dq)));
	}

	// initializes gradient image
	private void initGradient() {
		/*
		Using sobel
		for gx convolutes the following matrix
		
		|-1 0 1|
		Gx = |-2 0 2|
		|-1 0 1|
		*/
		for (int i = 1; i < rows - 1; ++i) {
			for (int j = 1; j < columns - 1; ++j) {
				gradientRows[i][j] = -1 * (imagePixels[i - 1][j - 1]) + 1 *
					(imagePixels[i + 1][j - 1]) - 2 * (imagePixels[i - 1][j]) + 2 *
						(imagePixels[i + 1][j]) - 1 * (imagePixels[i - 1][j + 1]) + 1 *
							(imagePixels[i + 1][j + 1]);
			}
		}

		// for gy convolutes the following matrix
		//
		// |-1 -2 -1|
		// Gy = | 0 0 0|
		// |+1 +2 +1|
		//
		for (int i = 1; i < rows - 1; ++i) {
			for (int j = 1; j < columns - 1; ++j) {
				gradientColumns[i][j] = -1 * (imagePixels[i - 1][j - 1]) + 1 *
					(imagePixels[i - 1][j + 1]) - 2 * (imagePixels[i][j - 1]) + 2 *
						(imagePixels[i][j + 1]) - 1 * (imagePixels[i + 1][j - 1]) + 1 *
							(imagePixels[i + 1][j + 1]);
			}
		}
		for (int i = 1; i < rows - 1; i++) {
			for (int j = 1; j < columns - 1; j++) {
				gradientr[i][j] = Math.sqrt(gradientRows[i][j] * gradientRows[i][j] +
					gradientColumns[i][j] * gradientColumns[i][j]);
			}
		}

		final double grMax = stream(gradientr).flatMapToDouble(Arrays::stream).max()
			.orElse(Double.NEGATIVE_INFINITY);
		for (int i = 0; i < gradientr.length; ++i) {
			for (int j = 0; j < gradientr[i].length; ++j) {
				gradientr[i][j] = 1.0 - gradientr[i][j] / grMax;
			}
		}
	}

	/*initializes laplacian image zero-crossings. Marks zero-crossings with 0, otherwise the value is 1*/
	private void initLaplacian() {

		// Using finite differences convolute
		// @formatter:off
        final double[][] laplacianKernel = {
		        { 0, 1, 0 },
                { 1, -4, 1 },
                { 0, 1, 0 }
		};
        // @formatter:on
		for (int i = 1; i < rows - 1; i++) {
			for (int j = 1; j < columns - 1; j++) {
				laplacian[i][j] = 0;
				for (int j2 = -1; j2 <= 1; ++j2) {
					for (int i2 = -1; i2 <= 1; ++i2) {
						laplacian[i][j] += imagePixels[i + i2][j + j2] *
							laplacianKernel[i2 + 1][j2 + 1];
					}
				}
			}
		}

		/*Search for zero crossing to binarize the result*/
		double[][] tempLap = new double[rows][columns];
		for (final double[] row : tempLap) {
			Arrays.fill(row, 1.0);
		}
		/*Check pixel neighbourhoods for zero-crossings*/
		final int[][] neighbourhood = new int[8][2]; // 8 connected neighbourhood
		for (int i = 1; i < rows - 1; i++) {
			for (int j = 1; j < columns - 1; j++) {
				tempLap[i][j] = 1;
				if (laplacian[i][j] == 0) {
					// No need to check neighbours
					tempLap[i][j] = 0;
				}
				else {
					// Check 8-connected neighbour
					for (int th = 0; th < 8; ++th) {
						neighbourhood[th][0] = i + LAPLACIAN_NEIGHBOUR_INDICES[th][0];
						neighbourhood[th][1] = j + LAPLACIAN_NEIGHBOUR_INDICES[th][1];
					}
					tempLap = checkNeighbours(tempLap, neighbourhood, new int[] { i, j });
				}
			}
		}
		/*OverWrite Laplacian*/
		for (int i = 0; i < rows; i++) {
			System.arraycopy(tempLap[i], 0, laplacian[i], 0, columns);
		}
	}

	/*updates Costs and Paths for a given point
		calculated over 8 directions N, NE, E, SE, S, SW, W, NW
		@param r target pixel x-coordinate
		@param c target pixel y-coordinate
	*/
	private void updateCosts(final int r, final int c, final double cost) {

		pixelCosts.poll();
		final int[][] neighbourhood = new int[8][2]; // 8 connected neighbourhood
		final int[][] neighbourIndices = { { -1, -1 }, { -1, 0 }, { -1, 1 }, { 0,
			1 }, { 1, 1 }, { 1, 0 }, { 1, -1 }, { 0, -1 } };
		// Check 8-connected neighbour
		for (int th = 0; th < 8; ++th) {
			neighbourhood[th][0] = r + neighbourIndices[th][0];
			neighbourhood[th][1] = c + neighbourIndices[th][1];
		}
		for (final int[] coordinates : neighbourhood) {
			final int x = coordinates[0];
			final int y = coordinates[1];
			if (x < 0 || x >= rows || y < 0 || y >= columns) {
				continue;
			}
			final int[] fromCoords = { r, c };
			final int[] pixelCoords = { x, y };
			final double distance = cost + edgeCost(r, c, x, y);
			pixelCosts.add(new PixelNode(pixelCoords, distance, fromCoords));
		}
		visited[r][c] = true;
	}

}
