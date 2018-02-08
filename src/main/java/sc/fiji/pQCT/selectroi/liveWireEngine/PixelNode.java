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

/**
 * Store pixel nodes in a Priority Queue so that Dijkstra can run on O(n log n)
 * The interface Comparable is required so that the Java class PriorityQueue
 * could be used The code is licensed under GPL 3.0 or newer
 */
class PixelNode implements Comparable<PixelNode> {

	private final int[] myIndex;
	private final double myDistance;
	private final int[] whereFrom;

	/**
	 * Constructor
	 * 
	 * @param index the index of the node
	 * @param distance the cost to the node
	 * @param whereFrom from which node we got to this node from
	 */
	PixelNode(final int[] index, final double distance, final int[] whereFrom) {
		myIndex = index;
		myDistance = distance;
		this.whereFrom = whereFrom;
	}

	@Override
	public int compareTo(final PixelNode other) {
		return Double.compare(myDistance, other.getDistance());
	}

	@Override
	public boolean equals(final Object o) {
		return o instanceof PixelNode && compareTo((PixelNode) o) == 0;
	}

	double getDistance() {
		return myDistance;
	}

	int[] getIndex() {
		return myIndex;
	}

	int[] getWhereFrom() {
		return whereFrom;
	}
}
