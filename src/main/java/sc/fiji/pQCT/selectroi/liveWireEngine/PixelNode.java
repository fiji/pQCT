
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

	public int compareTo(final PixelNode other) {
		return Double.compare(myDistance, other.getDistance());
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
