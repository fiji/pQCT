/*Taken from IvusSnakes ImageJ plugin http://ivussnakes.sourceforge.net/
	Licensed with the GPL 3.0 or above
*/

package sc.fiji.pQCT.selectroi.liveWireEngine;

// TODO Replace with a library class
class Vector2d {

	private double x;
    private double y;

	private Vector2d() {
		x = 0.0;
		y = 0.0;
	}

	Vector2d(final double x, final double y) {
		this.x = x;
		this.y = y;
	}

    private void setX(final double x) {
		this.x = x;
	}

	private void setY(final double y) {
		this.y = y;
	}

	private double getX() {
		return x;
	}

	private double getY() {
		return y;
	}

	// returns the modulus of the vector
    private double mod() {
		return Math.sqrt(x * x + y * y);
	}

	// returns this - other
	Vector2d sub(final Vector2d other) {
		final Vector2d ans = new Vector2d();
		ans.setX(getX() - other.getX());
		ans.setY(getY() - other.getY());
		return ans;
	}

	// returns a unit vector (versor) from this
	Vector2d getUnit() {
		final Vector2d ans = new Vector2d();
		if ((Math.abs(x) < 0.1) && (Math.abs(y) < 0.1)) {
			// if the vector is null
			// we'll return the direction (1,0) for it
			// this is because of a problem that happens when
			// we need to calculate the direction cost
			ans.setX(1.0);
			ans.setY(0.0);

		}
		else {
			ans.setX(getX() / mod());
			ans.setY(getY() / mod());
		}
		return ans;
	}

	double dotProduct(final Vector2d other) {
		// find angle between vectors
		return (getX() * other.getX() + getY() * other.getY());
	}
}
