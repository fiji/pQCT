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
