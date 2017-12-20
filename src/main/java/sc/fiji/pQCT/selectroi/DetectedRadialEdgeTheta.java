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

package sc.fiji.pQCT.selectroi;

public class DetectedRadialEdgeTheta implements
	Comparable<DetectedRadialEdgeTheta>
{

	private final double radius;
	int index;

	DetectedRadialEdgeTheta(final double radius, final int index) {
		this.radius = radius;
		this.index = index;
	}

	@Override
	public int compareTo(final DetectedRadialEdgeTheta o) {
		return Double.compare(radius, o.radius);
	}

	@Override
	public boolean equals(final Object o) {
		return o instanceof DetectedRadialEdgeTheta && compareTo(
			(DetectedRadialEdgeTheta) o) == 0;
	}
}
