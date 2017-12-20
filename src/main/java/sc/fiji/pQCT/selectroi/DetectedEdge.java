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

import java.util.Vector;

public class DetectedEdge implements Comparable<DetectedEdge> {

	public final Vector<Integer> iit; // indexes for x-coordinates
	public final Vector<Integer> jiit; // indexes for y-coordinates
	public final int area;
	public final int length;

	DetectedEdge(final Vector<Integer> iit, final Vector<Integer> jiit, final int area) {
		this.iit = iit;
		this.jiit = jiit;
		length = iit.size();
		this.area = area;
	}

	@Override
	public boolean equals(final Object o) {
		return o instanceof DetectedEdge && compareTo((DetectedEdge) o) == 0;
	}

	@Override
	public int compareTo(final DetectedEdge o) {
	    return Double.compare(area, o.area);
	}
}
