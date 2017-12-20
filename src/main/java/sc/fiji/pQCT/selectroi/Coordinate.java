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

// TODO Replace with a library class (same as Vector2d)
class Coordinate {

	double ii;
	double jj;

	Coordinate(final double ii, final double jj) {
		this.ii = ii;
		this.jj = jj;
	}

	Coordinate subtract(final Coordinate a) {
		return new Coordinate(ii - a.ii, jj - a.jj);
	}

	double maxVal() {
		return Math.max(Math.abs(ii), Math.abs(jj));
	}
}
