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

import java.awt.Color;
import java.awt.Component;
import java.awt.Polygon;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.ScrollbarWithLabel;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import sc.fiji.pQCT.selectroi.liveWireEngine.LiveWireCosts;

/**
 * LiveWire ImageJ plug-in modified from ivus snakes
 * (http://ivussnakes.sourceforge.net/) ImageJ plugin Changed the implementation
 * back to the one suggested in Barret &amp; Mortensen 1997. Interactive
 * live-wire boundary extraction. Medical Image Analysis (1996/7) volume 1,
 * number 4, pp 331-341.
 */
public class LiveWirePlugin implements PlugIn, MouseListener,
	MouseMotionListener, KeyListener, AdjustmentListener, MouseWheelListener
{

	private ImageCanvas canvas;
	private ImagePlus imp;
	private ImageWindow imw;
	private ScrollbarWithLabel stackScrollbar;
	private PolygonRoi roi;
	private Polygon polygon;
	private ArrayList<Polygon> polygons;
	private RoiManager rMan;
	private Overlay over;
	private int width;
	private int height;
	private int currentSlice;
	private int depth = -1;
	private LiveWireCosts lwc;

	@Override
	public void adjustmentValueChanged(final AdjustmentEvent e) {
		if (currentSlice != e.getValue()) {
			final int previousSlice = currentSlice;
			currentSlice = e.getValue();
			// Finalize ROI in the previous imp
			imp.setSlice(previousSlice);
			imp.setPosition(previousSlice);
			finalizeRoi();
			imp.setSlice(currentSlice);
			imp.setPosition(currentSlice);
			initLW();
		}
	}

	@Override
	public void keyPressed(final KeyEvent e) {
		if (e.getExtendedKeyCode() == KeyEvent.getExtendedKeyCodeForChar(
			KeyEvent.VK_Q) || e.getKeyChar() == 'q')
		{
			// Shut down the plug-in
			canvas.removeMouseListener(this);
			canvas.removeMouseMotionListener(this);
			canvas.removeKeyListener(this);
			stackScrollbar.removeAdjustmentListener(this);
			imw.removeMouseWheelListener(this);
		}
	}

	@Override
	public void keyReleased(final KeyEvent e) {}

	@Override
	public void keyTyped(final KeyEvent e) {}

	@Override
	public void mouseClicked(final MouseEvent e) {}

	@Override
	public void mouseDragged(final MouseEvent e) {}

	@Override
	public void mouseEntered(final MouseEvent e) {}

	@Override
	public void mouseExited(final MouseEvent e) {}

	@Override
	public void mouseMoved(final MouseEvent e) {
		if (polygon.npoints <= 0) {
			return;
		}
		// Visualize the segment to be added in real-time
		final int screenX = e.getX();
		final int screenY = e.getY();
		final int x = canvas.offScreenX(screenX);
		final int y = canvas.offScreenY(screenY);
		final int[] pX;
		final int[] pY;

		if ((e.getModifiersEx() & InputEvent.SHIFT_MASK) != 0 || (e
			.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0)
		{
			// If shift is pressed, visualize adding a straight line
			pX = new int[polygon.npoints + 1];
			pY = new int[polygon.npoints + 1];
			for (int i = 0; i < polygon.npoints; ++i) {
				pX[i] = polygon.xpoints[i];
				pY[i] = polygon.ypoints[i];
			}
			pX[polygon.npoints] = x;
			pY[polygon.npoints] = y;
		}
		else {
			// Visualize adding livewire segment
			int[][] fromSeedToCursor;
			while ((fromSeedToCursor = lwc.returnPath(x, y)) == null) {}
			pX = new int[polygon.npoints + fromSeedToCursor.length];
			pY = new int[polygon.npoints + fromSeedToCursor.length];
			for (int i = 0; i < polygon.npoints; ++i) {
				pX[i] = polygon.xpoints[i];
				pY[i] = polygon.ypoints[i];
			}
			for (int i = 0; i < fromSeedToCursor.length; ++i) {
				pX[polygon.npoints + i] = fromSeedToCursor[i][0];
				pY[polygon.npoints + i] = fromSeedToCursor[i][1];
			}
		}
		// Add the ROI
		imp.setRoi(new PolygonRoi(pX, pY, pX.length, Roi.POLYLINE), true);
	}

	@Override
	public void mousePressed(final MouseEvent e) {
		if (e.getClickCount() > 1) {
			finalizeRoi();
			init();
		}
	}

	@Override
	public void mouseReleased(final MouseEvent e) {
		if (e.getClickCount() >= 2) {
			// Ignore second and further clicks of a double click
			return;
		}
		final int screenX = e.getX();
		final int screenY = e.getY();
		final int x = canvas.offScreenX(screenX);
		final int y = canvas.offScreenY(screenY);
		// Backpedal polygon to the previous one if control is pressed
		if (!polygons.isEmpty() && ((e.getModifiersEx() &
			InputEvent.CTRL_MASK) != 0 || (e.getModifiersEx() &
				InputEvent.CTRL_DOWN_MASK) != 0))
		{
			// Get the previous polygon
			final Polygon tempP = polygons.get(polygons.size() - 1);
			final int[] pX = new int[tempP.npoints];
			final int[] pY = new int[tempP.npoints];
			for (int i = 0; i < tempP.npoints; ++i) {
				pX[i] = tempP.xpoints[i];
				pY[i] = tempP.ypoints[i];
			}
			polygon = new Polygon(pX, pY, pX.length);
			polygons.remove(polygons.size() - 1); /*Remove the previous polygon*/
			roi = new PolygonRoi(polygon, Roi.POLYLINE);
			imp.setRoi(roi, true);
			lwc.setSeed(pX[pX.length - 1], pY[pX.length - 1]);
		}
		else {
			// Add a new segment to the polygon
			if (polygon.npoints > 0) {
				// Store a copy of the previous polygon
				final int[] tX = new int[polygon.npoints];
				final int[] tY = new int[polygon.npoints];
				for (int i = 0; i < polygon.npoints; ++i) {
					tX[i] = polygon.xpoints[i];
					tY[i] = polygon.ypoints[i];
				}
				// Store the previous polygon
				polygons.add(new Polygon(tX, tY, tX.length));
				// If shift is pressed, add a straight line
				if ((e.getModifiersEx() & InputEvent.SHIFT_MASK) != 0 || (e
					.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0)
				{
					// Add a straight line
					polygon.addPoint(x, y);
				}
				else {
					// Add a livewire segment
					int[][] fromSeedToCursor;
					while ((fromSeedToCursor = lwc.returnPath(x, y)) == null) {}
					final int[] pX = new int[polygon.npoints + fromSeedToCursor.length];
					final int[] pY = new int[polygon.npoints + fromSeedToCursor.length];
					for (int i = 0; i < polygon.npoints; ++i) {
						pX[i] = polygon.xpoints[i];
						pY[i] = polygon.ypoints[i];
					}
					for (int i = 0; i < fromSeedToCursor.length; ++i) {
						pX[polygon.npoints + i] = fromSeedToCursor[i][0];
						pY[polygon.npoints + i] = fromSeedToCursor[i][1];
					}
					polygon = new Polygon(pX, pY, pX.length);
				}
				// Get, and set the ROI
				roi = new PolygonRoi(polygon, Roi.POLYLINE);
				imp.setRoi(roi, true);
			}
			else {
				polygon.addPoint(x, y);
				lwc.setSeed(x, y);
			}
			lwc.setSeed(x, y);
		}
	}

	@Override
	public void mouseWheelMoved(final MouseWheelEvent e) {
		final int rotation = e.getWheelRotation();
		final int rotatedSlice = currentSlice + rotation;
		if (rotatedSlice <= 0 || rotatedSlice > depth ||
			currentSlice == rotatedSlice)
		{
			return;
		}
		final int previousSlice = currentSlice;
		currentSlice += rotation;
		// Finalize ROI in the previous imp
		imp.setSlice(previousSlice);
		imp.setPosition(previousSlice);
		finalizeRoi();
		imp.setSlice(currentSlice);
		imp.setPosition(currentSlice);
		initLW();
	}

	@Override
	public void run(final String arg) {
		IJ.log("Started liveWire");
		imp = WindowManager.getCurrentImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}
		imw = WindowManager.getCurrentWindow();
		canvas = imw.getCanvas();

		if (imp.getImageStackSize() > 1) {
			depth = WindowManager.getCurrentImage().getImageStackSize();
			final Component[] components = imw.getComponents();
			for (final Component component : components) {
				if (component instanceof ScrollbarWithLabel) {
					stackScrollbar = ((ScrollbarWithLabel) component);
					stackScrollbar.addAdjustmentListener(this);
					imw.addMouseWheelListener(this);
					currentSlice = stackScrollbar.getValue();
					break;
				}
			}
		}

		width = imp.getWidth();
		height = imp.getHeight();

		IJ.log("Init LW " + width + " h " + height);
		initLW();

		// Pop up Roi Manager
		rMan = RoiManager.getInstance();
		if (RoiManager.getInstance() == null) {
			rMan = new RoiManager();
		}

		IJ.log("Add listeners");
		canvas.addMouseListener(this);
		canvas.addMouseMotionListener(this);
		canvas.addKeyListener(this);
	}

	private void finalizeRoi() {
		if (polygon.npoints <= 2) {
			return;
		}
		// Do not remove the last point, simply connect last point, and initial
		// point
		polygon.addPoint(polygon.xpoints[0], polygon.ypoints[0]);
		/*Create the ROI*/
		roi = new PolygonRoi(polygon, Roi.POLYGON);
		// Set roi color to differentiate ROIs from each other
		final int colorInd = over.size();
		final double angle = 2.0 * Math.PI * colorInd / 10.0;
		//@formatter:off
		final double[] colors = {
				0.5 + 0.5 * Math.sin(2.0 * Math.PI * (colorInd - 5) / 10.0), // R
				0.5 + 0.5 * Math.cos(angle), // G
				0.5 + 0.5 * Math.sin(angle) // B
		};
		//@formatter:on
		roi.setStrokeColor(new Color((float) colors[0], (float) colors[1],
			(float) colors[2]));
		// Add the roi to an overlay, and set the overlay active
		imp.setRoi(roi, true);
		over.add(roi);
		// Add the segmented area to the roiManager
		rMan.addRoi(roi);
	}

	/** Used to reset livewire when switching to another slice in a stack */
	private void initLW() {
		final double[][] pixels = new double[width][height];
		final short[] tempPointer = (short[]) imp.getProcessor().getPixels();
		for (int r = 0; r < height; ++r) {
			for (int c = 0; c < width; ++c) {
				pixels[c][r] = tempPointer[c + r * width];
			}
		}
		lwc = new LiveWireCosts(pixels);
		init();
	}

	/**
	 * Used to reset the polygon, and polygon list used to keep current polygon,
	 * and the history of the current polygon
	 */
	protected void init() {
		polygons = new ArrayList<>();
		polygon = new Polygon();

		over = imp.getOverlay();
		if (imp.getOverlay() == null) {
			over = new Overlay();
			imp.setOverlay(over);
		}
	}
}
