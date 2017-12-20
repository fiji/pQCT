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

	ImageJ density distribution analysis plugin
    Copyright (C) 2011 Timo Rantalainen
*/

package sc.fiji.pQCT;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.plugin.filter.Info;
import ij.text.TextPanel;

import java.util.StringTokenizer;

public class Export_Header implements PlugIn {
	String imageInfo;
	String site;
	double percent;
	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null)
			return;
		if (imp.getType() != ImagePlus.GRAY16){
			IJ.error("Distribution analysis expects 16-bit greyscale data");
			return;
		}
		imageInfo = new Info().getImageInfo(imp,imp.getChannelProcessor());
		String imageName;
		if (getInfoProperty(imageInfo,"File Name")!= null){
			imageName = getInfoProperty(imageInfo,"File Name");
		}else{
			if(imp.getImageStackSize() == 1){
				imageName = imp.getTitle();
				imageInfo+="File Name:"+imageName+"\n";
			}else{
				imageName = imageInfo.substring(0,imageInfo.indexOf("\n"));
				imageInfo+="File Name:"+imageName+"\n";
			}
		}
		
		//Try to get CT file additional header info
		String filePath = new String("");
		if (getInfoProperty(imageInfo,"File Path")!= null){
			filePath = getInfoProperty(imageInfo,"File Path");
		}
		String cName = "C"+imageName.substring(1,imageName.length());
		//IJ.log(imageInfo);
		
		//IJ.log(filePath+"\\"+cName);
		try{
			
			CTHeaderReader ct = new CTHeaderReader(filePath+cName);
			ct.read();
			//IJ.log(ct.site+" "+ct.percent);
			site = ct.site;
			percent = ct.percent;
		}catch(Exception err){
			IJ.log("CT header "+err.toString());
			site = new String("Unknown");
			percent = Double.NaN;
		}
		
		
		TextPanel textPanel = IJ.getTextPanel();
		if (textPanel == null) {textPanel = new TextPanel();}
		if (textPanel.getLineCount() == 0){writeHeader(textPanel);}
		
		String results = "";
		results = printResults(results, imp);
		textPanel.appendLine(results);
		textPanel.updateDisplay();
	}

	
	void writeHeader(TextPanel textPanel){
		String[] propertyNames = {"File Name","File Path","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","Object Length","Measurement Mask","Percent Length"};
		String headings = "";
		for (int i = 0;i<propertyNames.length;++i){
			headings+=propertyNames[i]+"\t";
		}
		textPanel.setColumnHeadings(headings);
	}

	String printResults(String results, ImagePlus imp){
		String[] propertyNames = {"File Name","File Path","Patient's Name","Patient ID","Patient's Birth Date","Acquisition Date","Pixel Spacing","ObjLen"};
		if (imp != null){
			if (getInfoProperty(imageInfo,"File Name")!= null){
				results+=getInfoProperty(imageInfo,"File Name")+"\t";
			}else{
				if(imp.getImageStackSize() == 1){
					results+=getInfoProperty(imageInfo,"Title")+"\t";
				}else{
					results+=imageInfo.substring(0,imageInfo.indexOf("\n"))+"\t";
				}
			}
			for (int i = 1;i<propertyNames.length;++i){
				results+=getInfoProperty(imageInfo,propertyNames[i])+"\t";
			}
			results+=site+"\t";
			results+=percent+"\t";
		}
		return results;
	}
	
	String getInfoProperty(String properties,String propertyToGet){
		String toTokenize = properties;
		StringTokenizer st = new StringTokenizer(toTokenize,"\n");
		String currentToken = null;
		while (st.hasMoreTokens() ) {
			currentToken = st.nextToken();
			if (currentToken.indexOf(propertyToGet) != -1){break;}
		}
		if (currentToken.indexOf(propertyToGet) != -1){
			String[] returnValue = currentToken.split(":",2);
			/*
			StringTokenizer st2 = new StringTokenizer(currentToken,":");
			String token2 = null;
			while (st2.hasMoreTokens()){
				token2 = st2.nextToken();
			}
			return token2.trim();
			*/
			return returnValue[1].trim();
		}
		return null;
	}
}
