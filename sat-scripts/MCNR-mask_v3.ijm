//## Copyright (C) 2016 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
//##
//## This program is free software: you can redistribute it and/or modify
//## it under the terms of the GNU General Public License as published by
//## the Free Software Foundation, either version 3 of the License, or
//## (at your option) any later version.
//##
//## This program is distributed in the hope that it will be useful,
//## but WITHOUT ANY WARRANTY; without even the implied warranty of
//## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//## GNU General Public License for more details.
//##
//## You should have received a copy of the GNU General Public License
//## along with this program.  If not, see <http://www.gnu.org/licenses/>.

dir1 = getDirectory("Choose Directory ");
//dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) 
{
    showProgress(i+1, list.length);
    filename = dir1 + list[i];
    string=replace(list[i],"(.tif)",""); 


    if (endsWith(filename, "MCN_ALN.dv")) 
    {

      run("Bio-Formats Importer", "open="+filename+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

      getDimensions(w,h,chan,sli,fra);

      MCNmask=getTitle();

      MCNmask=replace(MCNmask,"_MCN_ALN.dv","_ALN_MCNR-mask");
  
  //remove overlays if present (any overlay will stop the scale function):
      run("Remove Overlay");
  //make the same size
      run("Scale...","x=2 y=2 z=1.0 width=&w height=&h interpolation=Bicubic average create title=scaled");  
  //run("Threshold...");
      setThreshold(5.0000, 100.0000);
      run("NaN Background", "stack");
  //make it a binary mask (8bit)
      run("Convert to Mask", "method=Default background=Dark");



      //saves as a new tiff in dir1
      saveAs("TIFF",dir1+MCNmask);

//      closes all open images
      while (nImages>0) 
          { 
          selectImage(nImages); 
          close(); 
          }
    }
}

print("All files have been cleaned :D")
