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

dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) 
{
    showProgress(i+1, list.length);
    filename = dir1 + list[i];
    string=replace(list[i],"(.dv)",""); 


    if (endsWith(filename, "SIR.dv")) 
    {

      run("Bio-Formats Importer", "open="+filename+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

      getDimensions(w,h,chan,sli,fra);

      SIR=getTitle();

      RAW=replace(SIR,"_SIR.dv",".dv");
  //runs the SIMcheck 16bit mode threshold (16bit):
      run("Threshold and 16-bit Conversion", "auto-scale");
      THR=replace(SIR, "SIR.dv","SIR_THR");
      saveAs("TIFF",dir2+THR);     //saves as a new tiff in dir2
      THR=getTitle();
      close(); 

  //Open raw .dv
     run("Bio-Formats Importer", "open="+dir1+RAW+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
  //remove negatives if any
      run("Min...", "value=0 stack");
  //make the SIMcheck modulation contrast (32bit)
      run("Modulation Contrast Map", "calculate_mcnr_from_raw_data=&RAW camera_bit_depth=15 or,_specify_mcnr_stack=None reconstructed_data_stack=&SIR");
      
      MCM=getTitle();
      MCMtxt=replace(MCM, "_MCM","_MCNR-log.txt");
      saveAs("TIFF",dir2+MCM);

      MCN=replace(RAW, ".dv","_MCN");
      close(); 
      selectImage(RAW);
      close();
      selectImage(MCN);

  //saves as a new tiff in dir2
        saveAs("TIFF",dir2+MCN);

  //run("Threshold...");
  //    setThreshold(5.0000, 100.0000);
  //    run("NaN Background", "stack");
  //make it a binary mask (8bit)
  //    run("Convert to Mask", "method=Default background=Dark");

  //remove overlays if present (any overlay will stop the scale function):
  //    run("Remove Overlay");
  //make the same size
  //    run("Scale...","x=2 y=2 z=1.0 width=&w height=&h interpolation=Bicubic average create title=scaled");
 
      //saves as a new tiff in dir2
  //    saveAs("TIFF",dir2+MCNmask);

//      closes all open images
      while (nImages>0) 
          { 
          selectImage(nImages); 
          close(); 
          }

//      Saves then closes log
    if (isOpen("Log")) {
      selectWindow("Log");
      saveAs("text",dir2+MCMtxt);
      run("Close");
      }

    }

}

print("All files have been cleaned :D")