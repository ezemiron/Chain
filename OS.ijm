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


    if (endsWith(filename, "SIR_THR_ALN.tif")) 
    {

      open(filename);
	SIR=getTitle();
	THR=SIR;
   //   THR=replace(SIR, "(.tif)","_OS-ijm");
      saveAs("TIFF",dir2+THR);     //saves as a new tiff in dir2
    // saveAs("TIFF",dir2);     //saves as a new tiff in dir2
      close(); 


//      closes all open images
 //     while (nImages>0) 
   //       { 
     //     selectImage(nImages); 
       //   close(); 
         // }

//      Saves then closes log
  //  if (isOpen("Log")) {
    //  selectWindow("Log");
      //saveAs("text",dir2+MCMtxt);
      //run("Close");
     // }

    }

}

print("All files have been saved :D")



