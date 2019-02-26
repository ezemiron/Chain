
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory (Different from Source directory)");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(true);
for (i=0; i<list.length; i++) {
    showProgress(i+1, list.length);
    filename = dir1 + list[i];
    string1=replace(list[i],"(.tif)","_c2-H2B"); 
    string2=replace(list[i],"(.tif)","_c2-DAPI"); 
    
    
    if (endsWith(filename, ".tif")) {
    
    
    open(filename);
        filetitle=getTitle();
        run("Split Channels");
          
          selectWindow("C1-"+filetitle);
          titleC1=getTitle();
          
          selectWindow("C2-"+filetitle);
          titleC2=getTitle();
          
          selectWindow("C3-"+filetitle);
          titleC3=getTitle();

          run("Merge Channels...", 
          "c1=&titleC1 c2=&titleC2 create keep");
	  
	 // savename = replace(filetitle, ".tif","_c2-H2B");
	  saveAs("TIFF",dir2+string1);
	 // close();

	  run("Merge Channels...", 
          "c1=&titleC1 c2=&titleC3 create keep");
	  
	 // savename = replace(filetitle, ".tif","_c2-H2B");
	  saveAs("TIFF",dir2+string2);

//      closes all open images
      while (nImages>0) 
          { 
          selectImage(nImages); 
          close(); 
          }
  }
}


