
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) {
    showProgress(i+1, list.length);
    filename = dir1 + list[i];
    string=replace(list[i],"(.tif)",".tif"); 
    
    
    if (endsWith(filename, ".tif")) {
    
    
    open(filename);
        filetitle=getTitle();
        run("Duplicate...", "duplicate");
	selectWindow(filetitle);
	run("Duplicate...", "duplicate");
          
          selectWindow(filetitle);
          titleC1=getTitle();
          
          selectWindow(replace(filetitle,"(.tif)","-1.tif"));
	  titleC2=getTitle();
          //titleC2=replace(filetitle,"(.tif)","-1.tif");
          
          selectWindow(replace(filetitle,"(.tif)","-2.tif"));
	  titleC3=getTitle();
          //titleC3=replace(filetitle,"(.tif)","-2.tif");

          run("Merge Channels...", 
          "c1=&titleC1 c2=&titleC2 c3=&titleC3 create");
    
    saveAs("TIFF",dir2+string);     //saves as a new tiff in dir2
    close();
  }
}


