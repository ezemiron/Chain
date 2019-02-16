//Get list of open images for the dialog menu
n = nImages;
if (n==0){
  Dialog.create("Select images");
  Dialog.addMessage("Please select 16-bit Thresholded and MCNR images");
  Dialog.show();
  input_thresh = File.openDialog("Select 16-bit Thresholded data");
  input_MCNR = File.openDialog("Select MCNR image");
  run("Bio-Formats Importer", "open="+input+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
  n = nImages;
}
else {
  n = nImages;
  list = newArray(n);
  for (i=1; i<=n; i++){
    selectImage(i);
    list[i-1] = getTitle;
  }
  
//Create dialog menu of open images
Dialog.create("Select raw and reconstructed data");
Dialog.addChoice("Thresholded 16-bit reconstruction data", list);
Dialog.addChoice("MCNR image for mask", list);
Dialog.addNumber("Gaussian blur sigma = ", 1.0);
Dialog.addNumber("MCN Threshold value = ", 5.0);
Dialog.show();
  
th = Dialog.getChoice();
MC = Dialog.getChoice();
sigma = Dialog.getNumber();
min_thresh = Dialog.getNumber();

thresh = filename(th);
selectWindow(th);
rename(thresh);
MCNR = filename(MC);
selectWindow(MC);
rename(MCNR);

print(min_thresh);
MCF(thresh, MCNR, min_thresh);
run("Grays");
contrast_stack(thresh+"_MCF");

function MCF (THR, MCN, min_thresh){
	selectWindow(MCN);
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Scale...", "x=2 y=2 z=1.0 width=" + width + " height=" + height + "depth=" + slices + "interpolation=Bicubic average process create title="+MCN+"Scale");
    selectWindow(MCN+"Scale");
    getMinAndMax(min,max);
    setBatchMode(true);
    for (i=1; i<=nSlices; i++) {
       setSlice(i);
       getMinAndMax(min, max);
       setThreshold(min_thresh, max);
       run("Create Mask");
       run("Copy");
       selectWindow(MCN+"Scale");
       run("Paste");
       setSlice(1);
       resetThreshold;}
    setBatchMode(false);
    run("16-bit");
    run("Divide...", "value=255 stack");
    rename("Mask "+THR);
    imageCalculator("Multiply create stack", "Mask "+THR, THR);
    selectWindow("Result of Mask "+THR);
    rename(THR+"_MCF");
    run("Gaussian Blur...", "sigma="+sigma+" stack");
}

function filename (name){
	dot = indexOf(name, ".");
    if (dot>=0) name = substring(name, 0, dot);
    else if (dot == NaN) name = name;
    return name;
}

function substack(thresh, MCNR, array){
  for (i=0; i<channels; i++){
    if (array[i] == 1){
	  selectWindow(thresh);
	  Stack.getDimensions(width, height, channels, slices, frames);
	  run("Make Substack...", "channels="+i+" slices=1-"+slices);
      resetMinAndMax;
      rename("C"+i+"-"+thresh);
      selectWindow(MCNR);
      run("Make Substack...", "channels="+i+" slices=1-"+slices);
      resetMinAndMax;
      rename("C"+i+"-"+MCNR);
}}

function contrast_stack(stack){
  Stack.getDimensions(width, height, channels, slices, frames);
  setSlice(round(slices/2));
  run("Enhance Contrast", "saturated=0.35");
}
  