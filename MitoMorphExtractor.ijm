// This macro gets information about the mitochondria in each cell in each image and outputs a csv file
// To make these measurements per cell it needs corresponding mask images contained in the same folder (AutoCBfolder)
// the outputfile is a .txt file which can be read into R as a csv file

run("Fresh Start");
dir = getDir("Choose directory where all autoCB pictures are"); //SO NOT MAX Z
flist = getFileList(dir);

// Create a file to store the logdata in
resultsFile = File.open(dir+"MitoMorphResults_sigma03_new.txt");
resultsPath = dir+"MitoMorphResults_sigma03_new.txt";



setBatchMode(true);
setOption("CopyHeaders", false); 

for (i=0; i<(flist.length-1); i++){
	
	real_i = i;
	// Extract filenames from folder
	suff = substring(flist[i], lengthOf(flist[i])-8);
	while (matches(suff, "_[0][0-9]+.tif")==0 && i<(flist.length-1)){
		i++;
		suff = substring(flist[i], lengthOf(flist[i])-8);
	}
	
	pre = replace(flist[i], suff, "");				// extracts the shared namepart of the images
	fname = pre+suff; 								// creates a variable for the image
	fmask = replace(fname, ".tif", "mask.tiff"); 	// creates a variable for the corresponding mask
	print(fname);

	// Open file
	open(dir+"/"+fname);
	
	// Duplicate the mitochondrial channel (C4) and apply blur 
	// CHANGE NUMBER (duplicate channel= ... ) IF DIFFERENT
	selectWindow(fname);
	run("Duplicate...", "duplicate channels=4");
	C4 = getImageID();
	close("\\Others");
	// Even though the tubeness filter also does a blurring, the results somehow look best when I apply the blurring here first
	run("Gaussian Blur 3D...", "x=2 y=2 z=1");
	
	print(nSlices);
	// Analyze tubeness of the mitochondria and threshold the result
	tube = replace(fname, ".tif", "_tubenew.tif"); 
	if (File.exists(dir+tube) == 0){
		run("Tubeness", "sigma=0.3 use");
		save(dir+tube);
	}
	else {
		open(dir+tube);
	}
	setAutoThreshold("Moments dark no-reset stack");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Minimum background=Dark black");
	thresholdedTubeness = getImageID();
	mitoThresholded = replace(fname, ".tif", "_ThreshMito.tif");
	save(dir+mitoThresholded);

	// Skeletonize the thresholded image 
	run("Skeletonize (2D/3D)");
	skel = replace(fname, ".tif", "_skel.tif"); 
	save(dir+skel);

	// Convert original image to 8bit because that is what analyze skeleton needs later on
	selectImage(C4);
	run("8-bit");
	eightbit = replace(fname, ".tif", "_8bit.tif");		
	save(dir+eightbit);
	open(skel);
	
	// Convert mask to ROIS
	roiManager("reset");
	
	open(fmask);
	run("Duplicate...", " "); // duplicate to extract the rois
	setThreshold(1, 255);
	run("Analyze Particles...", "size=1000-Infinity add"); // this adds them to roi manager
	
	
	// Measure skeletons in each ROI 
	// The ROIS generated by cellprofiler are numbered by the intensity in the grayscale "mask"
	// measure all intensities add to array
	// sort array
	selectImage(fmask);
	roinrs = roiManager("count");
	arr = newArray(roinrs);
	for (r=0; r<roinrs; r++){
		roiManager("select", r);
		run("Measure");
		String.copyResults;
		intens = String.paste;
		intens = split(intens); 
		intensity = parseInt(intens[2]);
		print(intensity);
		run("Clear Results");
		arr[r]=intensity;
		if(r==roinrs){
			print("ERROR roinrs");
			break;
		}
	}
	sorted = Array.sort(arr);
	// loop again over the rois and over the sorted array
	selectImage(fmask);
	for (n=0; n<roinrs; n++){
		s = 0;
		roiManager("select", n);
		run("Measure");
		String.copyResults;
		Intens = String.paste;
		Intens = split(Intens); 
		Intensity = parseInt(Intens[2]);
		print("LOWER");
		print(Intensity);
		while (Intensity != sorted[s]){
			s += 1;
			if(s==roinrs){ 
			s=0;
			break;
			}
		}
		changeValues((Intensity-1),(Intensity+1), (s+1));
		run("Clear Results");
	}
	Roi.remove
	fmask2 = replace(fname, ".tif", "mask2.tiff"); 	// creates a variable for the corresponding mask
	save(dir+fmask2);


	for (j=0; j<roinrs; j++){
		
		// select roi and retrieve the numbering from CP
		open(dir+fmask2);
		roiManager("select", j);
		run("Measure");
		String.copyResults;
		intens = String.paste;
		intens = split(intens); 
		CP_number = parseInt(intens[2]);
		if (CP_number >= 253){
			CP_number=roinrs;
		}
		print(CP_number);
		
		unique_CP_Cell_ID = fname+"__Cell"+CP_number+"__"; 
		run("Clear Results");
		
		// Create a single mask with values 0 and 1
		newImage("currentROI-"+j+"", "8-bit black", 2048,2048,1);
		roiManager("select", j);
		run("Fill", "slice");
		changeValues(1,255, 1);
		
		// Multiply stack with this single mask 
		selectImage(skel);
		imageCalculator("multiply stack create", skel, "currentROI-"+j+"");
		singleSkel = getImageID();
	
		// Analyze skeleton in ROI
		open(eightbit);
		selectImage(singleSkel);
		run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] calculate show display original_image="+eightbit+"");


		// Get headers into results file
		if(real_i==0 && j==0){ 
		resHeadings = String.getResultsHeadings;
		resHeadings = replace(resHeadings, "# \\b" , "nr_"); 			// replace # with Nr-
		resHeadings = replace(resHeadings, "\\b \\b","_");				// replace spaces between two word boundaries with _
		resHeadings = "Mito_ID"+resHeadings;
		resHeadings = split(resHeadings, " \t\n\r"); 
		nr_columns = resHeadings.length;
		resHeadings = String.join(resHeadings);
		File.append(resHeadings+"\n", resultsPath);					
		}
		
		// Export mitochondrial measurements
		String.copyResults;
		res = String.paste;
		res = replace(res, "\\b \\b","_");				// replace spaces in headings with _
		res = split(res, " \t\n\r"); 
		
		for (k=0; k<(res.length-nr_columns); k+=nr_columns){
			// add a unique identifier to each mito ID so that we know which image and cell its in
			mito_ID = res[k];
			res[k] = unique_CP_Cell_ID+"Mito"+mito_ID;
			// add results line by line to the results file
			final_col = res[k+nr_columns];
			res[k+(nr_columns-1)] = final_col+"\n";
			row = Array.slice(res, k, (k+(nr_columns-1)));
			File.append(String.join(row), resultsPath);
		};

		// Remove the results from results window
		run("Clear Results");
		// Close all images except skel
		selectImage(skel);
		close("\\Others");
	};
	
	
	roiManager("delete");	// removes the rois of this image
	close("*"); 			// closes all open images
};


exit;






	