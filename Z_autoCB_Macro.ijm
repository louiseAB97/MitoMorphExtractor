// Louise Buijs 
// This macro processes a folder with vsi files obtained by the Olympus confocal spinning disc microscope
// i) creates a folder with max intensity Z-projections
// ii) creates a folder with auto contrasted and brightness tiffs 


run("Fresh Start");
dir = getDir("Choose directory where all pictures from your condition are");
flist = getFileList(dir);

File.makeDirectory(dir + "MaxZprojections");
Zfolder = dir + "MaxZprojections/";

File.makeDirectory(dir + "AutoCBtiffs");
ACBfolder = dir + "AutoCBtiffs/";

setBatchMode(true);

// color acquisition options
// *** To add other options: 
// 1. click on record in fiji 
// 2. Plugins>BioFormats>BioFormatsImporter,  and import an image using your own settings
// 3. Look in the recorder and copy the complete string after "open=[filename] "
BGRF = "color_mode=Custom rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1 series_0_channel_0_red=0 series_0_channel_0_green=0 series_0_channel_0_blue=255 series_0_channel_1_red=0 series_0_channel_1_green=255 series_0_channel_1_blue=0 series_0_channel_2_red=255 series_0_channel_2_green=0 series_0_channel_2_blue=0 series_0_channel_3_red=255 series_0_channel_3_green=0 series_0_channel_3_blue=255";
GRFB = "color_mode=Custom rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1 series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=255 series_0_channel_1_green=0 series_0_channel_1_blue=0 series_0_channel_2_red=255 series_0_channel_2_green=0 series_0_channel_2_blue=255 series_0_channel_3_red=0 series_0_channel_3_green=0 series_0_channel_3_blue=255";


for (i=0; i<flist.length; i++){
	if ( endsWith(flist[i], ".vsi") ){
		pStr = dir + flist[i];
		// Change depending on your channelsettings ***
		bfSuffix = BGRF;
		bfStr = "open=["+pStr+"] "+bfSuffix;

		run("Bio-Formats Importer", bfStr);

		
		// Apply autocontrast
		for (c=0; c<5; c++){
			Stack.setChannel(c);
			run("Enhance Contrast", "saturated=0.35");
		}
		
		// Save autocontrast tiffs in separate folder
		ACB = ACBfolder + replace(flist[i], ".vsi", ".tif");
		saveAs("Tiff", ACB);


		// Save Z projections in separate folder 
		run("Z Project...", "projection=[Max Intensity]");
		sStr = Zfolder + replace(flist[i], ".vsi", ".tif");
		saveAs("Tiff", sStr);
		
		run("Close All");

	}
	
}








	