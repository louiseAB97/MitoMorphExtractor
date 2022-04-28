# MitoMorphExtractor
Detailed description of how to use the complete pipeline



# Image acquisition prerequisites
- Keep acquisition settings and number of slices constant
- How I found the correct range of slices was by looking at the mitochondrial channel and first setting the bottom slice at at where the slice where the mitochondrial signal becomes faint after it being sharp and maximized in terms of spread, then I moved upwards (image a pyramid shape) and set the top slice at the a priori determined number of slices. I found that with cos7 cells 16 is usually good.
- IMPORTANT make sure the images are labeled as follows:  "DAPI_greenfluorophore_redfluorophor_farredfluorophor_###".   (where ### represents digits)
- you can put all experimental conditions in one folder


# Creating cell masks
- It is advisable preprocess the raw images with the Z_autoCB_Macro and obtain the cell masks using (an adapted version) of the CellProfiler pipeline provided in this repository (detailed instructions can be found within the pipeline itself)
In short: This pipeline uses the MaxZprojections to build a cell mask and creates measurements of intensity for single cells.
- If NOT using my pipeline: make sure the cellmasks are labeled as such: "DAPI_greenfluorophore_redfluorophor_farredfluorophor_###mask"
- IMPORTANT make sure the mask images are stored in the same folder as the images to be analyzed with the MitoMorphExtractor


# Running the MitoMorphExtractor
- Start up FIJI
- Open the MitoMorphExtractor.ijm file
- By default the pipeline assumes mitochondria are in channel 4, if this is not the case you need change the number in the code manually (just change the number)
- click on run and select the folder with your images and masks in there



# Parsing the output 
- I provided some example code in R to process the output file but feel free to use your own


# Notes on runtime
- Currently the most timeconsuming step is running the tubeness filter, this is done per image so making sure there is multiple cells within one image to save time
- Running a folder with 20 images takes about an hour on my macbook pro M1 
- Alternatively an implementation using CLIJx (GPU based) is provided in this repository. However there was a subtle difference in the output of the CLIJ tubeness and the ImageJ tubeness and I found the latter to be working better and this has solely to do with finding the right sigma, so feel free to tweak this
- Also you need a strong computer for the CLIJ version but it could be potentially powerful.
