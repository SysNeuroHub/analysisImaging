# Script to perform conversion of images from stereo to DMD
Stereo2DMD: wrapper to convert images defined in stereotaxic coordinates to images suited for DMD projection
(showCCFindividual): create DMD images for individual cortical areas (from Kim 2023 Neuron)

# Function to run during imaging experiment
imupdatepair(referenceImage, imageFolder): compare reference image with a latest image in imageFolder

# Scripts to create projection images in Stereotaxic coordinates:
showAllenCCFBregmaLambda_patches: square patches on CCF 
(showAllenCCFBregmaLambda_circles: circles on CCF)
(showAllenCCFBreagmaLambda: square patches on CCF (obsolete))
(showAllenCCFBreagmaLambda_star: stars on CCF)
showNatural: natural movies
showStar: create a png with a giant star and a cross hair for a reference image (star_800x500.png)

# Functions to convert Stereotaxic images to individual brains
applyStereo2DMD(images, bregma, MmPerPixel_img, mrangle, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir):  transforms images defined in stereotaxic coordinates to DMD space (500x800 pix)
saveEveryImages(imageStack, saveDir): save every images of a 3D image stack using screen2png. Resulting images are suited for storing in polyScan GUI
exportPng4DMD(savename, fig): export a png image with the original image resolution in fig using screen2png (called in saveEveryImages)

# Utility functions to convert from stereotaxic coordinate to DMD
registerStereo2CCF(bregma, MmPerPixel, path2Brain_template, usFactor): computes 2D transformation from stereotaxi coordinate to CCF looked from above, with upsampling factor usFactor (called in applyStereo2DMD)
getSurfaceData(mr_brain, direction, threshold, smoothSigma): Extracts a smooth surface map and intensity values from a 3D volume, using only voxel-based operations (no isosurface). (called in registerStereo2CCF)
Stereo2T2(imageStereo, usFactor, Vusinfo, tform3, surfDepth, mrangle): projectionmaps images defined in stereotaxic coordinates into CCF then to T2 space (called in applyStreo2DMD)
paintSurfaceToVolume(surfDepth, surfData, volSize): create 3D volume assigning surface data to detected surface voxels. (called in Stereo2T2)
