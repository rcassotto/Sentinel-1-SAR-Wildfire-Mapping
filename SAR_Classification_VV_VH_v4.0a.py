#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : {Ryan Cassotto}
# Created Date: 2025/02/21
# version ='4.0a'
#
## Version Notes: This corrects the Sigma0 normalization process by:
    ## 1) checking for overlap between scenes (to ensure normalization for same regions when multiple orbital passes are used)
    ## 2) Uses the same reference image date for VH as for VV. 
# ---------------------------------------------------------------------------
""" Sentinel 1 Classification change detection module """  
# ---------------------------------------------------------------------------
from image_proc_module import Image_proc
import Cluster
import rasterio
import glob
import numpy as np
import os
from rasterio.plot import reshape_as_raster
import argparse,ast
import time
from FIREDpy_ARIA_module_v1 import ARIA
import shutil
import inspect
import re
from termcolor import colored

# ---------------------------------------------------------------------------
# making the output directory
def mk_outdirectory(outpath_full):
    if not os.path.exists(outpath_full):
        print("Making output directory: ", outpath_full)
        os.makedirs(outpath_full)
    return


def convert_noData2NaNs(in_image):
   ### Remove no data values (i.e. values in which Sigma0 data == 0) 
    band_shape = np.shape(in_image)
    b1_flat = in_image.flatten()
    i_nan = np.argwhere(b1_flat == 0)  
    b1_flat[i_nan] = np.nan
    b1_raster = b1_flat.reshape(band_shape)
    return b1_raster
    
def write_luminosity_2_geotiff(tiff_outname, image_data, luminosity):
    profile_1band = image_data.profile.copy()  # copy geotiff meta data from input file   
    lum_out_raster = luminosity
    with rasterio.open(tiff_outname, 'w', **profile_1band) as dst:
        dst.write(lum_out_raster,1)
    return
    


def write_rgb_2_geotiff(tiff_outname, image_data, rgb):
    profile_3band = image_data.profile.copy()  # copy geotiff meta data from input file   
    profile_3band.update(count = 3)  # Change count from 1-band (input file) to 3-bands (output RGB image)
    rgb_out_raster = reshape_as_raster(rgb)  # reshape to output correctly as a raster (necessary to create RGB-style image)
    with rasterio.open(tiff_outname, 'w', **profile_3band) as dst:
        dst.write(rgb_out_raster)
    return


def Get_image_boundaries(vh_sigma_file): 
    image_data = rasterio.open(vh_sigma_file)
    lon = (image_data.bounds.left, image_data.bounds.right)
    lat = (image_data.bounds.top, image_data.bounds.bottom)
    minLat = np.min(lat)
    maxLat = np.max(lat)
    minLon = np.min(lon)
    maxLon = np.max(lon)
    geographic_bounds = [minLon, maxLon, minLat, maxLat]
    # image_data.close
    return geographic_bounds, image_data


def normalize_sigma0_images(sigma0_images, sigma0_dir, sigma_ref_img_fname): 

    ## Read in reference image
    print("Reading in Reference Luminance Image: ", sigma_ref_img_fname)
    ref_image_path_and_fname = os.path.join(sigma0_dir, sigma_ref_img_fname)
    ref_image = Image_proc(ref_image_path_and_fname)
    ref_image.read()
    
    
    ## Get image corner coordinates
    image_data = rasterio.open(ref_image_path_and_fname)
    lon = (image_data.bounds.left, image_data.bounds.right)
    lat = (image_data.bounds.top, image_data.bounds.bottom)
    minLat = np.min(lat)
    maxLat = np.max(lat)
    minLon = np.min(lon)
    maxLon = np.max(lon)
    geographic_bounds = [minLon, maxLon, minLat, maxLat]
    profile = image_data.profile.copy()  # copy geotiff meta data from input file   

    print("Geographic Bounds: ", geographic_bounds)
    
    ## Loop over subsequent luminance images
    for n in range(0, len(sigma0_images)):
        nFile_2Norm = sigma0_images[n]
        
        if nFile_2Norm == ref_image_path_and_fname:
            print(nFile_2Norm, "is the reference image. Skipping this file.")
            ref_norm_name_and_path = ref_image_path_and_fname.replace('.tif','_norm.tif')
            shutil.copy(ref_image_path_and_fname, ref_norm_name_and_path)
            continue
        else:  ## Perform normalization
            ## Read in secondary image
            print("Normalizing iteration", n, "of ", len(sigma0_images))
            print(sigma0_images[n])
          
            ## Next up: call to process aria
            sec_path, sec_fname = os.path.split(nFile_2Norm)
            
            png_outfilename_and_path = sigma0_dir + '/'+  sigma_ref_img_fname + "_" + sec_fname + ".png"
            check_file = os.path.isfile(png_outfilename_and_path)
            if check_file == False:
            
                sec_image = Image_proc(nFile_2Norm)
                sec_image.read()
                # print("Secondary Image:" , type(sec_image))
      
                aria = ARIA(geographic_bounds)
                print(png_outfilename_and_path)
                t_threshold=0.1
                aria_data = aria.process_ARIA(ref_image, [sec_image], png_outfilename_and_path, t=t_threshold, map_type=aria.simple_map, file_prefix='', show_hists=True)
                ### t is the threshold for coh change. If differences are greater than t, they will be highlighted, if not 0 is returned.
        
                #### Rev 4a update: save the normalized secondary image, as a geotiff
                for inmap in aria_data:
                    band = inmap.map # In Rev 4a, this is the normalized secondary coherence image
                    t_band_outname = nFile_2Norm.replace('.tif','_norm.tif')
                
                    print("Writing Normalized Secondary image ", t_band_outname)
                    with rasterio.open(t_band_outname, 'w', **profile) as dst:
                        dst.write(band,1)
                
                    ## Write normalized_tBand meta data file
                    write_metafile(t_band_outname, sigma0_dir, sigma_ref_img_fname, sec_path, sec_fname, sigma0_dir)
                    print("")
            else:
                print(nFile_2Norm, " already normalized. Moving on....")
                pass
                
    # ref_norm_name_and_path = ref_image_path_and_fname.replace('.tif','_norm.tif')
    # shutil.copy(ref_image_path_and_fname, ref_norm_name_and_path)
    return


def write_metafile(t_band_outname, infile_0_path, infile_0_name, infile_1_path, infile_1_name, outdirectory):
        Line1 = ( "Normalized Sigma0 Ouput Path + Directory: " + t_band_outname + "\n")
        Line2 = ( "Reference Image filename: " + infile_0_name + "\n")
        Line3 = ( "Refefence Image Path: " + infile_0_path + "\n")
        Line4 = ( "Secondary Image filename: " + infile_1_name + "\n")
        Line5 = ( "Secondary Image Path: " + infile_1_path + "\n")
        
        metafile_outname = t_band_outname.replace('.tif','_log.txt') # infile basename
        print('Writing Log file: ', metafile_outname)
        meta_file = open (metafile_outname, "w")
        meta_file.write(Line1)
        meta_file.write(Line2)
        meta_file.write(Line3)
        meta_file.write(Line4)
        meta_file.write(Line5)
        meta_file.close()
        return


def Run_classification_algorithm(Sigma_files,out_dir, use_norm_images):
    
    ncount=0
    for in_sigma_file in Sigma_files:
        # ---------------------------------------------------------------------------
        ## construct output directory name, make output directory     
        granule = in_sigma_file.split('/')[-1] # granule is a string with the original Sigma0 filename (e.g. S1A_IW_GRDH_1SDV_20150921T232856_20150921T232918_007820_00AE3B_E36F_Sigma0_VV)
        mk_outdirectory(out_dir)  # make output directory#
        
        ncount+=1
        time_start=time.time()
        print("Processing Classification Scene: ", ncount, " of ", len(Sigma_files), ":", granule)
     
        # ---------------------------------------------------------------------------
        # Read in images
        VV_image = Image_proc(in_sigma_file)
        VV_image.read()
        # print("VV Image type:", type(VV_image))       
       
        vh_sigma_file=in_sigma_file.replace("VV", "VH")
        VH_image = Image_proc(vh_sigma_file)
        VH_image.read()
        
        # ---------------------------------------------------------------------------
        ## geographic boundaries and data on one of the image files
        geographic_bounds, image_data = Get_image_boundaries(vh_sigma_file)
        
 
        # ---------------------------------------------------------------------------
        # Create RGB Image        
        sar_class = Cluster.cluster_image(geographic_bounds)
        rgb = sar_class.construct_rgb_composite(VH_image, VV_image)
        
        # ---------------------------------------------------------------------------
        ## Convert RGB to luminosity
        ### account for no data values of 1 in bands 1 and 2 
        b1 = rgb[:,:,0]
        b2 = rgb[:,:,1]
        b1_noData_filt = np.where(b1 == 1., np.nan, b1) # convert noData Values of 1 to 0
        b2_noData_filt = np.where(b2 == 1., np.nan, b2) # convert noData Values of 1 to 0
        luminosity = b1_noData_filt * 0.2126 + b2_noData_filt * 0.7152 + rgb[:,:,2]*0.0722
        
        luminosity_FiltNoData = np.where(luminosity > 3E38, 0, luminosity) # convert no data value (really high values) to 0
        
    
        # ---------------------------------------------------------------------------        
        ## Create luminosity and 3-band polarimetric 'false-color' output filenames
        if use_norm_images == "true":
            lum_tiff_outname = os.path.join(out_dir,granule.replace('_Sigma0_VV_norm.tif','_RGBLumin.tif'))
            falseColor_outname = os.path.join(out_dir,granule.replace('_Sigma0_VV_norm.tif','_FalseColorClass.tif'))

        else:
            lum_tiff_outname = os.path.join(out_dir,granule.replace('_Sigma0_VV.tif','_RGBLumin.tif'))
            falseColor_outname = os.path.join(out_dir,granule.replace('_Sigma0_VV.tif','_FalseColorClass.tif'))
   
        # ---------------------------------------------------------------------------        
        ## Write single-band luminosity and 3-band polarimetric 'false-color' outputs as geotiffs
        write_luminosity_2_geotiff(lum_tiff_outname, image_data, luminosity_FiltNoData)
        write_rgb_2_geotiff(falseColor_outname, image_data, rgb)            
        
        
        time_stop=time.time()
        print("Done. Elapsed Time: ", time_stop - time_start, " seconds."); print("")        
    return

def determine_reference_sigma0_image(sigma0_images):
    print("Determing Reference Image for Normalization")
    min_cube = np.zeros(len(sigma0_images))
    max_cube = np.zeros(len(sigma0_images))
    mean_cube = np.zeros(len(sigma0_images))
    
    for n in range(0, len(sigma0_images)):           
        print("Reading in Reference Image Coherence file: ", sigma0_images[n])
    
        ## Read in image
        lum_in = rasterio.open(sigma0_images[n])
        lum_image = lum_in.read(1)
        lum_image[lum_image==0] = np.nan # convert zeros to nans
        lum_image[lum_image>1E38] = np.nan ## check for images with other no Data values
        
    
        mean_cube[n] = np.nanmean(lum_image)
        max_cube[n] = np.nanmax(lum_image)
        min_cube[n] = np.nanmin(lum_image)
    
    print("mean values:" , mean_cube)
    print("max values:" , max_cube)
    print("min values:" , min_cube)
    
    print("Max (mean) value:", np.max(mean_cube))
    imax = int(np.argwhere(mean_cube == np.max(mean_cube)))
    ref_image_fname_and_path=sigma0_images[imax]  
    ref_image_path, sigma0_ref_img_fname = os.path.split(ref_image_fname_and_path)
    print("ref_image_fname:", sigma0_ref_img_fname)
    
    return sigma0_ref_img_fname


## Parse arguments from input file
def  parse_input_file(p):
    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)        

    outdirectory = args['out_dir']
    sigma0_dir = args['sigma_file_loc']
    sigma0_ref_img_fname_VV = args['sigma_ref_image']
    use_norm_images = args['use_norm_images']

    return outdirectory, sigma0_dir, sigma0_ref_img_fname_VV, use_norm_images

## Get list of variables containing a string pattern (e.g. Sigma_files_VV_track_*)
def get_matching_local_variables(pattern):
    frame = inspect.currentframe().f_back
    local_vars = frame.f_locals
    matching_vars = [name for name in local_vars if re.search(pattern, str(name))]
    return matching_vars

def find_matching_globals(pattern):
  matching_globals = {}
  for name, value in globals().items():
    if re.match(pattern, name):
      matching_globals[name] = value
  return matching_globals

## dynamically create variable using global()
def create_variable(name, value):
    globals()[name] = value


## check for consistent overlap over region
def check_for_consistent_overlap(Sigma_files_VV):
## Check for consistent overlap in scenes / subregions (this accounts for voids 
## in adjacent path/frames generated by cropping/merging scenes in previous steps)
##  Check overlap for scenes in Sigma0 directory. This parses scenes from adjacent orbital tracks (i.e. path/frames) to maintain
## consistency during the normalization process. Only one polarity is necessary for this check. 
## This approach uses no data values and creates mask for each scene, then sums non-nan pixels

    nCount=0 # initiate counter
    while len(Sigma_files_VV) > 1:
        
        ## Select first image as the reference image (to compare all others to)
        in_image_ref_and_path = Sigma_files_VV[0]
        Sigma_files_VV.remove(in_image_ref_and_path) # remove ref image from list
        Sigma_files_VV_tmp = [in_image_ref_and_path] ## create a list for current coincident path/frame pairings, staring with first

        ## Read in reference image, create a binary mask of valid pixels and sum number of valid pixels       
        in_image_ref_in = rasterio.open(in_image_ref_and_path)
        in_image_ref = in_image_ref_in.read(1)
        print("REf image:", in_image_ref_and_path); print("")
        in_image_ref_masked = np.where(in_image_ref != in_image_ref_in.nodata, 1, 0) # create a binary mask of in/valid pixels for reference image
        in_image_ref_in.close

        i_valid_ref = np.sum(in_image_ref_masked) ## sum number of valid pixels in ref image
        
        ## Loop over all other 'secondary' images to find overlapping scenes
        for n in range(0, len(Sigma_files_VV)):
            sec_img_and_path = Sigma_files_VV[n]
            print(sec_img_and_path)
            in_image_sec_in = rasterio.open(sec_img_and_path)
            in_image_sec = in_image_sec_in.read(1)
            in_image_sec_masked = np.where(in_image_sec != in_image_sec_in.nodata, 1, 0) # create a mask of valid/invalid pixels
            in_image_sec_in.close
            
            ## Sum ref and sec masks to find overlapping pixels
            images_summed = in_image_ref_masked + in_image_sec_masked  # sum reference and secondary images
            overlap_test = np.sum (np.where(images_summed == 2, 1, 0) ) # sum number of overlapping pixels
            prcnt_overlap = overlap_test / i_valid_ref * 100 # Calculate percent overlap
            print("Percent overlap: ", prcnt_overlap)
            print("")
            
            ## save filename and path into track one path
            if prcnt_overlap > 90:
                Sigma_files_VV_tmp.append(sec_img_and_path)
        
        ## Create a list of files with coincident path/frames; i.e. copy 'Sigma_files_VV_tmp' to unique variable name (Sigma_files_VV_track_#)
        list_name = f"Sigma_files_VV_track_{nCount}"
        create_variable(list_name, Sigma_files_VV_tmp)  
    
        ## remove duplicate files and repeat
        for n in range(1, len(Sigma_files_VV_tmp)):
            Sigma_files_VV.remove(Sigma_files_VV_tmp[n])
    
        ## Reset / increment key variables for next iteration
        nCount+=1
        del Sigma_files_VV_tmp
    
    ## Find global variables with Sigma_files_VV_track_#
    unique_sigma_file_tracks = find_matching_globals("Sigma_files_VV_track_*")
    print("")
    
    return unique_sigma_file_tracks


if __name__ ==  "__main__":
    # ---------------------------------------------------------------------------
    #Implement Arg parse to pull relevant input parameters from input.txt file
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',type=argparse.FileType('r'))
    p = parser.parse_args()
    (outdirectory, sigma0_dir, sigma0_ref_img_fname_VV, use_norm_images) = parse_input_file(p) 

    ## check whether initial sigma0_ref_img_fname_VV is specified
    if len(sigma0_ref_img_fname_VV)==0:
        initial_sigma0_ref_flag=1
    else:
        initial_sigma0_ref_flag=0
    print(colored(("Initial ref flag :......", initial_sigma0_ref_flag), 'green', attrs=['bold']))
    # print(colored(("Normalizing subregion ", (n+1) , " of ", len(unique_sigma_file_tracks)),'blue', attrs=['bold']))


    ## Get list of input files    
    Sigma_files_VV = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV.tif'))   
    Sigma_files_VH = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VH.tif'))   
    

    # ---------------------------------------------------------------------------
    ##  Check overlap for scenes in Sigma0 directory.     
    unique_sigma_file_tracks = check_for_consistent_overlap(Sigma_files_VV)


    ## Iterate over each list seperately to normalize different track/frames
    for n in range(0, len(unique_sigma_file_tracks)):
        print(colored('Normalizing subregion ' + str(n+1) + ' of ' +  str(len(unique_sigma_file_tracks)),'grey', attrs=['bold']))
        
        list_name = f"Sigma_files_VV_track_{n}"
        print("list name:", list_name)  # unique list_name
        in_Sigma_files_VV = unique_sigma_file_tracks[list_name]

        
        ### Make an equivalent list for VH files
        in_Sigma_files_VH=[]
        for nn in range(0, len(in_Sigma_files_VV)):    
            in_Sigma_files_VH.append(in_Sigma_files_VV[nn].replace('VV','VH'))
        print("")

        
    
        # ---------------------------------------------------------------------------
        ## Identify reference image, if one is not specified 
        #### Need to be careful, especially for big areas of interest. This assumes all scenes have the same overlap
        
        # if (len(sigma0_ref_img_fname_VV)==0 and initial_sigma0_ref_flag==1):
        if initial_sigma0_ref_flag==1:

            ### Loop to read in all images and find image with greatest mean; set this as reference iamge. 
            ## VV
            sigma0_ref_img_fname = determine_reference_sigma0_image(in_Sigma_files_VV)
            sigma0_ref_img_fname_VV = sigma0_ref_img_fname
            print(colored('Normalizing scene subset to Sigma0 image: ' + sigma0_ref_img_fname_VV,'blue')); print("")
          
            del sigma0_ref_img_fname
       
            ## VH
            sigma0_ref_img_fname_VH = sigma0_ref_img_fname_VV.replace('VV','VH') # RKC mod, 20240221 to ensure same image dates are used for reference
        
        else:
            sigma0_ref_img_fname_VH = sigma0_ref_img_fname_VV.replace('VV', 'VH')
            
           
        # ---------------------------------------------------------------------------
        ## New section to normalize ALL Sigma0 images. Skips this step if normalized Sigma0 images already exist.
        if use_norm_images == 'true':       
          ### Normalize all Sigma0 - VV images to the specified or determined reference image
          print("Normaling all Sigma0 - VV images to: ", sigma0_ref_img_fname_VV)
          normalize_sigma0_images(in_Sigma_files_VV, sigma0_dir, sigma0_ref_img_fname_VV)
          print("Done normalizing secondary VV images."); print("")
    
          ### Normalize all Sigma0 - VH images to the specified or determined reference image
          print("Now normaling all Sigma0 - VH images to: ", sigma0_ref_img_fname_VH)
          normalize_sigma0_images(in_Sigma_files_VH, sigma0_dir, sigma0_ref_img_fname_VH)
          print("Done normalizing secondary VH images."); print("")
    
          ### Create new list of *_VV_norm.tif
          Sigma_files = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV_norm.tif'))
        else:
          Sigma_files = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV.tif'))
        
        
        print("Sigma0 Files:", Sigma_files)
        
        # ---------------------------------------------------------------------------
        ## Create FalseColor Polarimetric Images and Luminance Images
        Run_classification_algorithm(Sigma_files,outdirectory, use_norm_images)  
    
    