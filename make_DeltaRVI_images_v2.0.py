#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : {Ryan Cassotto}
# Created Date: 2025/01/15
# version ='2.0'
#
# version 2.0 makes the following changes
#   - removes seconds from the image acquisition time stamp between coherence and Classification/Sigma0 image.
# ---------------------------------------------------------------------------
""" Python Module to calculate luminance change images based on coherence change images (i.e. coh change image dates) """
# ---------------------------------------------------------------------------
from image_proc_module import Image_proc
# from FIREDpy_coherence_module.FIREDpy_ARIA_module_v1 import ARIA
from FIREDpy_ARIA_module_v1 import ARIA

#import Cluster
import rasterio
import glob
import numpy as np
import os
#from rasterio.plot import reshape_as_raster
import argparse
import ast
import shutil
import socket


from collections import Counter
from datetime import datetime


## Hard Set parameters
lum_tDelta_cutoff_secs = 60
hostname = socket.gethostname()


# ---------------------------------------------------------------------------
# making the output directory
def mk_outdirectory(outpath_full):
    if not os.path.exists(outpath_full):
        print("Making output directory: ", outpath_full)
        os.makedirs(outpath_full)
    return


def filter_non_unique(img_date_sort):
    non_unique_dates = [item for item, count in Counter(
        img_date_sort).items() if count > 1]
    print(non_unique_dates)
    return non_unique_dates



def read_coh_image_dates(in_coh_file_and_path):
    in_coh_path, in_coh_fname = os.path.split(in_coh_file_and_path)
    fname_split = in_coh_fname.split( "_" )
    ref_image_date = fname_split[0]
    image1_date = fname_split[1]
    image2_date = fname_split[2]         
    return ref_image_date, image1_date, image2_date


def read_lumin_dates(luminance_images):
   # Get dates of classification images, save as a list
    indates_lumin = []
    for n in range(0, len(luminance_images)):
        infile_path, infile_name = os.path.split(luminance_images[n])
        infile_split = infile_name.split("_")
        infile_date = infile_split[4]
        indates_lumin.append(infile_date)
        # print(infile_date)

    return indates_lumin


def difference_luminance_images(lum1_path_and_name, lum2_path_and_name, in1_fname, in2_fname, output_directory, indates_lumin, i_lumin1, i_lumin2):
    print("Reading in First Luminance file: ", in1_fname)
    ref_in = rasterio.open(lum1_path_and_name)
    ref_image = ref_in.read(1)
    image_profile = ref_in.profile.copy()

    print("Reading in Second Luminance file: ", in2_fname)
    sec_in = rasterio.open(lum2_path_and_name)
    sec_image = sec_in.read(1)

    ### New 20240614 - Turn no data values to nans
    ref_image[ref_image>1E38] = np.nan ## check for images with other no Data values
    sec_image[sec_image>1E38] = np.nan ## check for images with other no Data values


    print("Differencing luminance images: ", in2_fname, " from ", in1_fname)
    # deltaLum = np.subtract(ref_image, sec_image)  ## original; looks for decreases in coherence only (ignores coherence increase)
    deltaLum = np.subtract(sec_image, ref_image)  ## secondary - reference. the lower the number (more negative) indicates greater change, greater fire loss)
    diff_eq = 'sec_image - ref_image'

  ## Write out differenced image
    tmpName = indates_lumin[i_lumin1] + "_" + indates_lumin[i_lumin2] + "_" + lum_suffix
    deltaLum_outfilename = tmpName.replace("rvi","DeltaRVI")
    print("DeltaLum outfilename:", deltaLum_outfilename)
    write_out_geotiff(deltaLum_outfilename, output_directory, image_profile, deltaLum)
    
    ## Write deltaLum log file
    deltaLum_logfile = deltaLum_outfilename.replace('.tif','_log.txt')
    deltaLum_logfile_and_path = os.path.join(output_directory, deltaLum_logfile)
    write_ImageDelta_metafile(deltaLum_logfile_and_path, in1_path, in1_fname, in2_path, in2_fname, diff_eq)
    return



def write_out_geotiff(outfilename, output_directory, image_profile, deltaLum):
    deltaLum_product_outfilename_and_path = os.path.join(output_directory, outfilename)
    print("Writing Delta Luminance Results to: ", deltaLum_product_outfilename_and_path)
    with rasterio.open(deltaLum_product_outfilename_and_path, 'w', **image_profile) as dst:
        dst.write(deltaLum,1)
    return

def write_metafile(t_band_outname, infile_0_path, infile_0_name, infile_1_path, infile_1_name, outdirectory):
        Line1 = ( "Normalized Lum Image Ouput Path + Directory: " + t_band_outname + "\n")
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

def normalize_luminance_images(luminance_images, lum_dir, lum_ref_img_fname): 
    ## Read in reference image
    print("Reading in Reference Luminance Image: ", lum_ref_img_fname)
    ref_image_path_and_fname = os.path.join(lum_dir, lum_ref_img_fname)
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
    for n in range(0, len(luminance_images)):
        nFile_2Norm = luminance_images[n]
        
        if nFile_2Norm == ref_image_path_and_fname:
            print(nFile_2Norm, "is the reference image. Skipping this file.")
            continue
        else:  ## Perform normalization
            ## Read in secondary image
            print("Normalizing iteration", n, "of ", len(luminance_images))
            print(luminance_images[n])
          
            ## Next up: call to process aria
            sec_path, sec_fname = os.path.split(nFile_2Norm)
            
            png_outfilename_and_path = lum_dir + '/'+  lum_ref_img_fname + "_" + sec_fname + ".png"
            check_file = os.path.isfile(png_outfilename_and_path)
            if check_file == False:
            
                sec_image = Image_proc(nFile_2Norm)
                sec_image.read()
                print("Secondary Image:" , type(sec_image))
                
                # sec_image_in = rasterio.open(nFile_2Norm)
                # sec_image = sec_image_in.read(1)
                # sec_image[sec_image>1E38] = np.nan ## UPDATE: 20240614; check for images with other no Data values
                # sec_image[sec_image==0] = np.nan # convert zeros to nans


                ## Next up: call to process aria
                #sec_path, sec_fname = os.path.split(nFile_2Norm)

                aria = ARIA(geographic_bounds)
           # png_outfilename_and_path = lum_dir + '/'+  lum_ref_img_fname + "_" + sec_fname + ".png"
                print(png_outfilename_and_path)
                t_threshold=0.1
                aria_data = aria.process_ARIA(ref_image, [sec_image], png_outfilename_and_path, t=t_threshold, map_type=aria.simple_map, file_prefix='', show_hists=True)
                ### t is the threshold for coh change. If differences are greater than t, they will be highlighted, if not 0 is returned.
        
                #### Rev 4a update: save the normalized secondary image, as a geotiff
                for inmap in aria_data:
                    band = inmap.map # In Rev 4a, this is the normalized secondary coherence image
                    print(type(band),len(band), np.size(band))
                    t_band_outname = nFile_2Norm.replace('.tif','_norm.tif')
                
                    print("Writing Normalized Secondary image ", t_band_outname)
                    with rasterio.open(t_band_outname, 'w', **profile) as dst:
                        dst.write(band,1)
                
                    ## Write normalized_tBand meta data file
                    write_metafile(t_band_outname, lum_dir, lum_ref_img_fname, sec_path, sec_fname, lum_dir)
                    print("")
            else:
                print(nFile_2Norm, " already normalized. Moving on....")
                pass
                
    ref_norm_name_and_path = ref_image_path_and_fname.replace('.tif','_norm.tif')
    shutil.copy(ref_image_path_and_fname, ref_norm_name_and_path)
    return

def write_ImageDelta_metafile(deltaLum_logfile_and_path, in1_path, in1_fname, in2_path, in2_fname, diff_eq):
    Line1 = ( "Delta Luminance log file and path: " + deltaLum_logfile_and_path + "\n")
    Line2 = ( "Lum Image 1 Input filename (minuend): " + in1_fname + "\n")
    Line3 = ( "Lum Image 1 Input Path: " + in1_path + "\n")
    Line4 = ( "Lum Image 2 Input filename (subtrehend): " + in2_fname + "\n")
    Line5 = ( "Lum Image 2 Input Path: " + in2_path + "\n")
    Line6 = ( "DeltaLum Equation: " + diff_eq + "\n")    
    Line7 = ( "Processing Machine: " + hostname + "\n")      
                         
    print('Writing Log file: ', deltaLum_logfile_and_path)
    meta_file = open (deltaLum_logfile_and_path, "w")
    meta_file.write(Line1)
    meta_file.write(Line2)
    meta_file.write(Line3)
    meta_file.write(Line4)
    meta_file.write(Line5)
    meta_file.write(Line6)
    meta_file.write(Line7)    
    meta_file.close()
    return
     

    

def determine_reference_luminance_image(luminance_images):
    print("Determing Reference Image for Normalization")
    min_cube = np.zeros(len(luminance_images))
    max_cube = np.zeros(len(luminance_images))
    mean_cube = np.zeros(len(luminance_images))
    
    for n in range(0, len(luminance_images)):           
        print("Reading in Reference Image Coherence file: ", luminance_images[n])
    
        ## Read in image
        lum_in = rasterio.open(luminance_images[n])
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
    print("index for reference image:", imax, type(imax))
    ref_image_fname_and_path=luminance_images[imax]  
    ref_image_path, lum_ref_img_fname = os.path.split(ref_image_fname_and_path)
    print("ref_image_fname:", lum_ref_img_fname)
    
    return lum_ref_img_fname

def find_closest_time(datearray, ref_date):
    closest_datetime = min(datearray, key=lambda x: abs(x - ref_date))
    return closest_datetime

def parse_input_data(p):
    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)
    
    lum_dir = args['rvi_file_loc']  
    output_directory = args['out_dir']
    lum_suffix = args['rvi_suffix']
    coh_change_dir = args['coh_change_file_loc']
    coh_change_suffix = args['coh_change_suffix']
    lum_ref_img_fname = args['rvi_ref_image']
    use_norm_images = args['use_norm_images']
    return lum_dir, output_directory, lum_suffix, coh_change_dir, coh_change_suffix, lum_ref_img_fname, use_norm_images


if __name__ == "__main__":
    # ---------------------------------------------------------------------------
    # Implement Arg parse to pull relevant input parameters from input.txt file
    # Use Argparge to get information from command line
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=argparse.FileType('r'))
    p = parser.parse_args()
    (lum_dir, output_directory, lum_suffix, coh_change_dir, coh_change_suffix, lum_ref_img_fname, use_norm_images) = parse_input_data(p)
    
    # make cropped2ROI output directory
    mk_outdirectory(output_directory)

    # Get list of coh Change input files (e.g. 20210212_20210224_20210308_cohChange.tif)
    infile_coh_str = "*" + coh_change_suffix
    coh_images = glob.glob(os.path.join(coh_change_dir, infile_coh_str))

    # Get list of luminance images (e.g. S1B_IW_GRDH_1SDV_20210401_RGBLumin_100m.tif)
    infile_lum_str = "*" + lum_suffix
    luminance_images = glob.glob(os.path.join(lum_dir, infile_lum_str)) # list of luminance images (path and filename)


    ## New section to normalize ALL luminance images to the user-specified reference image. 
    if use_norm_images == 'true':    
        
        ## Determine luminance reference image, if not specified by user in input file
        if len(lum_ref_img_fname)==0:
            ### Loop to read in all images and find image with greatest mean; set this as reference iamge. 
            lum_ref_img_fname = determine_reference_luminance_image(luminance_images)
            print("Luminance reference image:", lum_ref_img_fname)    
        
        ### Normalize all luminance images to the specified or determined reference image
        print("Normaling all luminance images to: ", lum_ref_img_fname)
        normalize_luminance_images(luminance_images, lum_dir, lum_ref_img_fname)
        print("Done normalizing secondary images."); print("")


    ## Create new list of luminance images, based on whether images are normalized or not. 
    if use_norm_images == 'true':    
        ## Calculate lumDelta images based on coherence images (i.e. loop over coh images)
        infile_lum_str_norm = infile_lum_str.replace('.tif','_norm.tif')
        print("Now calculating delta Luminance on normalized luminace images: ", infile_lum_str_norm)
        luminance_images_norm = glob.glob(os.path.join(lum_dir, infile_lum_str_norm)) # list of luminance images (path and filename)
    else:
        print("Now calculating delta Luminance on ORIGINAL luminace images: ")
        luminance_images_norm = luminance_images # copy of the orginal luminance images (path and filename)
   
    
    
    print("Normalized Luminance Image List:")
    for n in range(0, len(luminance_images_norm)):
        print(luminance_images_norm[n])
    print("Original Luminance Image List:")
    for n in range(0, len(luminance_images)):
        print(luminance_images[n])
    
    ## Read dates of normalized luminance images
    indates_lumin = read_lumin_dates(luminance_images_norm)
    print("Luminance Dates:", indates_lumin);
    print("luminance dates type:", type(indates_lumin)); print("")
    
    
    for n in range(0, len(coh_images)):
        print(n)
        
        ## read image dates for current coh change image 
        ref_image_date, image1_date, image2_date = read_coh_image_dates(coh_images[n]) # this includes filename and path
        print("infile:", coh_images[n])
        print("Coh ref date:", ref_image_date) #, type(ref_image_date))
        print("Coh image1 date:", image1_date) #, type(image1_date))
        print("Coh image2 date:", image2_date) #, type(image2_date))
        
        ## Convert luminance date and time strings to Python datetime format
        indates_lumin_datetime=[]
        for i in indates_lumin:
            tmp_timestamp = datetime.fromisoformat(i) #, format).date()
            indates_lumin_datetime.append(tmp_timestamp)   # list of luminance date and time strings in python datetime format
        
        print("indates lumin datime:" ,indates_lumin_datetime)
        print(type(indates_lumin_datetime))

        ## Convert coherence date and time strings to Python datetime format and find closest datetimes
        image1_datetime = datetime.fromisoformat(image1_date)
        image2_datetime = datetime.fromisoformat(image2_date)
        
        closest_datetime_1 = find_closest_time(indates_lumin_datetime, image1_datetime)
        print("Closest datetime to coherence ", image1_datetime, "is ", closest_datetime_1)            

        closest_datetime_2 = find_closest_time(indates_lumin_datetime, image2_datetime)
        print("Closest datetime to coherence ", image1_datetime, "is ", closest_datetime_2)            

        ## Create indices for luminance images with the closest datetime
        i_lumin1 = [i for i in range(len(indates_lumin_datetime)) if indates_lumin_datetime[i] == closest_datetime_1]
        i_lumin2 = [i for i in range(len(indates_lumin_datetime)) if indates_lumin_datetime[i] == closest_datetime_2]

        # ## Diag
        #  print("i_lumin1: ", i_lumin1, i_lumin1[0], type(i_lumin1))
        #  print(indates_lumin_datetime[i_lumin1[0]])
        #  print(indates_lumin_datetime[i_lumin2[0]])
        
        ## Need to verify time difference is less than x-minutes, proceed if it is, move on if it does not. 
        timeDelta_1 = indates_lumin_datetime[i_lumin1[0]] -  image1_datetime
        timeDelta_2 = indates_lumin_datetime[i_lumin2[0]] -  image2_datetime

        print("Time Delta Image 1 in seconds: ", timeDelta_1.total_seconds())
        print("Time Delta Image 2 in seconds: ", timeDelta_2.total_seconds())
                
        if len(i_lumin1) < 1 or len(i_lumin2) < 1:
            if len(i_lumin1) < 1:
                print("Luminance image cannot be found for: ", image1_date)
                
            if len(i_lumin2) < 1:
                print("Luminance image cannot be found for: ", image2_date)
            continue
        
        elif timeDelta_1.total_seconds() > lum_tDelta_cutoff_secs or timeDelta_2.total_seconds() > lum_tDelta_cutoff_secs:
            if timeDelta_1.total_seconds() > lum_tDelta_cutoff_secs:
                    print("Total time difference between coh and luminance of ", timeDelta_1.total_seconds(), " exceeds " , lum_tDelta_cutoff_secs , " seconds.")
                    print("Skipping this iteration.")
            if timeDelta_2.total_seconds() > lum_tDelta_cutoff_secs:
                print("Total time difference between coh and luminance of ", timeDelta_2.total_seconds(), " exceeds " , lum_tDelta_cutoff_secs , " seconds.")
                print("Skipping this iteration.")
            continue

        else:
            print("Found both luminance images")       
            i_lumin1 = int(i_lumin1[0]) # This assumes only 1 luminance image for this date
            i_lumin2= int(i_lumin2[0]) # This assumes only 1 luminance image for this date

            print("Normalized image 1:", indates_lumin[i_lumin1], ':', luminance_images_norm[i_lumin1])
            print("Normalized image 2:", indates_lumin[i_lumin2], ':', luminance_images_norm[i_lumin2])
            
            in1_path, in1_fname = os.path.split(luminance_images_norm[i_lumin1])
            in2_path, in2_fname = os.path.split(luminance_images_norm[i_lumin2])
            
            ## load and difference the images (image1 minus image2)
            # deltaLum, image_profile = difference_luminance_images(luminance_images_norm[i_lumin1], luminance_images_norm[i_lumin2], in1_fname, in2_fname)
            difference_luminance_images(luminance_images_norm[i_lumin1], luminance_images_norm[i_lumin2], in1_fname, in2_fname, output_directory, indates_lumin, i_lumin1, i_lumin2)
            
             
        print("")
        
        
    
    

