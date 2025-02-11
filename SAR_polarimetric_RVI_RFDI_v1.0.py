#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : {Ryan Cassotto}
# Created Date: 2025/02/03
# version ='1.0'
# 
# version 1.0 normalizes the sigma0 images first, then performs the calculations on the normalized sigma0 images
# ---------------------------------------------------------------------------
""" Sentinel 1 SAR polarimetric index calculations """  
# ---------------------------------------------------------------------------
from image_proc_module import Image_proc
#import Cluster
# from FIREDpy_coherence_module.FIREDpy_ARIA_module_v1 import ARIA
from FIREDpy_ARIA_module_v1 import ARIA

import rasterio
import glob
import numpy as np
import os
from rasterio.plot import reshape_as_raster
import argparse,ast
import time
import shutil

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
    i_nan = np.argwhere(b1_flat == 0)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
    b1_flat[i_nan] = np.nan
    b1_raster = b1_flat.reshape(band_shape)
    return b1_raster
    
def write_geotiff_single_band(tiff_outname, image_data, luminosity):
    profile_1band = image_data.profile.copy()  # copy geotiff meta data from input file   
    # profile_3band.update(count = 3)  # Change count from 1-band (input file) to 3-bands (output RGB image)
    # lum_out_raster = reshape_as_raster(luminosity)  # reshape to output correctly as a raster (necessary to create RGB-style image)
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
    return geographic_bounds, image_data




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
    print("index for reference image:", imax, type(imax))
    ref_image_fname_and_path=sigma0_images[imax]  
    ref_image_path, sigma0_ref_img_fname = os.path.split(ref_image_fname_and_path)
    print("ref_image_fname:", sigma0_ref_img_fname)
    
    return sigma0_ref_img_fname

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
                print("Secondary Image:" , type(sec_image))
                
                # sec_image_in = rasterio.open(nFile_2Norm)
                # sec_image = sec_image_in.read(1)
                # sec_image[sec_image>1E38] = np.nan ## UPDATE: 20240614; check for images with other no Data values
                # sec_image[sec_image==0] = np.nan # convert zeros to nans


                ## Next up: call to process aria
                #sec_path, sec_fname = os.path.split(nFile_2Norm)

                aria = ARIA(geographic_bounds)
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
                    write_metafile(t_band_outname, sigma0_dir, sigma_ref_img_fname, sec_path, sec_fname, sigma0_dir)
                    print("")
            else:
                print(nFile_2Norm, " already normalized. Moving on....")
                pass
                
    ref_norm_name_and_path = ref_image_path_and_fname.replace('.tif','_norm.tif')
    shutil.copy(ref_image_path_and_fname, ref_norm_name_and_path)
    return




def calculate_SAR_indices(Sigma_files,rvi_outdirectory, rfdi_outdirectory, use_norm_images):
    
    ncount=0
    for in_sigma_file in Sigma_files:
        # ---------------------------------------------------------------------------
        ## construct output directory name, make output directory     
        granule = in_sigma_file.split('/')[-1] # granule is a string with the original Sigma0 filename (e.g. S1A_IW_GRDH_1SDV_20150921T232856_20150921T232918_007820_00AE3B_E36F_Sigma0_VV)
        mk_outdirectory(rvi_outdirectory)  # make output directory#
        mk_outdirectory(rfdi_outdirectory)  # make output directory#
        
        
        ncount+=1
        time_start=time.time()
        print("Processing SAR Indexed Scenes: ", ncount, " of ", len(Sigma_files), ":", granule)
     
        # # ---------------------------------------------------------------------------
        # ## Read in images
        # VV_image = Image_proc(in_sigma_file)
        # VV_image.read()
        # print("VV Image type:", type(VV_image))       
       
        # vh_sigma_file=in_sigma_file.replace("VV", "VH")
        # VH_image = Image_proc(vh_sigma_file)
        # VH_image.read()
        
        
        # ---------------------------------------------------------------------------
        ### Read in images (UPDATE)
        in_sigma = rasterio.open(in_sigma_file)
        vv_band = in_sigma.read(1)
        print("in_sigma type:", type(in_sigma))
        VV_image = convert_noData2NaNs(vv_band)
        
     
        vh_sigma_file=in_sigma_file.replace("VV", "VH")
        in_sigma_vh = rasterio.open(vh_sigma_file)
        vh_band = in_sigma_vh.read(1)
        VH_image = convert_noData2NaNs(vh_band)
        
        
        # ---------------------------------------------------------------------------
        ## geographic boundaries and data on one of the image files
        geographic_bounds, image_data = Get_image_boundaries(vh_sigma_file)
        
        
        # Calculate RVI and output as a single band geotiff
   # RVI = 4VH / (VV + VH)
        rvi = (4 * VH_image) / (VV_image + VH_image)
       
        if use_norm_images == 'true': 
            rvi_tiff_outname = os.path.join(rvi_outdirectory,granule.replace('_Sigma0_VV_norm.tif','_rvi.tif'))
        else:
            rvi_tiff_outname = os.path.join(rvi_outdirectory,granule.replace('_Sigma0_VV.tif','_rvi.tif'))
        
        write_geotiff_single_band(rvi_tiff_outname, image_data, rvi)

        
        
        ## Calculate RFDI and output as a single band geotiff
####	RFDI = ( VV – VH ) / (VV + VH)
        rfdi = (VV_image - VH_image) / (VV_image + VH_image)

        if use_norm_images == 'true': 
            rfdi_tiff_outname = os.path.join(rfdi_outdirectory,granule.replace('_Sigma0_VV_norm.tif','_rfdi.tif'))
        else:
            rfdi_tiff_outname = os.path.join(rfdi_outdirectory,granule.replace('_Sigma0_VV.tif','_rfdi.tif'))
        
        write_geotiff_single_band(rfdi_tiff_outname, image_data, rfdi)      
        
        
        time_stop=time.time()
        print("Done. Elapsed Time: ", time_stop - time_start, " seconds."); print("")
    return

## Parse arguments from input file
def  parse_input_file(p):
    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)        

    sigma0_dir = args['sigma_file_loc']
    sigma0_ref_img_fname_VV = args['sigma_ref_image']
    use_norm_images = args['use_norm_images']    
    rvi_outdirectory = args['rvi_out_dir']
    rfdi_outdirectory = args['rfdi_out_dir']
    return sigma0_dir, sigma0_ref_img_fname_VV, use_norm_images, rvi_outdirectory, rfdi_outdirectory


if __name__ ==  "__main__":
    # ---------------------------------------------------------------------------
    ## Implement argparse to pull relevant input parameters from input.txt file
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',type=argparse.FileType('r'))
    p = parser.parse_args()
    (sigma0_dir, sigma0_ref_img_fname_VV, use_norm_images, rvi_outdirectory, rfdi_outdirectory) = parse_input_file(p)

    
    ## Get list of input files    
    Sigma_files_VV = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV.tif'))   
    Sigma_files_VH = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VH.tif'))   

    # ---------------------------------------------------------------------------
    ## Identify reference image, if one is not specified 
    if len(sigma0_ref_img_fname_VV)==0:
        ### Loop to read in all images and find image with greatest mean; set this as reference iamge. 
        ## VV
        sigma0_ref_img_fname = determine_reference_sigma0_image(Sigma_files_VV)
        sigma0_ref_img_fname_VV = sigma0_ref_img_fname
        print("Sigma0 reference image:", sigma0_ref_img_fname_VV)    
        del sigma0_ref_img_fname
    
        ## VH
        sigma0_ref_img_fname = determine_reference_sigma0_image(Sigma_files_VH)
        sigma0_ref_img_fname_VH = sigma0_ref_img_fname 
        print("Sigma0 reference image:", sigma0_ref_img_fname_VH)    
        del sigma0_ref_img_fname
    
    else:
        sigma0_ref_img_fname_VH = sigma0_ref_img_fname_VV.replace('VV', 'VH')
      
        
      
    # ---------------------------------------------------------------------------
    ## New section to normalize ALL Sigma0 images. 
    if use_norm_images == 'true':       
        ### Normalize all Sigma0 - VV images to the specified or determined reference image
        print("Normaling all Sigma0 - VV images to: ", sigma0_ref_img_fname_VV)
        normalize_sigma0_images(Sigma_files_VV, sigma0_dir, sigma0_ref_img_fname_VV)
        print("Done normalizing secondary VV images."); print("")

        ### Normalize all Sigma0 - VH images to the specified or determined reference image
        print("Now normaling all Sigma0 - VH images to: ", sigma0_ref_img_fname_VH)
        normalize_sigma0_images(Sigma_files_VH, sigma0_dir, sigma0_ref_img_fname_VH)
        print("Done normalizing secondary VH images."); print("")

        ### Create new list of *_VV_norm.tif
        Sigma_files = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV_norm.tif'))
    else:
        Sigma_files = glob.glob(os.path.join(sigma0_dir, '*Sigma0_VV.tif'))


    
    # ---------------------------------------------------------------------------
    calculate_SAR_indices(Sigma_files,rvi_outdirectory, rfdi_outdirectory, use_norm_images)  

