#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : {Ryan Cassotto}
# Created Date: 2025/02/07
# version ='3.0'
# ---------------------------------------------------------------------------
""" Python Module to combine coh change, lum change, dRVI, and dRFDI """
# ---------------------------------------------------------------------------

import rasterio
import glob
import numpy as np
import os
import argparse
import ast
from scipy.ndimage import gaussian_filter

import matplotlib.pyplot as plt
from scipy.optimize import fsolve

### Rev 03: updates the filtering, adds RVI, RFDI into the mix
### Rev 02: replaces OTSU bimodal selection with basic statistics mean and standar deviation to find the threshold values. 
### Rev 01: incorporates auto selection of thresholding through OTSU bimodal techniques



# ---------------------------------------------------------------------------
nSigma_coh=1
nSigma_lum=1
nSigma_dRVI=1
nSigma_dRFDI=1

def mk_outdirectory(outpath_full):
    if not os.path.exists(outpath_full):
        print("Making output directory: ", outpath_full)
        os.makedirs(outpath_full)
    return


def read_delLuminance_date_pairs(DelLum_images_list):
   # Get dates of classification images, save as a list
    indate_pairs_lum = [] # initialize indates_list
    for n in range(0, len(DelLum_images_list)):
        infile_path, infile_name = os.path.split(DelLum_images_list[n])
        i_date = infile_name.split("_")
        infile_date_str = i_date[0] + "_" + i_date[1]
        indate_pairs_lum.append(infile_date_str)
    return indate_pairs_lum 


def read_coh_change_date_pairs(Coh_change_image_list):
   # Get dates of classification images, save as a list
    indate_pairs_cohChange = [] # initialize indates_list
    for n in range(0, len(Coh_change_image_list)):
        infile_path, infile_name = os.path.split(Coh_change_image_list[n])
        i_date = infile_name.split("_")
        infile_date_str = i_date[1] + "_" + i_date[2]
        indate_pairs_cohChange.append(infile_date_str)
    return indate_pairs_cohChange 



def normalize_array_and_bin(subset_array_flat, N):
        subset_array_norm = (subset_array_flat - min(subset_array_flat))/(max(subset_array_flat) - min(subset_array_flat))
        bins = np.linspace(0, 1, N) # spacing of 256 discrete points between 0 and 1
        bin_count = np.histogram(subset_array_norm, bins)[0] # returns the count in each bin
        M = len(subset_array_norm)
        return bin_count, M # returns t=[0 1 2 ... 255], histogram, and number of elements in subarray


def count_peaks(x, y, cs, dcs):
    eps_tol=1e-16
    peak_cnt = 0
    valley_cnt = 0
    peaks = []
    valleys = []
    peak_heights = []
    for i in range(1,len(y)):
        if np.abs(y[i-1]) < eps_tol: continue
        else:
            d1 = dcs(x[i-1])
            d2 = dcs(x[i])
            sign1 = d1 / np.abs(d1)
            sign2 = d2 / np.abs(d2)
            test = sign1*sign2
            if test == -1.:
                if sign1 == 1.:
                    peak_cnt += 1
                    peaks.append(fsolve(dcs, x[i]))
                    peak_heights.append(y[i])
                else:
                    valley_cnt += 1
                    valleys.append(fsolve(dcs, x[i])[0])
    return peak_cnt, valley_cnt, peaks, valleys, peak_heights

def gaussian_convolution(y):
    gs = [0.2261, 0.5478, 0.2261] # from Uddin paper
    pad_y = np.zeros(len(y)+2)
    pad_y[1:-1] = y
    new_y = np.zeros(len(y))
    for k in range(1,len(pad_y)-1):
        t1 = pad_y[k-1] * gs[0] # term 1
        t2 = pad_y[k] * gs[1] # term 2
        t3 = pad_y[k+1] * gs[2] # term 3
        new_y[k-1] = t1+t2+t3
    return new_y


    
def plot_histograms_with_thresholds(coh_image, delLum_image, coh_threshold, lum_threshold, lum_inimage_path_and_name, outdirectory):
    nbins=128
    
    ## Prep coh image data
    coh_image = np.nan_to_num(coh_image)
    img_flat_coh = coh_image.flatten()
    img_flat_coh = img_flat_coh[img_flat_coh != 0]
    # y, bins = np.histogram(img_flat_coh, bins=nbins)

    ## Prep lum image data
    delLum_image = np.nan_to_num(delLum_image)
    img_flat_lum = delLum_image.flatten()
    img_flat_lum = img_flat_lum[img_flat_lum != 0]
    # y, bins = np.histogram(img_flat_lum, bins=nbins)
  
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
  
    ## subplot(221); coh change map
    coh_map = ax[0,0].imshow(coh_image, cmap='RdGy_r',vmin=-0.5, vmax=0.5)
    ax[0,0].set_title('Coh Change')
    fig.colorbar(coh_map, ax=ax[0,0], location='bottom')
    
    ### subplot(222) 
    ax[0,1].set_title('Coh Change Histogram ')   
    ax[0,1].hist(img_flat_coh, bins=nbins, color='skyblue', edgecolor='black')
    ax[0,1].set_title('Histogram w/Threshold (red): {:.2f}'.format(coh_threshold)) 
    ax[0,1].set_xlim(-0.5, 0.5)
    yLims = ax[0,1].get_ylim()
    yLine = np.linspace(yLims[0], yLims[1],50)
    xLine = coh_threshold * np.ones(len(yLine))
    ax[0,1].plot(xLine, yLine, color='red')


    
    ### subplot(223)
    lum_map = ax[1,0].imshow(delLum_image, cmap='PuOr',vmin=-0.25, vmax=0.25)
    ax[1,0].set_title('Delta Lum')
    fig.colorbar(lum_map, ax=ax[1,0], location='bottom')

    ### subplot(224)
    #ax[1,1]
    ax[1,1].hist(img_flat_lum, bins=nbins, color='magenta', edgecolor='black')
    ax[1,1].set_title('Histogram w//Threshold: {:.2f}'.format(lum_threshold)) 
    ax[1,1].set_xlim(-0.25, 0.25)
    yLims = ax[1,1].get_ylim()
    yLine = np.linspace(yLims[0], yLims[1],50)
    xLine = lum_threshold * np.ones(len(yLine))
    ax[1,1].plot(xLine, yLine, color='green')
    
    
    
    ### Save histogram and threshold figures
    in_image_path, in_image_name = os.path.split(lum_inimage_path_and_name)
    out_pngfilename = in_image_name.replace('.tif','_cohChange_thresh.png')  
    outfilename_and_path = os.path.join( outdirectory, out_pngfilename)
    supTitleStr= out_pngfilename
    print(supTitleStr)
    print("Saving output file: ", outfilename_and_path)
    plt.savefig(outfilename_and_path)
       
    return

def determine_cohChange_threshold(coh_image):
    print('Coherence Threshold not specified. Setting threshold to ', str(nSigma_coh), 'Std. Dev. of Scene Values')
    coh_image_mean = np.nanmean(coh_image)
    coh_image_sigma = np.nanstd(coh_image)
    coh_threshold = coh_image_mean - nSigma_coh * coh_image_sigma
    print("Coherence Threshold", str(nSigma_coh), "Sigma or ", str(coh_threshold))
    return coh_threshold

def determine_LumChange_threshold(delLum_image):
     print('Delta Luminance Threshold not specified. Setting threshold to ', str(nSigma_lum), 'Std. Dev. of Scene Values')
     lum_image_mean = np.nanmean(delLum_image)
     lum_image_sigma = np.nanstd(delLum_image)
     lum_threshold = lum_image_mean - nSigma_lum * lum_image_sigma
     print("Luminance Threshold", str(nSigma_lum), "Sigma or ", str(lum_threshold))
     return lum_threshold

def determine_dRVI_threshold(dRVI_image):
     print('dRVI Threshold not specified. Setting threshold to ', str(nSigma_lum), 'Std. Dev. of Scene Values')
     dRVI_image_mean = np.nanmean(dRVI_image)
     dRVI_image_sigma = np.nanstd(dRVI_image)
     dRVI_threshold = dRVI_image_mean - nSigma_dRVI * dRVI_image_sigma
     print("dRVI Threshold", str(nSigma_dRVI), "Sigma or ", str(dRVI_threshold))
     return dRVI_threshold

def determine_dRFDI_threshold(dRFDI_image):
     print('dRFDI Threshold not specified. Setting threshold to ', str(nSigma_lum), 'Std. Dev. of Scene Values')
     dRFDI_image_mean = np.nanmean(dRFDI_image)
     dRFDI_image_sigma = np.nanstd(dRFDI_image)
     dRFDI_threshold = dRFDI_image_mean + nSigma_dRFDI * dRFDI_image_sigma
     print("dRFDI Threshold", str(nSigma_dRFDI), "Sigma or ", str(dRFDI_threshold))
     return dRFDI_threshold
 

def apply_lpf (image2filter, image_fname_and_path, image_profile):
    print("Applying LPF to: ", image_fname_and_path)
    ## Filter deltaLum image with Gaussian Low Pass ***** NEW 20240501 ******
    ##### NOTE: radius produces a kerner that is 2*radius + 1 ALONG EACH AXIS, 
    lpf_image = gaussian_filter(image2filter, sigma=2, radius=2)
    in_image_path, in_image_name = os.path.split(image_fname_and_path)
    filt_image_outname = in_image_name.replace('.tif','_filt.tif')       
    write_out_geotiff(filt_image_outname, in_image_path, image_profile, lpf_image) # Write out filtered DeltaLum Image
    return lpf_image, filt_image_outname

def apply_hpf (image2filter, image_fname_and_path, image_profile):
    print("Applying HPF to: ", image_fname_and_path)
    ##### apply HPF by inverting the data, applying a lpf and then inverting it again. 
    ##### NOTE: radius produces a kerner that is 2*radius + 1 ALONG EACH AXIS, 
    image2filter_inverted = image2filter*-1 # invert the data
    tmp_image = gaussian_filter(image2filter_inverted, sigma=2, radius=2)
    hpf_image = tmp_image*-1 
    in_image_path, in_image_name = os.path.split(image_fname_and_path)
    filt_image_outname = in_image_name.replace('.tif','_filt.tif')       
    write_out_geotiff(filt_image_outname, in_image_path, image_profile, hpf_image) # Write out filtered DeltaLum Image
    return hpf_image, filt_image_outname


def read_cohChange_image(in_coh_image_and_path):
    ## Read in coherence change image
    print("Reading in Coh Change Image: ", in_coh_image_and_path)
    in_coh = rasterio.open(in_coh_image_and_path)
    coh_image = in_coh.read(1)
    image_profile = in_coh.profile.copy()  # For writing out geotiffs
    return coh_image, image_profile

def read_deltaLum_image(in_lum_image_and_path):
    print("Reading in Luminance Change Image: ", in_lum_image_and_path)
    delLum_in = rasterio.open(in_lum_image_and_path)
    delLum_image = delLum_in.read(1)
    return delLum_image

def read_dRVI_image(in_dRVI_image_and_path):
    print("Reading in dRVI Change Image: ", in_dRVI_image_and_path)
    dRVI_in = rasterio.open(in_dRVI_image_and_path)
    dRVI_image = dRVI_in.read(1)
    return dRVI_image

def read_dRFDI_image(in_dRFDI_image_and_path):
    print("Reading in dRFDI Change Image: ", in_dRFDI_image_and_path)
    dRFDI_in = rasterio.open(in_dRFDI_image_and_path)
    dRFDI_image = dRFDI_in.read(1)
    return dRFDI_image

   
# def combine_sar_images(in_coh_image_and_path, in_lum_image_and_path, coh_threshold, lum_threshold, outdirectory, indate_pair_string, use_filtered_Delta_imgs, use_filtered_Coh_imgs, save_threshold_figures): 
def combine_sar_images(in_coh_image_and_path, in_lum_image_and_path, in_dRVI_image_and_path, in_dRFDI_image_and_path, coh_threshold, lum_threshold, dRVI_threshold, \
                       dRFDI_threshold, outdirectory, indate_pair_string, use_filtered_Delta_imgs, use_filtered_Coh_imgs, save_threshold_figures):
   
    print('filtered flag:', use_filtered_Delta_imgs)

    ## ------------------------------------------------------------------------------
    ## Read in coherence change image; Determine threshold, if undefined
    (coh_image, image_profile) = read_cohChange_image(in_coh_image_and_path)
    if not coh_threshold:
       (coh_threshold) = determine_cohChange_threshold(coh_image)

    ## ------------------------------------------------------------------------------
    ## Read in lum change image; Determine threshold, if undefined
    delLum_image = read_deltaLum_image(in_lum_image_and_path)  
    if not lum_threshold:
       (lum_threshold) = determine_LumChange_threshold(delLum_image)

    ## ------------------------------------------------------------------------------
    ## Read in dRVI change image; Determine threshold, if undefined
    dRVI_image = read_dRVI_image(in_dRVI_image_and_path)  
    if not dRVI_threshold:
       (dRVI_threshold) = determine_dRVI_threshold(dRVI_image)
    

    ## ------------------------------------------------------------------------------
    ## Read in dRFDI change image; Determine threshold, if undefined
    dRFDI_image = read_dRFDI_image(in_dRFDI_image_and_path)  
    if not dRFDI_threshold:
       (dRFDI_threshold) = determine_dRFDI_threshold(dRFDI_image)
    
    ## ------------------------------------------------------------------------------
    ## plot histograms
    if save_threshold_figures==True:
        # (x, y, y_smooth, cs) = prep_histogram_plots(delLum_image)
        # plot_histograms_with_thresholds(delLum_image, x, y, y_smooth, lum_threshold, cs, in_lum_image_and_path, outdirectory)
        plot_histograms_with_thresholds(coh_image, delLum_image, coh_threshold, lum_threshold, in_lum_image_and_path, outdirectory)


    ## ------------------------------------------------------------------------------
    ## Apply filtering - Coherence
    ## Treat Coh seperately since we it is already filtered during the coherence generation.
    if use_filtered_Coh_imgs == "True":
        (lpf_image, filt_image_outname) = apply_lpf(coh_image, in_coh_image_and_path, image_profile)
        coh_fire = np.where(lpf_image < coh_threshold, 1., 0.).astype(np.int8) # Using LPF cohDelta
        print("Lum fire sum:", np.sum(coh_fire))
        ### Check for empty rasters or those with little data; return uniform array with 0.01
        if np.sum(coh_fire)<30:
            print(filt_image_outname, "is empty!!!!!!!!!");             coh_fire = np.ones(np.shape(coh_image))*0.01;         del lpf_image
    else:
        print("Using original (unfiltered) DeltaCoh images.")
        coh_fire = np.where(coh_image > coh_threshold, 1., 0.).astype(np.int8)

        
        
        
    ## ------------------------------------------------------------------------------
    ## Apply filtering - Radar indices and luminance; polarimetric combinations 
    if use_filtered_Delta_imgs == "True":

        ### apply LPF filter to deltaLum
        (lpf_image, filt_image_outname) = apply_lpf(delLum_image, in_lum_image_and_path, image_profile)
        lum_fire = np.where(lpf_image < lum_threshold, 1., 0.).astype(np.int8) # Using LPF deltaLum
        print("Lum fire sum:", np.sum(lum_fire))
        ### Check for empty rasters or those with little data; return uniform array with 0.01
        if np.sum(lum_fire)<30:
            print(filt_image_outname, "is empty!!!!!!!!!");             lum_fire = np.ones(np.shape(delLum_image))*0.01;         del lpf_image
        
        
        ### apply LFP filter to dRVI
        (lpf_image, filt_image_outname) = apply_lpf(dRVI_image, in_dRVI_image_and_path, image_profile)
        dRVI_fire = np.where(lpf_image < dRVI_threshold, 1., 0.).astype(np.int8) # Using LPF dRVI
        print("dRVI fire sum:", np.sum(dRVI_fire))
        ### Check for empty rasters or those with little data; return uniform array with 0.01
        if np.sum(dRVI_fire)<30:
            print(filt_image_outname, "is empty!!!!!!!!!");             dRVI_fire = np.ones(np.shape(dRVI_image))*0.01;         del lpf_image
        
        
        ### apply HFP filter to dRVI
        (hpf_image, filt_image_outname) = apply_hpf(dRFDI_image, in_dRFDI_image_and_path, image_profile)
        dRFDI_fire = np.where(hpf_image > dRFDI_threshold, 1., 0.).astype(np.int8) # Using LPF deltaLum
        print("dRFDI fire sum:", np.sum(dRFDI_fire))
        ### Check for empty rasters or those with little data; return uniform array with 0.01
        if np.sum(dRFDI_fire)<30:
            print(filt_image_outname, "is empty!!!!!!!!!");             dRFDI_fire = np.ones(np.shape(dRFDI_image))*0.01;         del hpf_image
        
            
      
      
    else: # create binary image based on non-filtered deltaLum, deltaRFDI, deltaRVI
        print("Using original (unfiltered) DeltaLum images.")
        lum_fire = np.where(delLum_image < lum_threshold, 1., 0.).astype(np.int8)

        print("Using original (unfiltered) dRVI images.")
        dRVI_fire = np.where(dRVI_image < lum_threshold, 1., 0.).astype(np.int8)

        print("Using original (unfiltered) dRDFI images.")
        dRFDI_fire = np.where(dRFDI_image > dRFDI_threshold, 1., 0.).astype(np.int8)


    ## ------------------------------------------------------------------------------
    ## write out indiv geotiffs (coh burn pixel and lum burn pixels)
    coh_fire_name = "FIREDpy_" + indate_pair_string + "_coh_fire.tif"
    lum_fire_name = "FIREDpy_" + indate_pair_string + "_lum_fire.tif"
    dRVI_fire_name = "FIREDpy_" + indate_pair_string + "_dRVI_fire.tif"
    dRFDI_fire_name = "FIREDpy_" + indate_pair_string + "_dRFDI_fire.tif"    
    write_out_geotiff(coh_fire_name, outdirectory, image_profile, coh_fire) # Write out coh fire pixels
    write_out_geotiff(lum_fire_name, outdirectory, image_profile, lum_fire) # Write out lum fire pixels
    write_out_geotiff(dRVI_fire_name, outdirectory, image_profile, dRVI_fire) # Write out lum fire pixels
    write_out_geotiff(dRFDI_fire_name, outdirectory, image_profile, dRFDI_fire) # Write out lum fire pixels
    
    
    ## ------------------------------------------------------------------------------
    ## combine rasters - keep all ppixels (diagnostic)==
    tmp_arr = np.add(coh_fire, lum_fire) # original     
    tmp_product_name = "FIREDpy_sar_combined_ALLPixels_" + indate_pair_string + ".tif"
    write_out_geotiff(tmp_product_name, outdirectory, image_profile, tmp_arr) # Write out combined - all pixels (unMasked)


    ### combine rasters, mask values > 1
    sar_product = np.where(tmp_arr > 1, 1., 0.).astype(np.int8)
    print("sar product sum:", np.sum(sar_product))
    ### Test for nearly empty array (10 was insufficient, 25 seems okay for 100 m)
    if np.sum(sar_product)<30:
        print(indate_pair_string, " combined SAR is empty!!!!!!!!!")
        sar_product = np.ones(np.shape(coh_image))*0.01
        #sar_product = np.zeros(np.shape(coh_image))     #### All Zeroes doesn't work
    
    sar_product_name = "FIREDpy_sar_combined_" + indate_pair_string + ".tif"
    write_out_geotiff(sar_product_name, outdirectory, image_profile, sar_product) # Write out out combined
   
   ## Write log file
    sar_combined_logname = sar_product_name.replace('.tif','_log.txt')
    sar_combined_logname_and_path = os.path.join(outdirectory,sar_combined_logname)
    write_sar_combination_log_file(sar_combined_logname_and_path, in_coh_image_and_path, in_lum_image_and_path, coh_threshold, lum_threshold)
   
    del coh_image, delLum_image
   
    return


def write_out_geotiff(outfilename, output_directory, image_profile, raster_2_write):
    output_product_outfilename_and_path = os.path.join(output_directory, outfilename)
    print("Writing output raster: ", output_product_outfilename_and_path)
    with rasterio.open(output_product_outfilename_and_path, 'w', **image_profile) as dst:
        dst.write(raster_2_write,1)
    return

def write_sar_combination_log_file(sar_combined_logname_and_path, in_coh_image_and_path, in_lum_image_and_path, coh_threshold, lum_threshold):
    lum_path, lum_fname = os.path.split(in_lum_image_and_path)
    coh_path, coh_fname = os.path.split(in_coh_image_and_path)

    Line1 = ( "Combined SAR Image Ouput Path + Directory: " + sar_combined_logname_and_path + "\n")
    Line2 = ( "Coh Image filename: " + coh_fname + "\n")
    Line3 = ( "Coh Image Path: " + coh_path + "\n")
    Line4 = ( "Lum Image filename: " + lum_fname + "\n")
    Line5 = ( "Lum Image Path: " + lum_path + "\n")
    Line6 = ( "Lum Threshold: " + str(lum_threshold) +"\n")
    Line7 = ( "Coh Threshold: " + str(coh_threshold) +"\n")        

    print('Writing Log file: ', sar_combined_logname_and_path)
    meta_file = open (sar_combined_logname_and_path, "w")
    meta_file.write(Line1)
    meta_file.write(Line2)
    meta_file.write(Line3)
    meta_file.write(Line4)
    meta_file.write(Line5)
    meta_file.write(Line6)
    meta_file.write(Line7)
    meta_file.close()
    return

def Generate_accumulated_fire_prod(image_product_list, outdirectory, img_type):
    image_data = rasterio.open(image_product_list[0])
    image_cube = np.zeros((image_data.height, image_data.width, len(image_product_list)), dtype=np.float32)
    image_profile = image_data.profile.copy()  # copy geotiff meta data from input file   

    print("")
    print("Generating accumulated burn area from " + img_type + " change.")
    for n in range(0, len(image_product_list)):
        print("Reading in:", image_product_list[n])

        in_prod = rasterio.open(image_product_list[n])
        image_prod = in_prod.read(1)

        ## add current image to data cube
        image_cube[:,:,n] = image_prod
    
    ### Sum over data cube
    cube_sum = np.sum(image_cube, axis=2, dtype=np.float32)
    accum_prod = np.where(cube_sum > 0.7, 1., 0.).astype(np.int8)
    
    accum_outname="FIREDpy_TotalBurnArea_" + img_type +".tif"
    write_out_geotiff(accum_outname, outdirectory, image_profile, accum_prod) # Write out coh fire pixels    
    
    write_accum_sar_prod_metafile(outdirectory, accum_outname, image_product_list)
    return
    

def write_accum_sar_prod_metafile(outdirectory, geotiff_outname, input_file_list):
    metafile_outname = geotiff_outname.replace('.tif','_Meta.txt') # infile basename
    metafile_outname_and_path = os.path.join(outdirectory, metafile_outname)
    print('Writing Log file: ', metafile_outname_and_path); print("")

    ## Make a list of output lines for metafile        
    lines = []
    lines.append("Meta File Name: " + metafile_outname_and_path)
    lines.append("Geotiff Output File: " + outdirectory + "/" + geotiff_outname)
    lines.append("")
    lines.append("Input Files:")
    for n in range(0,len(input_file_list)):
        #print(input_file_list[n])
        lines.append(input_file_list[n])            

    ## Write metafile
    meta_file = open (metafile_outname_and_path, "w")
    for line in lines:
            meta_file.write(line + "\n")
    meta_file.close()
    return

def parse_input_data(p):
    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)
    
    DeltaLum_path = args['DeltaLum_path']  
    outdirectory = args['out_dir']
    DeltaLum_suffix = args['DeltaLum_suffix']
    coh_change_path = args['coh_change_path']
    coh_change_suffix = args['coh_change_suffix']
    dRVI_path = args['dRVI_path']
    dRFDI_path = args['dRFDI_path']
    dRVI_suffix = args['dRVI_suffix']
    dRFDI_suffix = args['dRFDI_suffix']
    coh_threshold = args['coh_threshold']
    lum_threshold = args['lum_threshold']
    dRVI_threshold = args['dRVI_threshold']
    dRFDI_threshold = args['dRFDI_threshold']
    use_filtered_Delta_imgs = args['filter_Delta_images']
    use_filtered_Coh_imgs = args['filter_Coh_images']
        
    return DeltaLum_path, outdirectory, DeltaLum_suffix, coh_change_path, coh_change_suffix, dRVI_path, dRFDI_path, \
        dRVI_suffix, dRFDI_suffix, coh_threshold, lum_threshold, dRVI_threshold, dRFDI_threshold, use_filtered_Delta_imgs, use_filtered_Coh_imgs 



if __name__ == "__main__":
    # ---------------------------------------------------------------------------
    # Implement Arg parse to pull relevant input parameters from input.txt file
    # Use Argparge to get information from command line
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=argparse.FileType('r'))
    p = parser.parse_args()
    (DeltaLum_path, outdirectory, DeltaLum_suffix, coh_change_path, coh_change_suffix, dRVI_path, dRFDI_path, dRVI_suffix, \
     dRFDI_suffix, coh_threshold, lum_threshold, dRVI_threshold, dRFDI_threshold, use_filtered_Delta_imgs, use_filtered_Coh_imgs) = parse_input_data(p)

    ## ------------------------------------------------------------------------------
    # make output directory
    mk_outdirectory(outdirectory)
    

    ## ------------------------------------------------------------------------------
    ### Check for empty threshold strings. Initiatlize with "None" if empty, read in value, if not. 
    if not coh_threshold:
        coh_threshold=[]
    else:
        coh_threshold = float(coh_threshold)
    
    if not lum_threshold:
        lum_threshold=[]
    else:
        lum_threshold = float(lum_threshold)

    if not dRVI_threshold:
        dRVI_threshold=[]
    else:
        dRVI_threshold = float(dRVI_threshold)
    
    if not dRFDI_threshold:
         dRFDI_threshold=[]
    else:
         dRFDI_threshold = float(dRFDI_threshold)
        

    ## ------------------------------------------------------------------------------
    ## Get list of Delta Luminance input files 
    DelLum_images = glob.glob(os.path.join(DeltaLum_path, '*' + DeltaLum_suffix))
    for n in range(0, len(DelLum_images)): 
        print(DelLum_images[n])
    print("")

    ### dRVI
    dRVI_images = glob.glob(os.path.join(dRVI_path, '*' + dRVI_suffix))
    for n in range(0, len(dRVI_images)): 
        print(dRVI_images[n])
    print("")

    ### dRFDI 
    dRFDI_images = glob.glob(os.path.join(dRFDI_path, '*' + dRFDI_suffix))
    for n in range(0, len(dRFDI_images)): 
        print(dRFDI_images[n])
    print("")


    ## ------------------------------------------------------------------------------
    ## Get Luminance date pair strings
    indate_pairs_lum = read_delLuminance_date_pairs(DelLum_images)
    # print(indate_pairs_lum); print("")
    indate_pairs_dRVI = read_delLuminance_date_pairs(dRVI_images)
    indate_pairs_dRFDI = read_delLuminance_date_pairs(dRFDI_images)

    ## ------------------------------------------------------------------------------
    ## Get list of Coherence Change images
    Coh_change_images = glob.glob(os.path.join(coh_change_path, '*' + coh_change_suffix))
 
    ## ------------------------------------------------------------------------------
      ## Get coh Change date pair strings
    indate_pairs_cohChange = read_coh_change_date_pairs(Coh_change_images)
    # print("Indate Pairs Coh Change: ", indate_pairs_cohChange); print("")     
    print("Coh Threshold:", coh_threshold)
    print("Lum Threshold:", lum_threshold); print("")
    

    ## Combine SAR images - Loop over coh Change image
    for n in range(0,len(Coh_change_images)):

        in_coh_image_and_path = Coh_change_images[n] # coh Change file and path
        
        i_lumin1 = [i for i in range(len(indate_pairs_lum)) if indate_pairs_lum[i] == indate_pairs_cohChange[n]]  # index into indate_pairs_lum for current pair
        in_lum_image_and_path = DelLum_images[ (int(i_lumin1[0])) ]   # DelLum current filename and path
        
        i_dRVI1 = [i for i in range(len(indate_pairs_dRVI)) if indate_pairs_dRVI[i] == indate_pairs_cohChange[n]]   # index into indate_pairs_dRVI for current pair
        in_dRVI_image_and_path = dRVI_images[ (int(i_dRVI1[0])) ]  # dRVI current filename and path
        
        i_dRFDI1 = [i for i in range(len(indate_pairs_dRFDI)) if indate_pairs_dRFDI[i] == indate_pairs_cohChange[n]]   # index into indate_pairs_dRFDI for current pair
        in_dRFDI_image_and_path  = dRFDI_images[ (int(i_dRFDI1[0])) ]  # dRFDI current filename and path
        
        ## Sanity check
        print("coh change image:", in_coh_image_and_path)
        print("lum change image:", in_lum_image_and_path)
        print("dRVI change image:", in_dRVI_image_and_path)
        print("dRFDI change image:", in_dRFDI_image_and_path)
        print("")
  
        #combine_sar_images(in_coh_image_and_path, in_lum_image_and_path, coh_threshold, lum_threshold, outdirectory, indate_pairs_cohChange[n], use_filtered_Delta_imgs, use_filtered_Coh_imgs)
        # combine_sar_images(in_coh_image_and_path, in_lum_image_and_path, coh_threshold, lum_threshold, outdirectory, indate_pairs_cohChange[n], use_filtered_Delta_imgs, use_filtered_Coh_imgs, save_threshold_figures=True)
        combine_sar_images(in_coh_image_and_path, in_lum_image_and_path, in_dRVI_image_and_path, in_dRFDI_image_and_path, coh_threshold, lum_threshold, dRVI_threshold, \
                           dRFDI_threshold, outdirectory, indate_pairs_cohChange[n], use_filtered_Delta_imgs, use_filtered_Coh_imgs, save_threshold_figures=True)
        
        
        
        
        
        
        
        
        
        
    # ## Sum Coh Change images (accumulated fire perimeter)
    # img_type="coh"
    # coh_product_list = glob.glob(os.path.join(outdirectory, "*_" + img_type + "_fire.tif"))
    # Generate_accumulated_fire_prod(coh_product_list, outdirectory, img_type)


    #  ## Sum Lum Change images (accumulated fire perimeter)
    # img_type="lum"
    # lum_product_list = glob.glob(os.path.join(outdirectory, "*_" + img_type + "_fire.tif")) 
    # Generate_accumulated_fire_prod(lum_product_list, outdirectory, img_type)

    #  ## Sum Combined Change images (accumulated fire perimeter)
    # img_type="FIREDpy_sar_combined_2"
    # combined_product_list = glob.glob(os.path.join(outdirectory, img_type + "*" + ".tif")) 
    # Generate_accumulated_fire_prod(combined_product_list, outdirectory, img_type)










