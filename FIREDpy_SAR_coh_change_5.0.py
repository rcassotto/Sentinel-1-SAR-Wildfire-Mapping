#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
## Version 5
##  - built from 4a_draft
##  - cleans up from 4a
##  - picks first image as the reference, if none is specified. 
##  - auto selects sequential scenes (e.g. day1_day12 and day12_day24) instead of assuming all scenes in a given directory are already paired. This
##          simplifies the code, cleans up the directories and eliminates the need for parsing of scenes with different overlapping areas. 



## Version 4a_draft.py
Created on Thu Mar 17 15:22:31 2022
Modified Nov 7, 2023 to only look at sequences of coh images. 

@author: rcassotto


"""

from image_proc_module import Image_proc
from FIREDpy_ARIA_module_v1 import ARIA
import numpy as np
import argparse, ast
import glob
import os
import rasterio
#from datetime import datetime, timedelta



# make the output directory
def mk_outdirectory(outpath_full):
    if not os.path.exists(outpath_full):
        print("Making output directory: ", outpath_full)
        os.makedirs(outpath_full)
    return


def normalize_secondary_images(infile_0_path, infile_0_name, infile_1_path, infile_1_name, outdirectory):
       ### Create outfilename
    fname0_split = infile_0_name.split( "_" ) # coh filename 0
    fname1_split = infile_1_name.split( "_" ) # coh filename 1
    sec_0_datestr = fname1_split[0] 
    sec_1_datestr = fname1_split[1]
    t_band_outname = ( outdirectory + '/' + sec_0_datestr + '_' + sec_1_datestr + '_normalized_tBand.tif' )
    check_file = os.path.isfile(t_band_outname)
    if check_file == False:

        ## Read in images
        print("Reading in Reference Image Coherence file: ", infile_0_name)
        ref_image = Image_proc(infile_0_path + '/' + infile_0_name)
        ref_image.read()
        
        print("Reading in Secondary Image Coherence file: ", infile_1_name)
        sec_image = Image_proc(infile_1_path + '/' + infile_1_name)
        sec_image.read()
        print("Ref Image type:", type(ref_image))
        print("Sec Image type:", type(sec_image))
             
        ## Get image corner coordinates
        image_data = rasterio.open(infile_0_path + '/' + infile_0_name)
        lon = (image_data.bounds.left, image_data.bounds.right)
        lat = (image_data.bounds.top, image_data.bounds.bottom)
        minLat = np.min(lat)
        maxLat = np.max(lat)
        minLon = np.min(lon)
        maxLon = np.max(lon)
        geographic_bounds = [minLon, maxLon, minLat, maxLat]
        print("Geographic Bounds: ", geographic_bounds)
        
          
        ## Next up: call to process aria
        aria = ARIA(geographic_bounds)
        png_outfilename_and_path = outdirectory + '/'+  infile_0_name + "_" + infile_1_name + ".png"
        print(png_outfilename_and_path)
        t_threshold=0.2
        aria_data = aria.process_ARIA(ref_image, [sec_image], png_outfilename_and_path, t=t_threshold, map_type=aria.simple_map, file_prefix='', show_hists=True)
        ### t is the threshold for coh change. If differences are greater than t, they will be highlighted, if not 0 is returned.
        #### aria_data
        
        #### Rev 4a update: save t-band, the normalized secondary image, as a geotiff
        for inmap in aria_data:
            band = inmap.map # In Rev 4a, this is the normalized secondary coherence image
            print(type(band),len(band), np.size(band))

            profile = image_data.profile.copy()  # copy geotiff meta data from input file   
            print("Writing Normalized Secondary image ", t_band_outname, " to ", outdirectory)

            with rasterio.open(t_band_outname, 'w', **profile) as dst:
                dst.write(band,1)
                
            ## Write normalized_tBand meta data file
            write_metafile(t_band_outname, infile_0_path, infile_0_name, infile_1_path, infile_1_name, outdirectory)
                
    else:
        print("Normalized Coh Image exists. Moving on")
        print("Skipping: ", t_band_outname); print(" ")
        pass
    
def write_metafile(t_band_outname, infile_0_path, infile_0_name, infile_1_path, infile_1_name, outdirectory):
    Line1 = ( "Normalized Coh Image Ouput Path + Directory: " + t_band_outname + "\n")
    Line2 = ( "Coh Image 1 Input filename: " + infile_0_name + "\n")
    Line3 = ( "Coh Image 1 Input Path: " + infile_0_path + "\n")
    Line4 = ( "Coh Image 2 Input filename: " + infile_1_name + "\n")
    Line5 = ( "Coh Image 2 Input Path: " + infile_1_path + "\n")
    

    ##metafile_outname = ( outpath_full + '/' + in_dir_name + '_meta.txt' )
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


def write_coh_change_metafile(in_coh_path_1, in_coh_fname_1, in_coh_path_2, in_coh_fname_2, diff_eq, tiff_outfilename_and_path):       
    Line1 = ( "Coh Change Output File + Directory: " + tiff_outfilename_and_path + "\n")
    Line2 = ( "Coh Image 1 Input filename: " + in_coh_fname_1 + "\n")
    Line3 = ( "Coh Image 1 Input Path: " + in_coh_path_1 + "\n")
    Line4 = ( "Coh Image 2 Input filename: " + in_coh_fname_2 + "\n")
    Line5 = ( "Coh Image 2 Input Path: " + in_coh_path_2 + "\n")
    Line6 = ( "Coh Change Equation: " + diff_eq + "\n")

    ##metafile_outname = ( outpath_full + '/' + in_dir_name + '_meta.txt' )
    metafile_outname = tiff_outfilename_and_path.replace('.tif','_log.txt') # infile basename
    print('Writing Log file: ', metafile_outname)
    meta_file = open (metafile_outname, "w")
    meta_file.write(Line1)
    meta_file.write(Line2)
    meta_file.write(Line3)
    meta_file.write(Line4)
    meta_file.write(Line5)
    meta_file.write(Line6)
    meta_file.close()
    return        

def create_coh_img_pair_lists(Coh_files, tif_suffix):
    #### Create a list of coh image pairs to difference
    coh_img_pr1_fname_and_path=[]
    coh_img_pr2_fname_and_path=[]
    for nIterDate in range(0,len(Coh_files)-1):

        ### Get Date/Time string of secondary image
         in_coh_path_0, in_coh_fname_0 = os.path.split(Coh_files[nIterDate])
         in_coh_path_1, in_coh_fname_1 = os.path.split(Coh_files[nIterDate+1])
         
         coh_img_pr1_fname_and_path.append(Coh_files[nIterDate])
         coh_img_pr2_fname_and_path.append(Coh_files[nIterDate+1])

    return coh_img_pr1_fname_and_path, coh_img_pr2_fname_and_path
    



def Difference_coh_image_pairs(in_coh_path_1, in_coh_fname_1, in_coh_path_2, in_coh_fname_2, outdirectory):
    ### Create outfilename
    fname0_split = in_coh_fname_1.split( "_" )
    fname1_split = in_coh_fname_2.split( "_" )
    ref_0_datestr = fname0_split[0]
    ref_1_datestr = fname0_split[1]
    sec_0_datestr = fname1_split[0] 
    sec_1_datestr = fname1_split[1]
    
    if ref_1_datestr == sec_0_datestr:  ## only difference sequential images; the assumption is that pairs in directory are all 12-day and sequential
        tiff_outfilename_and_path = outdirectory + '/' + ref_0_datestr + '_' + sec_0_datestr + '_' + sec_1_datestr + '_cohChange.tif' 
        check_file = os.path.isfile(tiff_outfilename_and_path)
        
        if check_file == False:  # check if file exists, skip it if it does
            ## Read in images
            print("Reading in First Image Coherence file: ", in_coh_fname_1)
            ref_in = rasterio.open(in_coh_path_1 + '/' + in_coh_fname_1)
            ref_image = ref_in.read(1)
            
            print("Reading in Secondary Image Coherence file: ", in_coh_fname_2)
            sec_in = rasterio.open(in_coh_path_2 + '/' + in_coh_fname_2)
            sec_image = sec_in.read(1)
            
              ## Get image corner coordinates - this could be cleaned up
            (lat_arr,lon_arr,pixel_size_lat_m,pixel_size_lon_m,pixel_size_lat_deg,pixel_size_lon_deg,image_data) = get_image_geo_info(in_coh_path_1 + '/' + in_coh_fname_1)
          
            print("Differencing coh images: ", in_coh_fname_2, " from ", in_coh_fname_1)
            # diff_map = np.subtract(ref_image, sec_image)  ## original; 
            # diff_eq = 'ref_image - sec_image'
            diff_map = np.subtract(sec_image, ref_image)  ## update 20250206; Now coh loss will appear negative. 
            diff_eq = 'sec_image - ref_image'

            #### Write outfile - single band geotiff
            write_out_geotiff(image_data, diff_map, tiff_outfilename_and_path)
            
            #### Write coh change metafile (for accounting and sanity check)
            write_coh_change_metafile(in_coh_path_1, in_coh_fname_1, in_coh_path_2, in_coh_fname_2, diff_eq, tiff_outfilename_and_path)
            
            
        else:
            print("Coh Change file exists. Moving on")
            print("Skipping: ", tiff_outfilename_and_path); print(" ")
            pass
    else:
        print("Coherence image pairs are not sequential")
        pass

## Get image georeference information (e.g. lat, lon, pixel size, etc)
def get_image_geo_info(in_sigma_file):
        image_data = rasterio.open(in_sigma_file)
        lon = (image_data.bounds.left, image_data.bounds.right)
        lat = (image_data.bounds.top, image_data.bounds.bottom)
        minLat = np.min(lat)
        maxLat = np.max(lat)
        minLon = np.min(lon)
        maxLon = np.max(lon)
        medianLat = np.median(lat)
 
        gt = image_data.transform
        pixel_size_lon_deg = gt[0]
        pixel_size_lat_deg = gt[4]
        pixel_size_lat_m = gt[4]*111e3
        pixel_size_lon_m = gt[0]*111e3*np.cos(np.radians(medianLat))
        lon_arr = np.arange(minLon+0.5*pixel_size_lon_deg, maxLon, pixel_size_lon_deg, dtype=float)# lon array with 0.5 chip offset
        lat_arr = np.arange(maxLat+0.5*pixel_size_lat_deg, minLat, pixel_size_lat_deg, dtype=float)
        
        return lat_arr,lon_arr,pixel_size_lat_m,pixel_size_lon_m,pixel_size_lat_deg,pixel_size_lon_deg,image_data


def write_out_geotiff(image_data, diff_map, tiff_outfilename_and_path):  
        print('Saving output Geotiff: ', tiff_outfilename_and_path)
        profile = image_data.profile.copy()  # copy geotiff meta data from input file   
        with rasterio.open(tiff_outfilename_and_path, 'w', **profile) as dst:
           dst.write(diff_map,1)
           
def parse_input_file_data(parser):
    ### Function to parse input data
    parser.add_argument('filename',type=argparse.FileType('r'))
    p = parser.parse_args()
    
    ## Read in path to directory containing coherence images
    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)
        
    ref_image_fname = args['ref_image']
    coh_directory = args['coh_file_loc']
    tif_suffix = args['tif_suffix']
    outdirectory = args['out_dir']
    return ref_image_fname, coh_directory, tif_suffix, outdirectory

def define_ref_image(Coh_files):
  print("No Reference image specified for normalization.")
  ref_image_fname_and_path = Coh_files[0]
  ref_image_path, ref_image_fname = os.path.split(ref_image_fname_and_path)
  print("Using first image as the Reference image: ", ref_image_fname)
  return ref_image_fname

def parse_ref_and_sec_datestrings(coh_filenames):
   ref_coh_datestr = [] # init ref coh datestr list
   sec_coh_datestr = [] # init secondary coh datestr list
   for n in range(0, len(coh_filenames)):
       pStr = coh_filenames[n].split('_')
       ref_coh_datestr.append(pStr[0])    # reference datestr in each coherence pair
       sec_coh_datestr.append(pStr[1])    # secondary datestr in each coherence pair
   return ref_coh_datestr, sec_coh_datestr

def create_coh_delta_image_lists(sec_coh_datestr, coh_filenames):
   base_path, base_name = os.path.split(Coh_files[0])
   coh_img_pr1_fname_and_path = [] # init
   coh_img_pr2_fname_and_path = [] # init
   for n in range(0, len(sec_coh_datestr)):
       try:
           ii = ref_coh_datestr.index(sec_coh_datestr[n])  # ii is an index into ref_coh_datestr for the secondary image 
           coh_img_pr1_fname_and_path.append(base_path + '/' + coh_filenames[n]) # First image pair (e.g. day1_day12_coherence_VV.tif)
           coh_img_pr2_fname_and_path.append(base_path + '/' + coh_filenames[ii]) # Second image pair (e.g. day12_day24_coherence_VV.tif)
           print("Calculating Coh Difference on file pairs: ", coh_filenames[n], coh_filenames[ii])
       except ValueError:
           print("No reference image with date", sec_coh_datestr[n] , "found. Moving on.")
   return coh_img_pr1_fname_and_path, coh_img_pr2_fname_and_path


def create_ref_image_list(Coh_files, ref_image_fname, coh_directory):
    ref_list = []
    for nCopy in range(0,len(Coh_files)):
          ref_list.append(coh_directory + '/' + ref_image_fname)
    return ref_list


def normalize_images(Coh_files, ref_list, outdirectory):  
    for n in range(0,len(Coh_files)):
        print("Normalizing secondary images to reference image. Iteration ", n)
        print("Image1: ", ref_list[n])
        print("Image2: ", Coh_files[n])
        in_coh_path_1, in_coh_fname_1 = os.path.split(ref_list[n])
        in_coh_path_2, in_coh_fname_2 = os.path.split(Coh_files[n])
        normalize_secondary_images(in_coh_path_1, in_coh_fname_1, in_coh_path_2, in_coh_fname_2, outdirectory)
        print(" ")
    return
     

def Perform_coh_diff(coh_img_pr2_fname_and_path, coh_img_pr1_fname_and_path, outdirectory, tif_suffix):
    for n in range(0,len(coh_img_pr2_fname_and_path)):
        print("Processing Coh Image Difference Pair: ", n)
        in_coh_path_1, in_coh_fname_1 = os.path.split(coh_img_pr1_fname_and_path[n])
        in_coh_path_2, in_coh_fname_2 = os.path.split(coh_img_pr2_fname_and_path[n])
    
        if in_coh_fname_1 != ref_image_fname:
            in_coh_path_1 = outdirectory
            in_coh_fname_1 = in_coh_fname_1.replace(tif_suffix,'normalized_tBand.tif')
          
        in_coh_path_2 = outdirectory
        in_coh_fname_2 = in_coh_fname_2.replace(tif_suffix,'normalized_tBand.tif')
        Difference_coh_image_pairs(in_coh_path_1, in_coh_fname_1, in_coh_path_2, in_coh_fname_2, outdirectory)
        print(" ")
    return

if __name__ ==  "__main__":

    #### Parse input file information        
    parser = argparse.ArgumentParser()
    ref_image_fname, coh_directory, tif_suffix, outdirectory = parse_input_file_data(parser)    
    
    
    ### Make outdirectory, if it doesn't exist
    mk_outdirectory(outdirectory)
            
    ## Get a List of coherence images and paths
    Coh_files_unsorted = glob.glob(os.path.join(coh_directory, '*' + tif_suffix))   
    Coh_files = sorted(Coh_files_unsorted)  ## sort files. Should be chronologically given filename nomenclature (e.g. 20200827_20200908_coherence_VV.tif)
          
    ## Identify reference image from input file; Else pick the first image in the set
    if len(ref_image_fname)==0:
        ref_image_fname = define_ref_image(Coh_files)
        
    ### parse Reference and Secondary date strings to identify sequential scene pairs (day1_day2_day3) 
    coh_filenames = [os.path.basename(path) for path in Coh_files] # Create a list of filenames only, without full paths 
    (ref_coh_datestr, sec_coh_datestr) = parse_ref_and_sec_datestrings(coh_filenames)
   
    ### Create two lists of coherence images to then difference     
    (coh_img_pr1_fname_and_path, coh_img_pr2_fname_and_path) = create_coh_delta_image_lists(sec_coh_datestr, coh_filenames)
    
    ### Create a list of reference images (copies of reference image), the same length as secondary list, 
    ref_list = create_ref_image_list(Coh_files, ref_image_fname, coh_directory)
    
    ### Noramlize images to reference; Run Coh difference algorithm (i.e. aria); loop over coh_img_pr lists
    normalize_images(Coh_files, ref_list, outdirectory)
        
    ### Difference the coherence pairs
    print("length of coh_img_pr2:", len(coh_img_pr2_fname_and_path))
    Perform_coh_diff(coh_img_pr2_fname_and_path, coh_img_pr1_fname_and_path, outdirectory, tif_suffix)
   
    
  
    # print("length coh_filenames: ", len(coh_filenames))
    # print("length coh datestr1:", len(ref_coh_datestr))
    # print("length coh datestr2", len(sec_coh_datestr))
    # print("length coh coh_img1:", len(coh_img_pr1_fname_and_path))
    # print("length of coh_img_pr2:", len(coh_img_pr2_fname_and_path))




