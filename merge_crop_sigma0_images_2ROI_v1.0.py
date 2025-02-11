#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : {Ryan Cassotto}
# Created Date: 2024/02/01
# version ='1.0'
# ---------------------------------------------------------------------------
""" Python Module to merge adjacent Sentinel-1 classification scenes and crop to ROI"""
# ---------------------------------------------------------------------------
# from image_proc_module import Image_proc

import glob
import numpy as np
import os
import argparse
import ast
import sys
# import time


from collections import Counter
import json
from osgeo import gdal

### import gdal_merge and base path to it
import subprocess
gm = os.path.join('/usr/bin', 'gdal_merge.py')


# ---------------------------------------------------------------------------


def mk_outdirectory(outpath_full):
    if not os.path.exists(outpath_full):
        print("Making output directory: ", outpath_full)
        os.makedirs(outpath_full)
    return


def read_roi_gis_file(roi_path_and_filename, fire_buffer_prct):
    with open(roi_path_and_filename) as f:
        d = json.load(f)  # load input file (*.json) into dictionary, d

        # key to finding first, last date
        sd = d['features'][0]['properties']['ig_date']
        ed = d['features'][0]['properties']['last_date']
        fid = d['features'][0]['properties']['fid_1']
        climate = d['features'][0]['properties']['main_clim']
        print("start: ", sd)
        print("stop: ", ed)
        print("fid: ", fid)
        print("Main Climate: ", climate)

        # This is list of size 1xn with all the coordinates, including all points from multiple polygons
        kk = d['features'][0]['geometry']['coordinates']

        x_arr = []
        y_arr = []  # initialize x,y arrays
        # Loop over list of coordinates to extract vertices for multipolygon, if the ROI includes multipolygon
        for n in range(0, len(kk)):
            #print("Reading Polygon number: ", n)
            # inList is an list of coordinates embedded 3 layers deep; hence the [:][:][n]
            inList = kk[:][:][n]
            inList_length = len(inList[0][:])
            # print(inList[0][:][:])

            # Parse inList into x,y coordinates
            in_x_arr = []
            in_y_arr = []
            # separate into 1 x n arrays for each x,y coordiante
            for nn in range(0, inList_length):
                # Global array of all points, not just current polygon
                x_arr.append(inList[0][nn][0])
                y_arr.append(inList[0][nn][1])

                # local array of polygon points for current polygon
                in_x_arr.append(inList[0][nn][0])
                in_y_arr.append(inList[0][nn][1])

            print("Polygon ", n, "Extents: ", min(in_x_arr), max(
                in_x_arr), min(in_y_arr), max(in_y_arr), "degrees")

            # Sanity check for multiPolygon points
            poly_out = np.transpose(np.array((in_x_arr, in_y_arr)))
            csv_outname = roi_path_and_filename.replace(
                '.json', ('_Poly_' + str(n) + '.csv'))
            np.savetxt(csv_outname, poly_out, delimiter=',')

        # Confirmed! These are the max/min for all polygons within the multipolygon coordinates
        minx_gis = min(x_arr)
        maxx_gis = max(x_arr)
        miny_gis = min(y_arr)
        maxy_gis = max(y_arr)
        print("Min/Max Extents:", minx_gis, maxx_gis, miny_gis, maxy_gis)
        delta_lon = maxx_gis - minx_gis
        delta_lat = maxy_gis - miny_gis
        print("ROI Latitude Extents, Height:", miny_gis,
              maxy_gis, "(", delta_lat, ") deg")
        print("ROI Longitude Extents, Width:", minx_gis,
              maxx_gis, "(", delta_lon, ") deg")

        # Add buffer
        if miny_gis > maxy_gis:
            print("Coordinates are in Southern Hemisphere.")
            tmp_miny = miny_gis
            miny_gis = maxy_gis
            maxy_gis = tmp_miny
            del tmp_miny

        minx = minx_gis - delta_lon * int(fire_buffer_prct)/100
        maxx = maxx_gis + delta_lon * int(fire_buffer_prct)/100
        miny = miny_gis - delta_lat * int(fire_buffer_prct)/100
        maxy = maxy_gis + delta_lat * int(fire_buffer_prct)/100

        # sanity check - output csv file with roi points
        new_roi_out = ((minx, miny), (maxx, miny),
                       (maxx, maxy), ((minx, maxy)))
        new_roi_csvname = roi_path_and_filename.replace(
            '.json', ('_roi_wBuff.csv'))
        np.savetxt(new_roi_csvname, new_roi_out, delimiter=',')

        return minx, maxx, miny, maxy


def read_infile_dates(classification_image_list):
   # Get dates of classification images, save as a list
    indates = []
    for n in range(0, len(classification_image_list)):
        infile_path, infile_name = os.path.split(classification_image_list[n])
        i_date = infile_name.split("_")
        infile_date = i_date[4].split('T')[0]
        indates.append(infile_date)
        # print(infile_date)

    ### Convert the list to a numpy array, sort and save indices of sorted datevals
    indate_arr = np.array(indates)
    # This is an index into the 'classification_image_list' for image dates, sorted chronologically
    i_indate_sort_arr = np.argsort(indate_arr)
    img_date_sort_arr = indate_arr[i_indate_sort_arr]
    img_date_sort = img_date_sort_arr.tolist()
    i_indate_sort = i_indate_sort_arr.tolist()

    return img_date_sort, i_indate_sort


def filter_non_unique(img_date_sort):
    multi_image_dates = [item for item, count in Counter(img_date_sort).items() if count > 1]
    return multi_image_dates

def filter_single_image_dates(img_date_sort):
    single_image_dates = [item for item, count in Counter(img_date_sort).items() if count == 1]
    return single_image_dates

### Parse critical basefile info from first infilename
def read_basefile_info(sigma0_image_0):
    inPath, inName = os.path.split(sigma0_image_0)
    inName_parse = inName.split('_')
    sigma0_base = inName_parse[0] + '_' + inName_parse[1] + '_' + inName_parse[2] + '_' + inName_parse[3] +'_'
    base_polarity = inName_parse[-1].split('.')[0] # polarity of input file (e.g VV or HH)
    print("Base Name prefix:", sigma0_base)
    print("infile string:", infile_str)
    print("Infile polarity:", base_polarity)
    return sigma0_base, base_polarity
    
     
    
### Get complimentary SAR polarity. Assumes a fixed set of antenna polarity configurations for Sentinel-1
def  get_complimentary_polarity(base_polarity):
       if base_polarity == 'VV':
           comp_polarity = 'VH'
           
       elif base_polarity == 'HH':
           comp_polarity = 'HV'
       
       elif base_polarity =='VH':
           comp_polarity = 'VV'
           
       elif base_polarity =='HV':
           comp_polarity = 'HH'
           
       return comp_polarity


def crop_images_2_roi_wGdal(inimage_sort, out_tifname, minx, maxx, miny, maxy):
    ds = gdal.Open(inimage_sort)
    ds = gdal.Translate(out_tifname, ds, projWin=[minx, maxy, maxx, miny], noData=0)
    ds = None
    return 



def initiate_log_file(log_output_fname, cropped_outdirectory):
    stdoutOrigin = sys.stdout
    log_output_fname_and_path = os.path.join(cropped_outdirectory,log_output_fname)
    sys.stdout = open(log_output_fname_and_path, "w")
    print("Starting: merge_crop_sigma0_images_2ROI.py")
    return stdoutOrigin

def close_log_file(stdoutOrigin):
    sys.stdout.close()
    sys.stdout = stdoutOrigin
    return


if __name__ == "__main__":
    # ---------------------------------------------------------------------------
    # Implement Arg parse to pull relevant input parameters from input.txt file
    # Use Argparge to get information from command line
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=argparse.FileType('r'))
    p = parser.parse_args()

    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)
        
    # make cropped2ROI output 
    cropped_outdirectory = args['out_dir']
    mk_outdirectory(cropped_outdirectory) 
        
    ## Initiate log file
    log_output_fname = "merge_crop_sigma0_files.txt"
    stdoutOrigin = initiate_log_file(log_output_fname, cropped_outdirectory)
    
    # Read in fire ROI, add buffer from user-defined value
    fire_buffer_prct = args['fire_perimeter_buffer_prct']
    roi_path_and_filename = os.path.join( args['fire_roi_gis_path'], args['fire_roi_gis_file'])
    [minx, maxx, miny, maxy] = read_roi_gis_file( roi_path_and_filename, fire_buffer_prct)  # Read in GIS file (e.g. fid_11862.json)
 

    # Get list of input files (i.e. sigma0 images)
    infile_str = "*" + args['input_image_suffix']
    sigma0_path = args['sigma0_file_loc']
    sigma0_images = glob.glob(os.path.join(sigma0_path, infile_str))  
    
    ## Read important base file info
    sigma0_base, base_polarity = read_basefile_info(sigma0_images[0])
    comp_polarity = get_complimentary_polarity(base_polarity)

    # Read dates from filenames, sort and save indices into the original (unsorted) list
    img_date_sort, i_indate_sort = read_infile_dates(sigma0_images)

    # Creat a new list of classification files, sorted chronologically
    inimage_sort = []
    for n in range(0, len(sigma0_images)):
        inimage_sort.append(sigma0_images[i_indate_sort[n]])

    # Find images acquired on the same date (adjacent scenes),
    print("Total Number of Sigma0 images covering ROI:", len(sigma0_images))
    multi_image_dates = filter_non_unique(img_date_sort)
    img_date_sort_uniq = list(np.unique(img_date_sort))
    print("Number of Scenes after merging adjacent same-day images:", len(img_date_sort_uniq))
    print("")    
    
    single_image_dates = filter_single_image_dates(img_date_sort) # list of dates with single image acquisitions
    
    # Loop over all input images
    print("img date sort uniq:", len(img_date_sort_uniq), img_date_sort_uniq); print("")
    print("multi image dates:", multi_image_dates); print("")
    print("dates with single images:", single_image_dates);    print("")
    
  
    # Key variables so far
    # 1) inimage_sort: list of sigma0 images, sorted chronologically
    # 2) img_date_sort: list of dates (YYYYMMDD) for scenes in inimage_sort list.
    # 3) i_indate_sort: (perhaps less important) indices into the original sigma0_images list of files (NOT chronologically sorted)
    # 4) single_image_dates: list of dates with single acquisitions only (i.e. non-merge list)
    # 5) multi_image_dates: list of dates with multiple acquisitions/images (i.e.merge list)
    
    
    for n in range(0, len(img_date_sort_uniq)):
        print("")
        print("Processing:", img_date_sort_uniq[n], "    Day ", n+1, "of", len(img_date_sort_uniq))


        if img_date_sort_uniq[n] in multi_image_dates:
            
            ## index into img_date_sort for ducpliate days
            idate = [i for i, x in enumerate(img_date_sort) if x == img_date_sort_uniq[n]]
            
            num_scenes = len(idate)
            print("There are ", num_scenes, "SAR scenes for:", img_date_sort_uniq[n])
           
            merged_path_and_file_0 = os.path.join(sigma0_path, 'merged_0.tif')
            merged_path_and_file_1 = os.path.join(sigma0_path, 'merged_1.tif')
                        
            cmd_0 = "python3 " + gm + " -o " + merged_path_and_file_0 + " -n 0"  ## VV or HH image
            cmd_1 = "python3 " + gm + " -o " + merged_path_and_file_1 + " -n 0" ## VH or HV image
            
            ##### Create a string of Sigma0 files to merge for the current day
            for nDuplicate in range(0, len(idate)):
                print(inimage_sort[idate[nDuplicate]])
                ## VV or HH image
                infile_path_and_name = os.path.join(inimage_sort[idate[nDuplicate]])
                cmd_0 = cmd_0 + " " + infile_path_and_name              
                
                ## VH or HV image
                comp_filename = (inimage_sort[idate[nDuplicate]]).replace(base_polarity, comp_polarity)
                cmd_1 = cmd_1 + " " + comp_filename     
                
            ## Execute merge commands with gdal_merge
            print("Merging files to ", merged_path_and_file_0)
            print(cmd_0); print("")
            subprocess.call(cmd_0, shell=True)
            
            print("Merging files to ", merged_path_and_file_1)
            print(cmd_1); print("")
            subprocess.call(cmd_1,shell=True)
         
            ## Crop input image to ROI and output with _roi.tif suffix (VV or HH)
            print("Cropping", merged_path_and_file_0, "to", minx, max, miny, maxy)
            out_tifname = sigma0_base + img_date_sort_uniq[n] +'_' + args['input_image_suffix']
            out_tifname_and_path = os.path.join(cropped_outdirectory, out_tifname)
            print("Creating:", out_tifname_and_path); print("")
            crop_images_2_roi_wGdal(merged_path_and_file_0, out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)

            ## Crop input image to ROI and output with _roi.tif suffix (VH or HV)
            print("Cropping", merged_path_and_file_1, "to", minx, maxx, miny, maxy)
            comp_out_tifname_and_path = out_tifname_and_path.replace(base_polarity, comp_polarity)         
            print("Creating:", comp_out_tifname_and_path); print("")
            crop_images_2_roi_wGdal(merged_path_and_file_1, comp_out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)
          
            ### Delete the merged Sigma0 images. 
            os.remove(merged_path_and_file_0)
            os.remove(merged_path_and_file_1)
                    
            del out_tifname, out_tifname_and_path, comp_out_tifname_and_path, comp_filename

             

        else:  # for single image days
            print(img_date_sort_uniq[n], "has only 1 SAR scene.") 
            idate = [i for i, x in enumerate(img_date_sort) if x == img_date_sort_uniq[n]]
            print(inimage_sort[idate[0]])
            print("")
            
            ### Image 1 (VV or HH)
            print("Cropping", inimage_sort[idate[0]]  ,"to :", minx, maxx, miny, maxy)
            out_tifname = sigma0_base + img_date_sort_uniq[n] +'_' + args['input_image_suffix']
            out_tifname_and_path = os.path.join(cropped_outdirectory, out_tifname)
            print("Creating:", out_tifname_and_path); print("")
            crop_images_2_roi_wGdal(inimage_sort[idate[0]], out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)
  
            ### Image 2 - Complimentary image (VH or HV)
            comp_infilename_and_path=inimage_sort[idate[0]].replace(base_polarity, comp_polarity)
            print("Cropping", comp_infilename_and_path  ,"to :", minx, maxx, miny, maxy)
            comp_out_tifname_and_path = out_tifname_and_path.replace(base_polarity, comp_polarity)
            print("Creating:", comp_out_tifname_and_path); print("")
            crop_images_2_roi_wGdal(comp_infilename_and_path, comp_out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)

            del out_tifname, out_tifname_and_path, comp_out_tifname_and_path


    ### Close log file
    close_log_file(stdoutOrigin)
   











