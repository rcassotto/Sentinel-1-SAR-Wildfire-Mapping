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
from datetime import datetime


from collections import Counter
import json
from osgeo import gdal

### import gdal_merge and base path to it
import subprocess
#gm = os.path.join('/usr/bin', 'gdal_merge.py')  # cires research nodes
gm = os.path.join('/home/rcassotto/.conda/envs/FIREDpy/bin', 'gdal_merge.py')  # Skye


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


def read_infile_dates(coh_image_list):
   # Get dates of coherence images, save as a list
    indates_1 = []
    indates_2 = []
    indate_list = []
    for n in range(0, len(coh_image_list)):
        infile_path, infile_name = os.path.split(coh_image_list[n])
        i_date = infile_name.split("_")
#        infile_date = i_date[4].split('T')[0]
#        indates.append(infile_date)
        infile_date1 = i_date[0].split('T')[0]
        infile_date2 = i_date[1].split('T')[0]        
        indates_1.append(infile_date1)
        indates_2.append(infile_date2)
        indate_list.append(infile_date1 + '_' + infile_date2)
#        print("infilename: ", infile_name)
#        print("date1:", infile_date1)
#        print("date2:", infile_date2)
        
        

    ### Convert the list to a numpy array, sort and save indices of sorted datevals
#    indate_1_arr = np.array(indates_1)
#    indate_2_arr = np.array(indates_2)
    
#    print(indate_list)
    
    # This is an index into the 'coh_image_list' for image dates, sorted chronologically
#    i_indate_sort_arr = np.argsort(indate_arr)  # create an index of indates sorted chronologically
#    img_date_sort_arr = indate_arr[i_indate_sort_arr] # sort date array chronologically
#    img_date_sort = img_date_sort_arr.tolist()  # create a list of img dates sorted chronologically
#    i_indate_sort = i_indate_sort_arr.tolist() # create a list of the indices of the img_dates in chronological order

    return indate_list #img_date_sort, i_indate_sort


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
    coh_base = inName_parse[0] + '_' + inName_parse[1] + '_' + inName_parse[2] + '_' + inName_parse[3] +'_'
    base_polarity = inName_parse[-1].split('.')[0] # polarity of input file (e.g VV or HH)
    print("Base Name prefix:", coh_base)
    print("infile string:", infile_str)
#    print("Infile polarity:", base_polarity)
    return coh_base, base_polarity
    

def crop_images_2_roi_wGdal(inimage_sort, out_tifname, minx, maxx, miny, maxy):
    ds = gdal.Open(inimage_sort)
    ds = gdal.Translate(out_tifname, ds, projWin=[minx, maxy, maxx, miny], noData=0)
    ds = None
    return 



def initiate_log_file(log_output_fname, cropped_outdirectory):
    stdoutOrigin = sys.stdout
    log_output_fname_and_path = os.path.join(cropped_outdirectory,log_output_fname)
    sys.stdout = open(log_output_fname_and_path, "w")
    print("Starting: merge_crop_coh_images_2ROI.py")
    return stdoutOrigin

def close_log_file(stdoutOrigin):
    sys.stdout.close()
    sys.stdout = stdoutOrigin
    return


def write_merged_cropped_imaged_metafile(out_text_file_and_path, coh_images_sort, idate, cropped_outdirectory, out_tifname):
    out_text=[]
    out_text.append("Outfilename: " + out_tifname)
    out_text.append("Outfilepath: " + cropped_outdirectory)
    for nn in range(0, len(idate)):
        out_text.append("Input file(s): " + coh_images_sort[idate[nn]])
    out_text.append("Python Script: " + sys.argv[0])
    out_text.append("Creation Date: " + str(datetime.now()))
#    print(out_text)
   
    f = open(out_text_file_and_path,"w")
    for line in out_text:
        f.writelines(str(line)+"\n")
        print(str(line))
    f.close()
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
#    log_output_fname = "merge_crop_coh_files.txt"
#    stdoutOrigin = initiate_log_file(log_output_fname, cropped_outdirectory)
    
    # Read in fire ROI, add buffer from user-defined value
    fire_buffer_prct = args['fire_perimeter_buffer_prct']
    roi_path_and_filename = os.path.join( args['fire_roi_gis_path'], args['fire_roi_gis_file'])
    [minx, maxx, miny, maxy] = read_roi_gis_file( roi_path_and_filename, fire_buffer_prct)  # Read in GIS file (e.g. fid_11862.json)

    # Get list of input files (i.e. sigma0 images)
    infile_str = "*" + args['input_image_suffix']
    coh_path = args['coh_file_loc']
    coh_images = glob.glob(os.path.join(coh_path, infile_str))  
    
    ## Read important base file info
    coh_base, base_polarity = read_basefile_info(coh_images[0])
#    comp_polarity = get_complimentary_polarity(base_polarity)

    # Read dates from filenames, sort and save indices into the original (unsorted) list
    indate_list = read_infile_dates(coh_images) # list of primary and secondary image dates (no times); e.g. 20210301_20210313
    indate_arr = np.array(indate_list)
    coh_images_arr = np.array(coh_images)
    i_indate_sort_arr = np.argsort(np.array(indate_list))  # index of indates sorted chronologically
    indate_sort = indate_arr[i_indate_sort_arr] # sort date array chronologically
    coh_images_sort_arr = coh_images_arr[i_indate_sort_arr]
    img_date_sort = indate_sort.tolist()  # create a list of img dates sorted chronologically
    i_indate_sort = i_indate_sort_arr.tolist() # create a list of the indices of the img_dates in chronological order
    coh_images_sort = coh_images_sort_arr.tolist() # list of coh input files sorted chronologically
    
         
    print("Total Number of Coh images covering ROI:", len(coh_images_sort))
    single_image_dates = filter_single_image_dates(indate_sort) # list of primary and secondary dates (no times; e.g. 20210301_20210313) with single image acquisitions
    print("Number of single coh image days:", len(single_image_dates))
    
    ### Loop for single scene days; This is finding the correct image and associated filename, use this to iterate over single image days.
    for n in range(0, len(single_image_dates)):
        print(single_image_dates[n], "has only 1 SAR scene.") 
        idate = img_date_sort.index(single_image_dates[n])
        print(coh_images_sort[idate])
        print("")
        
        ### Image 1 (VV or HH)
        print("Cropping", coh_images_sort[idate]  ,"to :", minx, maxx, miny, maxy)
        out_tifname = single_image_dates[n] +'_' + args['input_image_suffix']
        out_tifname_and_path = os.path.join(cropped_outdirectory, out_tifname)
        print("Creating:", out_tifname_and_path); print("")
        crop_images_2_roi_wGdal(coh_images_sort[idate], out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)

        ### Write out log file for each image; turn into function. 
        out_text_file = out_tifname.replace('.tif','.txt')
        out_text_file_and_path = os.path.join(cropped_outdirectory, out_text_file)
        idate = [idate] # convert integer value to a list 
        write_merged_cropped_imaged_metafile(out_text_file_and_path, coh_images_sort, idate, cropped_outdirectory, out_tifname)
    
        del out_tifname, out_tifname_and_path, out_text_file_and_path, out_text_file



    # Loop for days with multiple scenes; Find images acquired on the same date (adjacent scenes),
    multi_image_dates = filter_non_unique(indate_sort) # list of primary and seconary image dates (e.g. 20210301_20210313) with multiple images
    print("Number of multi-image scenes:", len(multi_image_dates))
#    
    for n in range(0, len(multi_image_dates)):
        ## index into img_date_sort for ducpliate days
        idate = [i for i, x in enumerate(img_date_sort) if x == multi_image_dates[n]]
        print("There are ", len(idate), "SAR scenes for:", multi_image_dates[n])
        merged_path_and_file_0 = os.path.join(cropped_outdirectory, 'merged_0.tif')
        cmd_0 = "python3 " + gm + " -o " + merged_path_and_file_0 + " -n 0"  ## VV or HH image
            
        ##### Create a string of coh files to merge for the current day
        for nDuplicate in range(0, len(idate)):
            cmd_0 = cmd_0 + " " + coh_images_sort[idate[nDuplicate]]              
                
        ## Execute merge commands with gdal_merge
        print("Merging files to ", merged_path_and_file_0)
        print(cmd_0); print("")
        subprocess.call(cmd_0, shell=True)
           

        ## Crop input image to ROI and output with _roi.tif suffix (VV or HH)
        print("Cropping", merged_path_and_file_0, "to", minx, maxx, miny, maxy)
        out_tifname = multi_image_dates[n] +'_' + args['input_image_suffix']
        out_tifname_and_path = os.path.join(cropped_outdirectory, out_tifname)
        print("Creating:", out_tifname_and_path); print("")
        crop_images_2_roi_wGdal(merged_path_and_file_0, out_tifname_and_path, minx, maxx, miny, maxy)  # crop images to ROI (w/Buffer)

          ### Write out log file for each image; turn into function. 
        out_text_file = out_tifname.replace('.tif','.txt')
        out_text_file_and_path = os.path.join(cropped_outdirectory, out_text_file)
        write_merged_cropped_imaged_metafile(out_text_file_and_path, coh_images_sort, idate, cropped_outdirectory, out_tifname)
         
        
        ### Delete the merged coh temp file. 
        os.remove(merged_path_and_file_0)
        del out_tifname, out_tifname_and_path, out_text_file_and_path, out_text_file

             
#

#
#
#    ### Close log file
#    close_log_file(stdoutOrigin)
#   











