# SAR-Wildfire-Mapping
Python based tools to map the evolution of wildfire burn areas. The tools use both SLC-derived coherence images and GRD-derived polarimeteric observations to produce burn area maps. A diagram of the full workflow is provided above in the repository as a png file. Segments of it are provided below to detail the steps.   

## Overview
The workflow contains 6 major steps with some requiring several smaller steps. Future iterations will consolidate the workflow into fewer steps. For now, the workflow consists of:
  1) Download SLC and GRD data for regions of interest (ROI) from the Alaska SAR Facility (ASF).
  2) Pre-process SLC data to Coherence images using ESA's SNAP tool.
  3) Generate coherence change images.
  4) Pre-process GRD files to Radiometrically Terrain Corrected (RTC) Sigma0 backscatter images.
  5) Generate polarimetric change data.
  6) Generate binary fire products and combine results.


## Initial Setup 
Download or clone the repository to your local machine. 

It is highly recommended to create a dedicated python environment to operate in. Many python package tools are available; I use anaconda (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Once a dedicated environment is created and activated, install the necessary python dependencies (see below).

Once the python environment is configured, install additional non-python dependencies before running the workflow. 

### Python Dependencies
asf_search<br> fiona<br> gdal<br> geopandas<br> json<br> matplotlib<br> numpy<br> pandas<br> rasterio<br> scipy<br> shapely<br> simplekml<br> skimage<br>

### Non-Python Dependencies
Two additional dependecies are required to run the workflow; both are open sourced. 
- gdal.
- The SeNtinel Application Platform (SNAP).
- EarthData account (it's free!)

  #### _GDAL_
  The Geospatial Data Abstraction Library is called within the python framework. **_If you installed gdal with python in the   evironment above, you can copy the full path to the installed binary for gdal_translate and paste it where indicated in the python scripts below. You do not have to install gdal seperately_.** 

  #### _SNAP_
  ESA's SNAP program is required to pre-process the SLC and GRD data to coherence and sigma0 images, respectively. Follow the instructions on ESA's website to install SNAP (https://step.esa.int/main/download/snap-download/). Once installed, note the full path of the gpt tool. You will need to enter it in the necessary scripts identified below.
  
  _Note: A SNAP update may be necessary, even after a fresh install. <br> E.g._
       **_/local_path_to_snap/snap/bin/snap --nosplash --nogui --modules --update-all_**


## Detailed workflow

As mentioned, a png file of the full workflow is provide in the repository above. The following provides the individual steps of the workflow. 

### Step 1 - Download SAR data from ASF
<br>
For the purposes of this demonstration, a geopackage is provided (NIFC_2024_TX_Windy_Deuce.gkpg). The python script _FIREDpy_query_asf_v2.0.py_ will read an input file (_FIREDpy_query_ASF_input.txt_), parse the input argurments and download the data to the user specified output directory. As shown in the figure below, the user should modify the following inputs in the file prior to running the script. 

  - _gpkg_path_: full path where the geopackage file is stored.
  - _gpkg_file_: the geopackage filename (e.g. NIFC_2024_TX_Windy_Deuce.gkpg).
  - _SAR_file_type_: 'SLC' for coherence products, 'GRD' for polarimeteric Sigma0 files.
  - _SAR_download_path_: full path where the API should download SAR files to.
  - _OrbitalBuffer_days_: Number of days to add before the start date and after the end date for the desired time period of interest.
<br>

![FIREDpy-SAR Detection_zoom_step1](https://github.com/user-attachments/assets/b793ad49-adf6-4923-8bcf-0b096ecf739e)

<br>

  If no geopackge is available for your ROI, you can specify a wkt polygon and start/end dates in the API. An example for such case is provided here >> https://github.com/rcassotto/Sentinel-1-GRD-to-RTC-Pre-Processing/tree/main/ASF_API.

This step requires an EarthData account and credentials. Be sure to add your Earth Data account credentials to your bashrc file by executing the following commands in a terminal:      
         **_export ASF_API_PASS=&lt;password>_**   
         **_export ASF_API_USER=&lt;username>_**  
   _Note &lt;username> and &lt;password> should be replaced with the Earth Data account credentials without angled brackets but with single quotes (‘ ‘)._

To execute script:
1) Activate your dedicated python environment, if not already in the environment.
2) initiate script: **_python3 FIREDpy_query_asf_v2.0.py NIFC_2024_TX_Windy_Deuce.gkpg_**
   
<br><br>
### Step 2 - Pre-Process SLC data to Coherence Images
<br>
Once the data download is complete, continue to step 2: pre-process SLC to coherence images. This step should be completed before or in parallel with pre-process GRD to Sigma0 images (step 4). It must be completed before step 5 as the workflow relies on coherence change images (step 3) to generate polarimetric differenced images. 
<br>

![FIREDpy-SAR Detection_zoom_step2](https://github.com/user-attachments/assets/f25e28ba-a318-4c00-8505-d4d234bcf83a)

<br>

This step will create coherence images for all SLC files in the user specified directory. Consequently, this step can take several hours or days to complete, depending on the output resolution and number of input files. As with the data download step, users should modify an input text file (e.g. FIREDpy_process_coh_input_asc.txt) with their associated input values.  The fire_roi_polygon and fire_roi_path input arguments are not used in the current version, but will be used in future revisions. The sys_index_var should be 0 for the initial run.  If failures persist midway through batch processing coherence images, this value should be changed to reflect the next iteration of coherence image. For example, if 10 coherence images successfully completed and the 11th image failed, sys_index_var should be set to 11 to continue the batch process with the 11th iteration of coherence image pairs. The command below will generate coherence and intensity images for each coherence pair. 

Make the following changes prior to running this script for the first time:
  1) Open the script _Process_SLC2Coh_wSNAP_v1.0.py_ with a python editor.
  2) Perform a search and replace for the following fields
      - _base_snap_exe_no_aux_: full path for the gpt binary
      - _base_snap_exe_: full path for the gpt binary
      - _gdal_exe_: full path to gdal_translate 
      - _aux_poe_: full path for Sentinel-1 precise orbits
      - _res_orb_: full path for Sentinel-1 restituted orbits
      - _workflow_dir_: full path where the files in this repository were downloaded or cloned to.

Orbit Files: are critical to processing Sentinel-1 SLC or GRD data. As written, the script will try  to automatically download Sentinel-1 orbit files to the user's dedicated orbit directories. However, the location of the Copernicus host changed in 2023. As of Feb 2025, you can obtain Sentinel-1 orbit files through one of two methods below. Be sure to download the necessary orbit files prior to executing this script. 
    - A clever pythonized workaround from Scott Stanie (https://github.com/scottstanie/sentineleof/tree/master).
    - Through AWS: https://registry.opendata.aws/s1-orbits/.

Once correct pathways are updated in your python script, the input file is amended for your area of interest, and orbit files are downloaded, execute the python script with the following command:
  **_python3 Process_SLC2Coh_wSNAP_v1.0.py FIREDpy_process_coh_input_asc.txt_**

The coherence files will be located in a _Processed_data_ directory, with each image further located in subfolders organized by the date and times of the primary and secondary image pairs. A shell script (_move_tif_files2_coherence_tif_copies_dir.sh_) is provided to move these files into a new centralized directory called _coherence_tif_copies_ with the _Processed_data_ directory (e.g. _/Processed_data/coherence_tif_copies/_); it also removes the extraneous data generated during the pre-processing. 

The script will need to be executed from within the _Processed_data_ directory. 

<br><br><br>


### Step 3 - Generate Coherence Change Images

Two different scripts are used to generate coherence change images: 

  - _merge_crop_coherence_images_2ROI_v1.0.py_: crops the Sigma0 scenes to the geometry defined in the geopackage. It applies a user defined buffer to be added to the region as a percentage of the total bounding box of the geometry within the geopackage. When two or more scenes appear have the same date (e.g. adjacent Sentinel-1 paths), the script will merge these scenes into a single.  The script outputs a geotiff and a metafile with each having a date string prefix (e.g. YYYYMMDD_coherence.tif).
    
  - _FIREDpy_SAR_coh_change_rev05.py_: This script normalizes the stack of coherence images to a single image to reduce differences in coherence scenes between dates. It outputs the normalized coherence image, the coherence change image, log files, and a histogram of the normalization process.

    
   #### Initial Setup
  Make the following changes prior to running the scripts for the first time:
    1) Open the script _merge_crop_coherence_images_2ROI_v1.0.py_ with a python editor.
    2) Replace the path for gdal_merge.py with the full path on your local machine. The line to search for beging with "gm = ". Save and close the script in the python editor.
  
  #### Executing Step 3
  1) Use a text editor to make the following changes to the input file to merge and crop the coherence images (e.g. _CA_Monument_merge_crop_coh_input_des.txt_)
      - _coh_file_loc_: full path to the coherence image files from Step 2.
      - _out_dir_: full path for the user defined merged/cropped coherence output images.
      - _input_image_suffix_: ending suffix for the input files (e.g. 'coherence_VV.tif')
      - _fire_roi_gis_file_: filename for the geopackage defining your region of interest.
      - _fire_roi_gis_path_: full path for the geopackage defining your region of interest.
      - _fire_perimeter_buffer_prct_: numerical value representing a buffer as a percent of the image size based on the geometric bounding box defined in the geopackage.
  2) Execute the merge/crop step: **_python3 merge_crop_coherence_images_2ROI_v1.0.py CA_Monument_merge_crop_coh_input_des.txt_**
  3) Use a text editor to make the following changes to the input file to merge and crop the coherence images (e.g. _CA_Monument_FIREDpy_coh_change_input_des.txt_)
      - _coh_file_loc_: full path to the merged/cropped images above. 
      - _out_dir_: full path for coh change output files. 
      - _tif_suffix_: suffix for the input files (e.g. 'coherence_VV.tif').
      - _ref_image_: reference image for the normalization process. If not specified, the program will auto-select.
  4) Excecute the Coh Change script: **_python3 FIREDpy_SAR_coh_change_rev05.py CA_Monument_FIREDpy_coh_change_input_des.txt_**
        
The figure below illustrates the individual steps outlined above. 

![FIREDpy-SAR Detection_zoom_step3](https://github.com/user-attachments/assets/083d0fb7-8586-4a60-9536-cf25c1c3c6ce)


<br><br><br>
### Step 4 - Pre-process GRD files to RTC Sigma0 Images

The _RTC_V3.py_ script will use SNAP to pre-process GRD to RTC Sigma0 geotiffs. Specifically, it applies an orbit correction, removes border noise, calibrates the data, applies a speckle filter, terrain corrects and generates Sigma0 geotiffs; it also removes ancillary products generated during the pre-processing steps. 

  #### Initial Setup
  Make the following changes prior to running this script for the first time:
    1) Open the script _RTC_V3.py_ with a python editor.
    2) Perform a search and replace for the following fields
        - "usr/local/bin/gdal_translate" with the location of gdal_translate on your machine.
        - baseSNAP: replace the current path with the location of where the SNAP gpt binary is located.
  
  #### Executing Step 4
  1) Use a text editor to open the input file: _rtc_sample_inputs.txt_.
  2) Amend the inputs for your configuation
       - _DEM_: DEM path and filename (optional); the program will default to SRTM if not specified.
       - _grd_file_loc_: full path of the GRD zip files.
       - _output_dir_: full path of the output directory for the output files.
       - _pixsize_: desired output pixel size in meters.
  3) Initiate the python script: **_python3 RTC_V3.py rtc_sample_inputs.txt_**

  The program will generate several output files as shown in the figure below. 
        
![FIREDpy-SAR Detection_zoom_step4](https://github.com/user-attachments/assets/1e42a47d-0fd7-4a86-b051-09a41caecb1d)


<br><br><br>
### Step 5 - Generate Polarimetric Change Data

The generation of polarimetric change images includes several individual steps as outlined in the figure below. 
<br>

![FIREDpy-SAR Detection_zoom_step5](https://github.com/user-attachments/assets/4000ded6-a3a5-47b8-b206-0827b684766a)
<br><br>

  #### Merge/Crop Sigma0 Images
  _merge_crop_sigma0_images_2ROI_v1.0.py_: This script is similar to the merge/crop coherence script, but does so for the Sigma0 images. Like the merge/crop coherence module, the full path to gdal_merge will need to be updated in the python   script to reflect the location on your local machine.

  It utilizes an input file (e.g. _merge_sigma0_input.txt_) with the following inputs:
      - _sigma0_file_loc_: full path to the Sigma0 image files from Step 4.
      - _out_dir_: full path for the user defined merged/cropped sigma0 output images.
      - _input_image_suffix_: suffix for the input files (e.g. 'Sigma0_VV.tif'; the code will automatically search for VH counterparts)
      - _fire_roi_gis_file_: filename for the geopackage defining your region of interest.
      - _fire_roi_gis_path_: full path for the geopackage defining your region of interest.
      - _fire_perimeter_buffer_prct_: numerical value representing a buffer as a percent of the image size based on the geometric bounding box defined in the geopackage.

  Command to execute: **_python3 merge_crop_sigma0_images_2ROI_v1.0.py merge_sigma0_input.txt_**

  #### Make polarimetric images
  _SAR_Classification_VV_VH_v4.0.py_ creates a 3-band false color composite comprised of VV, VH, and VH/VV.  It also creates a single-band luminance image based on the false color composite. 

  The command utilizes an input file (e.g. _WindyDeuce_class_asc_input.txt_) with the following inputs:
      - _sigma0_file_loc_: full path to the cropped/merged Sigma0 image files. 
      - _out_dir_: full path for false color and luminance output images.
      - _sigma_ref_image_": user defined reference image for normalization. If none is specified, the script will auto select. 
      - _use_norm_images_: a boolean argument to specify whether the Sigma0 images should first be normalized. 

  _SAR_polarimetric_RVI_RFDI_v1.0.py_: creates radar vegetation indexed (RVI) images and radar forest degradation indexed (RFDI) images using the polarimetric Sigma0 files. It uses input files with similar input requirements as outlined in figure above. 
  
Once the appropriate input data has been modified for your data set, the commands can be executed by calling the python scripts with input files. For example, 
  **_python3 SAR_Classification_VV_VH_rev4.0.py WindyDeuce_class_asc_input.txt_**
  **_python3 SAR_polarimetric_RVI_RFDI_v1.0.py WindyDeuce_rvi_rfdi_input_asc.txt_**


#### Difference polarimetric images 
_make_lumDelta_images_v2.0.py_ differences the luminance data derived above. It pairs images based on the coherence change image dates; that is, it creates a DeltaLumin geotiff for each coherence change image.  

It utilizes an input file (e.g. _FIREDpy_input_lumDelta_ref_file.txt_) with the following inputs:
      - _lum_file_loc_: full path to the luminance image files generated above.
      - _out_dir_: full path for DeltaLumin output files.
      - _lumin_suffix_: suffix for the input files (e.g. 'RGBLumin.tif')
      - _coh_change_file_loc_: full path for the coherence change images.
      - _coh_change_suffix_": suffix of coherence change images (e.g. 'cohChange.tif')
      - _lum_ref_image_: filename for a reference luminance image for normalization.
      - _use_norm_images_: boolean input for normalizing luminance images prior to differencing. 

The _make_DeltaRVI_images_v2.0.py_ and _make_DeltaRFDI_images_v2.0.py_ scripts are very similar to _make_lumDelta_images_2.0.py_.  Their input files and values therein are also similar. The full input arguements are provided in the figure depicting Step 5 above. 

Once the appropriate input data has been modified for your data set, the commands can be executed by calling the python scripts with input files. For example, 
  **_python3 make_lumDelta_images_v2.0.py FIREDpy_input_lumDelta_ref_file.txt_**
  **_python3 make_DeltaRVI_images_v2.0.py FIREDpy_input_DeltaRVI_ref_file.txt_**
  **_python3 make_DeltaRFDI_images_v2.0.py FIREDpy_input_DeltaRFDI_ref_file.txt_**

<br><br><br>
### Step 6 - Generate Binary Fire Products and Combine Results

The final step in the workflow is to generate binary output products and combine the results using the script _combine_S1_lumChange_cohChange_products.py_. 

The input file contains several input arguments, including:
<br> &nbsp;      - _DeltaLum_path_: full path to the luminance change files.
<br> &nbsp;        - _coh_change_path_: full path to the coherence change files.
<br> &nbsp;        - _dRVI_path_: full path to the RVI change files.
<br> &nbsp;        - _dRFDI_path_: full path to the RFDI change files. 
<br> &nbsp;        - _DeltaLum_suffix_: suffix for the luminance change files (e.g. 'DeltaLumin.tif').
<br> &nbsp;        - _coh_change_suffix_: suffix for the coherence change files (eg. 'cohChange.tif').
 <br> &nbsp;       - _dRVI_suffix_: suffix for the RVI change files (e.g. 'DeltaRVI.tif').
<br> &nbsp;        - _dRFDI_suffix_: suffix for the RFDI change files (e.g. 'DeltaRFDI.tif').
<br> &nbsp;        - _out_dir_: full path for output products. 
<br> &nbsp;        - _coh_threshold_: coherence change cutoff; program will autoselect if none specified.
<br> &nbsp;        - _lum_threshold_: luminance change cutoff; program will autoselect if none specified.
<br> &nbsp;        - _dRVI_threshold_: dRVI change cutoff; program will autoselect if none specified.
<br> &nbsp;        - _dRFDI_threshold_: dRFDI change cutoff; program will autoselect if none specified.
<br> &nbsp;        - _filter_Delta_images_: boolean argument for spatially filtering the polarimetric images (e.g. DeltaLumin.tif, DeltaRVI.tif, DeltaRFDI.tif)
<br> &nbsp;        - _filter_Coh_images_: boolean argument for spatially filtering coherence change images. 

Once the appropriate input arguments have been modified for your data set, the command can be executed by calling the python script with an input file. For example, <br>
   &nbsp; **_python3 combine_S1_lumChange_cohChange_products.py FIREDpy_combine_SAR_products_WindyDeuce_asc.txt_**
<br><br>
![FIREDpy-SAR Detection_zoom_step6](https://github.com/user-attachments/assets/1bed0a59-f4f5-4b3a-832a-bfec1bdf2534)

<br><br>
