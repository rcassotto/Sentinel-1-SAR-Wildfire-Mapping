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
  The Geospatial Data Abstraction Library is called within the python framework. **_If you installed gdal with python in the   evironment above, you can copy and paste the full path to the gdal_translate binary where indicated in the python scripts   below. You do not have to install gdal seperately_.** 

  #### _SNAP_
  ESA's SNAP program is required to pre-process the SLC and GRD data to coherence and sigma0 images, respectively. Follow the instructions on ESA's website to install SNAP (https://step.esa.int/main/download/snap-download/). Once installed, note the full path of the gpt tool. You will need to enter it in the necessary scripts identified below.
  
  Note: A SNAP update may be necessary, even after a fresh install. 
    E.g. /local_path_to_snap/snap/bin/snap --nosplash --nogui --modules --update-all


## Detailed workflow

As mentioned, a png file of the full workflow is provide in the repository above. The following provides the individual steps of the workflow. 

### Step 1 - Download SAR data from ASF
<br>
For the purposes of this demonstration, a geopackage is provided (NIFC_2024_TX_Windy_Deuce.gkpg). The python script _FIREDpy_query_asf_v2.0.py_ will read an input file (_FIREDpy_query_ASF_input.txt_), parse the input argurments and download the data to the user specified output directory. As shown in the figure below, the user should modify the following inputs in the file prior to running the script. 

  - **gpkg_path**: full path where the geopackage file is stored.
  - **gpkg_file**: the geopackage filename (e.g. NIFC_2024_TX_Windy_Deuce.gkpg).
  - **SAR_file_type**: 'SLC' for coherence products, 'GRD' for polarimeteric Sigma0 files.
  - **SAR_download_path**: full path where the API should download SAR files to.
  - **OrbitalBuffer_days**: Number of days to add before the start date and after the end date for the desired time period of interest.
<br><br>

![FIREDpy-SAR Detection_zoom_step1](https://github.com/user-attachments/assets/b793ad49-adf6-4923-8bcf-0b096ecf739e)

<br><br>

  If no geopackge is available for your ROI, you can specify a wkt polygon and start/end dates in the API. An example for such case is provided here >> https://github.com/rcassotto/Sentinel-1-GRD-to-RTC-Pre-Processing/tree/main/ASF_API.

This step requires an EarthData account and credentials. Be sure to add your Earth Data account credentials to your bashrc file by executing the following commands in a terminal:      
         **_export ASF_API_PASS=&lt;password>_**   
         **_export ASF_API_USER=&lt;username>_**  
   _Note &lt;username> and &lt;password> should be replaced with the Earth Data account credentials without angled brackets but with single quotes (‘ ‘)._

To execute script:
1) Activate your dedicated python environment, if not already in the environment.
2) initiate script: **_python3 FIREDpy_query_asf_v2.0.py NIFC_2024_TX_Windy_Deuce.gkpg_**
     

<br><br><br>

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



![FIREDpy-SAR Detection_zoom_step3](https://github.com/user-attachments/assets/083d0fb7-8586-4a60-9536-cf25c1c3c6ce)


<br><br><br>
### Step 4 - Pre-process GRD files to RTC Sigma0 Images


![FIREDpy-SAR Detection_zoom_step4](https://github.com/user-attachments/assets/1e42a47d-0fd7-4a86-b051-09a41caecb1d)

<br><br><br>
### Step 5 - Generate Polarimetric Change Data

![FIREDpy-SAR Detection_zoom_step5](https://github.com/user-attachments/assets/4000ded6-a3a5-47b8-b206-0827b684766a)

<br><br><br>
### Step 6 - Generate Binary Fire Products and Combine Results

![FIREDpy-SAR Detection_zoom_step6](https://github.com/user-attachments/assets/1bed0a59-f4f5-4b3a-832a-bfec1bdf2534)


<br><br>
