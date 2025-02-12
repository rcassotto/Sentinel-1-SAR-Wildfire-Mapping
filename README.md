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
It is highly recommended to create a dedicated python environment to operate in. Many python package tools are available; I use anaconda (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Once a dedicated environment is created and activated, install the necessary python dependencies (see below). 



## Python Dependencies
asf_search<br> fiona<br> gdal<br> geopandas<br> json<br> matplotlib<br> numpy<br> pandas<br> rasterio<br> scipy<br> shapely<br> simplekml<br> skimage<br>

## Detailed workflow


![FIREDpy-SAR Detection_zoom_step1](https://github.com/user-attachments/assets/b793ad49-adf6-4923-8bcf-0b096ecf739e)


![FIREDpy-SAR Detection_zoom_step2](https://github.com/user-attachments/assets/f25e28ba-a318-4c00-8505-d4d234bcf83a)


![FIREDpy-SAR Detection_zoom_step3](https://github.com/user-attachments/assets/083d0fb7-8586-4a60-9536-cf25c1c3c6ce)


![FIREDpy-SAR Detection_zoom_step4](https://github.com/user-attachments/assets/1e42a47d-0fd7-4a86-b051-09a41caecb1d)


![FIREDpy-SAR Detection_zoom_step5](https://github.com/user-attachments/assets/4000ded6-a3a5-47b8-b206-0827b684766a)




![FIREDpy-SAR Detection_zoom_step6](https://github.com/user-attachments/assets/1bed0a59-f4f5-4b3a-832a-bfec1bdf2534)



## Download SLC or GRD data from ASF
Use ASF's vertex or API tools to download desired data.


