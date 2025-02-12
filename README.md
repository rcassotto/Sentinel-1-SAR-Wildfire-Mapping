# SAR-Wildfire-Mapping
Python based tools to map the evolution of wildfire burn areas. The tools use both SLC-derived coherence images and GRD-derived polarimeteric observations to produce burn area maps. A diagram of the full workflow is provided in the repository as a png file and included here for illustration.  

![FIREDpy-SAR Detection](https://github.com/user-attachments/assets/7df51b63-6817-4e90-9cf3-78dd15abd26e)

## Overview
The workflow contains 6 major steps with some requiring several smaller steps. Future iterations will consolidate the workflow into fewer steps. For now, the workflow consists of:
  1) Download SLC and GRD data for regions of interest (ROI) from the Alaska SAR Facility (ASF).
  2) Pre-process SLC data to Coherence images using ESA's SNAP tool.
  3) Generate coherence change images.
  4) Pre-process GRD files to Radiometrically Terrain Corrected (RTC) Sigma0 backscatter images.
  5) Generate polarimetric change data.
  6) Generate binary fire products and combine results.


## Initial Setup 
It is highly recommended to create a dedicated python environment to operate in. I use anaconda (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), but many tools are available. 

Once a dedicated environment is created and activated, install the necessary python dependencies (see below). 



## Python Dependencies
asf_search
fiona
gdal
geopandas
json
matplotlib
numpy
pandas
rasterio
scipy
shapely
simplekml
skimage


![FIREDpy-SAR Detection_zoom](https://github.com/user-attachments/assets/24a38cf5-5bbc-4071-9e8e-582626ec9bf0)




## Detailed workflow
## Download SLC or GRD data from ASF
Use ASF's vertex or API tools to download desired data.


