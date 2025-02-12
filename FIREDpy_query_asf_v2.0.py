## ASF's searchAPI solution with python wrapper

### rev02 only processes for a single, specified fid*json file
###         - also allows users to specify output directory and number of days to buffer fire observation start/end period

### rev01 adds buffered days at the start (sd) and end (ed) to get images before/after the fire; it also loops over all fid_json files in specified directory



import os
import asf_search as asf
import sys
from datetime import datetime
from datetime import timedelta
import argparse,ast
import geopandas as gpd

#infile_path='/Users/Ryan_old/Documents/1_Research/Projects/NASA_FIREDpy/ROIs/parsed_jsons'
#infile_path='/data2/NASA_FIREDpy/ROIs/parsed_jsons'
#infilename='fid_29.json'
#infile_path_and_filename=os.path.join(infile_path, infilename)



def make_output_dir(sar_download_path, sar_file_type, fid, climate, orbitDirection):
    ## Check outdir and create if not there 
#    outfile_path = infile_path.replace("site_geopackages","GRD_files")
    if not os.path.exists(sar_download_path):
        os.mkdir(sar_download_path)         

    sar_type_directory = os.path.join(sar_download_path, sar_file_type)
    if not os.path.exists(sar_type_directory):
        print("Making Parent Directory: ", sar_type_directory)
        os.mkdir(sar_type_directory)


    zip_file_directory = os.path.join(sar_download_path, sar_file_type,"fid_" + str(fid) + "_" + climate)
    if not os.path.exists(zip_file_directory):
        print("Making Parent Directory: ", zip_file_directory)
        os.mkdir(zip_file_directory)

    output_directory = os.path.join(sar_download_path, sar_file_type,"fid_" + str(fid) + "_" + climate, orbitDirection)
    if not os.path.exists(output_directory):
        print("Making Output Directory: ", output_directory)
        os.mkdir(output_directory)
    return output_directory

def parse_data_from_gpkg(data):
   ### Parse pertinent data
   geometry = data['geometry']
   minlon = geometry.total_bounds[0]
   minlat = geometry.total_bounds[1]
   maxlon = geometry.total_bounds[2]
   maxlat = geometry.total_bounds[3]
   
   fid = data['fid_1'][0]
   sd = data['ig_date'][0]
   ed = data['last_date'][0]
   climate = data['main_clim'][0]
   lc_name = data['lc_name'][0]
   print("start: ", sd) ; print("stop: ", ed); print("fid: ", fid); print("Main Climate: ", climate); print("Land Cover Name: ", lc_name)
   return minlon, maxlon, minlat, maxlat, fid, sd, ed, climate, lc_name
    

def open_log_file(fid, climate, output_directory, orbitDirection):
   stdoutOrigin = sys.stdout
   log_output_fname = "fid_" + str(fid) + "_" + climate + orbitDirection + "_asf_log.txt"
   log_output_fname_and_path = os.path.join(output_directory,log_output_fname)
   sys.stdout = open(log_output_fname_and_path, "w")
   return stdoutOrigin


def close_log_file(stdoutOrigin):
   sys.stdout.close()
   sys.stdout = stdoutOrigin
   return
   
def pull_asf_data(wkt_string, start_date, end_date, output_directory, orbitDirection):
   ### Extract USERNAME and PASSWORD from env
   pw = os.environ.get('ASF_API_PASS')
   un = os.environ.get('ASF_API_USER')
        
   # ###############################################################################
   #Query ASF and download data
   results = asf.geo_search(platform=[asf.PLATFORM.SENTINEL1], flightDirection=orbitDirection, intersectsWith=wkt_string, start=start_date, end=end_date, beamMode="IW", processingLevel="GRD_HD")
   session = asf.ASFSession()
   session.auth_with_creds(un, pw)
   print(results)  ## prints results to screen (or log file, if so specified)
   results.download(path=output_directory, session=session) # download results to a path


def download_asf_data_from_gpkg(infile_path_and_filename, sar_file_type, sar_download_path, nDay_orbital_buffer):
    
    ### Read in individual .gpkg file
    print('Reading input file: ', infile_path_and_filename)
    data = gpd.read_file(infile_path_and_filename)

     ### Parse pertinent data
    (minlon, maxlon, minlat, maxlat, fid, sd, ed, climate, lc_name) = parse_data_from_gpkg(data)
     
    ### Add buffer to start/end dates
    start_date = datetime.strftime(datetime.strptime(sd,"%Y-%m-%d") + timedelta(days=-nDay_orbital_buffer),"%Y-%m-%d")  # adds -48 days (3-4 orbits) to start date
    end_date = datetime.strftime(datetime.strptime(ed, "%Y-%m-%d") + timedelta(days=nDay_orbital_buffer),"%Y-%m-%d") # adds 48 days (3-4 orbits) to start date
     
    
    ## Get input path (for output directory)
    # infile_path = os.path.dirname(infile_path_and_filename)
    
    
     ### Set up ASF keyword variables common to both flightDirections (e.g. dates, ROI)   
    print("start: ", sd) ; print("stop: ", ed); print("fid: ", fid); print("Main Climate: ", climate)

    wkt_string = "polygon((" + str(maxlon) + " " + str(maxlat) + "," + str(maxlon) + " " + str(minlat) + ","  \
    + str(minlon) + " " + str(minlat) + "," + str(minlon) + " " + str(maxlat) + "," + str(maxlon) + " " + str(maxlat) + "))"
    print(wkt_string)
    
            
            
            ### Ascending
    orbitDirection = 'ASCENDING'
    output_directory_asc = make_output_dir(sar_download_path, sar_file_type, fid, climate, orbitDirection)  # create output directory
    stdoutOrigin = open_log_file(fid, climate, output_directory_asc, orbitDirection)    ### Setup log file
    pull_asf_data(wkt_string, start_date, end_date, output_directory_asc, orbitDirection)  ### Query ASF - Ascending
    close_log_file(stdoutOrigin)   ### Close log file
     
     
     ### Descending
    orbitDirection = 'DESCENDING'
    output_directory_des = make_output_dir(sar_download_path, sar_file_type, fid, climate, orbitDirection)  # create output directory
    stdoutOrigin = open_log_file(fid, climate, output_directory_des, orbitDirection)    ### Setup log file
    pull_asf_data(wkt_string, start_date, end_date, output_directory_des, orbitDirection)  ### Query ASF - Descending
    close_log_file(stdoutOrigin)   ### Close log file
     
        
    
    
    
    
    
if __name__ ==  "__main__":
    
    ### Read input parameter file (e.g. FIREDpy_query_ASF_input.txt)
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',type=argparse.FileType('r'))
    p = parser.parse_args()

    with p.filename as file:
        contents = file.read()
        args = ast.literal_eval(contents)
       
    ### Save input parameters into variables. 
   # json_dir = args['json_path'] # get location of SLC *.zip files from input.txt file
    #json_fname = args['json_file']
    gpkg_dir = args['gpkg_path']
    gpkg_fname = args['gpkg_file']
    
    sar_file_type = args['SAR_file_type']
    sar_download_path = args['SAR_download_path']
    nDay_orbital_buffer = args['OrbitalBuffer_days']  # Number of days to add to start/end dates of fire
    
    if len(nDay_orbital_buffer)==0:
        nDay_orbital_buffer = 48  # Default
    else:
        nDay_orbital_buffer = int(nDay_orbital_buffer)  # Convert nDay input 'string' to integer (for ASF API)
    

    infile_path_and_filename=os.path.join(gpkg_dir, gpkg_fname)

    download_asf_data_from_gpkg(infile_path_and_filename, sar_file_type, sar_download_path, nDay_orbital_buffer)
    
    
