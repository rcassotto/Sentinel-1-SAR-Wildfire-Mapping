from image_proc_module import Image_proc
import numpy as np
import matplotlib.pyplot as plt

class ARIA(Image_proc):

    def __init__(self, c_bounds):
        self.c_bounds = c_bounds

    """ Helper functions """

    # code adapted from https://github.com/mapbox/rio-hist
    def histogram_match(self, s_band, r_band, s_shape):
        s_band = np.ndarray.flatten(s_band) # flatten collapses array into one dimension
        r_band = np.ndarray.flatten(r_band)

        ### sort unique values of secondary band (and reference bad too)
        s_value, s_idx, s_count = np.unique(s_band, return_inverse=True, return_counts=True) 
                #### s_value: sorted unique values
                #### s_idx: indices to reconstruct original array from the unique array
                #### s_count: number of times each unique value appears in original array
        r_value, r_count = np.unique(r_band, return_counts=True)
       
        #### calculate the quantiles of each image by summing counts and dividing by the size of (or number of elements in) original image
        s_quantile = np.cumsum(s_count).astype(np.float64) / s_band.size
        r_quantile = np.cumsum(r_count).astype(np.float64) / r_band.size

        #### interpolate between the quantiles (historgrams) of each image, return 
        ##### y = np.interp(x,xp,fp) 
            #### y (output): interpolated values
            #### x: x-coords at which to evaluate interpolated values. 
            #### xp: x-coords of the data points. 
            #### fp: y-coordinates of data points.
        interp = np.interp(s_quantile, r_quantile, r_value) # interp is the interpolated values of s_quantile (secondary image histogram) evaluated at the ref image
        t_band = interp[s_idx] # t_band is the interpolated values of the secondary image, mapped back to the original secondary input coordinates
        t_band = t_band.reshape((s_shape[0], s_shape[1]))  # reshape this back to a 2D numpy array

        return t_band

    def plot_histograms(self, r_band, s_band, t_band, outfilename_and_path, n=256):
        r_band = np.ndarray.flatten(r_band)
        s_band = np.ndarray.flatten(s_band)
        t_band = np.ndarray.flatten(t_band)
        r_idxs = np.nonzero(r_band)[0]
        s_idxs = np.nonzero(s_band)[0]
        t_idxs = np.nonzero(t_band)[0]
        r_band = r_band[r_idxs]
        s_band = s_band[s_idxs]
        t_band = t_band[t_idxs]

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10,10))

        ax[0].hist(r_band, bins=n, density=True, alpha=0.5, color='r', label='Reference Histogram')
        ax[0].hist(s_band, bins=n, density=True, alpha=0.5, color='b', label='Secondary Histogram')
        ax[0].legend()
        ax[0].set_title("Before Matching - {} bins".format(n))

        # y0_lim=list(ax[0].axes.get_ylim())
        ax[0].set_xlim([-0.1, 1.1])  # RKC change: 20250204
        y0_lim=ax[0].axes.get_ylim()
        x0_lim=ax[0].axes.get_xlim()
        # y0_lim[1] = y0_lim[1] + 1
        # ax[0].axes.set_ylim(tuple(y0_lim))
        # print("YLIM: ", y0_lim, type(y0_lim))
        # print(y0_lim[1]*0.9)

        ### 20231130: RKC edit to print the number of elements less than, equal to, and greater than 0
        r_string = "Num r_band elements <0, =0, >0:", len(np.where(r_band < 0)[0]), len(np.where(r_band == 0)[0]), len(np.where(r_band > 0)[0])
        s_string = "Num s_band elements <0, =0, >0:", len(np.where(s_band < 0)[0]), len(np.where(s_band == 0)[0]), len(np.where(s_band > 0)[0])
        t_string = "Num t_band elements <0, =0, >0:", len(np.where(t_band < 0)[0]), len(np.where(t_band == 0)[0]),  len(np.where(t_band > 0)[0])
        ax[0].text(x0_lim[0]+0.2, y0_lim[1]*0.9, r_string, color = 'red')
        ax[0].text(x0_lim[0]+0.2, y0_lim[1]*0.8, s_string, color = 'blue')
        ######
        
        ax[1].hist(r_band, bins=n, density=True, alpha=0.5, color='r', label='Reference Histogram')
        ax[1].hist(t_band, bins=n, density=True, alpha=0.5, color='b', label='Target Histogram')
        ax[1].legend()
        ax[1].set_title("After Matching - {} bins".format(n))
        ax[1].set_xlim([-0.1, 1.1])  # RKC change: 20250204

        ##### RKC 20231130 edit
        x1_lim = ax[1].axes.get_xlim()
        y1_lim = ax[1].axes.get_ylim()
        # print("YLIM 1: ", y1_lim, type(y1_lim))
        ax[1].text(x1_lim[0]+0.1, y1_lim[1]*0.9, t_string)
        ########
        
#        plt.show()
        outfilename_and_path
        plt.savefig(outfilename_and_path)
        
        
        
        
        
        

    """ Some functions for making different kinds of coherence difference maps """

    def simple_map(self, arr, t=0):
        np.where(arr > t, arr, 0.)
        return arr 

    def binary_map(self, arr, t=0):
        flood = np.where(arr > t, 1., 0.).astype(np.int8)
        return flood

    def purple_map(self, arr, t=0): # assume causalty constraint
        rgba = np.zeros((arr.shape[0], arr.shape[1], 4)).astype(np.float32)
        pos = np.where(arr > 0.0)
        x = pos[0]
        y = pos[1]
        for i in range(len(x)):
            value = arr[x[i], y[i]]
            if value >= t:
                rgba[x[i], y[i], 3] = 1.0 # set completely opaque
                rgba[x[i], y[i], 0] = 0.5
                rgba[x[i], y[i], 2] = 1.0
        return rgba

    """ Driver to process using the ARIA method """

    # REQUIRED PARAMETERS
    # reference - Image object for the reference image
    # secondarys - list of Image objects for any and all secondary images
    # t - threshold used to determine whether loss of coherence is high enough to write to a map
    # map_type - choose a function to generate a numpy array that represents a colored or binary map for ARIA results

    # OPTIONAL PARAMETERS
    # file_prefix - prefix to all files being written to disk
    # show_hists - display plots for histogram matching

    # RETURNS
    # save_data - a list of SaveData objects
    def process_ARIA(self, reference, secondarys, png_outfilename_and_path, t, map_type, file_prefix='', show_hists=False):
        save_data= []
        r_str = reference.path.split('/')[-1] # grab file name
        r_str = r_str.split('.')[0] # remove file extension
        bounds = reference.c2b(self.c_bounds)
        for i in range(len(secondarys)):
            s_str = secondarys[i].path.split('/')[-1]
            s_str = s_str.split('.')[0]
            filename = file_prefix + '_' + r_str + '_to_' + s_str # not including file extension

            s_band, geo_bounds = self.crop(secondarys[i], bounds)
            r_band, _ = self.crop(reference, bounds)
            
            ## Convert Coh values == 0 to NaNs; RKC edit20231130
            r_band[r_band <= 0] = 'nan'
            s_band[s_band <= 0] = 'nan'
            
            t_band = self.histogram_match(s_band, r_band, s_band.shape)
            
            ### RKC addition, debug
            print(filename)

            if show_hists:
#                self.plot_histograms(r_band, s_band, t_band)
                self.plot_histograms(r_band, s_band, t_band, png_outfilename_and_path)
#                print('hi')
#            diff = np.subtract(r_band, t_band)  ## original; looks for decreases in coherence only (ignores coherence increase)
#            diff = np.absolute(np.subtract(r_band, t_band))     ## RKC mod to look for an increase in coherence (when secondary image has higher coherence than ref)
            
#            coh_map = map_type(diff, t)   
#            save_point = self.create_save(coh_map, filename, geo_bounds, reference.raster)
#            save_data.append(save_point)
##            
            #### Return t-band, the secondary image normalized to the reference image    
            save_point = self.create_save(t_band, filename, geo_bounds, reference.raster)
            save_data.append(save_point)
            
            

        return save_data
