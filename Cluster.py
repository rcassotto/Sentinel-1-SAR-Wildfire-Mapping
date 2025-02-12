#from source.Image import Image
from image_proc_module import Image_proc
import numpy as np
#import rasterio
import gc
import matplotlib.pyplot as plt


class cluster_image(Image_proc):

    def __init__(self, c_bounds):    # Initial function that gets run with the class
        self.c_bounds = c_bounds

    """ Create RGB composite """

    # both args are Image objects
    def construct_rgb_composite(self, VV_img, VH_img):
        # cropping
        bounds = VH_img.c2b(self.c_bounds)
        VH_arr, _ = VH_img.crop(VH_img, bounds)
        VH = Image_proc(None)
        #VH.band = VH_arr

        bounds = VV_img.c2b(self.c_bounds)
        VV_arr, _ = VV_img.crop(VV_img, bounds)
        VV = Image_proc(None)
        #VV.band = VV_arr
        #print("VV band type:", type(VV), type(VV.band))
        
        
#         ##RKC mod - 20240208
        # band_shape = np.shape(VV_arr)
        # b1_flat = VV_arr.flatten()
        # i_nan = np.argwhere(b1_flat == 0)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
        # print("VV inan", i_nan)
        # b1_flat[i_nan] = 0
        # b1_raster = b1_flat.reshape(band_shape)
        # VV.band = b1_raster
        # del b1_flat, b1_raster, i_nan, band_shape
        
        # band_shape = np.shape(VH_arr)
        # b1_flat = VH_arr.flatten()
        # i_nan = np.argwhere(b1_flat == 0)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
        # print("VH inan", i_nan)
        # b1_flat[i_nan] = 0
        # b1_raster = b1_flat.reshape(band_shape)
        # VH.band = b1_raster
        # del b1_flat, b1_raster, 
# ####

###### new mod  -- This is much improved!!!
        # VH_arr = np.where(VH_arr == 0., 0.1, VH_arr) # handle divide by zero # original
        # VV_arr = np.where(VV_arr == 0., 0.1, VV_arr) # handle divide by zero # original
        # VH.band = VH_arr
        # VV.band = VV_arr
#####
        ###### this changes zeros to nans; RKC mod 20240209
        VH_arr = np.where(VH_arr == 0., np.nan, VH_arr) # handle divide by zero # original
        VV_arr = np.where(VV_arr == 0., np.nan, VV_arr) # handle divide by zero # original
        ##### This masks the nans
        VH_mask = np.ma.array(VH_arr, mask=np.isnan(VH_arr))
        VV_mask = np.ma.array(VV_arr, mask=np.isnan(VV_arr))

        VH.band = VH_mask # VH data as numpy array that ignores NaNs
        VV.band = VV_mask # VV data as numpy array that ignores NaNs

###

        # convert linear to dB
        VH.convert_to_dB()
        VV.convert_to_dB()

        # quantile clipping
        VH.quantile_clip(upper_quantile=0.99)
        VV.quantile_clip(upper_quantile=0.99)

        # map to interval [0,1]
        VH.map_to_interval(0, 1)
        VV.map_to_interval(0, 1)

        # make the ratio band
        ratio_band = np.divide(VH.band, VV.band)
        #ratio_band = np.where(VV.band == 0., 0.1, ratio_band) # handle divide by zero # original
      
         
      ### RKC mod 20240209
        #ratio_band = np.where(VH_arr == 0, np.nan, ratio_band) # handle divide by zero # original
        # band_shape = np.shape(ratio_band)
        # rb_1d = ratio_band.flatten();
        # i_nan = np.argwhere(rb_1d == 1)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
        # print("ratio inan", i_nan)
        # rb_1d[i_nan] = np.nan
        # ratio_band = rb_1d.reshape(band_shape)
        
        
        ratio_mask = np.ma.array(ratio_band, mask=np.isnan(ratio_band))
        # ratio_mask = np.ma.array(ratio_band, mask=(ratio_band==1))
        del ratio_band
        ratio_band = ratio_mask
        ####
        
        ## Add code to remask or re-nan no data values AFTER the normalization. 20240212
        band_shape = np.shape(ratio_band)
        rb_1d = ratio_band.flatten();
        i_nan = np.argwhere(rb_1d == np.isnan)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
        print("ratio inan", i_nan)
        # rb_1d[i_nan] = np.nan
        ratio_band = rb_1d.reshape(band_shape)
        
        
        
        #ratio_band = np.where(VH_band == 0., 0., ratio_band)
        ratio = Image_proc(None)
        ratio.band = ratio_band
        ratio.map_to_interval(0, 1)
        
        ## Add code to remask or re-nan no data values AFTER the normalization. 20240212
        # band_shape = np.shape(ratio_band)
        rb_1d = ratio_band.flatten();
        # i_nan = np.argwhere(rb_1d == 1)  # identify pixels common in bands 1 and 2 that are equal to 1 (no data)
        # print("ratio inan", i_nan)
        rb_1d[i_nan] = np.nan
        ratio_band = rb_1d.reshape(band_shape)
        ratio.band = ratio_band
               
        

        # return a numpy array because Image objects are only tested for single channel images
        rgb = np.zeros((ratio.band.shape[0], ratio.band.shape[1], 3))
        rgb[:,:,0] = VH.band
        rgb[:,:,1] = VV.band
        rgb[:,:,2] = ratio.band

        return rgb

    def plot_rgb_histograms(self, rgb, n=256):
        vh = np.ndarray.flatten(rgb[:,:,0])
        vv = np.ndarray.flatten(rgb[:,:,1])
        rat = np.ndarray.flatten(rgb[:,:,2])
        vh_idxs = np.nonzero(vh)[0]
        vv_idxs = np.nonzero(vv)[0]
        r_idxs = np.nonzero(rat)[0]
        vh = vh[vh_idxs]
        vv = vv[vv_idxs]
        rat = rat[r_idxs]

        fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10,10))

        ax[0].hist(vv, bins=n, density=True, alpha=0.5, color='r')
        ax[0].set_title("VV (red)")
        ax[1].hist(vh, bins=n, density=True, alpha=0.5, color='g')
        ax[1].set_title("VH (green)")
        ax[2].hist(rat, bins=n, density=True, alpha=0.5, color='b')
        ax[2].set_title("VH/VV (blue)")

        ax[3].hist(vv, bins=n, density=True, alpha=0.5, color='r')
        ax[3].hist(vh, bins=n, density=True, alpha=0.5, color='g')
        ax[3].hist(rat, bins=n, density=True, alpha=0.5, color='b')
        ax[3].set_title("All together")

        plt.show()

    """ Driver to process ISODATA """

    # based on http://web.pdx.edu/~jduh/courses/Archive/geog481w07/Students/Vassilaros_ISODATA.pdf

    # img - band of an Image object
    def isodata(self, img, random_init=False, start_num=20, max_iter=10, std_threshold=0.1, merge_threshold=0.1):

        # initialize starting cluster centers
        if random_init: centers = [np.random.uniform(size=3) for _ in range(start_num)]
        else:
            rows = [np.random.choice(range(img.shape[0]), replace=False) for _ in range(start_num)]
            cols = [np.random.choice(range(img.shape[1]), replace=False) for _ in range(start_num)]
            centers = [img[row, col] for row,col in list(zip(rows, cols))]
        center_mats = []
        for center in centers:
            center_mat = np.ones((img.shape[0], img.shape[1], 3))
            center_mat[:,:,:] = center
            center_mats.append(center_mat)
        dist_mats = np.zeros((img.shape[0], img.shape[1], len(centers)))
        for k in range(len(centers)):
            center_mat = center_mats[k]
            dist_mat = np.linalg.norm(img-center_mat, axis=2) # Euclidean distance between each pixel and a certain cluster
            dist_mats[:,:,k] = dist_mat
        assignments = np.argmin(dist_mats, axis=2)
        cluster_ids = np.unique(assignments)
        # free memory
        del center_mats
        del dist_mats
        gc.collect()

        # iterate
        count = 0
        cont = True
        while cont:
            # perform merge test
            n = len(cluster_ids)
            neg_count = -1
            for i in range(n):
                for j in range(i+1,n):
                    id1 = cluster_ids[i]
                    id2 = cluster_ids[j]
                    if id1 > 0 and id2 > 0: # if the clusters have not already been merged
                        c1 = centers[id1]
                        c2 = centers[id2]
                        inter_dist = np.linalg.norm(c1-c2)
                        if inter_dist < merge_threshold: # merge clusters that have centers that are close together
                            print('Merging clusters %s and %s into cluster %s' % (id1, id2, neg_count))
                            cluster_ids[i] = neg_count
                            cluster_ids[j] = neg_count
                            assignments = np.where(assignments == id1, neg_count, assignments)
                            assignments = np.where(assignments == id2, neg_count, assignments)
                            neg_count -= 1

            # perform standard deviation test and re-define centers
            centers = []
            cluster_ids = np.unique(cluster_ids) # collapse duplicates from merged clusters
            for cluster_id in cluster_ids:
                rows, cols = np.where(assignments == cluster_id)
                pixels = img[rows, cols, :]
                std = np.std(pixels, axis=0, ddof=1)
                test_std = np.where(std > std_threshold, 1, 0)
                if np.sum(test_std) > 0:
                    # split the cluster
                    print('Splitting cluster number %s with std %s' % (cluster_id, np.max(std)))
                    faulty_dim = np.argmax(std) # choose the dimension with the highest std
                    dim_vec = np.array(pixels[:, faulty_dim])
                    cut_off = np.quantile(dim_vec, 0.5)
                    l = np.where(dim_vec <= cut_off)[0]
                    sub_pix1 = pixels[l,:]
                    u = np.where(dim_vec > cut_off)[0]
                    sub_pix2 = pixels[u,:]
                    centers.append(np.mean(sub_pix1, axis=0))
                    centers.append(np.mean(sub_pix2, axis=0))
                else:
                    centers.append(np.mean(pixels, axis=0))

            # re-assign the clusters
            center_mats = []
            for center in centers:
                center_mat = np.ones((img.shape[0], img.shape[1], 3))
                center_mat[:,:,:] = center
                center_mats.append(center_mat)
            dist_mats = np.zeros((img.shape[0], img.shape[1], len(centers)))
            for k in range(len(centers)):
                center_mat = center_mats[k]
                dist_mat = np.linalg.norm(img-center_mat, axis=2) # Euclidean distance between each pixel and a certain cluster
                dist_mats[:,:,k] = dist_mat
            assignments = np.argmin(dist_mats, axis=2)
            cluster_ids = np.unique(assignments)
            # free memory
            del center_mats
            del dist_mats
            gc.collect()

            # stop if maximum iterations is reached
            count += 1
            if count >= max_iter: cont = False

        return assignments


    def k_means(self, img, k, random_init=False, max_iter=20):
        # initialize starting cluster centers
        if random_init: centers = [np.random.uniform(size=3) for _ in range(k)]
        else:
            rows = [np.random.choice(range(img.shape[0]), replace=False) for _ in range(k)]
            cols = [np.random.choice(range(img.shape[1]), replace=False) for _ in range(k)]
            centers = [img[row, col] for row,col in list(zip(rows, cols))]
        center_mats = []
        for center in centers:
            center_mat = np.ones((img.shape[0], img.shape[1], 3))
            center_mat[:,:,:] = center
            center_mats.append(center_mat)
        dist_mats = np.zeros((img.shape[0], img.shape[1], len(centers)))
        for k in range(len(centers)):
            center_mat = center_mats[k]
            dist_mat = np.linalg.norm(img-center_mat, axis=2) # Euclidean distance between each pixel and a certain cluster
            dist_mats[:,:,k] = dist_mat
        assignments = np.argmin(dist_mats, axis=2)
        cluster_ids = np.unique(assignments)
        # free memory
        del center_mats
        del dist_mats
        gc.collect()

        count = 0
        while count < max_iter:
             # re-assign the clusters
            center_mats = []
            for center in centers:
                center_mat = np.ones((img.shape[0], img.shape[1], 3))
                center_mat[:,:,:] = center
                center_mats.append(center_mat)
            dist_mats = np.zeros((img.shape[0], img.shape[1], len(centers)))
            for k in range(len(centers)):
                center_mat = center_mats[k]
                dist_mat = np.linalg.norm(img-center_mat, axis=2) # Euclidean distance between each pixel and a certain cluster
                dist_mats[:,:,k] = dist_mat
            assignments = np.argmin(dist_mats, axis=2)
            cluster_ids = np.unique(assignments)
            # free memory
            del center_mats
            del dist_mats
            gc.collect()

            count += 1

        return assignments
