#%% Import necessary Libraries

import numpy as np
from osgeo import gdal
import glob, os
import logging
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.plot import show
from pathlib import Path
import cv2
from astropy.convolution import convolve
import warnings
import time
from scipy.spatial import cKDTree
import geopandas as gpd 
from shapely.geometry import box 
from matplotlib import pyplot

gdal.UseExceptions()

print('Libraries loaded')

#%% INPUT

# Unziped Download Folder from GFM:
input_dir = Path(r'C:/UrbanFloods/SriLanka/srilanka_2024-10-14T00_25_26_2024_11_20T13_39_01_033640/') 
output_dir = Path(r'C:/UrbanFloods/SriLanka/output-v2/')

dtm_path = 'C:/UrbanFloods/SriLanka/fabdem_srilanka.tif'

# Option to specify AOI if not specified, tiles are mosaiced
roi_path = 'C:/UrbanFloods/SriLanka/srilanka_roi.gpkg' 

likelihood = True # True, if likelihood layer should be used as flood layer
likelihood_threshold = 40 # in Percent

if likelihood == True:
    print('Likelihood Layer will be thresholded with ' + str(likelihood_threshold) + '%.')

dst_crs = 'none' # if no crs is specified, rasters will be reprojected into EPSG 3857

param_tiling      = False    # True to run FLEXTH on the tiled inputs
param_tile_inputs = True   # True to tile the inputs. If input is already tiled, select False. Relevant if "param_tiling = True"
param_tile_size   = 20000   # size of the squared tiles (pixels)

print('Inputs are specified') # Print what the inputs are

# For further details on the following parameters see https://doi.org/10.5194/nhess-24-2817-2024

# Water level estimation method (options: 'method_A', 'method_B'): 
param_WL_estimation_method = 'method_A'  


# Select output: Water depth ("WD") , Water level ("WL"), both ("WL_WD")
param_output_map = "WL_WD"


# Parameters
param_threshold_slope           =  0.1     # S_max : border pixels steeper than this (D_z/D_x) are not used to estimate water level
param_size_gaps_close           =  0.01    # A_g   : up to this size (in km2) gaps in the flood map will be closed

param_border_percentile         =  0.5     # P: assign water level based on the percentile P of the elevation of border cells (valid if Method B is selected)
param_max_number_neighbors      =  100     # N_max: number of border pixels used to compute water level at a single location
param_min_flood_size            =  10     # N_min: if the number of valid pixels along the border is less than this, WL is estimated based on the distribution of the elevation of the pixels inside the flooded area
param_inner_quantile            =  0.98    # P*: for  flooded areas that don't meet the criteria above, uses this percentile of the elevation of the pixels inside the flooded area to estimate water levels
param_inverse_dist_exp          =  1       # alpha: inverse distance weighting exponent used to interpolate WL inside flooded areas 


param_max_propagation_distance  =  10    # D_max: maximum propagation distance in km
param_distance_range            =  10    # A_1/2: flooded areas of this size (km2) reach half of the maximum propagation distance  

param_WD_star                   =  10    # WD*: dummy water depth in cm assigned if estimated WL<DTM in initially delineated flooded areas


#%% Mosaicing 
def mosaic_tiles(input_dir, pattern):
    id = pattern.split('_')[1].split('*')[0]
    q = os.path.join(input_dir, pattern)
    flood_rst = glob.glob(q)
    
    rst_lst = []
    
    for fp in flood_rst:
        src = rasterio.open(fp)
        rst_lst.append(src)
        
    mosaic, out_trans = merge(rst_lst)
    
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                     "transform": out_trans,
                     "crs": src.crs
                     }
                  )
    
    with rasterio.open(input_dir / id +'_mosaic.tif', "w", **out_meta) as dest:
        dest.write(mosaic)

# File Path Patterns
lyr_ids = ['*ENSEMBLE_UNCERTAINTY*.tif', '*ENSEMBLE_EXCLAYER*.tif', '*ENSEMBLE_OBSWATER*.tif', '*ENSEMBLE_FLOOD*.tif']

# Mosaicing
for i in range(len(lyr_ids)):
    mosaic_tiles(input_dir, lyr_ids[i])


#%% Reprojections

# check if mosaics are already in the folder

# Checkpoint for Reprojections - if they are even necessary.
if dst_crs == 'none':
    dst_crs = 'EPSG:3857'

# Threshold the Uncertainty layer
if likelihood == True:
    with rasterio.open(input_dir / 'UNCERTAINTY_mosaic.tif') as src:
        rst = src.read(1)
        likelihood_thr = (rst > likelihood_threshold).astype(np.uint8) 
        nodata = src.nodata
        likelihood_thr[rst == nodata] = nodata
    
        out_meta = src.meta.copy()
        out_meta.update(dtype=rasterio.uint8)

        with rasterio.open(input_dir / 'Likelihood_thr.tif', "w", **out_meta) as dst:
            dst.write(likelihood_thr, 1)
    
    input_flood_delineation = input_dir / 'Likelihood_thr.tif'
    print('Likelihood Layer successfully thresholded.')
    
else:
    input_flood_delineation = input_dir / 'FLOOD_mosaic.tif'

# reproject the flood raster
with rasterio.open(input_flood_delineation) as src:
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)
    kwargs = src.meta.copy()
    raster_bounds = src.bounds
    kwargs.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height
    })
    output_filename = "flood.tif" if roi_path is None else "flood_reproj.tif"
    with rasterio.open(input_dir / output_filename, 'w', **kwargs) as dst:
        for i in range(1, src.count + 1):
            reproject( 
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)
    
# Crop the flood raster to the ROI
if roi_path != None:
    print('Flood Layer is clipped to AOI...')
    roi = gpd.read_file(roi_path)
    roi_reproj = roi.to_crs(dst_crs) 
    
    with rasterio.open(input_dir / 'flood_reproj.tif') as src:
        raster_bounds = src.bounds
        raster_gdf = gpd.GeoDataFrame({"geometry": [box(*raster_bounds)]}, crs = dst_crs)
        roi_crp = gpd.clip(roi_reproj, raster_gdf) # has to be cropped to the extent of the input rasters
        
        # checkpoint if roi intersects the mosaic
        if roi_crp.empty:
            raise ValueError("The AOI is empty. Please check the input AOI file.")
            
        out_image, out_transform = mask(src, roi_crp.geometry, crop = True)
        
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
            })
        with rasterio.open(input_dir / 'flood.tif', "w", **out_meta) as dst:
            dst.write(out_image)

print('Flood Layer reprojected.')

#%%
flood_dataset = rasterio.open(input_dir / 'flood.tif')
flood_array   = flood_dataset.read(1)

transform = flood_dataset.transform
flood_array_nrows, flood_array_ncol =  flood_array.shape
Dx= transform[0]
Dy= transform[4]
minX= transform[2]
maxY= transform[5]
maxX= minX + Dx*flood_array_ncol
minY= maxY + Dy*flood_array_nrows

proj = flood_dataset.crs.wkt
outputBounds=[minX, minY, maxX, maxY]

# function to reproject all input rasters other than the input flood delineation
def reproj(input_path, output_path, proj, outputBounds, continuous_input = True):
    
    raster_dataset = rasterio.open(input_path)
    input_proj       = raster_dataset.crs.wkt
    output_proj      = proj

    if continuous_input == True:
        resampling_method = 'bilinear'# other options: 'average', 'cubic', 'lanczos' ... 
    elif continuous_input == False:
        resampling_method = 'mode'   # other options: 'nearest' ...
    
    gdal.Warp(str(output_path), 
              str(input_path), 
              outputBounds=outputBounds,
              cropToCutline    = True,
              outputBoundsSRS = output_proj, 
              warpMemoryLimit= 5000,
              srcSRS= input_proj, 
              dstSRS= output_proj, 
              xRes=Dx,
              yRes=Dy, 
              resampleAlg= resampling_method,
              targetAlignedPixels = False,
              creationOptions = ["COMPRESS=ZSTD", "BIGTIFF=YES", "TILED=YES"]    # compression options: ZSTD, DEFLATE, LZW
              )
    
    #flood_dataset.close()
    raster_dataset.close()

reproj(dtm_path, input_dir / 'dtm.tif', proj, outputBounds, continuous_input = True)
reproj(input_dir / 'EXCLAYER_mosaic.tif', input_dir / 'exclusion.tif', proj, outputBounds, continuous_input = False)
reproj(input_dir / 'OBSWATER_mosaic.tif', input_dir / 'obswater.tif', proj, outputBounds, continuous_input = False)

#%% FLEXTH

##############
#FUNCTIONS######
##################

#tiles and input_raster raster into output_dir
def tiling(input_raster, source_dir, param_tile_size):
    from rasterio.windows import Window
    
    file_name = os.path.splitext(os.path.basename(input_raster))[0]
    
    output_folder = source_dir / 'input_tiled'
    
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    
    # Open the GeoTIFF file
    with rasterio.open(input_raster) as src:
        # Get the raster size
        rastersize_width, rastersize_height = src.width, src.height
    
    
        # Calculate the number of tiles in x and y directions
        num_tiles_x = (rastersize_width  + param_tile_size - 1) // param_tile_size
        num_tiles_y = (rastersize_height + param_tile_size - 1) // param_tile_size
    
    
        # Loop through each tile and extract it from the original raster
        for i in range(num_tiles_x):
            for j in range(num_tiles_y):
                # Calculate the tile's bounding box
                x_off = i * param_tile_size
                y_off = j * param_tile_size
                x_size = min(param_tile_size, rastersize_width  - x_off)
                y_size = min(param_tile_size, rastersize_height - y_off)
    
                # Create a window object for the tile
                window = Window(x_off, y_off, x_size, y_size)
    
                # Read the tile from the original raster
                tile_data = src.read(1, window=window)
    
                # Create a new GeoTIFF file for the tile
                tile_profile = src.profile
                tile_profile.update({
                    'width': x_size,
                    'height': y_size,
                    'transform': rasterio.windows.transform(window, src.transform),
                    'compress': 'deflate'
                })
                with rasterio.open(output_folder/f'{file_name}_tile_{i+1}_{j+1}.tif', 'w', **tile_profile) as dst:
                    dst.write(tile_data, 1)
    
    


    
#  identifies the indices of neighboring cells in an array 
#  i,j : indices of the target location; 
#  n_row and n_col:  number or rows and columns of the  target array 
#  connectivity: type of connectivity to use (4 , 8)
 
def ij_neighbors(i,j,n_row,n_col,connectivity):  
    if connectivity == 4:
        neighbors=np.array([ [i-1,j], [i,j+1], [i+1,j],  [i,j-1]])    
    elif connectivity == 8:
        neighbors=np.array([[i-1,j-1],[i-1,j], [i-1,j+1], [i,j+1], [i+1,j+1], [i+1,j], [i+1,j-1], [i,j-1]])  
    elif connectivity == 24:
        neighbors=np.array([[i-2, j-2],[i-2, j-1], [i-2,j], [i-2,j+1], [i-2,j+2], [i-1, j-2], [i-1,j-1], [i-1,j], [i-1,j+1], [i-1,j+2], [i,j-2], [i,j-1], [i,j+1], [i,j+2], [i+1,j-2], [i+1,j-1], [i+1,j], [i+1,j+1],[i+1,j+2],[i+2,j-2],[i+2,j-1], [i+2,j], [i+2,j+1], [i+2,j+2]  ]     ) 
    else:
        raise ValueError(f"Unsupported connectivity value: {connectivity}")
    
        
    for indice in range(len(neighbors)):
        if neighbors[indice,0]> n_row-1 or neighbors[indice,0]<0 :
            neighbors[indice,0]=i
            neighbors[indice,1]=j
        if neighbors[indice,1]> n_col-1 or neighbors[indice,1]<0 :
            neighbors[indice,0]=i
            neighbors[indice,1]=j                     
            
    return neighbors.astype('uint32')    
   
    


# function that computes the weighted quantiles out of a PDF
def weighted_quantile(values, percentile, weights=None):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: 
    :computes one percentile across the columns of a 2d array
    :percentile should be in [0, 1]!
    :param values: numpy.array with data
    :param percentile: the value of the quantile
    :param sample_weight: array-like of the same length as `array`
    :return: numpy.array with computed percentile for each row of data.
    """
    values = np.array(values)
    percentile = np.array(percentile)
    if weights is None:
        weights= np.ones(np.shape(values))
    weights = np.array(weights)
    assert np.all(percentile >= 0) and np.all(percentile <= 1), \
        'percentile should be in [0, 1]'


    sorter = np.argsort(values)
    values = np.take_along_axis(values, sorter, axis = 1)   
    weights = np.take_along_axis(weights, sorter, axis = 1)      

    weighted_quantiles = np.cumsum(weights, axis = 1 ) - 0.5 * weights

    weighted_quantiles /= np.sum(weights, axis = 1)[:, np.newaxis] * np.ones ( np.shape(weighted_quantiles))
    
    
    return values [np.array( np.arange(0,np.shape(values)[0],1)  ), np.argmin( (weighted_quantiles - percentile  )**2, axis = 1) ] 


##################
##################
# main function ##
##################
##################

def flood_processing(flood_path, dtm_path, exclusion_path, obswater_path, permanent_water_path):
       
    flood_georeferenced  =   rasterio.open(  flood_path )
    flood                =   rasterio.open(  flood_path ).read(1).astype('uint8')
    dtm                  =   rasterio.open(  dtm_path   )
    dtm_nodata           =   dtm.nodata
    dtm                  =   dtm.read(1)
    
    dtm[dtm == dtm_nodata]  = np.nan
    
    
    
    
    #   checks if optional inputs are provided:
    # - exclusion mask : areas masked from flood mapping
    # - observed water : flood water + permanent waters
    # - permanent water: permanent and/or seasonal water bodies
    
    if os.path.isfile(exclusion_path):
        exclusion      =  rasterio.open(  exclusion_path ).read(1).astype('uint8')
    else:
        exclusion      = np.full_like(flood, 0)
        
        
    if os.path.isfile(obswater_path):
        obswater      =  rasterio.open(  obswater_path ).read(1).astype('uint8')
    else:
        obswater      =  np.copy(flood)
    
    
    if os.path.isfile(permanent_water_path):
        permanent_water      =  rasterio.open(  permanent_water_path ).read(1).astype('uint8')
    else:
        permanent_water      =  np.full_like(flood, 0)
    
    
    
    
    #geotranformation
    trans =  flood_georeferenced.transform
    trans*(0,0)  # goes from COL-ROW index to X-Y:  trans*(ncol,nrow)=(x_coordinate,y_coordinate)
    
    #extract the dimension of the raster
    size_equi7_tile = dtm.shape    
    n_row, n_col    = size_equi7_tile 
    
    # pixel size in m
    L = np.abs(np.array(trans*(0,1))[1]-np.array(trans*(0,0))[1])  
     
          
    #local slope
    slope_max =  np.gradient( dtm, L)
    slope_max =  np.sqrt(slope_max[0]**2 + slope_max[1]**2 )   
    
    slope_steep   =   np.array(slope_max >  param_threshold_slope).astype('uint8')
    
    
    
    #############################################################################
    ##connected components anaysis with cv2 - identifies contiguous flooded areas        
    #############################################################################
    
    #makes sure there are no values other than 1 and 0
    exclusion[(exclusion !=1) & (exclusion !=0)]                     = 1  
    flood[(flood !=1) & (flood !=0)]                                 = 0
    obswater[(obswater !=1) & (obswater !=0)]                        = 0
    permanent_water[(permanent_water !=1) & (permanent_water !=0)]   = 0
    
     
    
    #some kernels 
    kernel_1       = np.ones((3,3),np.uint8)
    kernel_1_cross = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype = np.uint8)
    kernel_2       = np.ones((5,5),np.uint8)
    
    

    #runs morphological closing to remove small holes and very irregular borders from flood and obswater
    flood    = cv2.morphologyEx(flood   ,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
    obswater = cv2.morphologyEx(obswater,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
     
    #dilate exclusion mask
    exclusion = cv2.dilate(exclusion  ,  kernel_1)
    
    
    

    #############################################
    #closes small holes in flood map with  cv2 ##
    #############################################
    print ('Closing small gaps in flood map...')
    
    threshold_n_pixels = (param_size_gaps_close*1e6 / L**2).astype(int)
    
    binary_map_complementary_flood               = np.logical_not(flood).astype('uint8')
    z_num_labels, z_labels, z_stats, z_centroids = cv2.connectedComponentsWithStats(binary_map_complementary_flood,connectivity=8) 

    for i in range(len(z_stats)):
        if np.abs(z_stats[i,4]) < threshold_n_pixels:
            flood    [ z_stats[i,1]: z_stats[i,1] + z_stats[i,3] , z_stats[i,0]:z_stats[i,0] + z_stats[i,2] ] [z_labels[ z_stats[i,1]: z_stats[i,1] + z_stats[i,3] , z_stats[i,0]:z_stats[i,0] + z_stats[i,2] ] ==i] = 1
            obswater [ z_stats[i,1]: z_stats[i,1] + z_stats[i,3] , z_stats[i,0]:z_stats[i,0] + z_stats[i,2] ] [z_labels[ z_stats[i,1]: z_stats[i,1] + z_stats[i,3] , z_stats[i,0]:z_stats[i,0] + z_stats[i,2] ] ==i] = 1
    
    
    
    
    
    
    if os.path.isfile(input_dir / 'permanent_water.tif'):
        permanent_water = cv2.morphologyEx(permanent_water,    cv2.MORPH_CLOSE, kernel_1_cross,  iterations = 2 )
        
    else:
        permanent_water = obswater - flood
        permanent_water[(permanent_water!=0) & (permanent_water!=1) ] = 0
    


    # !!!  opencv with connectivity = 4 has problems in some versions of the package better use connectivity = 8 !!!!
           
    # Applying cv2.connectedComponents() - possibly change the connectivity type 
    # z_stats contains: [xtop, ytop, xwidth, ywidth, #pixels] - first element of stats correspondes to background value
    
    z_num_labels, z_labels, z_stats, z_centroids = cv2.connectedComponentsWithStats(flood,connectivity=8)   
    
      
    flood_dilated          =   cv2.dilate(flood  ,  kernel_1)
    flood_eroded           =   cv2.erode (flood  ,  kernel_1_cross)
    
    exclusion_dilated      =   cv2.dilate(exclusion ,  kernel_1)

    permanent_water_dilated =   cv2.dilate(permanent_water ,  kernel_1)
    
    

    
    flood_border           =   flood_dilated    - flood_eroded 
    flood_border_inner     =   flood            - flood_eroded
    
    
    slope_steep_dilated      =   cv2.dilate(slope_steep  ,  kernel_1)
    

    exclusion_dilated      [ (exclusion       ==0)  & (exclusion_dilated       ==1)  & (flood==0) & (flood_dilated==1)] = 0
    permanent_water_dilated[ (permanent_water ==0)  & (permanent_water_dilated ==1)  & (flood==0) & (flood_dilated==1)] = 0
    slope_steep_dilated    [ (slope_steep      ==0) & (slope_steep_dilated      ==1) & (flood==0) & (flood_dilated==1)] = 0
    
    
    flood_border          =     flood_border *  np.logical_not( exclusion_dilated )   *  np.logical_not( permanent_water_dilated)   *    np.logical_not(slope_steep_dilated)  
    
    dtm_border_flood      =   dtm * flood_border
      
    
   
    
    ## AVERAGES THE DEM ALONG THE BORDER 
    dtm_border_flood_smooth= convolve ( dtm_border_flood, kernel_2, mask = flood_border == 0, preserve_nan = True )
    dtm_border_flood_smooth[flood_border==0] = np.nan     
    
    
    flood_border =   flood_border * flood_border_inner


    #######################################################
    ## WATER LEVEL ESTIMATION IN INITIALLY FLOODED AREAS ##
    #######################################################
    
    
    z_water_level = np.full(size_equi7_tile, np.nan, dtype='float32')
    
    ROW, COL   =  np.mgrid[0:n_row,0:n_col].astype('uint32')
    
    row_flood                          =   np.ma.masked_array(ROW,      z_labels==0).compressed()
    col_flood                          =   np.ma.masked_array(COL,      z_labels==0).compressed()
    flood_labels_compressed            =   np.ma.masked_array(z_labels, z_labels==0).compressed()
        
    row_flood_border                   =   np.ma.masked_array(ROW,                     flood_border==0).compressed()
    col_flood_border                   =   np.ma.masked_array(COL,                     flood_border==0).compressed()        
    flood_borders_labels_compressed    =   np.ma.masked_array(z_labels,                flood_border==0).compressed()
    dtm_border_flood_smooth_compressed =   np.ma.masked_array(dtm_border_flood_smooth, flood_border==0).compressed()

    
    #  assigns water levels to initial flooded areas based on the DTM values along the borders
    #  of flooded areas that are not contiguous with the exclusion mask

     
    print('Assigns water elevation to initial flooded areas')
    

    for i in range (1,len(z_stats)):

        temp_position_flood           =  np.stack( ( row_flood[flood_labels_compressed==i],col_flood[flood_labels_compressed==i]),axis = 1)
        temp_position_border          =  np.stack( ( row_flood_border[flood_borders_labels_compressed==i],col_flood_border[flood_borders_labels_compressed==i]),axis = 1)
        temp_dtm_border_flood_smooth  =  dtm_border_flood_smooth_compressed[flood_borders_labels_compressed==i]
        
        len_border  = len(temp_position_border)
 

        if  len_border < param_min_flood_size: 
                
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                z_water_level[temp_position_flood[:,0], temp_position_flood[:,1] ]  =  np.nanquantile(dtm[temp_position_flood[:,0], temp_position_flood[:,1] ], param_inner_quantile)
    
    
        else:   
    
            

            
            # Build a cKDTree for each set of points
            tree_BORDER = cKDTree(temp_position_border)
    
            param_workers_cKDTree =  8
            distances, indices = tree_BORDER.query(temp_position_flood, k=param_max_number_neighbors,workers = param_workers_cKDTree)

            distances = distances.astype('float32')
            distances[distances==0]  = 1
            indices   = indices.astype('uint32')


            if distances.ndim == 1:
                distances = np.expand_dims(distances, axis=1)
                indices   = np.expand_dims(indices, axis=1)



            if param_max_number_neighbors > len_border:
                indices   = indices  [:, :len_border]
                distances = distances[:, :len_border]


            
            #####################
            #####################
            ## choose the weights
            #####################
            #####################
                
                
            ################################
            # inverse distance weighting IDW
            ################################
                
            weights  =   1 / (distances ** param_inverse_dist_exp).astype('float32')

                
              
             ##############################
             ##############################
             ## choose WL estimation method
             ##############################
             ##############################


            if param_WL_estimation_method == 'method_A':
            #METHOD A : water level is interpoleated using IDW or EDW weights
                    
                
                dtm_temp_statistic    =  (np.sum((temp_dtm_border_flood_smooth[indices] * weights).astype('float32'), axis = 1 )  / np.sum( weights , axis = 1 ) ).astype('float32')
                
                
            elif param_WL_estimation_method == 'method_B':
                #METHOD B : water level is a distence-weighted QUANTILE of all closest dtm neighboring cells 
        
                    
                dtm_temp_statistic    =   weighted_quantile(temp_dtm_border_flood_smooth[indices], param_border_percentile, weights  =  weights)
                             
                    
    
    
                
            z_water_level[temp_position_flood[:,0], temp_position_flood[:,1]]    =    dtm_temp_statistic 
    
    
     
    

    ###############
    #################
    ##FLOOD EXPANSION
    ###################
    #####################
 
    if os.path.isfile(exclusion_path) and np.sum(flood)>0: 
     
        
        print('Augmenting flooded areas...')
        
    
        param_connectivity                 =  8         # <<INPUT connectivity used to spread flooded areas
    
         
     
        z_water_level_masked_compressed         =   np.ma.masked_array(z_water_level, flood_border_inner==0).compressed() 
        z_water_level_masked_compressed_initial =   np.ma.masked_array(z_water_level, flood_border_inner==0).compressed()
        row_flood_border                        =   np.ma.masked_array(ROW,           flood_border_inner==0).compressed()
        col_flood_border                        =   np.ma.masked_array(COL,           flood_border_inner==0).compressed()
        flood_labels_border                     =   np.ma.masked_array(z_labels,      flood_border_inner==0).compressed()
        
        FLOOD = np.array ( [flood_labels_border,  row_flood_border,  col_flood_border,  z_water_level_masked_compressed, z_water_level_masked_compressed_initial ]  ).T
        
        #sorts based on decreasing water elevation
        FLOOD_sort_descending= np.flipud(FLOOD[FLOOD[:, 3].argsort()])
        
        #size of initial flooded areas in m2
        dim_flooded_area = z_stats[:,4] * L**2
        dim_flooded_area[dim_flooded_area<0]=0  #  prevents potential negative values caused by overflow 
        
        
        ########################
        # max expansion distance  
        ########################
       

        ##EXPONENTIAL
        threshold_distance = param_max_propagation_distance * 1000   *  ( 1 -   pow(2,   -  dim_flooded_area   / (param_distance_range * 1e6  )  )        )
        
        
        #avoids singularities in case no flood expansion is set
        threshold_distance = threshold_distance + 1
        
        
            
        FLOOD_sort_descending_list= list (FLOOD_sort_descending)
          
        z_water_level_augmented = np.copy(z_water_level)
        
        
        z_labels_augmented      = np.copy(z_labels)
        temp_dist_from_border   = np.zeros(z_labels.shape, dtype = 'float16')
        
        
        
        i=0
        i_counter  =  len ( FLOOD_sort_descending_list )
        i_stop     =  len ( FLOOD_sort_descending_list )
        total_processed = 0
        
        
        
           
        while i < i_counter: 
            
            if i == i_stop:
                FLOOD_sort_descending = np.vstack(FLOOD_sort_descending_list)
                FLOOD_sort_descending = FLOOD_sort_descending[i_stop: ,:]
                FLOOD_sort_descending = np.flipud(FLOOD_sort_descending[FLOOD_sort_descending[:, 3].argsort()]) 
                FLOOD_sort_descending_list= list (FLOOD_sort_descending)
                i=0
                i_stop    = len(FLOOD_sort_descending)
                i_counter = len(FLOOD_sort_descending) 
        
        
            vicini =  ij_neighbors(FLOOD_sort_descending[i,1], FLOOD_sort_descending[i,2] , n_row, n_col, param_connectivity)  
             
         
        
            label     = FLOOD_sort_descending[i,0].astype('uint32')
            
            
            initial_wl= FLOOD_sort_descending[i,4]  
            
    
            temp_= ( L +  temp_dist_from_border[FLOOD_sort_descending[i,1].astype('uint32'), FLOOD_sort_descending[i,2].astype('uint32')]  ) / threshold_distance[label]
            
            wl = (  initial_wl  -  ( initial_wl  - dtm[vicini[:,0],vicini[:,1]]   )  * temp_  )     * np.heaviside(1-temp_ , 0 )   +    dtm[vicini[:,0],vicini[:,1]]  *     (  1- np.heaviside(1-temp_ , 0 )    )
           
            wl[wl>FLOOD_sort_descending[i,3]] = FLOOD_sort_descending[i,3]
            
            
                                                                                                                    
            condition =  (   
                
                         np.heaviside( (exclusion[vicini[:,0],vicini[:,1]] > 0)  + (obswater[vicini[:,0],vicini[:,1]] > 0) + (permanent_water[vicini[:,0],vicini[:,1]] > 0)   , 0).astype('bool')          # 1. must be in the exlusion or permanent water layer
    
                       * ( dtm[vicini[:,0],vicini[:,1]] <  (wl-0.01) ) # FLOOD_sort_descending[i,3] #  z_flooded_areas_augmented[i][2]                   # 2. dtm must be lower than the corresponding water level    
                     
                       * np.heaviside(  np.isnan( z_water_level_augmented[vicini[:,0],vicini[:,1]])         , 0  ).astype('bool')                        # 3. must be nodata (i.e. not alreay assigned)
    
                         )
        
           
            temp_dist_from_border[vicini[:,0][condition],vicini[:,1][condition]]    = L * np.array([1.414, 1, 1.414, 1, 1.414, 1, 1.414, 1])[condition]     + temp_dist_from_border[FLOOD_sort_descending[i,1].astype('uint32'), FLOOD_sort_descending[i,2].astype('uint32')]
            
            z_water_level_augmented[vicini[:,0][condition],vicini[:,1][condition]]  =  wl[condition]   
            
            FLOOD_sort_descending_list.append  (  np.array([ label*(condition[condition]) ,  vicini[:,0][condition],  vicini[:,1][condition]  ,  wl[condition] , initial_wl*(condition[condition])  ]).T  ) 
            
            z_labels_augmented[vicini[:,0][condition],vicini[:,1][condition]]  =  label 
            
            i_counter+= np.sum(condition)
                        
                                 
            i+=1 
            
            #print(i)
            total_processed +=1
        
    
        

        ##############################
        ## smoothing expanded WL ####
        ##############################
    
         
        print('Smoothing... \n')
               
        
        param_number_smoothings      = 20   #<<<<<<  INPUT
        
        
        z_water_level_augmented_initial = np.copy(z_water_level_augmented)
        z_water_level_augmented[ np.isnan(z_water_level_augmented_initial) ] = dtm[  np.isnan(z_water_level_augmented_initial)  ]
        
        
        target_i, target_j  =  np.where( (~np.isnan(z_water_level_augmented_initial)) & (flood == 0) )
        target_i = target_i.astype(np.int32) ; target_j = target_j.astype(np.int32)
        
        
        
        connectivity = 24      #<<<<<<  INPUT
        
        
        neighbors_i = np.zeros((target_i.shape[0], connectivity)).astype(np.uint32)
        neighbors_j = np.zeros((target_i.shape[0], connectivity)).astype(np.uint32)
        
          
        
        for i in range(target_i.shape[0]):
            
            neighbors = ij_neighbors(target_i[i], target_j[i], n_row, n_col, connectivity)
            neighbors_i[i, :] = neighbors[:, 0]
            neighbors_j[i, :] = neighbors[:, 1]
                    
        
        
        mask_1 = np.where(~np.isnan(z_water_level))

        
        for v in range(param_number_smoothings):
            z_water_level_augmented[target_i,target_j] = np.mean(z_water_level_augmented[neighbors_i,neighbors_j], axis = 1 )
            
            z_water_level_augmented[ mask_1]  = z_water_level[  mask_1  ]
        
        z_water_level_augmented[ np.isnan(z_water_level_augmented_initial) ] = np.nan
        

    else:
        z_water_level_augmented = np.copy(z_water_level)

    

     
     
    ######
    # WD #
    ######
    WD                        = ( ( z_water_level_augmented - dtm ) * 100  ).astype('float32')   
    
    WD[WD<0] = 0 
    WD[WD > 0]                = WD[WD > 0] + param_WD_star
    WD[(WD==0) & (flood>0)  ] = param_WD_star
    
    
    
    #dummy water depth assigned to permanent water bodies
    WD[permanent_water==1] = 9999 
    
    WD[np.isnan(WD)] = 0
    
    
    
    ######
    # WL #
    ######
 
    WL          = ( z_water_level_augmented   ).astype('float32')   
    
    WL[WD > 0]  = WL[WD > 0] + param_WD_star/100
    WL[WD==0 ]  = np.nan
    
    
    
    #dummy water level assigned to permanent water bodies
    WL[permanent_water==1] = 9999 
    
    
    
    
    ########
    ##SAVING
    ########
    
    print('Saving \n')
    
    prefix  = os.path.splitext(os.path.basename(flood_path))[0]
    
    if param_output_map == "WD" or  param_output_map == "WL_WD" :
    
        with rasterio.open(
            output_dir / f'_WD_{prefix}_{param_WL_estimation_method}_Smax_{param_threshold_slope}_Nmax_{param_max_number_neighbors}_a_{param_inverse_dist_exp}_Dmax_{param_max_propagation_distance}_A12_{param_distance_range}_gaps_{param_size_gaps_close}.tif',
            mode="w",
            driver="GTiff",
            compress='deflate',
            height=n_row,
            width=n_col,
            count=1,
            dtype=rasterio.uint16,
            crs=flood_georeferenced.crs,
            transform=flood_georeferenced.transform,
            nodata= 0
            ) as water_depth:
                water_depth.write(WD, 1)
                
                
    if param_output_map == "WL" or  param_output_map == "WL_WD" :
    
        with rasterio.open(
            output_dir / f'_WL_{prefix}_{param_WL_estimation_method}_Smax_{param_threshold_slope}_Nmax_{param_max_number_neighbors}_a_{param_inverse_dist_exp}_Dmax_{param_max_propagation_distance}_A12_{param_distance_range}_gaps_{param_size_gaps_close}.tif',
            mode="w",
            driver="GTiff",
            compress='deflate',
            height=n_row,
            width=n_col,
            count=1,
            dtype=rasterio.float32,
            crs=flood_georeferenced.crs,
            transform=flood_georeferenced.transform,
            ) as water_level:
                water_level.write(WL, 1)
    
   
    
    print('Finish \n')
      
      
    
    
    flood_georeferenced.close()
    
# Function to print and calculate statistics on processed flood map
def FLEXTH_stats():
    with (
    rasterio.open(
        #output_dir / f'_WL_flood_{param_WL_estimation_method}_Smax_{param_threshold_slope}_Nmax_{param_max_number_neighbors}_a_{param_inverse_dist_exp}_Dmax_{param_max_propagation_distance}_A12_{param_distance_range}_gaps_{param_size_gaps_close}.tif'
        output_dir / 'WD_method_A_Smax_0.2_Nmax_100_a_1_Dmax_10_A12_3_gaps_0.05.tif'
            ) as src,
    rasterio.open(
        input_dir / 'flood.tif') as src2
    ):
        flexth_rst = src.read(1)
        initial_flood = src2.read(1)

        pixel_width, pixel_height = src.res  
        pixel_area = abs(pixel_width * pixel_height) 
        
        mask = (flexth_rst > 10) & (flexth_rst != 9999)
        initial_flood_area = (initial_flood == 1)

        count_after = np.sum(mask)
        count_before = np.sum(initial_flood_area)
        
        # Calculate the total area
        total_area_after = np.round((count_after * pixel_area) / 1000000, decimals = 2)
        total_area_before = np.round((count_before * pixel_area) / 1000000, decimals = 2)
        
        print(f'Initial Flooded Area: {total_area_before}km2\nFlooded Area with FLEXTH: {total_area_after}km2')
        
        # Plotting
        fig, (ax_before, ax_after) = pyplot.subplots(1,2,figsize = (14, 7))
        show(src2, ax=ax_before)
        show(src, ax=ax_after)

# Function to Plot the Map of the flooded area
# Plot Histogramm of water depths
        
########        
## RUN ! 
########     
  
   
    
if __name__ == '__main__':
    
    print('\nFLEXTH\n')
    start = time.time()  
    
    if param_tiling == True: 
       
        
        #tiles all the inputs if param_tile_inputs is True 
        if param_tile_inputs == True:             
            print('Tiling the input rasters...')
            
            tiling(input_dir / 'flood.tif' , input_dir, param_tile_size)
            
            tiling(input_dir / 'dtm.tif'   , input_dir, param_tile_size)
            
            if os.path.isfile(input_dir / 'exclusion.tif'):
                tiling(input_dir / 'exclusion.tif' , input_dir, param_tile_size)
                
            if os.path.isfile(input_dir / 'obswater.tif'):
                tiling(input_dir / 'obswater.tif' , input_dir, param_tile_size)
                
            if os.path.isfile(input_dir / 'permanent_water.tif'):
                tiling(input_dir / 'permanent_water.tif' , input_dir, param_tile_size)
            



    
        input_list_flood = glob.glob( str(input_dir/"**/*flood_tile*.tif")  , recursive=True)
       
        for index, flood_path in enumerate(input_list_flood):
            
            print(f'Processing tile {index+1} out of {len(input_list_flood)}')
            
            dtm_path             =  flood_path.replace("flood_tile", "dtm_tile") 
            exclusion_path       =  flood_path.replace("flood_tile", "exclusion_tile") 
            obswater_path        =  flood_path.replace("flood_tile", "obswater_tile")  
            permanent_water_path =  flood_path.replace("flood_tile", "permanent_water_tile") 
    
    
            flood_processing(flood_path           = flood_path ,
                             dtm_path             = dtm_path  , 
                             exclusion_path       = exclusion_path ,
                             obswater_path        = obswater_path , 
                             permanent_water_path = permanent_water_path )   

       
        
        
    else:
        flood_processing(flood_path           = input_dir / 'flood.tif' ,
                         dtm_path             = input_dir / 'dtm.tif'   , 
                         exclusion_path       = input_dir / 'exclusion.tif' ,
                         obswater_path        = input_dir / 'obswater.tif' , 
                         permanent_water_path = input_dir / 'permanent_water.tif' )        
    
    
        
    
    
        

    end = time.time()
    print(f'Total processing time: {int((end - start)/60)} minutes')

