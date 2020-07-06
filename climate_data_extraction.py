# -*- coding: utf-8 -*-


import gdal, ogr, osr, gdalnumeric
import numpy as np
# from netCDF4 import Dataset
# from datetime import datetime
from matplotlib import pyplot as plt
import os
import ast
import fnmatch
import pandas as pd
from netCDF4 import num2date
import datetime
from pyproj import Proj, transform
from pyproj import Transformer

# there are several input raster datasets in slightly different format
# define read_function for each of them


class met():
    """
    /***************************************************************************
     Extract climate data at specific locations, given via an ESRI shapefile

     The shapefile needs an attribute of data type "numeric" to extract the data!
     In the example below, the "OBJECTID" attribute field is numeric. If you want to
     use e.g. RGI-code, you need to transform it first into something that is strictly numeric.
     That is because these values are burned into a raster (rasterization) and can thus not be
     of another format


     ***************************************************************************/


    """
    # TODO: Include a conversion from shapefile_ATTR_FIELD as 'str' to 'numeric' to allow for non-numeric att fields

    def __init__(self):
        self.mstation_fname = None
        self.mstation_variable = None
        self.mstation_pixels = None
        self.mstation_data = None

        # the gridded dataset
        self.era5_fname = None
        self.era5_dir = None
        self.era5_variable = None
        self.era5_subset = None  # should be provided as list [lonmin, lonmax, latmin, latmax]
        self.era5_data = None
        self.era5_ds = None
        self.era5_layers = None
        self.era5_res_x = None
        self.era5_res_y = None
        self.era5_lons = None
        self.era5_lats = None
        self.era5_data = None
        self.era5_time = None
        self.lon_grid = None
        self.lat_grid = None
        self.era5_nan_negative = -10e16
        self.era5_nan_positive = +10e16

        # output SpatialReference
        self.era5_geotransform = None
        self.era5_outSpatialRef = osr.SpatialReference()
        self.era5_outSpatialRef.ImportFromEPSG(4326)
        self.coordTrans = None
        # self.era5_latitudes = None
        # self.era5_longitudes = None

        # output file
        self.output_dir = None

        # the glacier mask in the input raster from polygon
        self.era5_mask_poly = None  # from the polygon body
        # self.era5_mask_outline = None  # from the polygon outline
        self.era5_mask_points = None  # from the point shapefile OR from the polygon outline
        self.gdal_dtype = gdal.GDT_Int32
        self.gdal_nodata_val = 0

        # the extracted values
        self.era5_time_series_extracted = None

        # shapefile
        self.shapefile_driver = ogr.GetDriverByName('ESRI Shapefile')
        self.shapefile = None
        self.shapefile_layer = None
        self.shapefile_spatialRef = None
        self.shapefile_authorityCode = None
        self.shapefile_spatialRef_EPSG = None
        self.shapefile_nFeat = None
        self.shapefile_class = None
        # the attribute field must be numeric - create a numeric lookup table for any possible non-numeric value here
        self.shapefile_ATTR_FIELD = None

    # def test_ATTR_FIELD_num(self):
    #     # TODO: Test if self.ATTR_FIELD is numeric.
    #     # If not, create a numeric lookup table.
    #

    def rasterize_poly(self):
        """
        Rasterizes shapefile to have same extent and resolution as input file

        """
        model_dataset = self.era5_ds
        shape_layer = self.shapefile.GetLayer()
        # mem_drv = gdal.GetDriverByName('MEM')
        mem_drv = gdal.GetDriverByName('MEM')
        mem_raster = mem_drv.Create(
            '',
            self.era5_res_x,
            self.era5_res_y,
            1,
            self.gdal_dtype
        )
        mem_raster.SetProjection(model_dataset.GetProjection())
        mem_raster.SetGeoTransform(self.era5_geotransform)
        mem_band = mem_raster.GetRasterBand(1)
        mem_band.Fill(self.gdal_nodata_val)
        mem_band.SetNoDataValue(self.gdal_nodata_val)

        # http://gdal.org/gdal__alg_8h.html#adfe5e5d287d6c184aab03acbfa567cb1
        # http://gis.stackexchange.com/questions/31568/gdal-rasterizelayer-doesnt-burn-all-polygons-to-raster
        err = gdal.RasterizeLayer(
            mem_raster,
            [1],
            shape_layer,
            None,
            None,
            [1],
            options=["ATTRIBUTE=%s"%self.shapefile_ATTR_FIELD]
        )
        assert err == gdal.CE_None
        return mem_raster.ReadAsArray()

    def read_ERA5(self):
        varstr = "_" + self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_file = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp][0]
        print("Reading %s"%m_file)
        self.era5_fname = self.era5_dir + m_file

        # with GDAL
        self.era5_ds = gdal.Open(self.era5_fname)
        self.era5_geotransform = self.era5_ds.GetGeoTransform()
        # lats and lons
        self.era5_layers = self.era5_ds.RasterCount
        self.era5_res_x = self.era5_ds.RasterXSize
        self.era5_res_y = self.era5_ds.RasterYSize

        # pixel center coordinates
        self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
            self.era5_res_x) * self.era5_geotransform[1]
        self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
            self.era5_res_y) * self.era5_geotransform[5]

        # pixel outer coords
        self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
        self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]
        # scale factor and offset
        scale_str = '%s#scale_factor'%str(self.era5_variable)
        offset_str = '%s#add_offset' % str(self.era5_variable)
        self.era5_scale_factor = float(self.era5_ds.GetMetadata_Dict()[scale_str])
        self.era5_add_offset = float(self.era5_ds.GetMetadata_Dict()[offset_str])
        time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
        time_calendar = self.era5_ds.GetMetadata_Dict()['time#calendar']
        time_units = self.era5_ds.GetMetadata_Dict()['time#units']

        era5_time = list(ast.literal_eval(time_string))
        self.era5_time = num2date(era5_time, time_units, time_calendar)
        self.era5_years = [x.year for x in self.era5_time]
        self.era5_months = [x.month for x in self.era5_time]
        self.era5_days = [x.day for x in self.era5_time]

        self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))
        for slice_i in range(self.era5_layers):
            self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(slice_i + 1).ReadAsArray() * \
                                            self.era5_scale_factor + self.era5_add_offset

        self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
        self.lon_grid = np.array(self.llgrid[0])
        self.lat_grid = np.array(self.llgrid[1])

    def read_CHELSAcru(self):
        print('Reading Chelsa dataset...')

        varstr = "_" + self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.tif') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith("CHELSA")]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = len(m_files)

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        for m_file in m_files:
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # lats and lons
                # self.era5_layers = self.era5_ds.RasterCount
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2+1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2+1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            fname_base = self.era5_fname.split('/')[-1]

            self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(1).ReadAsArray(xoff=int(lon_x1),
                                                                                      yoff=int(lat_y2),
                                                                                      win_xsize=self.era5_res_x,
                                                                                      win_ysize=self.era5_res_y)

            self.era5_time.append(datetime.datetime(int(fname_base.split('_')[3]), int(fname_base.split('_')[2]), 1))

            slice_i += 1
        print('Read in Chelsa data.')

    def read_CHELSAts(self):
        print('Reading Chelsa dataset...')

        varstr = "_" + self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.tif') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith("CHELSA")]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = len(m_files)

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        for m_file in m_files:
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # lats and lons
                # self.era5_layers = self.era5_ds.RasterCount
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                # TODO - subset - adjust and overwrite the previous statements
                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2+1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2+1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            fname_base = self.era5_fname.split('/')[-1]

            self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(1).ReadAsArray(xoff=int(lon_x1),
                                                                                      yoff=int(lat_y2),
                                                                                      win_xsize=self.era5_res_x,
                                                                                      win_ysize=self.era5_res_y)

            self.era5_time.append(datetime.datetime(int(fname_base.split('_')[2]), int(fname_base.split('_')[3]), 1))

            slice_i += 1

        print('Read in Chelsa data.')

    def read_HAR30_mon(self):
        print('Reading HAR30 (reprojected) dataset...')
        #
        varstr = "_" + self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith("har_")]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = None  # len(m_files) * 12

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        for m_file in m_files:

            # m_file = "har_d30km_m_2d_prcp_2001_WGS84.nc"
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)
            rasterband_n = self.era5_ds.RasterCount

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # lats and lons
                self.era5_layers = len(m_files) * rasterband_n
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                # lonlat grid
                self.llcgrid = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)

                # TODO - subset - adjust and overwrite the previous statements

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            # TIME not working after the R transform/reprojection
            # instead use the layer number as days and take the year from the filename
            # time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
            f_dates = None
            # YEAR
            f_year = int(m_file.split('_')[5])
            # MONTHS
            if self.era5_ds.RasterCount == 12:
                f_months = [ii + 1 for ii in range(12)]
                f_days = [1 for ii in range(12)]
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')
                f_dates = f_dates[f_dates.day == 1]
            elif self.era5_ds.RasterCount >= 365:
                print('assuming this is a daily dataset ...  results are wrong if this is not the case ...')
                # f_months = [ii + 1 for ii in range(12)]
                # f_year = 2002
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')


            time_string = f_dates
            # era5_time = []
            # era5_time.append([time_string_i for time_string_i in time_string])
            # np.asarray(era5_time)
            self.era5_time.append([time_string_i for time_string_i in time_string])

            # rasterband_cnt = 1
            for rb in range(rasterband_n):
                self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                               yoff=int(lat_y2),
                                                                                               win_xsize=self.era5_res_x,
                                                                                               win_ysize=self.era5_res_y)

                # self.era5_time.append(datetime.datetime(int(fname_base.split('_')[2]), int(fname_base.split('_')[1]), 1))
                slice_i += 1

        # set manually the threshold for no data values
        nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                          self.era5_data > self.era5_nan_positive)
        self.era5_data[nan_exceed_thresh] = np.nan

        # flatten the time list
        self.era5_time = [tst for tstl in self.era5_time for tst in tstl]

        print('Read in HAR 30 (reprojected) data.')

    def read_HAR30_day(self, v2=False):
        print('Reading HAR30 (reprojected) dataset...')
        #
        varstr = "_" + self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith("har_")]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]
        if v2:
            varstr = "_" + self.era5_variable + "_"
            m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
            m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
            m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
            m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
            # make sure not to read in temp-files
            m_files = [filename for filename in m_files if filename.startswith("HAR")]
            # make sure not to read hidden files
            m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = None  # len(m_files) * 12

        # subset to significantly reduce processing time
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        for m_file in m_files:

            # m_file = "har_d30km_m_2d_prcp_2001_WGS84.nc"
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)
            rasterband_n = self.era5_ds.RasterCount

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # lats and lons
                f_year_st = int(m_files[0].split('_')[5])
                f_year_en = int(m_files[-1].split('_')[5])
                dates_length = pd.date_range(start=str(f_year_st) + '-01-01', end=str(f_year_en) + '-12-31')
                self.era5_layers = len(dates_length)
                # self.era5_layers = len(m_files) * rasterband_n

                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                # lonlat grid
                self.llcgrid = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)

                # subset - adjust and overwrite the previous statements

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            # TIME - use the layer number as days and take the year from the filename
            f_dates = None
            # YEAR
            f_year = int(m_file.split('_')[5])
            # MONTHS
            if self.era5_ds.RasterCount == 12:
                f_months = [ii + 1 for ii in range(12)]
                f_days = [1 for ii in range(12)]
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')
                f_dates = f_dates[f_dates.day == 1]
            elif self.era5_ds.RasterCount >= 365:
                print('assuming this is a daily dataset ...  results are wrong if this is not the case ...')
                # f_months = [ii + 1 for ii in range(12)]
                # f_year = 2002
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')

            time_string = f_dates
            self.era5_time.append([time_string_i for time_string_i in time_string])

            for rb in range(rasterband_n):
                self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                               yoff=int(lat_y2),
                                                                                               win_xsize=self.era5_res_x,
                                                                                               win_ysize=self.era5_res_y)

                slice_i += 1

        # set manually the threshold for no data values
        nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                          self.era5_data > self.era5_nan_positive)
        self.era5_data[nan_exceed_thresh] = np.nan

        # flatten the time list
        self.era5_time = [tst for tstl in self.era5_time for tst in tstl]

        print('Read in HAR 30 (reprojected) data.')

    def read_APHR(self, v='1101'):
        print('Reading APHRODITE data ...')
        # v='1101'

        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_onlyfiles.sort()
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith("APHRO_")]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = None  # len(m_files) * 12

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        for m_file in m_files:

            # m_file = "har_d30km_m_2d_prcp_2001_WGS84.nc"
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL -- NOT WORKING - BANDS EMPTY (GetRasterBand())
            self.era5_ds = gdal.Open('NETCDF:"' + self.era5_fname + '":%s' % self.era5_variable)
            if not self.era5_ds:
                print('Something went wrong with opening the file. Check the variable name "e.g. with gdalinfo')
                return
            # self.era5_ds = gdal.Open(self.era5_fname)
            rasterband_n = self.era5_ds.RasterCount

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())

                # lats and lons
                if v == '1101':
                    f_year_st1 = m_files[0].split('_')[3]
                    f_year_st = int(f_year_st1.split('.')[1])

                    f_year_en1 = m_files[-1].split('_')[3]
                    f_year_en = int(f_year_en1.split('.')[1])
                elif v == '1101_EXR1':
                    f_year_st1 = m_files[0].split('_')[4]
                    f_year_st = int(f_year_st1.split('.')[1])

                    f_year_en1 = m_files[-1].split('_')[4]
                    f_year_en = int(f_year_en1.split('.')[1])
                elif v == '1808':
                    f_year_st1 = m_files[0].split('_')[4]
                    f_year_st = int(f_year_st1.split('.')[1])

                    f_year_en1 = m_files[-1].split('_')[4]
                    f_year_en = int(f_year_en1.split('.')[1])
                elif v == '1808R1':
                    f_year_st1 = m_files[0].split('_')[3]
                    f_year_st = int(f_year_st1.split('.')[1])

                    f_year_en1 = m_files[-1].split('_')[3]
                    f_year_en = int(f_year_en1.split('.')[1])
                elif v == '1901':
                    f_year_st1 = m_files[0].split('_')[3]
                    f_year_st = int(f_year_st1.split('.')[1])

                    f_year_en1 = m_files[-1].split('_')[3]
                    f_year_en = int(f_year_en1.split('.')[1])
                else:
                    print('Unspecified version APHRODITE dataset identifier (>>1101<< e.g.)')
                    return

                dates_length = pd.date_range(start=str(f_year_st) + '-01-01', end=str(f_year_en) + '-12-31')
                self.era5_layers = len(dates_length)

                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                # lonlat grid
                self.llcgrid = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)

                # - subset - adjust and overwrite the previous statements
                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            # TIME - use the layer number as days and take the year from the filename
            f_dates = None
            # YEAR
            if v == '1101':
                f_year_st1 = m_file.split('_')[3]
                f_year_st = int(f_year_st1.split('.')[1])
                f_year = int(f_year_st)
            elif v == '1101_EXR1':
                f_year_st1 = m_file.split('_')[4]
                f_year_st = int(f_year_st1.split('.')[1])
                f_year = int(f_year_st)
            elif v == '1808':
                f_year_st1 = m_file.split('_')[4]
                f_year_st = int(f_year_st1.split('.')[1])
                f_year = int(f_year_st)
            elif v == '1808R1':
                f_year_st1 = m_file.split('_')[3]
                f_year_st = int(f_year_st1.split('.')[1])
                f_year = int(f_year_st)
            elif v == '1901':
                f_year_st1 = m_file.split('_')[3]
                f_year_st = int(f_year_st1.split('.')[1])
                f_year = int(f_year_st)
            else:
                print('Unspecified version APHRODITE dataset identifier (>>1101<< e.g.)')
                return

            # MONTHS
            if self.era5_ds.RasterCount == 12:
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')
                f_dates = f_dates[f_dates.day == 1]
            elif self.era5_ds.RasterCount == 365 or self.era5_ds.RasterCount == 366:
                print('assuming this is a daily dataset ...  results are wrong if this is not the case ...')
                f_dates = pd.date_range(start=str(f_year) + '-01-01', end=str(f_year) + '-12-31')

            time_string = f_dates
            self.era5_time.append([time_string_i for time_string_i in time_string])

            for rb in range(rasterband_n):
                self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                               yoff=int(lat_y2),
                                                                                               win_xsize=self.era5_res_x,
                                                                                               win_ysize=self.era5_res_y)

                slice_i += 1
        # set manually the threshold for no data values
        nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative, self.era5_data > self.era5_nan_positive)
        self.era5_data[nan_exceed_thresh] = np.nan

        # flatten the time list
        self.era5_time = [tst for tstl in self.era5_time for tst in tstl]

        print('Read in APHRODITE data.')

    # def read_ERA5_hourly(self, subset = False):
    #     varstr = self.era5_variable + "_"
    #     m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
    #     m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
    #     m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
    #     m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
    #     m_files.sort()
    #
    #     m_files = [filename for filename in m_files if filename.startswith(self.era5_variable)]
    #     # make sure not to read hidden files
    #     m_files = [filename for filename in m_files if not filename.startswith(".")]
    #
    #     self.era5_layers = None  # len(m_files) * 12
    #
    #     # subset to significantly reduce processing time
    #     # self.era5_subset = [70., 73., 37., 41]
    #     if subset:
    #         if not self.era5_subset:
    #             print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
    #             return
    #         else:
    #             lon_min = self.era5_subset[0]
    #             lon_max = self.era5_subset[1]
    #             lat_min = self.era5_subset[2]
    #             lat_max = self.era5_subset[3]
    #
    #     slice_i = 0
    #     self.era5_geotransform = None
    #     self.era5_time = []
    #     for m_file in m_files:
    #
    #         # m_file = "har_d30km_m_2d_prcp_2001_WGS84.nc"
    #         print("Reading %s" % m_file)
    #         self.era5_fname = self.era5_dir + m_file
    #
    #         # with GDAL
    #         self.era5_ds = gdal.Open(self.era5_fname)
    #         rasterband_n = self.era5_ds.RasterCount
    #
    #         # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
    #         if not self.era5_geotransform:
    #             self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
    #
    #             # get times - some files might not start at specific time -  so read them all now in a loop to extract the time
    #             times_all = []
    #             for fname_tmp in m_files:
    #                 print('Extracting time code from %s' %fname_tmp)
    #                 tmp_ds = gdal.Open(self.era5_dir + fname_tmp)
    #                 time_string = tmp_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
    #                 time_calendar = tmp_ds.GetMetadata_Dict()['time#calendar']
    #                 time_units = tmp_ds.GetMetadata_Dict()['time#units']
    #
    #                 tmp_time = list(ast.literal_eval(time_string))
    #                 tmp_time_date = num2date(tmp_time, time_units, time_calendar)
    #                 times_all.append(tmp_time_date)
    #
    #             self.era5_time = [tst for tstl in times_all for tst in tstl]
    #             self.era5_layers = len(self.era5_time)
    #
    #             scale_str = '%s#scale_factor' % str(self.era5_variable)
    #             offset_str = '%s#add_offset' % str(self.era5_variable)
    #             self.era5_scale_factor = float(self.era5_ds.GetMetadata_Dict()[scale_str])
    #             self.era5_add_offset = float(self.era5_ds.GetMetadata_Dict()[offset_str])
    #
    #             # lats and lons
    #             self.era5_res_x = self.era5_ds.RasterXSize
    #             self.era5_res_y = self.era5_ds.RasterYSize
    #
    #             # pixel center coordinates
    #             self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
    #                 self.era5_res_x) * self.era5_geotransform[1]
    #             self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
    #                 self.era5_res_y) * self.era5_geotransform[5]
    #
    #             # lonlat grid
    #             self.llcgrid = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)
    #
    #             # subset - adjust and overwrite the previous statements
    #             if subset:
    #                 lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
    #                 lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
    #                 lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
    #                 lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))
    #
    #                 # overwrite lons_center and lats_center
    #                 self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
    #                 self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])
    #
    #                 # overwrite resx/resy
    #                 self.era5_res_x = len(self.era5_lons_center)
    #                 self.era5_res_y = len(self.era5_lats_center)
    #
    #                 # overwrite era5_geotransform
    #                 self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
    #                 self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])
    #
    #             # pixel outer coords
    #             self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
    #             self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]
    #
    #             self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
    #             self.lon_grid = np.array(self.llgrid[0])
    #             self.lat_grid = np.array(self.llgrid[1])
    #
    #             self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))
    #
    #         for rb in range(tmp_time_date):
    #             if subset:
    #                 self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray(xoff=int(lon_x1),
    #                                                                                                yoff=int(lat_y2),
    #                                                                                                win_xsize=self.era5_res_x,
    #                                                                                                win_ysize=self.era5_res_y) * \
    #                                         self.era5_scale_factor + self.era5_add_offset
    #             else:
    #                 self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray() * \
    #                                         self.era5_scale_factor + self.era5_add_offset
    #             slice_i += 1
    #
    #     # set manually the threshold for no data values
    #     nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
    #                                       self.era5_data > self.era5_nan_positive)
    #     self.era5_data[nan_exceed_thresh] = np.nan
    #     print('Read in ERA5 data.')

    def read_ERA5_hourly(self, subset=False):
        varstr = self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        m_files.sort()

        m_files = [filename for filename in m_files if filename.startswith(self.era5_variable)]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = None  # len(m_files) * 12

        # subset to significantly reduce processing time
        # self.era5_subset = [70., 73., 37., 41]
        # if subset:
        #     if not self.era5_subset:
        #         print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
        #         return
        #     else:
        #         lon_min = self.era5_subset[0]
        #         lon_max = self.era5_subset[1]
        #         lat_min = self.era5_subset[2]
        #         lat_max = self.era5_subset[3]

        rs_stack_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        raster_stack_list = []
        times_all = []
        n_rasterband_total = 0

        for m_file in m_files:

            # m_file = "har_d30km_m_2d_prcp_2001_WGS84.nc"
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)
            rasterband_n = self.era5_ds.RasterCount
            n_rasterband_total = n_rasterband_total + rasterband_n
            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())

                # get times - some files might not start at specific time -  so read them all now in a loop to extract the time
                time_calendar = self.era5_ds.GetMetadata_Dict()['time#calendar']
                time_units = self.era5_ds.GetMetadata_Dict()['time#units']

                scale_str = '%s#scale_factor' % str(self.era5_variable)
                offset_str = '%s#add_offset' % str(self.era5_variable)
                self.era5_scale_factor = float(self.era5_ds.GetMetadata_Dict()[scale_str])
                self.era5_add_offset = float(self.era5_ds.GetMetadata_Dict()[offset_str])

                # lats and lons
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                # lonlat grid
                self.llcgrid = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)

                # subset - adjust and overwrite the previous statements
                if subset:
                    lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                    lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                    lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                    lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                    # overwrite lons_center and lats_center
                    self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                    self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                    # overwrite resx/resy
                    self.era5_res_x = len(self.era5_lons_center)
                    self.era5_res_y = len(self.era5_lats_center)

                    # overwrite era5_geotransform
                    self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                    self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                # make individual array for each file and combine them only at the end into one final array
                # use a dictionary with the year name as keys - also store the time in this one and combine in the end


            time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
            tmp_time = list(ast.literal_eval(time_string))
            tmp_time_date = num2date(tmp_time, time_units, time_calendar)
            times_all.append(tmp_time_date)
            era5_data = np.zeros((self.era5_res_y, self.era5_res_x, len(tmp_time_date)))


            slice_i = 0
            for rb in range(len(tmp_time_date)):
                if subset:
                    era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                              yoff=int(lat_y2),
                                                                                              win_xsize=self.era5_res_x,
                                                                                              win_ysize=self.era5_res_y) * \
                                               self.era5_scale_factor + self.era5_add_offset
                else:
                    era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(rb + 1).ReadAsArray() * \
                                               self.era5_scale_factor + self.era5_add_offset
                slice_i += 1

            # add the layer stack of an individual file to the raster stack dictionary
            raster_stack_list.append(era5_data)

        self.era5_data = np.dstack(raster_stack_list)
        self.era5_time = [tst for tstl in times_all for tst in tstl]

        # set manually the threshold for no data values
        nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                          self.era5_data > self.era5_nan_positive)
        self.era5_data[nan_exceed_thresh] = np.nan
        print('Read in ERA5 data.')

    def read_CRUTS(self):
        print('Reading CRU-TS dataset...')

        varstr = self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        m_files = [filename for filename in m_files if filename.startswith(varstr)]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        self.era5_layers = len(m_files) * 12  # monthly dataset with 12 rasters per file

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        times_all = []
        for m_file in m_files:
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open(self.era5_fname)

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # lats and lons
                # self.era5_layers = self.era5_ds.RasterCount
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # time
                time_calendar = self.era5_ds.GetMetadata_Dict()['time#calendar']
                time_units = self.era5_ds.GetMetadata_Dict()['time#units']

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            fname_base = self.era5_fname.split('/')[-1]

            for lyr in range(12):
                self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(lyr + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                                yoff=int(lat_y2),
                                                                                                win_xsize=self.era5_res_x,
                                                                                                win_ysize=self.era5_res_y)
                slice_i += 1

            time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
            tmp_time = list(ast.literal_eval(time_string))
            tmp_time_date = num2date(tmp_time, time_units, time_calendar)
            # replace day with 1
            tmp_time_date = [tmp_time_date_i.replace(day=1) for tmp_time_date_i in tmp_time_date]
            times_all.append(tmp_time_date)

            self.era5_time = [tst for tstl in times_all for tst in tstl]

            nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                              self.era5_data > self.era5_nan_positive)
            self.era5_data[nan_exceed_thresh] = np.nan

        print('Read in CRU-TS data.')

    def read_GPCC(self):
        print('Reading GPCC dataset...')

        # varstr = self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.nc') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        m_files = [filename for filename in m_files if not filename.startswith(".")]

        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        self.era5_time = []
        times_all = []
        for m_file in m_files:
            print("Reading %s" % m_file)
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open('NETCDF:{0}:{1}'.format(self.era5_fname, self.era5_variable))

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())

                # Layers
                self.era5_layers = self.era5_ds.RasterCount

                # adjust if negative lons (-180 to 180)
                if self.era5_geotransform[0] < 0:
                    self.era5_geotransform[0] = self.era5_geotransform[0] + self.era5_geotransform[0] * (-1)
                # lats and lons
                # self.era5_layers = self.era5_ds.RasterCount
                self.era5_res_x = self.era5_ds.RasterXSize
                self.era5_res_y = self.era5_ds.RasterYSize

                # time
                time_calendar = self.era5_ds.GetMetadata_Dict()['time#calendar']
                time_units = self.era5_ds.GetMetadata_Dict()['time#units']

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y) * self.era5_geotransform[5]

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            # fname_base = self.era5_fname.split('/')[-1]

            for lyr in range(self.era5_layers):
                self.era5_data[:, :, slice_i] = self.era5_ds.GetRasterBand(lyr + 1).ReadAsArray(xoff=int(lon_x1),
                                                                                                yoff=int(lat_y2),
                                                                                                win_xsize=self.era5_res_x,
                                                                                                win_ysize=self.era5_res_y)
                slice_i += 1

            time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
            tmp_time = list(ast.literal_eval(time_string))
            tmp_time_date = num2date(tmp_time, time_units, time_calendar)
            times_all.append(tmp_time_date)

            self.era5_time = [tst for tstl in times_all for tst in tstl]

            nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                              self.era5_data > self.era5_nan_positive)
            self.era5_data[nan_exceed_thresh] = np.nan

        print('Read in GPCC data.')

    def read_GPM(self):
        print('Reading GPCC dataset...')

        # files
        # varstr = self.era5_variable + "_"
        m_onlyfiles = [f for f in os.listdir(self.era5_dir) if os.path.isfile(os.path.join(self.era5_dir, f))]
        m_split = [fnmatch.fnmatch(tmp_file, '*.HDF5') for tmp_file in m_onlyfiles]
        m_files = [m_onlyfiles[ii] for ii in np.arange(len(m_onlyfiles)) if m_split[ii]]
        # m_files = [m_files_tmp for m_files_tmp in m_files if varstr in m_files_tmp]
        # make sure not to read in temp-files
        # m_files = [filename for filename in m_files if filename.startswith(varstr)]
        # make sure not to read hidden files
        m_files = [filename for filename in m_files if not filename.startswith(".")]
        m_files.sort()

        self.era5_layers = len(m_files)

        # get date from filenames

        m_str = [m_files[ii].split('.')[4].split("-")[0] for ii in np.arange(len(m_files))]
        self.era5_time = [datetime.datetime.strptime(all_dates_i, '%Y%m%d') for all_dates_i in m_str]

        # subset to significantly reduce processing time
        # self.era5_subset = [71.47, 71.63, 39.5, 39.7]
        if not self.era5_subset:
            print('Please define a subset "metObject.era5_subset = [lonmin, lonmax, latmin, latmax]"')
            return
        else:
            lon_min = self.era5_subset[0]
            lon_max = self.era5_subset[1]
            lat_min = self.era5_subset[2]
            lat_max = self.era5_subset[3]

        slice_i = 0
        self.era5_geotransform = None
        # self.era5_time = []
        # times_all = []
        for m_file in m_files:
            print("Reading %s" % m_file)
            # m_file = m_files[0]
            self.era5_fname = self.era5_dir + m_file

            # with GDAL
            self.era5_ds = gdal.Open('HDF5:"{0}"://Grid/{1}'.format(self.era5_fname, self.era5_variable))

            # if the metaparameters and geospatial data are not assigned yet read it, otherwise skip this redundant step
            if not self.era5_geotransform:
                self.era5_geotransform = list(self.era5_ds.GetGeoTransform())
                # weird format: write info manually
                # also the file is transposed - for the subset lat/lon needs to be switched
                self.era5_geotransform = [-180, 0.1, 0.0, 90.0, 0.0, -0.1]
                # switch geotransform to 0-360 system at the end

                # SWITCHED X AND Y HERE
                self.era5_res_x0 = self.era5_ds.RasterYSize
                self.era5_res_y0 = self.era5_ds.RasterXSize

                # pixel center coordinates
                self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
                    self.era5_res_x0) * self.era5_geotransform[1]
                self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
                    self.era5_res_y0) * self.era5_geotransform[5]

                lon_x1 = np.argmin(np.abs(self.era5_lons_center - lon_min))
                lon_x2 = np.argmin(np.abs(self.era5_lons_center - lon_max))
                lat_y1 = np.argmin(np.abs(self.era5_lats_center - lat_max))
                lat_y2 = np.argmin(np.abs(self.era5_lats_center - lat_min))

                # overwrite lons_center and lats_center
                self.era5_lons_center = self.era5_lons_center[lon_x1:lon_x2 + 1] + (0.5 * self.era5_geotransform[1])
                self.era5_lats_center = self.era5_lats_center[lat_y1:lat_y2 + 1] + (0.5 * self.era5_geotransform[5])

                # overwrite resx/resy
                self.era5_res_x = len(self.era5_lons_center)
                self.era5_res_y = len(self.era5_lats_center)

                # overwrite era5_geotransform
                self.era5_geotransform[0] = self.era5_lons_center[0] - (0.5 * self.era5_geotransform[1])
                self.era5_geotransform[3] = self.era5_lats_center[0] - (0.5 * self.era5_geotransform[5])

                # pixel outer coords
                self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
                self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

                self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
                self.lon_grid = np.array(self.llgrid[0])
                self.lat_grid = np.array(self.llgrid[1])

                self.era5_data = np.zeros((self.era5_res_y, self.era5_res_x, self.era5_layers))

            self.era5_data[:, :, slice_i] = np.flipud(
                self.era5_ds.GetRasterBand(1).ReadAsArray(xoff=int(self.era5_res_y0 - lat_y2),
                                                          yoff=int(lon_x1),
                                                          win_xsize=self.era5_res_y,
                                                          win_ysize=self.era5_res_x).T)

            slice_i += 1

            # time_string = self.era5_ds.GetMetadata_Dict()['NETCDF_DIM_time_VALUES'][1:-1]
            # tmp_time = list(ast.literal_eval(time_string))
            # tmp_time_date = num2date(tmp_time, time_units, time_calendar)
            # replace day with 1
            # tmp_time_date = [tmp_time_date_i.replace(day=1) for tmp_time_date_i in tmp_time_date]
            # times_all.append(tmp_time_date)

            # self.era5_time = [tst for tstl in times_all for tst in tstl]

        nan_exceed_thresh = np.logical_or(self.era5_data < self.era5_nan_negative,
                                          self.era5_data > self.era5_nan_positive)
        self.era5_data[nan_exceed_thresh] = np.nan

        # fix the coordinates
        # self.era5_geotransform
        if self.era5_geotransform[0] < 0:
            self.era5_geotransform[0] = self.era5_geotransform[0] + 360

        self.era5_lons_center = self.era5_geotransform[0] + (0.5 * self.era5_geotransform[1]) + np.arange(
            self.era5_res_x) * self.era5_geotransform[1]
        self.era5_lats_center = self.era5_geotransform[3] + (0.5 * self.era5_geotransform[5]) + np.arange(
            self.era5_res_y) * self.era5_geotransform[5]

        # pixel outer coords
        self.era5_lons = self.era5_geotransform[0] + np.arange(self.era5_res_x + 1) * self.era5_geotransform[1]
        self.era5_lats = self.era5_geotransform[3] + np.arange(self.era5_res_y + 1) * self.era5_geotransform[5]

        self.llgrid = np.meshgrid(self.era5_lons, self.era5_lats, sparse=False)
        self.lon_grid = np.array(self.llgrid[0])
        self.lat_grid = np.array(self.llgrid[1])

        print('Read in GPM data.')

    def read_shape(self):
        ''''
        # differentiate between points and polygons
        # for points: take 4 x,y coordinates from raster that are closest
        # for polygons: take all raster cells that are overlayn; if only one raster cell --> take 4 closest
        '''

        if not os.path.isfile(self.mstation_fname):
            print('Shapefile not found. Stopping.')
            return

        self.shapefile = self.shapefile_driver.Open(self.mstation_fname)
        self.shapefile_layer = self.shapefile.GetLayer()
        self.shapefile_spatialRef = self.shapefile_layer.GetSpatialRef()
        self.shapefile_authorityCode = self.shapefile_spatialRef.GetAuthorityCode(None)
        print('EPSG code: %s'%int(self.shapefile_authorityCode))
        self.shapefile_spatialRef_EPSG = self.shapefile_spatialRef.ImportFromEPSG(int(self.shapefile_authorityCode))
        self.shapefile_nFeat = self.shapefile_layer.GetFeatureCount()

        # create the CoordinateTransformation
        self.coordTrans = osr.CoordinateTransformation(self.shapefile_spatialRef, self.era5_outSpatialRef)

        # test if polygon or point
        aa_feat = self.shapefile_layer.GetFeature(0)
        aa_class = aa_feat.geometry().GetGeometryName()
        self.shapefile_class = aa_feat.geometry().GetGeometryName()
        self.shapefile_class_type = None
        if aa_class == 'POINT' or aa_class == 'Point':
            self.shapefile_class_type = 'point'
            print(self.shapefile_class_type)

        elif aa_class == 'POLYGON' or aa_class == 'Polygon':
            self.shapefile_class_type = 'polygon'
            print(self.shapefile_class_type)
        else:
            print('Could not identify if point or polygon. Stopping!')


        # get points if POINT class
        if self.shapefile_class_type == 'point':
            # loop through the input features
            inFeature = self.shapefile_layer.GetNextFeature()
            self.latlon = []
            glacier_locations_xy = {}

            while inFeature:
                # get the input geometry
                geom = inFeature.GetGeometryRef()

                # get RGI OBJECTID (numeric)
                objectid = inFeature[self.shapefile_ATTR_FIELD]

                # # re-project the geometry
                # geom.Transform(self.coordTrans)
                p = geom.GetPoints()
                self.latlon.append(p)

                # get xy coordinates and save them into glacier_locations_xy together with OBJECTID
                x, y = zip(*p)
                xc = [(n - self.era5_geotransform[0]) / self.era5_geotransform[1] for n in x]
                yc = [(n - self.era5_geotransform[3]) / self.era5_geotransform[5] for n in y]
                yc_r = [np.floor(yc_i).astype('int') for yc_i in yc]
                xc_r = [np.floor(xc_i).astype('int') for xc_i in xc]

                # get unique pixels
                xyc_r = list(set(zip(yc_r, xc_r)))
                glacier_cells = xyc_r
                glacier_locations_xy[objectid] = {'glacier_cells': glacier_cells}
                inFeature = self.shapefile_layer.GetNextFeature()
            self.era5_mask_points = glacier_locations_xy

        # get points if POLYGON class
        elif self.shapefile_class_type == 'polygon':

            # This function will convert the rasterized clipper shapefile
            # to a mask for use within GDAL.
            self.era5_mask_poly = self.rasterize_poly()
            self.latlon = []
            glacier_locations_xy = {}
            n = self.shapefile_nFeat

            for ni in range(n):
                inFeature = self.shapefile_layer.GetFeature(ni)

                # get RGI OBJECTID (numeric)
                objectid = inFeature[self.shapefile_ATTR_FIELD]
                geom = inFeature.GetGeometryRef()
                points = geom.GetGeometryRef(0).GetPoints()

                # use points to check which raster cells have parts of the polygon within them
                self.latlon.append(points)
                x, y = zip(*points)

                xc = [(n - self.era5_geotransform[0]) / self.era5_geotransform[1] for n in x]
                yc = [(n - self.era5_geotransform[3]) / self.era5_geotransform[5] for n in y]
                yc_r = [np.floor(yc_i).astype('int') for yc_i in yc]
                xc_r = [np.floor(xc_i).astype('int') for xc_i in xc]

                xyc_r = list(set(zip(yc_r, xc_r)))
                glacier_cells = xyc_r
                glacier_locations_xy[objectid] = {'glacier_cells': glacier_cells}

            self.era5_mask_points = glacier_locations_xy
            # ^ glaciers cells must be added to the rasterized polygon file!
            # Combine mask_poly and mask_outline
            _add = []
            for k in self.era5_mask_points.keys():
                # ind_glac_id = np.where(self.era5_mask_poly == k)
                ind_glac_id = np.asarray(self.era5_mask_poly == k).nonzero()

                a = ind_glac_id[0]
                b = ind_glac_id[1]
                _add = np.asarray([a, b]).T
                [self.era5_mask_points[k]['glacier_cells'].append(_add_i) for _add_i in _add]

    def extract_values(self):
        '''
        Use the mask_points (polygon or point shapefile) coordinates to
        extract the values from the raster file

        # 1) go trough each key of self.era5_mask_points and extract the corresponding x,y matrix coordinates
        # 2) extract the time-series from the raster file
        '''
        if not self.era5_mask_points:
            print('Read in shapefile for extraction first')
            return

        self.era5_time_series_extracted = {}
        # EXTRACT
        print('Extracting time series for glacier IDs ...')
        for key in self.era5_mask_points.keys():
            if any(np.array(self.era5_mask_points[key]['glacier_cells']).flatten() < 0):
                print('It seems like the raster is smaller than the shapefile. Stopping process ...')
                return
            xy = self.era5_mask_points[key]['glacier_cells']
            xy = np.asarray(xy)

            ts_tmp = self.era5_data[xy[:, 0], xy[:, 1], :]
            lons = self.era5_lons_center[xy[:, 1]]
            lats = self.era5_lats_center[xy[:, 0]]
            # pyplot.scatter(lons, lats, c=ts_tmp[:,0])
            self.era5_time_series_extracted[key] = {"lons": lons, "lats": lats, "ts": ts_tmp}
        print('Done extracting time series.')
        
    def save_data(self):
        '''
        # 1) save into dictionary with the same keys as the self.era5_time_series_extracted
        '''
        if self.era5_time_series_extracted:
            for key in self.era5_time_series_extracted.keys():
                list1, list2 = self.era5_time_series_extracted[key]['lons'], self.era5_time_series_extracted[key][
                    'lats']
                head = np.array([str(a) + ';' +  str(b) for a, b in zip(list1, list2)])
                data = self.era5_time_series_extracted[key]['ts'].T
                df_out = pd.DataFrame(data, columns=head, index=self.era5_time)
                outfile = self.output_dir + self.era5_variable + '_' + str(key) + '.txt'
                df_out.to_csv(outfile)
                print('Output file %s written.'%outfile)
        else:
            print('Run the "extract_values()" command first.')

    def save_grid(self):
        '''
        save the raster stack instead of the extracted values in simple txt-file format
        :return: csv file output
        '''
        # use reshape to get each pixel ts as a single column
        reshape_ts = (self.era5_data.shape[0] * self.era5_data.shape[1], -1)
        reshape_ll = (self.era5_data.shape[0] * self.era5_data.shape[1])

        llgrid_center = np.meshgrid(self.era5_lons_center, self.era5_lats_center, sparse=False)
        lon_grid = np.array(llgrid_center[0])
        lat_grid = np.array(llgrid_center[1])

        data = self.era5_data.reshape(reshape_ts).T
        lons_all = lon_grid.reshape(reshape_ll)
        lats_all = lat_grid.reshape(reshape_ll)

        head = np.array([str(a) + ';' + str(b) for a, b in zip(lons_all, lats_all)])
        df_out = pd.DataFrame(data, columns=head, index=self.era5_time)
        outfile = self.output_dir + 'subset_' + self.era5_variable + '.txt'
        df_out.to_csv(outfile)
        print('Output file %s written.' % outfile)

# --- GPM 3IMERGM.06 --- #
melf = met()
melf.era5_dir = '/Volumes/DATA/GPM/GPM_3IMERGM.06/'
melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_GPM-3IMERG_mon/'
melf.era5_variable = 'precipitation'
melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
melf.era5_subset = [69, 73, 38, 41]

melf.read_GPM()
melf.shapefile_ATTR_FIELD = 'OBJECTID'
melf.read_shape()
melf.extract_values()
melf.save_data()
melf.save_grid()

# # --- GPCC --- #
# melf = met()
# melf.era5_dir = '/Volumes/DATA/GPCC/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_GPCC/'
# melf.era5_variable = 'precip'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [69, 73, 38, 41]
# melf.read_GPCC()
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- CRU TS 4.03 --- #
# melf = met()
# melf.era5_dir = '/Volumes/DATA/CRU/cru-ts-4.03/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_CRU-TS/'
# melf.era5_variable = 'pre'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [69, 73, 38, 41]
# melf.read_CRUTS()
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()

# # --- ERA5 daily
# melf = met()
# melf.era5_dir = '/Users/pohle/TMP/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/'
# melf.era5_variable = 'tp'
# # melf.mstation_fname = '/Volumes/EX1/data/glacier/marlene/test_2_glaciers_WGS84_TEST.shp'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# # melf.mstation_fname = '/Volumes/EX1/data/glacier/marlene/test_2_glaciers_p2.shp'
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_ERA5_hourly(subset=False)
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# --- HAR30 monthly timeseries --- #
# melf = met()
# melf.era5_dir = '/Volumes/EX1/data/HAR/d30km/m/2d_WGS84/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_mon/'
# melf.era5_variable = 'prcp'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# # read the dataset
# melf.read_HAR30_mon()
# # define the attribute field of the shapefile that is used to classify the different glaciers (must be numeric)
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# # read in the shapefile
# melf.read_shape()
# # extract the values where the shapefile covers the raster file
# melf.extract_values()
# # write the results into a text file
# melf.save_data()


# # --- HAR30 daily --- #
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/HAR/v14/d30km/d/2d_WGS84/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR30_day/'
# melf.era5_variable = 'prcp'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_HAR30_day()
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- HAR10 V2 daily --- #
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/HAR/v200/d10km/d/2d_WGS84/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_HAR10_V2_day/'
# melf.era5_variable = 'prcp'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_HAR30_day(v2=True)
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- CHELSA TS monthly --- #
# melf = met()
# melf.era5_dir = '/Volumes/EX1/data/CHELSA/timeseries/prec/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_chelsats_mon/'
# melf.era5_variable = 'prec'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [71.49, 71.61, 39.58, 39.67]
# melf.read_CHELSAts()
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- CHELSA CRUTS monthly --- #
# melf = met()
# melf.era5_dir = '/Volumes/EX1/data/CHELSA/timeseries20c/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_chelsacruts_mon/'
# melf.era5_variable = 'prec'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [71.49, 71.61, 39.58, 39.67]
# melf.read_CHELSAcru()
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- APHRODITE --- #
# # V1101 - precip
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/APHRODITE/V1101/'
# fname = 'APHRO_MA_025deg_V1101.1951.nc'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1101/'
# melf.era5_variable = 'precip'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_APHR(v='1101')
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- APHRODITE --- #
# # V1101_EXR1 - precip
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/APHRODITE/V1101_EXR1/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1101_EXR1/'
# melf.era5_variable = 'precip'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_APHR(v='1101_EXR1')
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()


# # --- APHRODITE --- #
# # V1808 - TAVE
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/APHRODITE/V1808/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1808/'
# melf.era5_variable = 'tave'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_APHR(v='1808')
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()

# # --- APHRODITE --- #
# # V1808R1 - precip
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/APHRODITE/V1801R1/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1808R1/'
# melf.era5_variable = 'precip'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_APHR(v='1808R1')
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()

# # --- APHRODITE --- #
# # V1901 - precip
# melf = met()
# melf.era5_dir = '/Volumes/DATA_FAT/APHRODITE/V1901/'
# melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_APHR_V1901/'
# melf.era5_variable = 'precip'
# melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
# melf.read_APHR(v='1901')
# melf.shapefile_ATTR_FIELD = 'OBJECTID'
# melf.read_shape()
# melf.extract_values()
# melf.save_data()

