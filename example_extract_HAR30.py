from climate_data_extraction import met


# --- ERA5 daily
melf = met()
melf.era5_dir = '/Volumes/DATA/ECMWF/ERA5/reanalysis-era5-single-levels/'
melf.output_dir = '/Users/pohle/PycharmProjects/data-download/marlene/output_ERA5_hourly/'
melf.era5_variable = 'tp'
melf.mstation_fname = "/Volumes/EX1/data/glacier/marlene/abramov_poly.shp"
# melf.era5_subset = [70., 73., 37., 41]
melf.read_ERA5_hourly(subset=False)
melf.shapefile_ATTR_FIELD = 'OBJECTID'
melf.read_shape()
melf.extract_values()
melf.save_data()



