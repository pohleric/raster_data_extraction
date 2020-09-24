# raster_data_extraction
Collection of wrapper functions to extract data from various gridded data products, including:

- High Asia Refined (HAR)
    - (daily + monthly)
        - V1.4 
        - V2.0

- APHRODITE
  - V1101
  - V1101-EXR1
  - V1808
  - V1808-R1
  - V1901

- GPCC

- CHEALSA
  - TS
  - CRU

- ERA5
  - hourly
  - monthly
  
- GRACE
    - RL06M JPL
    
- MODIS
    - MOD10CM
  
Data can be extracted using either ESRI shapefiles; from point features or polygon features (single or multi), or using a bounding box.
Extracted data will be stored as time series in csv-format. Each column represents a pixel's time series. The header will have the name of the coordinates.
