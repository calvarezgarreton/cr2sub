## CR2SUB: consolidated monthly groundwater level database for Chile

The CR2SUB database includes groundwater level (GWL) observations from 1,137 wells maintained by the Water Bureau of Chile (Dirección General de Aguas, DGA), downloaded from https://snia.mop.gob.cl/BNAConsultas/reportes, and processed into a homogenized monthly database for the full domain of Chile, for the period **1957-2025**. The database also includes metadata for each observations well and a range of topographical and hydroclimatic attributes computed from ancillary information.

The database has been consolidated within the framework of the Center for Climate and Resilience Research (CR2, ANID/FONDAP/1523A0002) and is part of the research project ANID/FONDECYT/11240924.

### cr2sub

**cr2sub_v1_mon.csv**<br>
 Consolidated monthly time series of raw groundwater level (GWL) observations, downloaded from the [DGA website](https://snia.mop.gob.cl/BNAConsultas/reportes) and stored in the _input_ folder. When more than one record was available for a given month, the values were averaged to obtain a single monthly estimate.

**cr2sub_v1_mon_clean.csv**<br>
 Monthly time series after applying an outlier-removal procedure. The code used for outlier detection and removal is available in the _scripts_ folder.

**cr2sub_v1_metadata.csv**<br>
Metadata of the observations wells, including the following information provided in the xsl files downloaded from DGA website (available in _input_ folder):
<small>
- cr2sub_id: unique observation well identifier in cr2sub database. The cr2sub_id is the same as the DGA well code, without the identification last digit.<br>
- dga_well_code: unique observation well identifier in DGA database, as reported in the xls downloaded file.<br>
- dga_well_name: observation well name in DGA database, as reported in the xls downloaded file.<br>
- dga_well_basin: basin in which the well is located, as reported by DGA in the xls downloaded file.<br>
- dga_well_subbasin: sub-basin in which the well is located, as reported by DGA in the xls downloaded file.<br>
- dga_well_lat: latitude of the well (_degrees, minutes, seconds_), as reported by DGA in the xls downloaded file.<br>
- dga_well_lon: longitude of the well (_degrees, minutes, seconds_), as reported by DGA in the xls downloaded file.<br>
- dga_well_utm_north: latitude of the well (_m UTM_), as reported by DGA in the xls downloaded file.<br>
- dga_well_utm_east: longitude of the well (_m UTM_), as reported by DGA in the xls downloaded file.<br>
- dga_well_elev: elevation of the well (_m a.s.l._), as reported by DGA in the xls downloaded file.<br>
</small>

**cr2sub_v1_attributes.csv**<br>
Attributes of the observations wells computed based on ancillary information. These attributes include metadata columns and the following computed characteristics:
<small>
- cr2sub_id: unique observation well identifier in cr2sub database. The cr2sub_id is the same as the DGA well code, without the identification last digit.<br>
- dga_well_code: unique observation well identifier in DGA database, as reported in the xls downloaded file.<br>
- dga_well_name: observation well name in DGA database, as reported in the xls downloaded file.<br>
- dga_well_basin: basin in which the well is located, as reported by DGA in the xls downloaded file.<br>
- dga_well_subbasin: sub-basin in which the well is located, as reported by DGA in the xls downloaded file.<br>
- dga_well_lat: latitude of the well (_degrees, minutes, seconds_), as reported by DGA in the xls downloaded file.<br>
- dga_well_lon: longitude of the well (_degrees, minutes, seconds_), as reported by DGA in the xls downloaded file.<br>
- dga_well_utm_north: latitude of the well (_m UTM_), as reported by DGA in the xls downloaded file.<br>
- dga_well_utm_east: longitude of the well (_m UTM_), as reported by DGA in the xls downloaded file.<br>
- dga_well_elev: elevation of the well (_m a.s.l._), as reported by DGA in the xls downloaded file.<br>
- cr2sub_lat: latitude of the well (_degrees_). The original DGA coordinate was corrected in some cases.<br> 
- cr2sub_lon: longitude of the well (_degrees_). The original DGA coordinate was corrected in some cases.<br>  
- cr2sub_utm_north_h19: latitude of the well (_m UTM H19_). The original DGA coordinate was corrected in some cases.<br>
- cr2sub_utm_south_h19: longitude of the well (_m UTM H19_). The original DGA coordinate was corrected in some cases.<br>
- cr2sub_elev: elevation of the well location (cr2sub_lat, cr2sub_lon) in _m a.s.l._, computed based on FABDEM (Forest And Building removed copernicus Digital
Elevation Model), downloaded from [University of Bristol repository](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn) and stored in the _input_ folder.<br>	
- cr2sub_slp: slope of the well location (cr2sub_lat, cr2sub_lon) in _%_ [**** rodrigo please chequear****], computed based on FABDEM.<br>	
- cr2sub_mean_gwl: mean depth of GWL (_m_), computed for all available records in _cr2sub_v1_mon.csv_.<br>
- cr2sub_sd_gwl: standard deviation of GWL (_m_), computed for all available records in _cr2sub_v1_mon.csv_<br>
- cr2sub_cv_gwl: coefficient of variation of GWL (_-_), computed as cr2sub_mean_gwl/cr2sub_sd_gwl.<br>
- cr2sub_mean_gwl_clean: mean depth of GWL without outliers (_m_), computed for all available data in _cr2sub_v1_mon_clean.csv_.<br>
- cr2sub_sd_gwl_clean: standard deviation of GWL without outliers (_m_), computed for all available data in _cr2sub_v1_mon_clean.csv_.<br>	
- cr2sub_cv_gwl_clean: coefficient of variation of GWL without outliers (_-_), computed as cr2sub_mean_gwl_clean/cr2sub_sd_gwl_clean.<br>	
- cr2sub_in_basin_camels: gauge_id from the smallest CAMELS-CL basin where the well is located. The basin was identified by intersecting the well location (cr2sub_lat, cr2sub_lon) with CAMELS-CL v2 polygons.<br>  
- cr2sub_camels_elev: mean elevation (_m a.s.l._) of the CAMELS-CL basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_camels_slp: mean slope (_%_) [**** rodrigo please chequear****] of the CAMELS-CL basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_camels_pr_yr: mean annual precipitation (_mm_) of the CAMELS-CL basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010<br>  
- cr2sub_camels_pet_yr: mean annual potential evapotranspiration (_mm_) of the CAMELS-CL basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010<br>  
- cr2sub_camels_aridity: aridity index of the CAMELS-CL basin where the well is located, computed as cr2sub_camels_pet_yr/cr2sub_camels_pr_yr<br>  
- cr2sub_camels_snowf: mean annual snow fraction (_-_) of the CAMELS-CL basin where the well is located [**** rodrigo please chequear****] [rodrigo aquí por favor especificar brevemente la metodología ---- computed based on CAMELS-CL v2 and CR2MET v2.5 for the period 1980-2010]<br>  
- cr2sub_in_basin_bna: id from the smallest BNA basin where the well is located. The basin was identified by intersecting the well location (cr2sub_lat, cr2sub_lon) with BNA polygons downloaded from [DGA website](https://dga.mop.gob.cl/mapoteca-digital/).<br>  
- cr2sub_bna_elev: mean elevation (_m a.s.l._) of the BNA basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_bna_slp: mean slope (_m a.s.l._) of the BNA basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_bna_pr_yr: mean annual precipitation (_mm_) of the BNA basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010<br>	
- cr2sub_bna_pet_yr: mean annual potential evapotranspiration (_mm_) of the BNA basin where the well is located, computed based CR2MET v2.5 for the period 1980-2010<br>  
- cr2sub_bna_aridity: aridity index of the BNA basin where the well is located, computed as cr2sub_bna_pet_yr/cr2sub_bna_pr_yr<br>  	
- cr2sub_bna_snowf: mean annual snow fraction (_-_) of the BNA basin where the well is located [**** rodrigo please chequear****] [rodrigo aquí por favor especificar brevemente la metodología ---- computed based on CR2MET v2.5 for the period 1980-2010]<br> 
- cr2sub_clsoilmap_awc_0_100cm: available water content (_mm_) in soil horizon 0-100 cm at well location (cr2sub_lat, cr2sub_lon), computed from [CLSoilMap](https://www.nature.com/articles/s41597-023-02536-x) database (see readme file in _input_ folder). <br> 
- cr2sub_clsoilmap_awc_100_200cm: available water content (_mm_) in soil horizon 100-200 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_ksat_0_100cm: hydraulic conductivity (_cm/day_) in soil horizon 0-100 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_ksat_100_200cm: hydraulic conductivity (_cm/day_) in soil horizon 100-200 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_bulkd_0_100cm: bulk density (_g/cm³_) in soil horizon 0-100 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_bulkd_100_200cm: bulk density (_g/cm³_)in soil horizon 100-200 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_clay_0_100cm: clay content (_%_) in soil horizon 0-100 cm at well location, computed from CLSoilMaps. <br>   
- cr2sub_clsoilmap_clay_100_200cm: clay content (_%_) in soil horizon 100-200 cm at well location, computed from CLSoilMaps. <br>  
- cr2sub_clsoilmap_sand_0_100cm: sand content (_%_) in soil horizon 0-100 cm at well location, computed from CLSoilMaps. <br>   
- cr2sub_clsoilmap_sand_100_200cm: sand content (_%_) in soil horizon 100-200 cm at well location, computed from CLSoilMaps. <br>  
</small>


### Scripts:
To support reproducible workflows for updating the database, this repository provides R scripts to download and process GWL observations available at DGA website. An initial bulk download procedure was done including the initial period of GWL records (1957–04-01) until 2025-03-31. Then, to exemplify how new updates can be processed, a new download was made on 2025-06-10 and merged into the initial bulk downloaded data. Subsequent updates can be performed following this procedure.

The directory is structured as follow:
<div style="font-size: 65%">

    IN/                         # Raw input data (xls files from DGA)
        base_2025-03-31/        # Historical data set spanning 1957–04-01 to 2025-03-31
        update_2025-06-10/      # Updated batch downloaded from DGA on 2025-06-10
    
    OUT/                        # Processed outputs
        base_2025-03-31/        # Historical data set spanning 1957–04-01 to 2025-03-31
          cr2sub_metadata.csv   # GWL observation well metadata
          cr2sub_mon.csv        # Processed monthly time series
          cr2sub_raw.csv        # Raw daily time series
       
        merge/                  # Merged data
          cr2sub_metadata.csv   # GWL observation well metadata
          cr2sub_mon.csv        # Processed monthly time series
          cr2sub_raw.csv        # Raw daily time series
           
        update_2025-06-10/      # Updated batch downloaded from DGA on 2025-06-10
          cr2sub_metadata.csv   # GWL observation well metadata
          cr2sub_mon.csv        # Processed monthly time series
          cr2sub_raw.csv        # Raw daily time series 
    
    scripts/                      
        00_functions.R          # Functions for reading, transforming, and converting data
        01_read_xls.R           # Script to process individual xls files
        02_merge_timeseries.R   # Script to merge time series
    
    README.md                   # Markdown documentation

</div>


##### 1) Download new data: 

The first step is to download the data from DGA website (https://snia.mop.gob.cl/BNAConsultas/reportes) and store it into the folder `IN/`. The folder `base_2025-03-31/` contains a bulk download of historical GWL data covering the period 1957–04-01 to 2025-03-31. Subsequent updates should be stored in folders prefixed with `update_` followed by the date of the update (e.g., `update_2025-06-10/`), each representing a new partial update from DGA. If the new downloaded spreadsheet contains fewer wells than the original dataset (`base_2025-03-31/`), the processed observations from unavailable wells are filled with NA in the next step.


##### 2) Extract and process DGA xls files: 

`01_read_xls.R`: reads xls files downloaded from the DGA website, extracts groundwater levels, cleans metadata, and exports raw and monthly time series per well.

For each update, the script `01_read_xls.R` should be run to generate three output files: 1) **cr2sub_metadata.csv**: table of well metadata. 2) **cr2sub_raw.csv**: records as reported by DGA. 3) **cr2sub_mon.csv**: monthly time series. If more than one observation was available during a month, the average is computed as monthly value.

The metadata (cr2sub_metadata.csv) includes the following information: 

- well_id: original DGA code
- well_name: original DGA name
- well_basin: BNA basin where the observation well is located
- well_subbasin: BNA subbasin where the observation well is located
- well_elev: elevation of the well location
- well_lat: latitude reported by DGA (epsg:4326) 
- well_lon: longitude reported by DGA (epsg:4326)
- well_north: latitude reported by DGA (UTM)
- well_east: longitude reported by DGA (UTM)


##### 3) Consolidate 

`02_merge_timeseries.R`: merges CSVs from available updates, resolves duplicated columns, and ensures one value per well per date, using the first non-`NA` when overlaps exist. The output dataset merges the base dataset with each subsequent update, ensuring that metadata and time series remain integrated and up to date.

