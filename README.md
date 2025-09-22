# CR2SUB: consolidated monthly groundwater level database for Chile

The **CR2SUB** database compiles groundwater level (GWL) observations from 1,137 wells maintained by the Water Bureau of Chile (Dirección General de Aguas, DGA), downloaded from [DGA hydrometeorological data portal](https://snia.mop.gob.cl/BNAConsultas/reportes), and processed into a homogenized monthly database for the period **1957-2025**. The database also includes metadata for each observation well and a range of topographical and hydroclimatic attributes computed from ancillary information.  

This repository is structured in three main folders, explained in detail below: 
- **cr2sub**: csv files with GWL timeseries and well attributes.
- **input**: data used for processing the database.
- **scripts**: scripts used to process the database.

All data and scripts are openly provided to ensure full reproducibility of results. The workflow also enables users to update the CR2SUB database whenever new DGA records are downloaded. 

**Developer**: The CR2SUB database has been consolidated within the framework of the Center for Climate and Resilience Research (CR2, ANID/FONDAP/1523A0002) and is part of the research project ANID/FONDECYT/11240924. Responsible: Camila Alvarez-Garreton. Colaborators: Rodrigo Marinao Rivas, Diego Dinamarca.

**License**: This dataset is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. 

**Access**: [CR2SUB - Zenodo](https://zenodo.org/records/xxxxxxx); [CR2SUB - Github](https://github.com/calvarezgarreton/cr2sub).  

**Citation**:  
Zenodo: Alvarez-Garreton, C. (2025). *CR2SUB: Consolidated monthly groundwater level database for Chile* [Data set]. Zenodo. [https://doi.org/10.xxxx/zenodo.xxxxxx](https://doi.org/10.xxxx/zenodo.xxxxxx).<br>
Article: Alvarez-Garreton, C., Boisier, J. P., Marinao, R., et al. (in prep./2025). *Recovery and memory of surface and groundwater systems after droughts in Chile*. *Hydrology and Earth System Sciences* (in prep.).  


## cr2sub

**cr2sub_v1.1_gwl_mon.csv**<br>
 Consolidated monthly time series of raw groundwater level (GWL) observations. Raw GLW observations were downloaded as xls spreadsheets from the [DGA website](https://snia.mop.gob.cl/BNAConsultas/reportes), stored in _input_ and processed through the codes provided in _scripts_. When more than one record was available for a given month, the values were averaged to obtain a single monthly estimate.

**cr2sub_v1.1_gwl_mon_clean.csv**<br>
 Monthly time series after applying an outlier-removal procedure. The code for outlier detection and removal is available in _scripts_.

**cr2sub_v1.1_attributes.csv**<br>
Metadata of the observation wells extracted from the xsl files downloaded from DGA website (_input_), and a range of climatic, topographical and pedological attributes computed from on ancillary information (available in _input_):
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
- cr2sub_elev: elevation of the well location (cr2sub_lat, cr2sub_lon) in _m a.s.l._, computed based on FABDEM (Forest And Building removed copernicus Digital Elevation Model), available from [University of Bristol repository](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn).<br>	
- cr2sub_slp: slope of the well location (cr2sub_lat, cr2sub_lon) in _%_ [**** rodrigo please chequear****], computed from FABDEM.<br>	
- cr2sub_mean_gwl: mean depth of GWL (_m_), computed for all available records in _cr2sub_v1.1_gwl_mon.csv_.<br>
- cr2sub_sd_gwl: standard deviation of GWL (_m_), computed for all available records in _cr2sub_v1.1_gwl_mon.csv_.<br>
- cr2sub_cv_gwl: coefficient of variation of GWL (_-_), computed as cr2sub_mean_gwl/cr2sub_sd_gwl.<br>
- cr2sub_mean_gwl_clean: mean depth of GWL without outliers (_m_), computed for all available data in _cr2sub_v1.1_gwl_mon_clean.csv_.<br>
- cr2sub_sd_gwl_clean: standard deviation of GWL without outliers (_m_), computed for all available data in _cr2sub_v1.1_gwl_mon_clean.csv_.<br>	
- cr2sub_cv_gwl_clean: coefficient of variation of GWL without outliers (_-_), computed as cr2sub_mean_gwl_clean/cr2sub_sd_gwl_clean.<br>	
- cr2sub_in_basin_camels: gauge_id from the smallest CAMELS-CL basin where the well is located. The basin was identified by intersecting the well location (cr2sub_lat, cr2sub_lon) with CAMELS-CL v2021 polygons and store as lookup table in _input_.<br>  
- cr2sub_camels_elev: mean elevation (_m a.s.l._) of the CAMELS-CL basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_camels_slp: mean slope (_%_) [**chequear**] of the CAMELS-CL basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_camels_pr_yr: mean annual precipitation (_mm_) of the CAMELS-CL basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010. CR2MET v2.5 daily data were downloaded from https://ftp.cr2.cl/browse/cr2met/v2.5 and processed into annual values (_input_).<br>  
- cr2sub_camels_pet_yr: mean annual potential evapotranspiration (_mm_) of the CAMELS-CL basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010.<br>  
- cr2sub_camels_aridity: aridity index of the CAMELS-CL basin where the well is located, computed as cr2sub_camels_pet_yr/cr2sub_camels_pr_yr.<br>  
- cr2sub_camels_snowf: snow fraction (_-_) of the CAMELS-CL basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010.<br>  
- cr2sub_in_basin_bna: id from the smallest BNA basin where the well is located. The basin was identified by intersecting the well location (cr2sub_lat, cr2sub_lon) with BNA polygons downloaded from [DGA spatial data repository](https://dga.mop.gob.cl/mapoteca-digital/) and store as lookup table in _input_.<br>  
- cr2sub_bna_elev: mean elevation (_m a.s.l._) of the BNA basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_bna_slp: mean slope (%) [**chequear**] of the BNA basin where the well is located, computed based on FABDEM.<br>  
- cr2sub_bna_pr_yr: mean annual precipitation (_mm_) of the BNA basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010<br>	
- cr2sub_bna_pet_yr: mean annual potential evapotranspiration (_mm_) of the BNA basin where the well is located, computed based CR2MET v2.5 for the period 1980-2010<br>  
- cr2sub_bna_aridity: aridity index of the BNA basin where the well is located, computed as cr2sub_bna_pet_yr/cr2sub_bna_pr_yr<br>  	
- cr2sub_bna_snowf: snow fraction (_-_) of the BNA basin where the well is located, computed based on CR2MET v2.5 for the period 1980-2010.<br> 
- cr2sub_clsoilmap_awc_0_100cm: available water content (_mm_) in soil horizon 0-100 cm at well location (cr2sub_lat, cr2sub_lon), computed from [CLSoilMaps database](https://www.nature.com/articles/s41597-023-02536-x), and stored in _input_. <br> 
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

## input
**BNA_basin_polygons**: shapefile with BNA basins boundaries downloaded from [DGA spatial data repository](https://dga.mop.gob.cl/mapoteca-digital/).<br>  

**camels_cl_basins**: shapefile with CAMELS-CL v2021 basins boundaries downloaded from [CAMELS-CL platform](https://www.camels.cr2.cl/).<br>  

**CLSoilmaps**: soil data computed from [CLSoilMap databse](https://www.nature.com/articles/s41597-023-02536-x). 

**cr2met**: precipitation, potential evpotranspiration and snowfraction data downloaded from https://ftp.cr2.cl/browse/cr2met/v2.5 and processed into annual values.

**dem**: elevation and slope data obtained from FABDEM (Forest And Building removed copernicus Digital
Elevation Model), available from [University of Bristol repository](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn).<br>

**DGA_GWL_observations**: raw GWL observations downloaded as xls spreadsheets from the [DGA website](https://snia.mop.gob.cl/BNAConsultas/reportes).
<small>
- dga_xls_1950-01-01_2025-03-31: a bulk download of GWL observations spreadsheets including the initial period of GWL records (1957–04-01) until 2025-03-31. These spreadsheets were downloaded manually for all available observations wells in Chile.
- dga_xls_update_2025-06-10 and dga_xls_update_2025-08-05: additional downloads made on 2025-06-10 amd 2025-08-05, respectively, that permits to exemplify how new updates can be processed and merged into the initial bulk downloaded data.
</small>

**other_data**:
<small>
- Chile_boundary_simplified_from_SIIT: shapefile with national boundaries downloaded from https://www.bcn.cl/siit/mapas_vectoriales. 
- cr2sub_v1.1_join_table_with_bna_basins: lookup table with the codes of the BNA basin(s) in which the CR2SUB observation wells are located.<br>  
- cr2sub_v1.1_join_table_with_camels_cl_v2021_basins: lookup table with the codes of the CAMELS-CL v2021 basin(s) in which the CR2SUB observation wells are located.<br>  

## scripts
This folder contains the scripts used to process the database. These can be run by users to update CR2SUB if new data from DGA is downloaded.
