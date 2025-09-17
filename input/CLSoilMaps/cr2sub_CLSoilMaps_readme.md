## Soil Database for cr2sub wells

File name: **cr2sub_CLSoilMaps_data.csv**

### General description

This file contains a database of soil variables associated with wells distributed across different regions of Chile. The data were extracted from the **CLSoilMaps** product [Dinamarca et al., 2023](https://www.nature.com/articles/s41597-023-02536-x), a national database of soil properties. The database can be downloaded from [CLSoilMaps zenodo repository](https://doi.org/10.5281/zenodo.7464210). 

### Content details

* **Data source**: CLSoilMaps (Dinamarca et al., 2023)
* **Spatial extraction unit**: average of the variable within a circular buffer of 500 meters around each well
* **File format**: CSV
* **Missing values (NA)**: indicate wells located outside the CLSoilMaps coverage, such as water bodies, glaciers, or urban areas

### Contained variables

The following soil properties were extracted for each well:

| Variable | Description             | Unit   |
| -------- | ----------------------- | ------ |
| `ksat`   | Hydraulic conductivity  | cm/day |
| `Clay`   | Clay content            | %      |
| `Sand`   | Sand content            | %      |
| `Bulkd`  | Bulk density            | g/cm³  |
| `AWC`    | Available water content | mm     |

These variables are reported for each of the following 6 soil horizons:

* 0–5 cm
* 5–15 cm
* 15–30 cm
* 30–60 cm
* 60–100 cm
* 100–200 cm

The column names in the file follow the structure:
`<variable>_<depth>`
For example: `ksat_0_5`, `clay_5_15`, `bulkd_30_60`, etc.

### Synthetic horizon: 0–100 cm

Average or cumulative values were calculated for a composite horizon of 0–100 cm. The calculation methods are:

#### Thickness-weighted average (for `ksat`, `Clay`, `Sand`, `Bulkd`):

For a property \$P\$, the aggregated value for the 0–100 cm profile was calculated as:

$$
P_{0-100} = \frac{P_1 d_1 + P_2 d_2 + P_3 d_3 + P_4 d_4 + P_5 d_5}{d_1 + d_2 + d_3 + d_4 + d_5}
$$

Where:

* \$P\_1\$ to \$P\_5\$ are the property values for each of the first five horizons
* \$d\_1\$ to \$d\_5\$ are the thicknesses (in cm) of each horizon:

  * \$d\_1 = 5\$ (0–5 cm)
  * \$d\_2 = 10\$ (5–15 cm)
  * \$d\_3 = 15\$ (15–30 cm)
  * \$d\_4 = 30\$ (30–60 cm)
  * \$d\_5 = 40\$ (60–100 cm)

#### Direct sum of layers (for `AWC`):

$$
awc_{0-100} = awc_{0-5} + awc_{5-15} + awc_{15-30} + awc_{30-60} + awc_{60-100}
$$
