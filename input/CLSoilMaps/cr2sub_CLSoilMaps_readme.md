# Soil Database for cr2sub wells

File name: `cr2sub_CLSoilMaps_data_{tag}_{version}.csv` (for example,
`cr2sub_CLSoilMaps_data_cr2sub_v1.2.csv`). The tag and version are configured
in `scripts/main_pipeline.py`.

## General description

This file contains a database of soil variables associated with wells across Chile. The data were extracted from **CLSoilMaps** (Dinamarca et al., 2023), a national database of soil properties.

## Content details

* **Data source**: CLSoilMaps (Dinamarca et al., 2023)
* **Spatial extraction unit**: average of the variable within a circular buffer of 500 meters around each well
* **File format**: CSV
* **Missing values (NA)**: indicate wells located outside the CLSoilMaps coverage, such as water bodies, glaciers, or urban areas

## Contained variables

The CSV contains ten CLSoilMaps attributes per well. Each value is the mean within a 500 m buffer around that well.

| Variable | Description             | Unit   |
| -------- | ----------------------- | ------ |
| `ksat`   | Hydraulic conductivity  | cm/day |
| `Clay`   | Clay content            | %      |
| `Sand`   | Sand content            | %      |
| `Bulkd`  | Bulk density            | g/cm³  |
| `AWC`    | Available water content | mm     |

The CSV reports each property for two depth ranges:

* 0–100 cm
* 100–200 cm

The identifier column is `cr2sub_id`; the other columns use the prefix
`cr2sub_clsoilmap_` followed by the property and depth, for example
`cr2sub_clsoilmap_awc_0_100cm` and
`cr2sub_clsoilmap_ksat_100_200cm`. The 0–100 cm values are calculated from
the five source layers from 0 to 100 cm; the 100–200 cm values come from the
100–200 cm source layer.

## Synthetic horizon: 0–100 cm

Average or cumulative values were calculated for a composite horizon of 0–100 cm. The calculation methods are:

### Thickness-weighted average (for `ksat`, `Clay`, `Sand`, `Bulkd`):

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

### Direct sum of layers (for `AWC`):

$$
awc_{0-100} = awc_{0-5} + awc_{5-15} + awc_{15-30} + awc_{30-60} + awc_{60-100}
$$

---


# Base de datos de suelos de pozos de Chile

Nombre archivo: `cr2sub_CLSoilMaps_data_{tag}_{version}.csv`

El pipeline define `{tag}` y `{version}` en `scripts/main_pipeline.py` (por
ejemplo, `cr2sub_CLSoilMaps_data_cr2sub_v1.2.csv`) para mantener cada versión
separada y facilitar su comparación.

## Descripción general

Este archivo contiene una base de datos de variables de suelo asociadas a pozos distribuidos en distintas regiones de Chile. Los datos fueron extraídos del producto **CLSoilMaps** (Dinamarca et al., 2023), una base de datos nacional de propiedades edáficas.

## Detalles del contenido

* Fuente de datos: CLSoilMaps (Dinamarca et al., 2023)
* Unidad espacial de extracción: promedio de la variable en un buffer circular de 500 metros alrededor de cada pozo
* Formato de archivo: CSV
* Valores faltantes (NA): indican pozos ubicados fuera de la cobertura de CLSoilMaps, como cuerpos de agua, glaciares o zonas urbanas

## Variables contenidas

El CSV contiene diez atributos CLSoilMaps por pozo. Cada valor corresponde al promedio dentro de un buffer de 500 m alrededor del pozo.

| Variable | Descripción              | Unidad |
| -------- | ------------------------ | ------ |
| `ksat`   | Conductividad hidráulica | cm/día |
| `Clay`   | Contenido de arcilla     | %      |
| `Sand`   | Contenido de arena       | %      |
| `Bulkd`  | Densidad aparente        | g/cm³  |
| `AWC`    | Humedad aprovechable     | mm     |

El CSV presenta cada propiedad para dos rangos de profundidad:

* 0–100 cm
* 100–200 cm

La columna identificadora es `cr2sub_id`; las demás columnas usan el prefijo
`cr2sub_clsoilmap_` seguido de la propiedad y la profundidad, por ejemplo
`cr2sub_clsoilmap_awc_0_100cm` y
`cr2sub_clsoilmap_ksat_100_200cm`. Los valores de 0–100 cm se calculan a partir
de las cinco capas de origen entre 0 y 100 cm; los valores de 100–200 cm
provienen de la capa de origen de 100–200 cm.

## Horizonte sintético: 0–100 cm

Se calcularon valores promedio o acumulados para un horizonte compuesto de 0–100 cm. Los métodos de cálculo son:

### Promedio ponderado por grosor de horizonte (para `ksat`, `Clay`, `Sand`, `Bulkd`):

Para una propiedad $P$, el valor agregado para el perfil 0–100 cm se calculó como:

P_0_100 = (P1 * d1 + P2 * d2 + P3 * d3 + P4 * d4 + P5 * d5) / (d1 + d2 + d3 + d4 + d5)

Donde:

    P1 a P5 son los valores de la propiedad en cada uno de los primeros cinco horizontes

    d1 a d5 son los espesores (en cm) de cada horizonte:

        d1 = 5 (0–5 cm)

        d2 = 10 (5–15 cm)

        d3 = 15 (15–30 cm)

        d4 = 30 (30–60 cm)

        d5 = 40 (60–100 cm)


### Suma directa de capas (para `AWC`):

awc_0_100 = awc_0_5 + awc_5_15 + awc_15_30 + awc_30_60 + awc_60_100
