# Evidencia de revisión — migración de input de niveles DGA

## Versiones comparadas

| | v1.1 (referencia) | v1.2 (esta rama) |
|---|---|---|
| Rama | `cr2sub_v1.1` | `cr2sub_v1.2` |
| Input DGA | 844 reportes `dga_xls_*` por pozo | export único `AM006T0034616_1.xlsx` |
| Carpeta de outputs | `versions/v1.1/` | `versions/v1.2/` |

Los CSV de cada carpeta llevan la versión en el nombre
(`cr2sub_v1.1_*` y `cr2sub_v1.2_*`), según la variable `version` de los
scripts de cada rama. Los resultados de
la comparación (v1.1 → v1.2) están en `comparison_v1.1_vs_v1.2/`.

**Objetivo:** reemplazar `input/DGA_GWL_obs_descarga/dga_xls_*` (844 reportes
individuales por pozo, descargados manualmente del portal DGA en 3 tandas)
por el archivo único `input/DGA_GWL_obs_descarga/AM006T0034616_1.xlsx`
(export tabular largo, 147.180 filas, 1 fila por pozo-fecha-valor).

**Scripts modificados:** [`scripts/01_consolidate.R`](../scripts/01_consolidate.R),
[`scripts/functions/functions_consolidate_cr2sub.R`](../scripts/functions/functions_consolidate_cr2sub.R)
(nueva función `deg2dec_compact` para el formato de coordenadas del nuevo
archivo). De paso se corrigió un bug preexistente, no relacionado a la
migración, en [`scripts/03_process_attributes.R`](../scripts/03_process_attributes.R)
(`cuencas_BNA.shp` → `Cuencas_BNA.shp`, mismatch de mayúsculas que impedía
correr el paso 3).

**Método:** se guardó una copia de los 3 outputs de la v1.1 (input
viejo) en `versions/v1.1/` antes de tocar nada, se corrió el
pipeline completo (`01_consolidate.R` → `02_remove_outliers.ipynb` →
`03_process_attributes.R`) con el input nuevo, y se comparó celda por celda
contra la v1.1 con [`scripts/compare_outputs.py`](scripts/compare_outputs.py)
(reproducible, ver sección "Cómo reproducir" abajo). Los resultados de esa
comparación quedan en `comparison_v1.1_vs_v1.2/`.

## 1. Cobertura de pozos — la brecha

| | Pozos v1.1 | Pozos v1.2 | En común |
|---|---|---|---|
| `gwl_mon.csv` / `gwl_mon_clean.csv` / `attributes.csv` | 1137 | 1133 | 1133 |

La brecha son **4 pozos** que existen en las descargas viejas pero no en el
export único nuevo (`comparison_v1.1_vs_v1.2/wells_missing_in_new_source.csv`):

| Código DGA | Nombre | Cuenca | Subcuenca |
|---|---|---|---|
| 02942002-5 | PIQUE AV-1 | QS. entre Q. La Negra y Q. Pan de Azúcar | Quebrada de Taltal |
| 02942008-4 | PIQUE AV - 42 | QS. entre Q. La Negra y Q. Pan de Azúcar | Quebrada de Taltal |
| 02942011-4 | ESSAN AV - 48 | QS. entre Q. La Negra y Q. Pan de Azúcar | Quebrada de Taltal |
| 02942023-8 | ESSAN AV - 48 | QS. entre Q. La Negra y Q. Pan de Azúcar | Quebrada de Taltal |

Los 4 son puntos de monitoreo privados (prefijo ESSAN/PIQUE AV, minería) en
la cuenca de Taltal, no estaciones públicas DGA — su ausencia del export del
portal es esperable y no se considera un defecto de la migración.

No hay pozos en el archivo nuevo que no existieran antes (0 pozos nuevos
"de más").

## 2. Comparación de metadatos (`attributes.csv`)

De los 1133 pozos en común, **17 valores de metadata** (en 13 pozos)
difieren en más de 0.01 entre el input viejo y el nuevo
(`comparison_v1.1_vs_v1.2/metadata_diffs.csv`):

| cr2sub_id | Nombre | Campo | v1.1 | v1.2 |
|---|---|---|---|---|
| 1700042 | SALAR BELLAVISTA | lat | -21.305833 | -20.734444 |
| 1700042 | SALAR BELLAVISTA | lon | -69.520278 | -69.618056 |
| 1700042 | SALAR BELLAVISTA | elev | 946 | 970 |
| 4400037 | PARCELA 48 CUATRO ESQUINAS | lat | -29.731944 | -29.918889 |
| 4400037 | PARCELA 48 CUATRO ESQUINAS | lon | -71.230278 | -71.269722 |
| 4308010 | POZO CALVARIO | lon | -70.594722 | -70.502778 |
| 8114002 | POZO COIHUECO en Internado Municipal | lon | -71.871667 | -71.838056 |
| 1700051/057/076/099-106 (JICA, 9 pozos) | — | elevación | — | -25 a -35 m |

Son correcciones/actualizaciones reales en la base de metadatos de DGA entre
la fecha en que se bajaron los reportes viejos y la fecha del export nuevo
(re-levantamientos de coordenadas/cota), no errores de parsing: el formato
de coordenadas del archivo nuevo (`DDDMMSS` compacto) es uniforme en los
1134 pozos y se validó contra ubicaciones conocidas antes de escribir el
parser. Quedan documentadas para que quien revise sepa que no son
sorpresa de la migración, pero no se "corrigieron" en ningún sentido — se
usa el valor que trae cada fuente.

Las métricas agregadas (`cr2sub_mean_gwl`, `cr2sub_sd_gwl`) cambian
levemente en ~110-120 pozos — esperado, porque ahora hay más meses de datos
por pozo (huecos rellenados + serie extendida).

## 3. Series que evidencian la extensión temporal

**Rango de fechas:**

| Archivo | v1.1 | v1.2 |
|---|---|---|
| `gwl_mon.csv` | 1957-04 → 2024-12 (813 meses) | 1957-04 → 2026-08 (833 meses) |
| `gwl_mon_clean.csv` | igual | igual |

Ejemplo concreto — pozo `5423019` (el de mayor cobertura reciente), serie
mensual completa (`cr2sub_v1.2_gwl_mon.csv`), generada con
[`scripts/plot_series_extension.py`](scripts/plot_series_extension.py):

![Serie de tiempo pozo 5423019: registro existente en azul, puntos agregados por la actualización en naranja](comparison_v1.1_vs_v1.2/series_extension_5423019.png)

En azul, el registro que ya existía con el input viejo: 344 meses entre
1974-04 y 2024-12 (los cortes de la línea son huecos reales de medición, no
artefactos del gráfico). En naranja, los 20 meses que agrega la
actualización, de 2025-01 a 2026-08, enganchados directo a la cola del
registro existente — no hay reescritura del histórico, solo extensión hacia
adelante. Todos los puntos del tramo en común entre ambas fuentes (hasta
2024-12) son idénticos; el detalle completo de celdas que antes eran `NA` y
ahora tienen valor está en `comparison_v1.1_vs_v1.2/gwl_mon_gap_filled.csv` (420 en
`gwl_mon.csv`, 840 en `gwl_mon_clean.csv`, repartidas en huecos intermedios
de otros pozos además de esta cola).

## 4. Cobertura mensual por región, 2020–2026

Número de pozos con dato en cada mes, por región (norte a sur), para v1.1
y v1.2 con la misma escala de color. La región de cada pozo se asigna por
ubicación (cruce espacial con los límites regionales; el polígono más
cercano si el punto cae fuera), no por el código DGA, cuyo prefijo es la
macrocuenca BNA y no la región. El panel superior (total mensual) superpone ambas
versiones: v1.2 en azul grueso y v1.1 en naranja fino. Generado con
[`scripts/plot_coverage_by_region.py`](scripts/plot_coverage_by_region.py).

![Cobertura mensual de pozos por región, v1.1 vs v1.2, 2020-2026](comparison_v1.1_vs_v1.2/cobertura_mensual_por_region_v1.1_vs_v1.2_2020-2026.png)

- Hasta dic-2024 las dos versiones son prácticamente idénticas; la línea
  naranja corre dentro de la azul.
- Desde 2025 solo v1.2 tiene datos (hasta ago-2026), concentrados en
  Valparaíso, Tarapacá, Coquimbo y O'Higgins. Las demás regiones no
  aportan meses nuevos.
- En pozos presentes en ambas versiones, v1.2 suma 600 pozos-mes en
  2020–2026 (557 desde 2025 y 43 de huecos rellenados antes) y pierde 1.
  Por región: Valparaíso 362, Tarapacá 127, Coquimbo 63, O'Higgins 43,
  Arica y Parinacota 3, Maule 2.
- Los otros 21 pozos-mes que v1.2 no tiene en el periodo corresponden a
  los 4 pozos privados de Taltal que no están en el export nuevo (sección 1).
- La cobertura cae desde ~2020 en ambas versiones y las regiones al sur de
  O'Higgins casi no tienen datos desde 2021–2023.

## 5. Consistencia de valores (no solo cobertura)

Comparando las **921.129 celdas** (pozo × fecha) presentes en ambas fuentes
en `gwl_mon.csv`:

| v1.1 → v1.2 | Cantidad |
|---|---|
| Valores idénticos o ambos `NA` | 920.708 |
| `NA` antes → con valor ahora (hueco rellenado) | 420 |
| Con valor antes → `NA` ahora (dato perdido) | 1 |
| **Valores contradictorios** (mismo pozo/fecha, valor distinto) | **0** |

El único dato "perdido" es el pozo `1700028` en 2024-04
(`comparison_v1.1_vs_v1.2/gwl_mon_lost.csv`). En `gwl_mon_clean.csv` (post
remoción de outliers) la cifra sube a 16 perdidos / 840 rellenados — es
esperable, porque el filtro de outliers recalcula con más contexto
temporal y puede marcar como outlier algo que antes, con una serie más
corta, no se evaluaba igual.

**Conclusión:** cero conflictos de valor en ~920 mil puntos comparados. El
archivo único es, en la práctica, un superset de las 3 carpetas viejas
(mismos valores donde coinciden, más puntos donde antes había huecos),
salvo los 4 pozos privados fuera del export del portal.

## 6. Atributos columna a columna (1133 pozos en común)

Generado con [`scripts/compare_attributes.py`](scripts/compare_attributes.py)
(`comparison_v1.1_vs_v1.2/attribute_column_diffs.csv`). De las 45 columnas de
`attributes.csv`, solo `cr2sub_id` es idéntica; los cambios se agrupan así:

| Grupo | Pozos que cambian | Causa |
|---|---|---|
| `dga_well_basin` / `subbasin` / `name` | 1133 / 712 / 45 | La DGA cambió el formato del texto (`RIO COPIAPO` → `R.  Copiapo`; `EL TRANSITO` → `POZO EL TRANSITO`). No es un cambio de ubicación |
| Coordenadas y cota `dga_*` / `cr2sub_*` | 14 a 33 según columna | Re-levantamientos en la base DGA (sección 2). Arrastran `cr2sub_elev` y `cr2sub_slp` (14 pozos) |
| `cr2sub_mean/sd/cv_gwl` (raw y clean) | 61 a 142 | Más meses de datos por pozo (secciones 3 y 5) |
| `cr2sub_*camels*` | 273 cambian de cuenca, 77 pasan de NA a tener cuenca, 1 pasa a NA | Cambio de CAMELS-CL v1 (516 cuencas) a v2 (860). Las variables climáticas y topográficas de la cuenca se recalculan sobre el nuevo polígono |
| `cr2sub_*bna*` | 9 | Nuevo criterio de asignación: cuenca que contiene al pozo, con buffer de 500 m solo si no cae en ninguna (antes el buffer se aplicaba a todos los pozos) |
| `cr2sub_clsoilmap_*` | 0 | Sin cambios; ver sección 8 |

Pozos sin cuenca CAMELS: 310 en v1.1 → 230 en v1.2. El pozo `5221007`
(PAPUDO Y ZAPALLAR, mismas coordenadas) tenía cuenca en v1 y ya no en v2,
porque el polígono cambió.

## 7. Pipeline

Corrida completa con el input nuevo (`main_pipeline.py`: 01 → 01b → 02 → 03).
Las salidas coinciden con las secciones 1 a 5; los CSV de `versions/v1.2/` son
copia de `cr2sub/`. La asignación a cuencas ahora se regenera con
`input/cr2sub_basin_lookup_tables/lookup_table_script.R` (intersección
vectorial, sin dependencias externas).

## 8. CLSoilMaps entre pozos y entre versiones

Scripts: [`scripts/validate_clsoilmaps_extraction.R`](scripts/validate_clsoilmaps_extraction.R)
(re-extrae el suelo de todos los pozos con los rasters originales) y
`compare_attributes.py`. Salidas: `clsoilmaps_*.csv` y
`clsoilmaps_stored_vs_reextracted.png` en `comparison_v1.1_vs_v1.2/`.

**v1.1 vs v1.2.** Los valores de suelo son idénticos en 1131 de los 1133
pozos comunes. Los otros 2 son los pozos nuevos (`05120007` y `05120011`) que
v1.1 tenía en NA y que ahora se extrajeron de los rasters. Pozos con suelo en
NA: 66 de 1133 (fuera de la cobertura de CLSoilMaps).

**Coherencia espacial.** Para cada pozo y su vecino más cercano a menos de
5 km (884 pares), la correlación del valor de suelo es 0,77 a 0,95 (por
ejemplo arena 0-100 cm: 0,95; arcilla 0-100 cm: 0,92), frente a ~0 para
pares al azar. Los dos pozos nuevos quedan en el rango de sus vecinos
(AWC 0-100 cm: 130 y 126 mm; vecinos entre 120 y 139 mm).

**Verificaciones físicas.** Sin valores negativos y densidad aparente entre
0,7 y 2,0 g/cm³ en todos los pozos. Hay 9 pozos con arcilla + arena > 100 %
en 100-200 cm; son promedios sobre un buffer de 500 m, no mediciones, y
conviene tenerlo presente.

**Re-extracción vs valores guardados** (1066 pozos con dato en ambos):
correlación 0,997 a 0,999; mediana de la diferencia relativa 0,06 % a 0,5 %.
Solo 106 a 598 pozos por variable son idénticos al milésimo, así que el
buffer no se construye exactamente igual que el original, pero la diferencia
es pequeña (mediana de AWC 0-100 cm: 0,1 mm sobre una mediana de 136 mm).

**Hallazgo importante.** Los tres pozos con mayor desvío son pozos cuyas
coordenadas cambiaron entre v1.1 y v1.2 (`4308010`, `1700042` y `8114002`,
sección 2). Su suelo guardado se calculó con la ubicación vieja y no se
recalculó, porque el paso 01b solo extrae pozos que faltan en el CSV. El
cuarto pozo que cambió de coordenadas, `4400037`, tiene valor guardado
(AWC 0-100 cm: 94,7 mm) pero la re-extracción da NA. Además `6011022` y
`4400017` están en NA en el valor guardado y con dato en la re-extracción
(136,7 y 124,4 mm). Conviene hacer que el paso 01b recalcule también los
pozos cuyas coordenadas cambian y revisar esos casos de borde de cobertura.

## Cómo reproducir

Todo junto (tras `main_pipeline.py`): `review_dga_xlsx_migration/scripts/regenerate_review.sh`. Por partes:

```bash
cd review_dga_xlsx_migration/scripts
python3 compare_outputs.py        # series y metadatos; solo stdlib
Rscript ../../review_dga_xlsx_migration/scripts/validate_clsoilmaps_extraction.R   # desde la raíz del repo
python3 compare_attributes.py     # atributos y CLSoilMaps; pandas, matplotlib
python3 plot_series_extension.py; python3 plot_coverage_by_region.py   # figuras; geopandas
```

Lee de `../versions/{v1.1,v1.2}/*.csv` y escribe a
`../comparison_v1.1_vs_v1.2/`. Para regenerar `versions/v1.2/` desde cero, correr el
pipeline (`scripts/main_pipeline.py`, requiere pandas para el paso 02)
y copiar los 3 CSV de `cr2sub/` ahí.
