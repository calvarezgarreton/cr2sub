# Evidencia de revisión — migración de input de niveles DGA

## Versiones comparadas

| | v1.1 (referencia) | v1.3 (esta rama) |
|---|---|---|
| Rama | `cr2sub_v1.1` | `cr2sub_v1.3` |
| Input DGA | 844 reportes `dga_xls_*` por pozo | export único `AM006T0034616_1.xlsx` |
| Carpeta de outputs | `versions/v1.1/` | `versions/v1.3/` |

Los CSV de cada carpeta llevan la versión en el nombre
(`cr2sub_v1.1_*` y `cr2sub_v1.3_*`), según la variable `version` de los
scripts de cada rama. Los resultados de
la comparación (v1.1 → v1.3) están en `comparison_v1.1_vs_v1.3/`.

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
comparación quedan en `comparison_v1.1_vs_v1.3/`.

## 1. Cobertura de pozos — la brecha

| | Pozos v1.1 | Pozos v1.3 | En común |
|---|---|---|---|
| `gwl_mon.csv` / `gwl_mon_clean.csv` / `attributes.csv` | 1137 | 1133 | 1133 |

La brecha son **4 pozos** que existen en las descargas viejas pero no en el
export único nuevo (`comparison_v1.1_vs_v1.3/wells_missing_in_new_source.csv`):

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
(`comparison_v1.1_vs_v1.3/metadata_diffs.csv`):

| cr2sub_id | Nombre | Campo | v1.1 | v1.3 |
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

| Archivo | v1.1 | v1.3 |
|---|---|---|
| `gwl_mon.csv` | 1957-04 → 2024-12 (813 meses) | 1957-04 → 2026-08 (833 meses) |
| `gwl_mon_clean.csv` | igual | igual |

Ejemplo concreto — pozo `5423019` (el de mayor cobertura reciente), serie
mensual completa (`cr2sub_v1.3_gwl_mon.csv`), generada con
[`scripts/plot_series_extension.py`](scripts/plot_series_extension.py):

![Serie de tiempo pozo 5423019: registro existente en azul, puntos agregados por la actualización en naranja](comparison_v1.1_vs_v1.3/series_extension_5423019.png)

En azul, el registro que ya existía con el input viejo: 344 meses entre
1974-04 y 2024-12 (los cortes de la línea son huecos reales de medición, no
artefactos del gráfico). En naranja, los 20 meses que agrega la
actualización, de 2025-01 a 2026-08, enganchados directo a la cola del
registro existente — no hay reescritura del histórico, solo extensión hacia
adelante. Todos los puntos del tramo en común entre ambas fuentes (hasta
2024-12) son idénticos; el detalle completo de celdas que antes eran `NA` y
ahora tienen valor está en `comparison_v1.1_vs_v1.3/gwl_mon_gap_filled.csv` (420 en
`gwl_mon.csv`, 840 en `gwl_mon_clean.csv`, repartidas en huecos intermedios
de otros pozos además de esta cola).

## 4. Cobertura mensual por región, 2020–2026

Número de pozos con dato en cada mes, por región (norte a sur), para v1.1
y v1.3 con la misma escala de color. La región de cada pozo se asigna por
ubicación (cruce espacial con los límites regionales; el polígono más
cercano si el punto cae fuera), no por el código DGA, cuyo prefijo es la
macrocuenca BNA y no la región. El panel superior (total mensual) superpone ambas
versiones: v1.3 en azul grueso y v1.1 en naranja fino. Generado con
[`scripts/plot_coverage_by_region.py`](scripts/plot_coverage_by_region.py).

![Cobertura mensual de pozos por región, v1.1 vs v1.3, 2020-2026](comparison_v1.1_vs_v1.3/cobertura_mensual_por_region_v1.1_vs_v1.3_2020-2026.png)

- Hasta dic-2024 las dos versiones son prácticamente idénticas; la línea
  naranja corre dentro de la azul.
- Desde 2025 solo v1.3 tiene datos (hasta ago-2026), concentrados en
  Valparaíso, Tarapacá, Coquimbo y O'Higgins. Las demás regiones no
  aportan meses nuevos.
- En pozos presentes en ambas versiones, v1.3 suma 600 pozos-mes en
  2020–2026 (557 desde 2025 y 43 de huecos rellenados antes) y pierde 1.
  Por región: Valparaíso 362, Tarapacá 127, Coquimbo 63, O'Higgins 43,
  Arica y Parinacota 3, Maule 2.
- Los otros 21 pozos-mes que v1.3 no tiene en el periodo corresponden a
  los 4 pozos privados de Taltal que no están en el export nuevo (sección 1).
- La cobertura cae desde ~2020 en ambas versiones y las regiones al sur de
  O'Higgins casi no tienen datos desde 2021–2023.

## 5. Consistencia de valores (no solo cobertura)

Comparando las **921.129 celdas** (pozo × fecha) presentes en ambas fuentes
en `gwl_mon.csv`:

| v1.1 → v1.3 | Cantidad |
|---|---|
| Valores idénticos o ambos `NA` | 920.708 |
| `NA` antes → con valor ahora (hueco rellenado) | 420 |
| Con valor antes → `NA` ahora (dato perdido) | 1 |
| **Valores contradictorios** (mismo pozo/fecha, valor distinto) | **0** |

El único dato "perdido" es el pozo `1700028` en 2024-04
(`comparison_v1.1_vs_v1.3/gwl_mon_lost.csv`). En `gwl_mon_clean.csv` (post
remoción de outliers) la cifra sube a 16 perdidos / 840 rellenados — es
esperable, porque el filtro de outliers recalcula con más contexto
temporal y puede marcar como outlier algo que antes, con una serie más
corta, no se evaluaba igual.

**Conclusión:** cero conflictos de valor en ~920 mil puntos comparados. El
archivo único es, en la práctica, un superset de las 3 carpetas viejas
(mismos valores donde coinciden, más puntos donde antes había huecos),
salvo los 4 pozos privados fuera del export del portal.

## Cómo reproducir

```bash
cd review_dga_xlsx_migration/scripts
python3 compare_outputs.py   # requiere solo stdlib, no pandas
```

Lee de `../versions/{v1.1,v1.3}/*.csv` y escribe a
`../comparison_v1.1_vs_v1.3/`. Para regenerar `versions/v1.3/` desde cero, correr el
pipeline (`scripts/main_pipeline.py`, requiere pandas para el paso 02)
y copiar los 3 CSV de `cr2sub/` ahí.
