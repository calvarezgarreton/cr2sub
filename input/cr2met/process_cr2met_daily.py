#!/usr/bin/env python3
"""Resume los archivos diarios de CR2MET en archivos mensuales y anuales."""

from datetime import datetime
import os
from pathlib import Path
import re
import sys

import netCDF4
import numpy as np
import xarray as xr


RAIZ_PROYECTO = Path(__file__).resolve().parents[2]
CARPETA_SALIDA = RAIZ_PROYECTO / "input" / "cr2met"

# Periodo inclusivo de procesamiento. Usa None para no limitar ese extremo.
FECHA_INICIO = "1960-01-01"  # Ejemplo: "1960-01-01"
FECHA_FIN = "2025-12-31"  # Ejemplo: "2026-07-31"

# Nombre de la variable en el archivo: tipo de agregación
PRODUCTOS = {
    "pr": "sum",
    "sf": "sum",
    "et0": "sum",
    "tmean": "mean",
}

PATRON_FECHA = re.compile(r"_(\d{4})_(\d{2})_005deg\.nc$")

# Hallazgo: et0 y tmean R1 de mayo de 2026 tienen variable 2D y tiempo único
# (2026-05-01); la redescarga fue idéntica. Se agrega el eje time al leerlos.


def carpeta_cr2met():
    """Lee de .env la carpeta local de los archivos CR2MET."""
    archivo_env = RAIZ_PROYECTO / ".env"
    if archivo_env.is_file():
        for linea in archivo_env.read_text(encoding="utf-8").splitlines():
            linea = linea.strip()
            if not linea or linea.startswith("#") or "=" not in linea:
                continue
            clave, valor = linea.split("=", 1)
            os.environ.setdefault(clave.strip(), valor.strip().strip("\"'"))

    ruta = os.environ.get("CR2MET_DATA_DIR")
    if not ruta:
        raise RuntimeError(
            "Define CR2MET_DATA_DIR en .env para ubicar los archivos CR2MET."
        )
    return Path(ruta).expanduser().resolve()


def mes_del_archivo(archivo):
    """Lee el año y el mes desde el nombre del archivo."""
    fecha = PATRON_FECHA.search(archivo.name)
    if fecha is None:
        raise ValueError(f"No se reconoce la fecha en {archivo.name}")
    return int(fecha.group(1)), int(fecha.group(2))


def archivos_de_mes(version, producto):
    """Busca los NetCDF mensuales que contienen datos diarios."""
    carpeta = carpeta_cr2met() / version / producto / f"{version}_day"
    return sorted(carpeta.glob("*_day_*_005deg.nc"), key=mes_del_archivo)


def meses_esperados(primer_mes, ultimo_mes):
    """Crea la lista de meses entre dos fechas, incluyendo ambas."""
    meses = []
    anio, mes = primer_mes
    while (anio, mes) <= ultimo_mes:
        meses.append((anio, mes))
        if mes == 12:
            anio += 1
            mes = 1
        else:
            mes += 1
    return meses


def leer_periodo():
    """Convierte las fechas configuradas al formato usado por Python."""
    try:
        inicio = (
            datetime.strptime(FECHA_INICIO, "%Y-%m-%d").date()
            if FECHA_INICIO
            else None
        )
        fin = (
            datetime.strptime(FECHA_FIN, "%Y-%m-%d").date()
            if FECHA_FIN
            else None
        )
    except ValueError as error:
        raise ValueError("Las fechas deben tener el formato AAAA-MM-DD.") from error

    if inicio and fin and inicio > fin:
        raise ValueError("La fecha de inicio debe ser anterior o igual a la fecha fin.")
    return inicio, fin


def etiqueta_anio(fecha):
    """Devuelve el año de una fecha para incluirlo en el nombre del archivo."""
    if fecha is None:
        return "sin_limite"
    return datetime.strptime(fecha, "%Y-%m-%d").strftime("%Y")


def archivo_en_periodo(archivo, inicio, fin):
    """Comprueba si el mes del archivo se cruza con el periodo solicitado."""
    anio, mes = mes_del_archivo(archivo)
    if inicio and (anio, mes) < (inicio.year, inicio.month):
        return False
    if fin and (anio, mes) > (fin.year, fin.month):
        return False
    return True


def crear_archivo_salida(ruta, nombre_variable, muestra, agregacion, atributos):
    """Crea un NetCDF vacío para ir escribiendo un resultado por periodo."""
    archivo = netCDF4.Dataset(ruta, "w", format="NETCDF4")
    archivo.createDimension("time", None)
    archivo.createDimension("lat", len(muestra["lat"]))
    archivo.createDimension("lon", len(muestra["lon"]))

    tiempo = archivo.createVariable("time", "f8", ("time",))
    tiempo.units = "days since 1900-01-01 00:00:00"
    tiempo.calendar = "standard"
    tiempo.long_name = "time"

    for coordenada in ("lat", "lon"):
        variable = archivo.createVariable(coordenada, "f8", (coordenada,))
        variable[:] = muestra[coordenada]
        variable.units = "degrees_north" if coordenada == "lat" else "degrees_east"

    datos = archivo.createVariable(
        nombre_variable,
        "f4",
        ("time", "lat", "lon"),
        zlib=True,
        complevel=4,
        fill_value=netCDF4.default_fillvals["f4"],
    )
    for atributo in ("long_name", "units"):
        if atributo in atributos:
            setattr(datos, atributo, atributos[atributo])
    if agregacion == "sum" and atributos.get("units") == "mm/day":
        datos.units = "mm"
    datos.cell_methods = f"time: {agregacion}"
    archivo.setncattr("aggregation", agregacion)
    archivo.setncattr("source", "CR2MET v2.5 y v2.5_R1")
    return archivo


def escribir_periodo(archivo, indice, fecha, valores):
    """Escribe una capa de datos y su fecha en el archivo NetCDF."""
    tiempo = archivo.variables["time"]
    tiempo[indice] = netCDF4.date2num(
        datetime(fecha[0], fecha[1], 1),
        units=tiempo.units,
        calendar=tiempo.calendar,
    )
    datos = archivo.variables[next(
        nombre for nombre in archivo.variables if nombre not in ("time", "lat", "lon")
    )]
    datos[indice, :, :] = np.ma.masked_invalid(valores.astype(np.float32))


def procesar_producto(
    producto, agregacion, fin_v25, archivos_v25, archivos_r1, inicio, fin
):
    """Genera un archivo mensual y otro anual para una variable."""
    archivos = list(archivos_v25)
    archivos.extend(
        archivo for archivo in archivos_r1 if mes_del_archivo(archivo) > fin_v25
    )
    archivos = [
        archivo for archivo in archivos if archivo_en_periodo(archivo, inicio, fin)
    ]
    archivos.sort(key=mes_del_archivo)

    if not archivos:
        raise ValueError("No se encontraron archivos diarios.")

    meses = [mes_del_archivo(archivo) for archivo in archivos]
    esperados = meses_esperados(meses[0], meses[-1])
    meses_encontrados = set(meses)
    faltantes = [mes for mes in esperados if mes not in meses_encontrados]
    if faltantes:
        ejemplo = ", ".join(f"{anio}-{mes:02d}" for anio, mes in faltantes[:5])
        raise ValueError(f"Faltan meses de datos: {ejemplo}")

    carpeta_salida_mensual = (
        carpeta_cr2met() / "v2.5_v2.5_R1_mon_yr"
    )
    carpeta_salida_mensual.mkdir(parents=True, exist_ok=True)
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)
    inicio_nombre = etiqueta_anio(FECHA_INICIO)
    fin_nombre = etiqueta_anio(FECHA_FIN)
    salida_mensual = carpeta_salida_mensual / (
        f"CR2MET_{producto}_v2.5_v2.5_R1_mon_{inicio_nombre}_{fin_nombre}.nc"
    )
    salida_anual = CARPETA_SALIDA / (
        f"CR2MET_{producto}_v2.5_v2.5_R1_yr_{inicio_nombre}_{fin_nombre}.nc"
    )
    temporal_mensual = salida_mensual.with_suffix(".tmp.nc")
    temporal_anual = salida_anual.with_suffix(".tmp.nc")

    archivo_mensual = None
    archivo_anual = None
    suma_anual = None
    dias_anuales = None
    anio_actual = None
    indice_mensual = 0
    indice_anual = 0
    muestra = None

    try:
        for ruta in archivos:
            anio, mes = mes_del_archivo(ruta)
            with xr.open_dataset(ruta) as datos:
                variable = datos[producto]
                if variable.dims == ("time", "lat", "lon"):
                    diario = np.asarray(variable.values, dtype=np.float64)
                elif (
                    variable.dims == ("lat", "lon")
                    and datos.sizes.get("time") == 1
                    and "time" in datos.variables
                ):
                    diario = np.asarray(variable.values, dtype=np.float64)[np.newaxis, ...]
                else:
                    raise ValueError(
                        f"{ruta.name}: dimensiones inesperadas para {producto}: "
                        f"{variable.dims}."
                    )

                latitud = np.asarray(datos["lat"].values)
                longitud = np.asarray(datos["lon"].values)
                atributos = dict(variable.attrs)
                fechas = np.asarray(datos["time"].values).astype("datetime64[D]")

            if len(fechas) != diario.shape[0]:
                raise ValueError(
                    f"{ruta.name}: hay {len(fechas)} fechas, pero "
                    f"{diario.shape[0]} capas de datos."
                )

            if inicio:
                fechas_inicio = np.datetime64(inicio.isoformat(), "D")
                dias_seleccionados = fechas >= fechas_inicio
            else:
                dias_seleccionados = np.ones(len(fechas), dtype=bool)
            if fin:
                fechas_fin = np.datetime64(fin.isoformat(), "D")
                dias_seleccionados &= fechas <= fechas_fin

            if not np.any(dias_seleccionados):
                continue
            diario = diario[dias_seleccionados]

            if archivo_mensual is None:
                muestra = {"lat": latitud, "lon": longitud}
                archivo_mensual = crear_archivo_salida(
                    temporal_mensual, producto, muestra, agregacion, atributos
                )
                archivo_anual = crear_archivo_salida(
                    temporal_anual, producto, muestra, agregacion, atributos
                )
            elif not (
                np.array_equal(latitud, muestra["lat"])
                and np.array_equal(longitud, muestra["lon"])
            ):
                raise ValueError(
                    f"{ruta.name}: la grilla lat/lon no coincide con los "
                    "archivos anteriores."
                )

            dias_validos = np.isfinite(diario)
            cantidad_dias = dias_validos.sum(axis=0)
            suma_dias = np.where(dias_validos, diario, 0).sum(axis=0)

            if agregacion == "sum":
                valor_mensual = np.where(cantidad_dias > 0, suma_dias, np.nan)
            else:
                valor_mensual = np.divide(
                    suma_dias,
                    cantidad_dias,
                    out=np.full_like(suma_dias, np.nan),
                    where=cantidad_dias > 0,
                )

            escribir_periodo(archivo_mensual, indice_mensual, (anio, mes), valor_mensual)
            indice_mensual += 1

            if anio_actual is not None and anio != anio_actual:
                valor_anual = calcular_valor(suma_anual, dias_anuales, agregacion)
                escribir_periodo(archivo_anual, indice_anual, (anio_actual, 1), valor_anual)
                indice_anual += 1
                suma_anual = None
                dias_anuales = None

            if suma_anual is None:
                suma_anual = np.zeros_like(suma_dias)
                dias_anuales = np.zeros_like(cantidad_dias)

            suma_anual += suma_dias
            dias_anuales += cantidad_dias
            anio_actual = anio
            print(f"Procesado: {producto} {anio}-{mes:02d}")

        valor_anual = calcular_valor(suma_anual, dias_anuales, agregacion)
        escribir_periodo(archivo_anual, indice_anual, (anio_actual, 1), valor_anual)
        archivo_mensual.close()
        archivo_anual.close()
        archivo_mensual = None
        archivo_anual = None
        temporal_mensual.replace(salida_mensual)
        temporal_anual.replace(salida_anual)
        print(f"Guardado: {salida_mensual}")
        print(f"Guardado: {salida_anual}")
    finally:
        if archivo_mensual is not None:
            archivo_mensual.close()
        if archivo_anual is not None:
            archivo_anual.close()
        temporal_mensual.unlink(missing_ok=True)
        temporal_anual.unlink(missing_ok=True)


def calcular_valor(suma, cantidad, agregacion):
    """Calcula el total o el promedio usando los días disponibles."""
    if agregacion == "sum":
        return np.where(cantidad > 0, suma, np.nan)
    return np.divide(
        suma,
        cantidad,
        out=np.full_like(suma, np.nan),
        where=cantidad > 0,
    )


def main():
    try:
        carpeta_salida_mensual = (
            carpeta_cr2met() / "v2.5_v2.5_R1_mon_yr"
        )
    except (OSError, RuntimeError) as error:
        print(f"Error de configuración: {error}", file=sys.stderr)
        return 1

    try:
        inicio, fin = leer_periodo()
    except ValueError as error:
        print(f"Periodo no válido: {error}", file=sys.stderr)
        return 1

    print(
        "Periodo solicitado: "
        f"{inicio.isoformat() if inicio else 'inicio disponible'} a "
        f"{fin.isoformat() if fin else 'fin disponible'}."
    )

    # La última fecha se busca en todos los archivos v2.5 disponibles.
    archivos_v25_por_producto = {
        producto: archivos_de_mes("v2.5", producto) for producto in PRODUCTOS
    }
    todos_los_archivos_v25 = [
        archivo
        for archivos in archivos_v25_por_producto.values()
        for archivo in archivos
    ]
    if not todos_los_archivos_v25:
        print("No se encontraron archivos diarios v2.5.", file=sys.stderr)
        return 1

    fin_v25 = max(mes_del_archivo(archivo) for archivo in todos_los_archivos_v25)
    print(f"Los datos v2.5 terminan en {fin_v25[0]}-{fin_v25[1]:02d}.")

    errores = []
    for producto, agregacion in PRODUCTOS.items():
        try:
            archivos_v25 = archivos_v25_por_producto[producto]
            if not archivos_v25:
                print(f"No se encontraron archivos v2.5 para {producto}; se omite.")
                continue
            archivos_r1 = archivos_de_mes("v2.5_R1", producto)
            procesar_producto(
                producto,
                agregacion,
                fin_v25,
                archivos_v25,
                archivos_r1,
                inicio,
                fin,
            )
        except (OSError, ValueError, KeyError) as error:
            print(f"Error en {producto}: {error}", file=sys.stderr)
            errores.append(producto)

    if errores:
        print(
            "\nNo se pudieron procesar estas variables: " + ", ".join(errores),
            file=sys.stderr,
        )
        return 1

    print(
        "\nProceso terminado. Mensuales: "
        f"{carpeta_salida_mensual}; anuales: {CARPETA_SALIDA}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
