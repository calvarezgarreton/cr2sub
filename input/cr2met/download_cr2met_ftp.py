#!/usr/bin/env python3
"""Descarga los NetCDF mensuales diarios de CR2MET."""

from html.parser import HTMLParser
import argparse
import os
from pathlib import Path
import re
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, urlopen


URL_BASE = "https://ftp.cr2.cl"
PRODUCTOS = ("pr", "et0", "tmean", "sf")
RAIZ_PROYECTO = Path(__file__).resolve().parents[2]


def carpeta_cr2met():
    """Lee de .env la carpeta local para los archivos CR2MET."""
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
            "Define CR2MET_DATA_DIR en .env para indicar dónde guardar CR2MET."
        )
    return Path(ruta).expanduser().resolve()


class LectorDeEnlaces(HTMLParser):
    """Guarda los enlaces encontrados en una página web."""

    def __init__(self):
        super().__init__()
        self.enlaces = []

    def handle_starttag(self, etiqueta, atributos):
        if etiqueta == "a":
            enlace = dict(atributos).get("href")
            if enlace:
                self.enlaces.append(enlace)


def leer_enlaces(url):
    """Abre una página del servidor y devuelve los enlaces que contiene."""
    solicitud = Request(url, headers={"User-Agent": "CR2MET downloader"})
    with urlopen(solicitud, timeout=120) as respuesta:
        pagina = respuesta.read().decode("utf-8")

    lector = LectorDeEnlaces()
    lector.feed(pagina)
    return lector.enlaces


def productos_publicados(version):
    """Devuelve los productos disponibles para una versión."""
    url = f"{URL_BASE}/browse/cr2met/{version}/"
    productos = set()

    for enlace in leer_enlaces(url):
        ruta = urlsplit(urljoin(URL_BASE, enlace)).path
        carpeta = ruta.rstrip("/").split("/")[-1]
        if carpeta in PRODUCTOS:
            productos.add(carpeta)

    return productos


def archivos_netcdf(version, producto):
    """Devuelve los nombres y direcciones de los NetCDF."""
    carpeta = f"{version}/{producto}/{version}_day"
    url = f"{URL_BASE}/browse/cr2met/{carpeta}/"
    prefijo_descarga = f"/download/cr2met/{carpeta}/"
    archivos = []

    for enlace in leer_enlaces(url):
        direccion = urljoin(URL_BASE, enlace)
        ruta = urlsplit(direccion).path
        nombre = ruta.split("/")[-1]

        # Los archivos diarios se identifican por "_day_" y terminan en ".nc".
        es_netcdf_diario = (
            ruta.startswith(prefijo_descarga)
            and "_day_" in nombre
            and nombre.endswith(".nc")
        )
        if es_netcdf_diario:
            archivos.append((nombre, direccion))

    return archivos


def mes_del_archivo(nombre):
    """Obtiene el año y mes que aparecen en el nombre del archivo."""
    fecha = re.search(r"_(\d{4})_(\d{2})_005deg\.nc$", nombre)
    if not fecha:
        raise ValueError(f"No se reconoce la fecha del archivo: {nombre}")
    return int(fecha.group(1)), int(fecha.group(2))


def descargar_archivo(url, destino):
    """Descarga un archivo, sin reemplazar uno que ya exista."""
    if destino.is_file() and destino.stat().st_size > 0:
        print(f"Ya existe, se omite: {destino}")
        return

    destino.parent.mkdir(parents=True, exist_ok=True)
    archivo_temporal = destino.with_name(destino.name + ".part")

    try:
        solicitud = Request(url, headers={"User-Agent": "CR2MET downloader"})
        with urlopen(solicitud, timeout=120) as respuesta:
            tamano_esperado = respuesta.headers.get("Content-Length")
            bytes_descargados = 0

            with archivo_temporal.open("wb") as archivo:
                while True:
                    bloque = respuesta.read(1024 * 1024)
                    if not bloque:
                        break
                    archivo.write(bloque)
                    bytes_descargados += len(bloque)

        if bytes_descargados == 0:
            raise OSError("El servidor entregó un archivo vacío.")
        if tamano_esperado and bytes_descargados != int(tamano_esperado):
            raise OSError("La descarga quedó incompleta.")

        archivo_temporal.replace(destino)
        print(f"Descargado: {destino}")
    finally:
        archivo_temporal.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="Descarga los NetCDF diarios de CR2MET."
    )
    parser.add_argument(
        "--version",
        choices=("v2.5", "v2.5_R1"),
        help="Versión que quieres descargar. Si no se indica, descarga ambas.",
    )
    parser.add_argument(
        "--producto",
        choices=PRODUCTOS,
        help="Variable que quieres descargar. Si no se indica, descarga todas.",
    )
    opciones = parser.parse_args()

    try:
        carpeta_base = carpeta_cr2met()
    except (OSError, RuntimeError) as error:
        print(f"Error de configuración: {error}", file=sys.stderr)
        return 1

    versiones = (opciones.version,) if opciones.version else ("v2.5", "v2.5_R1")
    productos = (opciones.producto,) if opciones.producto else PRODUCTOS
    errores = []

    for version in versiones:
        try:
            publicados = productos_publicados(version)
        except (HTTPError, URLError, OSError) as error:
            print(f"Error al consultar {version}: {error}", file=sys.stderr)
            errores.append(version)
            continue

        for producto in productos:
            if producto not in publicados:
                print(f"{version}/{producto} no está publicado; se omite.")
                continue

            try:
                archivos = archivos_netcdf(version, producto)
                if not archivos:
                    print(f"No hay archivos NetCDF para {version}/{producto}; se omite.")
                    continue

                if version == "v2.5_R1":
                    archivos_v25 = archivos_netcdf("v2.5", producto)
                    if not archivos_v25:
                        raise ValueError(
                            f"No se encontró la versión v2.5 de {producto} "
                            "para calcular desde qué mes descargar R1."
                        )
                    ultimo_mes_v25 = max(
                        mes_del_archivo(nombre) for nombre, _ in archivos_v25
                    )
                    archivos = [
                        (nombre, url)
                        for nombre, url in archivos
                        if mes_del_archivo(nombre) > ultimo_mes_v25
                    ]
                    print(
                        f"{version}/{producto}: se descargan los meses posteriores "
                        f"a {ultimo_mes_v25[0]}-{ultimo_mes_v25[1]:02d}."
                    )

                for nombre, url in archivos:
                    destino = (
                        carpeta_base / version / producto / f"{version}_day" / nombre
                    )
                    try:
                        descargar_archivo(url, destino)
                    except (HTTPError, URLError, OSError, ValueError) as error:
                        print(f"Error al descargar {nombre}: {error}", file=sys.stderr)
                        errores.append(f"{version}/{producto}/{nombre}")
            except HTTPError as error:
                if error.code == 404:
                    print(f"{version}/{producto} no está publicado; se omite.")
                else:
                    print(f"Error en {version}/{producto}: {error}", file=sys.stderr)
                    errores.append(f"{version}/{producto}")
            except (URLError, OSError, ValueError) as error:
                print(f"Error en {version}/{producto}: {error}", file=sys.stderr)
                errores.append(f"{version}/{producto}")

    if errores:
        print(f"\nHubo errores en {len(errores)} carpeta(s).", file=sys.stderr)
        return 1

    print("\nDescarga completada.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
