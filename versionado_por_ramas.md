# Versionado y releases de CR2SUB

## Convención

- `main` es la rama canónica y contiene la versión vigente del producto.
- Cada release se identifica con un tag anotado e inmutable, como `v1.2.0`.
- Las ramas de trabajo son temporales y se integran a `main` mediante pull request.
- Las ramas históricas `cr2sub_v1.1` y `cr2sub_v1.2` se conservan como referencia durante la transición; no son el lugar para desarrollar nuevas versiones.
- Los archivos de resultados cuyo nombre incluye la versión (por ejemplo, `cr2sub_v1.2_attributes.csv`) mantienen esa versión aunque estén en `main`.

## Flujo de trabajo

1. Actualizar `main` y crear una rama de trabajo: `git switch -c feature/descripcion main`.
2. Implementar y validar los cambios en esa rama.
3. Abrir un pull request hacia `main`; revisar los cambios y las validaciones antes de integrarlo.
4. Para una nueva release, actualizar los metadatos y nombres de salida de la versión, validar los archivos publicados y crear un tag anotado en el commit final de `main`, por ejemplo: `git tag -a v1.3.0 -m "CR2SUB v1.3.0"`.
5. Publicar el tag con `git push origin v1.3.0`. Los tags de release no se mueven ni se reutilizan.

## Elegir el tipo de incremento

- **Patch** (`v1.2.1`): corrección sin cambios incompatibles en los datos o su estructura.
- **Minor** (`v1.3.0`): actualización compatible de datos, atributos o funcionalidad.
- **Major** (`v2.0.0`): cambio incompatible en el esquema o interpretación del producto.

## Datos pesados

Los archivos de entrada y sus rutas se documentan en el README principal. Los datos grandes que no se versionan se configuran localmente mediante `.env`; nunca se incluyen rutas específicas de una máquina en Git. Los resultados pequeños y necesarios para reproducir o consumir una release pueden versionarse si su licencia y tamaño lo permiten.

## Referencias históricas

Los tags preservan snapshots exactos de cada release. Las ramas antiguas pueden mantenerse temporalmente para facilitar la comparación o reproducir resultados, pero las correcciones futuras deben hacerse desde `main` y publicarse mediante tags nuevos, no modificando los tags anteriores.
