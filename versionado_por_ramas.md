# Versionado por ramas

Estrategia para mantener varias versiones de un mismo producto en un solo repositorio: la versión se identifica por la **rama**, no por el nombre de una carpeta.

## Idea

- La estructura de carpetas es la misma en todas las versiones. Los nombres de carpetas y archivos no llevan sufijos como `_v1` o `_v2`.
- Cada versión vive en su propia rama de larga duración (por ejemplo `app_producto_v1`, `app_producto_v2`).
- `main` tiene la estructura común y la documentación general; no es una versión del producto.
- Para ver o ejecutar una versión, se cambia de rama; no se busca otra carpeta.

## Por qué

- **Rutas estables.** Scripts, Dockerfiles y documentación apuntan siempre a las mismas rutas, sin importar la versión.
- **Historial por versión.** Cada versión conserva su propio historial y se puede corregir sin tocar las demás.
- **Comparación directa.** Las diferencias entre versiones se ven con `git diff rama_v1..rama_v2`, en lugar de comparar carpetas duplicadas.
- **Sin copias.** No hay carpetas paralelas que se desincronicen.

## Estructura

```
repo/
├── app/            # misma ruta en todas las versiones
├── data_pipeline/  # scripts de preparación de datos
├── docs/
└── ...
```

| Rama | Contenido |
|---|---|
| `main` | Estructura común, documentación y utilidades compartidas |
| `app_producto_v1` | Versión 1 |
| `app_producto_v2` | Versión 2 |
| `app_producto_v2_<tema>` | Ramas de trabajo que salen de una versión |

## Cómo trabajar

1. **Cambiar de versión:** `git switch app_producto_v1`.
2. **Corregir algo en una versión:** se hace en su rama. Si aplica también a otra, se trae con `git cherry-pick <commit>` o con un merge, según el caso.
3. **Trabajo nuevo:** se crea una rama desde la versión correspondiente (`git switch -c app_producto_v2_tema app_producto_v2`) y se integra de vuelta a esa versión.
4. **Archivos que se mueven entre versiones:** se traen de una rama a otra de forma explícita, con `git checkout <rama> -- <ruta>` o con cherry-pick, y no con copias manuales.

## Cuidados

- **Datos grandes o no versionados.** Los datos que no están en git (insumos, bases, imágenes pesadas) pueden diferir entre versiones. Se debe documentar en cada rama qué datos necesita.
- **Compatibilidad de cambios compartidos.** Un cambio que se lleva de una versión a otra puede depender de datos o columnas que solo existen en la primera. Se revisa antes de integrarlo.
- **Documentar la convención.** Se deja escrito qué rama corresponde a qué versión (por ejemplo, en el README o en el documento de contexto), para que no dependa de memoria.
- **Borrar ramas con cuidado.** Las ramas de versión son la única copia de esa versión; no se eliminan ni se reescribe su historial.

## Cuándo conviene

- Las versiones comparten estructura y se parecen bastante.
- Se necesita seguir manteniendo versiones anteriores.
- Se quiere un historial claro de cómo cambió cada versión.

## Cuándo no

- Las versiones difieren tanto que casi no comparten código: conviene repositorios separados.
- Se necesita desplegar o leer varias versiones a la vez desde una misma copia de trabajo: en ese caso sirven los `git worktree`, una carpeta por rama, o carpetas con versión.
