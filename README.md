# STAC2Cube

STAC2Cube is a Python library for automated construction of multitemporal geospatial data cubes from collections available through INPE's SpatioTemporal Asset Catalog (STAC) service.

The library integrates STAC querying, spatial filtering, variable classification, spatial harmonization, overlap aggregation, provenance recording, visualization, and export in an xarray-based multidimensional structure.

## Features

- Query multitemporal collections from INPE's STAC service
- Filter scenes according to a minimum spatial coverage threshold
- Automatically classify continuous, categorical, and count variables
- Automatically select compatible resampling and overlap-aggregation methods
- Harmonize raster assets to a common projected grid
- Preserve source-item and processing provenance
- Organize results as an xarray `DataArray`
- Visualize single-band images and three-band compositions
- Export cubes to GeoTIFF and NetCDF

## Current scope

STAC2Cube is designed for multitemporal STAC collections whose selected assets are directly accessible as raster images readable through Rasterio/GDAL, such as GeoTIFF and JPEG 2000.

The current implementation supports continuous, categorical, and count variables. Collections whose data are available only as packaged archives or products without a temporal dimension may require additional preprocessing or a different workflow.

## Installation

STAC2Cube is currently available directly from the GitHub repository.

Clone the repository:

```bash
git clone https://github.com/matheus-m-rodrigues/STAC2Cube.git
cd STAC2Cube
```

Install the package:

```bash
pip install .
```

For development, an editable installation can be used:

```bash
pip install -e .
```

PyPI installation will be made available in a future release.

## Basic usage

```python
import stac2cube
```

Collections and assets are defined through a configuration dictionary. Logical band names used in the output cube are mapped to their corresponding STAC asset identifiers.

Example:

```python
settings = {
    "S2-16D-2": {
        "bands": {
            "red": "B04",
            "green": "B03",
            "blue": "B02",
            "nir": "B08"
        }
    }
}
```

A cube can then be generated with:

```python
cube = stac2cube.cube(
    bounds=[
        -46.10365788688668,
        -23.384428852848448,
        -45.67793767204293,
        -23.062614816171813
    ],
    start_date="2025-01-01",
    end_date="2025-02-01",
    treatment="mean",
    settings=settings,
    min_coverage=0.10
)
```

The resulting object is an xarray `DataArray` organized by:

```text
collection × time × band × y × x
```

## Sentinel-1 RTC example

STAC2Cube can also process multitemporal SAR raster collections such as Sentinel-1 RTC.

```python
settings_s1 = {
    "sentinel-1-rtc-1": {
        "bands": {
            "vh": "Gamma0_VH",
            "vv": "Gamma0_VV"
        }
    }
}

cube_s1 = stac2cube.cube(
    bounds=[
        -46.10365788688668,
        -23.384428852848448,
        -45.67793767204293,
        -23.062614816171813
    ],
    start_date="2025-01-01",
    end_date="2025-03-01",
    treatment="mean",
    settings=settings_s1,
    target_crs="EPSG:31983",
    target_spatialres=10,
    min_coverage=0.10
)
```

## Visualization

Three-band compositions can be displayed with:

```python
stac2cube.plot_cube(
    cube,
    bands=("red", "green", "blue")
)
```

Single-band images are also supported:

```python
stac2cube.plot_cube(
    cube_s1,
    bands="vh"
)
```

## Export

Cubes can be exported to GeoTIFF:

```python
stac2cube.export_cube(
    cube,
    export_format="geotiff"
)
```

to NetCDF:

```python
stac2cube.export_cube(
    cube,
    export_format="netcdf"
)
```

or to both formats:

```python
stac2cube.export_cube(
    cube,
    export_format="both"
)
```

## Processing metadata

STAC2Cube stores processing and provenance information in the resulting cube, including:

- target coordinate reference system
- affine transformation
- target spatial resolution
- spatial bounds
- minimum coverage threshold
- number of STAC items found and used
- source item identifiers
- resampling methods applied
- overlap-aggregation methods
- processing timestamps and duration

## Repository

Source code and development versions are available at:

https://github.com/matheus-m-rodrigues/STAC2Cube
