import pystac_client
import matplotlib.pyplot as plt
import rasterio
import numpy
import xarray as xr
import math
import os
import json
from pyproj import Transformer
from pyproj.crs import CRS
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from datetime import datetime, timezone
from affine import Affine
from shapely.geometry import box, shape
from pyproj import Geod
import warnings

def geometry_area(geometry, geod):
    if geometry.is_empty:
        return 0.0
    area, _ = geod.geometry_area_perimeter(geometry)
    return abs(area)
    
def normalize(array, brightness_factor=1):
    if not numpy.isfinite(array).any():
        return numpy.zeros_like(array, dtype='float32')
    array_min = numpy.nanpercentile(array,2)
    array_max = numpy.nanpercentile(array,98)
    if array_max == array_min:
        return numpy.zeros_like(array, dtype='float32')
    normalized = ((array-array_min)/(array_max-array_min))
    return numpy.clip(normalized * brightness_factor, 0,1)

def plot_cube(
    cube,
    bands=('red', 'green', 'blue'),
    collection=None,
    time=None,
    brightness_factor=1
):
    if isinstance(bands, str):
        bands = (bands,)

    if len(bands) not in (1, 3):
        raise ValueError(
            'Visualization requires either one band or exactly three bands'
        )

    available_bands = set(cube.coords['band'].values)

    missing_bands = [
        band for band in bands
        if band not in available_bands]

    if missing_bands:
        raise ValueError(
            f'Bands {missing_bands} were not found in the cube'
        )

    if collection is None:
        collections = cube.coords['collection'].values
    else:
        collections = [collection]

    if time is None:
        times = cube.coords['time'].values
    else:
        times = [time]

    spatial_bounds = cube.attrs['spatial_bounds']
    x_min, y_min, x_max, y_max = spatial_bounds

    for collection_name in collections:
        for date in times:
            if len(bands) == 1:
                cube_plot = cube.sel(collection=collection_name,time=date,band=bands[0])

                if numpy.isnan(cube_plot.values).all():
                    continue

                image = normalize(cube_plot.values,brightness_factor)

            else:
                cube_plot = cube.sel(
                    collection=collection_name,
                    time=date,
                    band=list(bands))

                if numpy.isnan(cube_plot.values).all():
                    continue

                cube_rgb = cube_plot.transpose('y', 'x', 'band').values

                image = normalize(cube_rgb,brightness_factor)
                
            map_width = x_max - x_min
            map_height = y_max - y_min
            proportion = map_width / map_height
            base = 8

            if proportion >= 1:
                figsize = (base, base / proportion)
            else:
                figsize = (base * proportion, base)
                
            fig, ax = plt.subplots(figsize=figsize)
            if len(bands) == 1:
                ax.imshow(image,extent=[x_min, x_max, y_min, y_max],origin='upper',cmap='gray')
                ax.set_title(
                    f'{collection_name}, {date}, {bands[0]}')
            else:
                ax.imshow(image,extent=[x_min, x_max, y_min, y_max],origin='upper')

                ax.set_title(
                    f'{collection_name}, {date}')

            ax.set_aspect('equal')
            ax.ticklabel_format(
                style='plain',
                axis='both',
                useOffset=False)

            ax.set_xlabel('Easting (m)')
            ax.set_ylabel('Northing (m)')

            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.show()
            plt.close(fig)

def variable_type(asset, asset_key=None):

    band_info = asset.extra_fields.get('eo:bands',[])

    description = ''
    common_name = ''
    name = ''
    center_wavelength = 0
    full_width_half_max = 0

    if len(band_info) > 0:
        band_meta = band_info[0]
        description = str(band_meta.get('description') or '').upper()
        common_name = str(band_meta.get('common_name') or '').upper()
        name = str(band_meta.get('name') or '').upper()
        center_wavelength = (band_meta.get('center_wavelength') or 0)
        full_width_half_max = (band_meta.get('full_width_half_max') or 0)
    asset_key_text = str(asset_key or '').upper()
    title = str(getattr(asset, 'title', None) or '').upper()

    metadata_text = (
        f'{description} '
        f'{common_name} '
        f'{name} '
        f'{asset_key_text} '
        f'{title}')

    categorical_terms = ('CLASSIFICATION','QUALITY','MASK','SCL','QA_PIXEL','ST_QA','CMASK')
    count_terms = ('COUNT','CLEAROB','TOTALOB')
    continuous_terms = ('GAMMA0','SIGMA0','BETA0','BACKSCATTER') 

    if any(term in metadata_text for term in categorical_terms):
        return 'categorical'

    if any(term in metadata_text for term in count_terms):
        return 'count'
        
    if any(term in metadata_text for term in continuous_terms):
        return 'continuous'

    if (center_wavelength > 0 or full_width_half_max > 0):
        return 'continuous'

    continuous_assets = {'COASTAL','BLUE','GREEN','RED','REDEDGE','REDEDGE1','REDEDGE2','REDEDGE3','NIR','NIR08','NIR09','SWIR','SWIR16','SWIR22','LWIR','LWIR11','LWIR12','PAN','AOT','WVP'}
    sentinel_spectral_assets = {'B01','B02','B03','B04','B05','B06','B07','B08','B8A','B09','B10','B11','B12'}
    if (asset_key_text in continuous_assets or asset_key_text in sentinel_spectral_assets):
        return 'continuous'

    return 'unknown'

def validate_bounds(bounds):
    if not isinstance(bounds, (list, tuple)):
        raise ValueError('bounds must be a bounding box or a sequence of coordinate pairs')
    if (len(bounds) == 4 and all(isinstance(value, (int, float, numpy.number)) for value in bounds)):
        xmin, ymin, xmax, ymax = bounds
    elif (len(bounds) > 0 and all(isinstance(point, (list, tuple)) and len(point) == 2 for point in bounds)):
        for point in bounds:
            x, y = point
            if not isinstance(x, (int, float, numpy.number)):
                raise ValueError('All longitude values must be numeric')
            if not isinstance(y, (int, float, numpy.number)):
                raise ValueError('All latitude values must be numeric')
        x_coords = [point[0] for point in bounds]
        y_coords = [point[1] for point in bounds]
        xmin = min(x_coords)
        ymin = min(y_coords)
        xmax = max(x_coords)
        ymax = max(y_coords)

    else:
        raise ValueError('bounds must be [xmin, ymin, xmax, ymax] or a sequence of (longitude, latitude) pairs')

    if not numpy.isfinite([xmin, ymin, xmax, ymax]).all():
        raise ValueError('bounds coordinates must be finite numbers')

    if xmin < -180 or xmax > 180:
        raise ValueError('Longitude values must be between -180 and 180 degrees')

    if ymin < -90 or ymax > 90:
        raise ValueError(
            'Latitude values must be between -90 and 90 degrees')

    if xmin >= xmax:
        raise ValueError('xmin must be smaller than xmax')

    if ymin >= ymax:
        raise ValueError('ymin must be smaller than ymax')

    return [float(xmin), float(ymin), float(xmax), float(ymax)]
    
def choose_resampling(asset, mode='auto', asset_key=None):
    mode = mode.strip().lower()
    if mode == 'nearest':
        return Resampling.nearest
    elif mode == 'bilinear':
        return Resampling.bilinear
    elif mode == 'auto':
        m_type = variable_type(asset,asset_key=asset_key)
        if m_type in ('categorical', 'count'):
            return Resampling.nearest
        elif m_type == 'continuous':
            return Resampling.bilinear
        else:
            return Resampling.nearest
    else:
        raise ValueError('Unsupported resampling format')
        
def choose_aggregation(asset,treatment,asset_key=None):
    a_type = variable_type(asset,asset_key=asset_key)
    treatment = treatment.strip().lower()
    if a_type == 'continuous':
        return treatment
    elif a_type == 'categorical':
        return 'mode'
    elif a_type == 'count':
        return 'max'
    else:
        raise ValueError(
            f'Unable to determine variable type '
            f'for asset {asset_key}')
  
def cube(bounds, start_date, end_date, treatment, settings, target_crs=None, target_spatialres=None, resampling='auto', min_coverage=0.0):
    if not isinstance(min_coverage, (int,float)):
        raise ValueError('min_coverage must be a number between 0 and 1')
    if not 0 <= min_coverage <= 1:
        raise ValueError('min_coverage must be between 0 and 1')
    treatment = treatment.strip().lower()
    valid_treatments = {'mean', 'median', 'min', 'max'}
    resampling = resampling.strip().lower()
    valid_resampling = {'auto', 'nearest', 'bilinear'}
    if resampling not in valid_resampling:
        raise ValueError(f'Unsupported resampling {resampling}.'
                         f'Choose one of: {sorted(valid_resampling)}')
    if treatment not in valid_treatments:
        raise ValueError(f'Unsupported treatment {treatment}.'
                         f'Choose one of: {sorted(valid_treatments)}')
    try: 
        start_date_dt = datetime.strptime(start_date, '%Y-%m-%d')
        end_date_dt = datetime.strptime(end_date, '%Y-%m-%d')
    except(TypeError, ValueError):
        raise ValueError('start_date and end_date must use the YYYY-MM-DD format')
    if start_date_dt > end_date_dt:
        raise ValueError('start_date cannot be later than end_date')
    date_range = f'{start_date}/{end_date}' 


    stac_endpoint = 'https://data.inpe.br/bdc/stac/v1/'
    service = pystac_client.Client.open(stac_endpoint) 
    processing_start_at = datetime.now(timezone.utc)
    id_collection_list = [c.id for c in service.get_collections()] 
    selected_collections = list(settings.keys())
    for collection in selected_collections:
        if collection not in id_collection_list:
            raise ValueError(f'Collection {collection} was not found in the STAC service')

    for collection in selected_collections:
        sett_collection = settings[collection]
        if not isinstance(sett_collection, dict):
            raise ValueError(f'Configuration for {collection} must be a dictionary')
        if 'bands' not in sett_collection:
            raise ValueError(f"Collection {collection} must contain a 'bands' mapping")
        if not isinstance(sett_collection['bands'], dict):
            raise ValueError(f"'bands' for collection {collection} must be a dictionary")
        if len(sett_collection['bands']) == 0:
            raise ValueError(f"The collection not have any band.")

    bounds = validate_bounds(bounds)


    
    item_search = service.search(bbox=bounds, datetime=date_range, collections=selected_collections) 
    items = list(item_search.items())
    number_items_found = len(items)
    if len(items) == 0:
        raise ValueError('No images were found for the selected collection, spatial extent, and date range')
    geod = Geod(ellps='WGS84')
    aoi_geometry = box(*bounds)
    aoi_area = geometry_area(aoi_geometry, geod)
    filtered_items = []
    coverage_by_item = {}
    for item in items:
        if item.geometry is None:
            continue
        item_geometry = shape(item.geometry)
        intersection = aoi_geometry.intersection(item_geometry)    
        intersection_area = geometry_area(intersection, geod)
        coverage = intersection_area/aoi_area
        coverage_by_item[item.id] = coverage
        if coverage >= min_coverage:
            filtered_items.append(item)
    items = filtered_items
    number_items_used = len(items)
    if len(items) == 0:
        raise ValueError('No images satisfy the minimum coverage requirement')
    for item in items:
        col_name = item.collection_id
        config_bands = settings[col_name]['bands']
        expected_assets = set(config_bands.values())
        available_asset = set(item.assets.keys())
        absent_assets = (expected_assets - available_asset)
        if absent_assets:
            raise ValueError(f'Assets {sorted(absent_assets)} were not found'
                             f'in item {item.id} from collection {col_name}')
    date = [i.datetime for i in items] 
    image_dates = [] 
    for i in date: 
        dat_string = str(i) 
        dat_string = dat_string[0:10] 
        image_dates.append(dat_string) 
    image_dates = sorted(list(set(image_dates))) 
    print(image_dates) 
    
    print(f'Available images: ({len(items)}): ', *items, sep='\n')
    projections_found = set() 
    resolutions_found = set() 
    
    for i in items: 
        name = i.collection_id 
        bands = settings[name]['bands'] 
        epsg = i.properties.get('proj:epsg') 
        if epsg: 
            projections_found.add(epsg)
        with rasterio.open(i.assets[list(bands.values())[0]].href) as src:
            resolucao = src.res[0]
            resolutions_found.add(resolucao)
            if not epsg:
                epsg_code = src.crs.to_epsg()
                if epsg_code:
                    projections_found.add(epsg_code)
                else:
                    projections_found.add(src.crs.to_string())
                    
    projections_found.discard(None) 
    if target_crs is None: 
        if len(projections_found) == 1: 
            target_crs = list(projections_found)[0]
            print()
            print(f'Target CRS: {target_crs}\n') 
        else: 
            raise ValueError(f'Multiple CRS were found: {projections_found}. Set target_crs explicitly')
            print()
            print(f'Target CRS: {target_crs}\n') 
        
    if target_spatialres is None: 
        if len(resolutions_found) == 1:
            target_spatialres = list(resolutions_found)[0]
            print()
            print(f'The size of pixel in the images is {target_spatialres}\n')
        else: 
            raise ValueError(f'Multiple spatial resolutions were found: {resolutions_found}. Set target_spatialres explicity')
    if not isinstance(target_spatialres, (int, float)) or target_spatialres <= 0:
        raise ValueError('target_spatialres must be a positive number')
        print(f'Target spatial resolution: {target_spatialres} m\n') 
            
    
    in_proj = CRS.from_epsg(4326) 
    out_proj = CRS.from_user_input(target_crs)
    if not out_proj.is_projected:
        raise ValueError('Target CRS must be a projected CRS.')
    for axis in out_proj.axis_info:
        if not math.isclose(axis.unit_conversion_factor, 1.0):
            raise ValueError('Target CRS must use metres as linear units.')

    
    tTransform = Transformer.from_crs(in_proj, out_proj, always_xy = True) 
    x1,y1,x2,y2 = bounds 
    x1_reproj, y1_reproj = tTransform.transform(x1,y1) 
    x2_reproj, y2_reproj = tTransform.transform(x2,y2) 

    x1_aligned = math.floor(x1_reproj/target_spatialres)*target_spatialres
    x2_aligned = math.ceil(x2_reproj/target_spatialres)*target_spatialres
    y1_aligned = math.floor(y1_reproj/target_spatialres)*target_spatialres
    y2_aligned = math.ceil(y2_reproj/target_spatialres)*target_spatialres
    
    dist_x = x2_aligned - x1_aligned 
    dist_y = y2_aligned - y1_aligned 

    width = math.ceil(dist_x/target_spatialres) 
    height = math.ceil(dist_y/target_spatialres) 
    cube_transform = from_origin(x1_aligned, y2_aligned, target_spatialres, target_spatialres)
    print()
    print(f'The coodinates in WGS84 are: {x1},{y1},{x2},{y2}\n') 
    print(f'The coordinates in {target_crs} are: {x1_aligned},{y1_aligned},{x2_aligned},{y2_aligned}\n') 
    print(f'Cube dimensions: {width} x {height} pixels\n') 

    columns = numpy.arange(width)
    lines = numpy.arange(height)
    xAxis, y_0 = rasterio.transform.xy(cube_transform, rows= numpy.zeros_like(columns),cols=columns,offset='center') 
    x_0, yAxis = rasterio.transform.xy(cube_transform, rows= lines,cols=numpy.zeros_like(lines),offset='center') 

    cube_bands_set = set()  
    for info_collection in settings.values(): 
        key_band_sensor = info_collection['bands'].keys() 
        for name in key_band_sensor: 
            cube_bands_set.add(name) 

    cube_bands = sorted(list(cube_bands_set)) 

    unique_collections = list(settings.keys()) 
    coords = {'collection': unique_collections,'time': image_dates, 'band': cube_bands, 'y': yAxis, 'x': xAxis} 

    matrix = numpy.full([len(unique_collections),len(image_dates), len(cube_bands), len(yAxis), len(xAxis)], numpy.nan, dtype='float32') 

    cube = xr.DataArray(data=matrix, coords=coords, dims=['collection','time','band','y', 'x']) 
    
    source_item = []
    for item in items: 
        source_item.append({'id':item.id, 'collection':item.collection_id, 'datetime':item.properties.get('datetime'),
                            'start_datetime':item.properties.get('start_datetime'),
                           'end_datetime':item.properties.get('end_datetime'), 'coverage':coverage_by_item[item.id]})


    counts = {}
    group_median = {}
    group_mode = {}
    resampling_applied = {}
    aggregation_applied = {}
    for item in items: 
        clean_date = str(item.datetime)[0:10] 
        col_name = item.collection_id 
        bands = settings[col_name]['bands'] 
        for x,y in bands.items():
            resampling_method = choose_resampling(item.assets[y], mode=resampling, asset_key=y)
            aggregation_method = choose_aggregation(item.assets[y], treatment, asset_key=y)
            resampling_key = f'{col_name}:{x}:{y}'
            resampling_applied[resampling_key] = resampling_method.name
            aggregation_applied[resampling_key] = aggregation_method
            key = (col_name, clean_date, x)
            if aggregation_method == 'median':
                if key not in group_median:
                    group_median[key] = []
                group_median[key].append((item.assets[y].href, resampling_method))
                continue
            elif aggregation_method == 'mean':
                if key not in counts:
                    counts[key] = numpy.zeros((height, width), dtype='uint16')
                counting = counts[key]
            elif aggregation_method == 'mode':
                if key not in group_mode:
                    group_mode[key]=[]
                group_mode[key].append((item.assets[y].href, resampling_method))
                continue
                    
            with rasterio.open(item.assets[y].href) as src: 
                with WarpedVRT(src, crs=out_proj, transform=cube_transform, width=width, height=height, resampling=resampling_method) as vrt: 
                    rRead = vrt.read(1, masked=True).astype('float32') 
                    nAnRead = rRead.filled(numpy.nan) 
                    cube_value = cube.loc[col_name,clean_date,x,:,:].values 

                    valid_cube = ~numpy.isnan(cube_value)
                    new_valid = ~numpy.isnan(nAnRead)

                    both_valid = valid_cube & new_valid
                    only_cube = valid_cube & ~new_valid
                    only_new = ~valid_cube & new_valid
                    both_nan = ~valid_cube & ~new_valid

                    fusion = numpy.full_like(nAnRead, numpy.nan, dtype='float32')
                    fusion[only_cube] = cube_value[only_cube]
                    fusion[only_new] = nAnRead[only_new]

                    if aggregation_method   == 'min': 
                        fusion[both_valid] = numpy.minimum(cube_value[both_valid], nAnRead[both_valid])
                    elif aggregation_method   == 'max':
                        fusion[both_valid] = numpy.maximum(cube_value[both_valid], nAnRead[both_valid])
                    elif aggregation_method   == 'mean':
                        counting[only_new] = 1
                        fusion[both_valid] = (cube_value[both_valid]*counting[both_valid]+nAnRead[both_valid])/(counting[both_valid]+1) 
                        counting[both_valid] += 1
                    cube.loc[col_name, clean_date, x, :, :] = fusion

    for key, hrefs in group_median.items():
        col_name, clean_date, x = key
        pile = numpy.empty((len(hrefs), height, width),dtype='float32')
        for index, (link, resampling_method) in enumerate(hrefs):
            with rasterio.open(link) as src:
                with WarpedVRT(src,crs=out_proj,transform=cube_transform,width=width,height=height,resampling=resampling_method) as vrt:
                    rRead = vrt.read(1,masked=True).astype('float32')
                    pile[index] = rRead.data
                    if numpy.ma.is_masked(rRead):
                        pile[index][rRead.mask] = numpy.nan
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=r"^All-NaN slice encountered$", category=RuntimeWarning)
            median = numpy.nanmedian(pile, axis=0, overwrite_input=True)
        cube.loc[col_name,clean_date,x,:,:] = median
        del pile, median, rRead

    for key, hrefs in group_mode.items():
        col_name, clean_date, x = key
        pile = numpy.empty((len(hrefs), height, width),dtype='float32')
        for index, (link, resampling_method) in enumerate(hrefs):
            with rasterio.open(link) as src:
                with WarpedVRT(src, crs=out_proj, transform=cube_transform, width=width, height=height, resampling=resampling_method) as vrt:
                    rRead = vrt.read(1,masked=True, out_dtype='float32')
                    pile[index] = rRead.data
                    if numpy.ma.is_masked(rRead):
                        pile[index][rRead.mask] = numpy.nan
        mode = numpy.full(pile.shape[1:],numpy.nan,dtype='float32')
        high_count = numpy.zeros(pile.shape[1:],dtype='uint16')
        for observation in pile:
            current_count = numpy.sum(pile == observation,axis=0,dtype='uint16')
            valid = ~numpy.isnan(observation)
            update = (valid & (current_count > high_count))
            mode[update] = observation[update]
            high_count[update] = (current_count[update])
        cube.loc[col_name,clean_date,x,:,:] = mode
        del (pile,mode,high_count,current_count,observation,valid,update,rRead)
        
    processing_end_at = datetime.now(timezone.utc)
    processing_time = processing_end_at - processing_start_at
    processing_time = processing_time.total_seconds()
    minutes, seconds = divmod(processing_time,60)
    processing_duration_text = f'{int(minutes)} min {seconds:.2f} s'
    epsg = out_proj.to_epsg()
    if epsg is not None:
        cube.attrs['crs'] = f'EPSG:{epsg}'
    else:
        cube.attrs['crs'] = out_proj.to_wkt()
    cube.attrs['transform'] = [cube_transform.a, cube_transform.b, cube_transform.c, cube_transform.d, cube_transform.e, cube_transform.f]
    cube.attrs['spatial_resolution'] = target_spatialres
    cube.attrs['spatial_bounds'] = [x1_aligned, y1_aligned, x2_aligned, y2_aligned]
    cube.attrs['treatment'] = treatment
    cube.attrs['resampling_requested'] = resampling
    cube.attrs['stac_endpoint'] = stac_endpoint
    cube.attrs['source_items'] = json.dumps(source_item)
    cube.attrs['resampling_applied'] = json.dumps(resampling_applied)
    cube.attrs['aggregation_applied'] = json.dumps(aggregation_applied)
    cube.attrs['processing_start_at'] = processing_start_at.isoformat()
    cube.attrs['processing_finished_at'] = processing_end_at.isoformat()
    cube.attrs['processing_duration'] = processing_duration_text
    cube.attrs['min_coverage'] = min_coverage
    cube.attrs['items_found'] = number_items_found
    cube.attrs['items_used'] = number_items_used
    
    return cube 

def export_geotiff(cube, output_folder='exported_cube'):
    os.makedirs(output_folder, exist_ok=True)
    source_items = json.loads(cube.attrs['source_items'])
    
    width = cube.sizes['x']
    height = cube.sizes['y']
    collections = cube.coords['collection'].values
    dates = cube.coords['time'].values
    bands = cube.coords['band'].values
    number_of_bands = len(bands)
    crs = cube.attrs['crs']
    transform = Affine(*cube.attrs['transform'])
    source_item = json.loads(cube.attrs['source_items'])
    resampling_applied = json.loads(cube.attrs['resampling_applied'])
    aggregation_applied = json.loads(cube.attrs['aggregation_applied'])
    for collection in collections:
        for date in dates:
            date_text = str(date)[:10]
            source_item_id = []
            for item in source_item:
                collection_item = item['collection']
                item_date = item['datetime'][0:10]
                if collection_item == str(collection) and item_date == date_text:
                    source_item_id.append(item['id'])
            data = cube.sel(collection=collection, time=date)
            if numpy.isnan(data.values).all():
                continue
            name_tiff = f'{output_folder}/{collection}_{date_text}.tif'
            profile = {'driver':'GTiff',
                       'dtype':'float32',
                       'nodata':numpy.nan,
                       'width':width,
                       'height':height,
                       'count':number_of_bands,
                       'crs':crs,
                       'transform':transform,
                       'compress':'lzw'}
            with rasterio.open(name_tiff,'w',**profile) as dst:
                dst.update_tags(collection_id=str(collection),
                                    acquisition_date=date_text,
                                    treatment=cube.attrs['treatment'],
                                    resampling_requested=cube.attrs['resampling_requested'],
                                    aggregation_applied = cube.attrs['aggregation_applied'],
                                    stac_endpoint=cube.attrs['stac_endpoint'],
                                    processing_start_at=cube.attrs['processing_start_at'],
                                    processing_end_at=cube.attrs['processing_finished_at'],
                                    processing_duration=cube.attrs['processing_duration'],
                                    source_item_ids=json.dumps(source_item_id))
                for x,y in enumerate(bands):
                    band_matrix = data.sel(band=y).values
                    layer_number = x+1
                    dst.write(band_matrix, layer_number)
                    dst.set_band_description(layer_number, str(y))
                    for key, method in resampling_applied.items():
                        collection_key, band_key, asset_key = key.split(':')
                        if collection_key == str(collection) and band_key == str(y):
                            aggregation_method = aggregation_applied[key]
                            dst.update_tags(layer_number, asset=asset_key, resampling_applied=method, aggregation_applied=aggregation_method)
                            break
                print(f'Saved successfully: {name_tiff}')
                    
def export_netcdf(cube,file_name='exported_cube/data_cube.nc'):
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    cube_export = cube.copy()
    resampling_applied = json.loads(cube.attrs['resampling_applied'])
    aggregation_applied = json.loads(cube.attrs['aggregation_applied'])
    collections = cube_export.coords['collection'].values
    bands = cube_export.coords['band'].values
    asset_matrix = numpy.empty((len(collections), len(bands)),dtype=object)
    resampling_matrix = numpy.empty((len(collections), len(bands)), dtype=object)
    aggregation_matrix = numpy.empty((len(collections), len(bands)),dtype=object)
    for key, method in resampling_applied.items():
        collection, band, asset = key.split(':')
        index_collection = numpy.where(collections==collection)[0][0]
        index_band = numpy.where(bands == band)[0][0]
        asset_matrix[index_collection, index_band]= asset
        resampling_matrix[index_collection, index_band]=method
        aggregation_matrix[index_collection, index_band] = aggregation_applied[key]
    asset_matrix = asset_matrix.astype(str)
    resampling_matrix = resampling_matrix.astype(str)
    aggregation_matrix = aggregation_matrix.astype(str)
    cube_export = cube_export.assign_coords(asset=(('collection','band'),asset_matrix),resampling_method=(('collection','band'),resampling_matrix), aggregation_method=(('collection','band'),aggregation_matrix))
        
    cube_export.name = 'data_cube'
    cube_export.to_netcdf(file_name)
    print(f'Saved successfully: {file_name}')


def export_cube(cube, export_format='both', output_dir='exported_cube', netcdf_name='data_cube.nc'):
    export_format = export_format.strip().lower()
    if export_format == 'geotiff':
        export_geotiff(cube, output_folder=output_dir)
    elif export_format == 'netcdf':
        file_name = os.path.join(output_dir, netcdf_name)
        export_netcdf(cube, file_name=file_name)
    elif export_format == 'both':
        export_geotiff(cube, output_folder=output_dir)
        file_name = os.path.join(output_dir, netcdf_name)
        export_netcdf(cube, file_name=file_name)
    else:
        raise ValueError("export_format must be 'geotiff', 'netcdf', or 'both")
