"""IPBES global coastal vulnerability calculation."""
import shutil
import time
import os
import math
import logging
import multiprocessing
import traceback
import sys

import numpy
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import rtree
import shapely
import shapely.wkb
import shapely.ops
import pygeoprocessing

import Task

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')

_TARGET_WORKSPACE = "ipbes_cv_workspace"

_GLOBAL_POLYGON_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\Land\global_polygon.shp"

_GLOBAL_WWIII_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\CoastalProtection\Input\WaveWatchIII.shp"

# The global bounding box to do the entire analysis
# This range was roughly picked to avoid the poles
# [min_lat, min_lng, max_lat, max_lng]
_GLOBAL_BOUNDING_BOX_WGS84 = [-180, -60, 180, 60]

# This is the lat/lng grid size to slice the runs into, annoying since it's
# lat/lng, but if you have a better idea lets hear it.
# The 3.0 degrees comes from the fact that UTM zones are 6 degrees wide so
# half of that plus some buffer should be good enough
_WGS84_GRID_SIZE = 3.0

# Once global grid is cut, it is reprojected into UTM with this underlying
# square grid cell size
_UTM_GRID_SIZE = 250

_GLOBAL_GRID_VECTOR_FILE_PATTERN = 'global_grid.shp'
_GLOBAL_FEATURE_INDEX_FILE_PATTERN = 'global_feature_index.dat'
_GRID_POINT_FILE_PATTERN = 'grid_points_%d.shp'
_WIND_EXPOSURE_POINT_FILE_PATTERN = 'rei_points_%d.shp'
_WORK_COMPLETE_TOKEN = 'work.complete_%s'
_WORK_COMPLETE_PATH = os.path.join(
    _TARGET_WORKSPACE, 'work_tokens')

def main():
    """Entry point."""
    work_queue = multiprocessing.JoinableQueue(multiprocessing.cpu_count())
    for _ in range(multiprocessing.cpu_count()):
        multiprocessing.Process(
            target=Task.worker, args=(work_queue,)).start()

    global_thread_map = {}
    global_task_lock = multiprocessing.Lock()

    if not os.path.exists(_TARGET_WORKSPACE):
        os.makedirs(_TARGET_WORKSPACE)
    if not os.path.exists(_WORK_COMPLETE_PATH):
        os.makedirs(_WORK_COMPLETE_PATH)

    landmass_bounding_rtree_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_FEATURE_INDEX_FILE_PATTERN)

    rtree_done_token = os.path.join(
        _WORK_COMPLETE_PATH, _WORK_COMPLETE_TOKEN % ('rtree_task'))
    build_rtree_task = Task.Task(
        _build_feature_bounding_box_rtree, (
            _GLOBAL_POLYGON_PATH, landmass_bounding_rtree_path,
            rtree_done_token),
        [rtree_done_token], [], global_thread_map, global_task_lock)

    global_grid_vector_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_GRID_VECTOR_FILE_PATTERN)
    global_grid_token = os.path.join(
        _WORK_COMPLETE_PATH, _WORK_COMPLETE_TOKEN % ('global_grid_task'))

    grid_edges_of_vector_task = Task.Task(
        _grid_edges_of_vector, (
            _GLOBAL_BOUNDING_BOX_WGS84, _GLOBAL_POLYGON_PATH,
            landmass_bounding_rtree_path, global_grid_vector_path,
            _WGS84_GRID_SIZE, global_grid_token),
        [global_grid_token], [build_rtree_task],
        global_thread_map, global_task_lock)

    grid_edges_of_vector_task()

    global_grid_vector = ogr.Open(global_grid_vector_path)
    global_grid_layer = global_grid_vector.GetLayer()
    grid_count = global_grid_layer.GetFeatureCount()

    smallest_feature_size = 2000
    max_fetch_distance = 60000

    for grid_id in [3]:
        LOGGER.info("Calculating grid %d of %d", grid_id, grid_count)
        grid_point_path = os.path.join(
            _TARGET_WORKSPACE, _GRID_POINT_FILE_PATTERN % (grid_id))
        temp_workspace = os.path.join(
            _TARGET_WORKSPACE, 'grid_%d' % grid_id)
        work_complete_token_path = os.path.join(
            _WORK_COMPLETE_PATH, _WORK_COMPLETE_TOKEN % ('grid_%d' % grid_id))
        create_shore_points_task = Task.Task(
            _create_shore_points, (
                global_grid_vector_path, grid_id, landmass_bounding_rtree_path,
                _GLOBAL_POLYGON_PATH, _GLOBAL_WWIII_PATH,
                smallest_feature_size, temp_workspace, grid_point_path,
                work_complete_token_path),
            [work_complete_token_path], [], global_thread_map,
            global_task_lock)

        temp_workspace = os.path.join(
            _TARGET_WORKSPACE, 'wind_exposure_%d' % grid_id)
        work_complete_token_path = os.path.join(
            _WORK_COMPLETE_PATH, _WORK_COMPLETE_TOKEN % (
                'wind_exposure_%d' % grid_id))
        target_wind_exposure_point_path = os.path.join(
            _TARGET_WORKSPACE, _WIND_EXPOSURE_POINT_FILE_PATTERN % grid_id)
        calculate_wind_exposure_task = Task.Task(
            _calculate_wind_exposure, (
                grid_point_path, landmass_bounding_rtree_path,
                _GLOBAL_POLYGON_PATH, temp_workspace, smallest_feature_size,
                max_fetch_distance, target_wind_exposure_point_path,
                work_complete_token_path),
            [work_complete_token_path], [create_shore_points_task],
            global_thread_map, global_task_lock)

    calculate_wind_exposure_task()

    for _ in range(multiprocessing.cpu_count()):
        work_queue.put('STOP')
    work_queue.join()


def _calculate_wind_exposure(
        base_shore_point_vector_path,
        landmass_bounding_rtree_path, landmass_vector_path, temp_workspace,
        smallest_feature_size, max_fetch_distance,
        target_fetch_point_vector_path, work_complete_token_path):
    """Calculate wind exposure for each shore point.

    Parameters:
        base_shore_point_vector_path (string): path to a point shapefile
            representing shore points that should be sampled for wind
            exposure.
        landmass_bounding_rtree_path (string): path to an rtree bounding box
            for the landmass polygons.
        landmass_vector_path (string): path to landmass polygon vetor.
        temp_workspace (string): path to a directory that can be created for
            temporary workspace files
        smallest_feature_size (float): smallest feature size to detect in
            meters.
        max_fetch_distance (float): maximum fetch distance for a ray in
            meters.
        target_fetch_point_vector_path (string): path to target point file,
            will be a copy of `base_shore_point_vector_path`'s geometry with
            an 'REI' (relative exposure index) field added.
        work_complete_token_path (string): path to file that once created
            signals that the function is complete.
    """
    if os.path.exists(temp_workspace):
        shutil.rmtree(temp_workspace)
    os.makedirs(temp_workspace)

    temp_utm_clipped_vector_path = os.path.join(
        temp_workspace, 'utm_clipped_landmass.shp')
    temp_fetch_rays_path = os.path.join(
        temp_workspace, 'fetch_rays.shp')

    # reproject base_shore_point_vector_path to utm coordinates
    base_shore_info = pygeoprocessing.get_vector_info(
        base_shore_point_vector_path)

    utm_spatial_reference = _get_utm_spatial_reference(
        base_shore_info['bounding_box'])
    base_spatial_reference = osr.SpatialReference()
    base_spatial_reference.ImportFromWkt(base_shore_info['projection'])

    pygeoprocessing.reproject_vector(
        base_shore_point_vector_path, utm_spatial_reference.ExportToWkt(),
        target_fetch_point_vector_path)

    utm_bounding_box = pygeoprocessing.get_vector_info(
        target_fetch_point_vector_path)['bounding_box']

    # extend bounding box for max fetch distance
    utm_bounding_box = [
        utm_bounding_box[0] - max_fetch_distance,
        utm_bounding_box[1] - max_fetch_distance,
        utm_bounding_box[2] + max_fetch_distance,
        utm_bounding_box[3] + max_fetch_distance]

    # get lat/lng bounding box of utm projected coordinates

    # get global polygon clip of that utm box
    # transform local box back to lat/lng -> global clipping box
    base_clipping_box = pygeoprocessing.transform_bounding_box(
        utm_bounding_box, utm_spatial_reference.ExportToWkt(),
        base_spatial_reference.ExportToWkt(), edge_samples=11)
    base_clipping_shapely = shapely.geometry.box(*base_clipping_box)

    landmass_vector_rtree = rtree.index.Index(
        os.path.splitext(landmass_bounding_rtree_path)[0])

    landmass_vector = ogr.Open(landmass_vector_path)
    landmass_layer = landmass_vector.GetLayer()

    # this will hold the clipped landmass geometry
    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
    temp_clipped_vector_path = os.path.join(
        temp_workspace, 'clipped_geometry_vector.shp')
    temp_clipped_vector = esri_shapefile_driver.CreateDataSource(
        temp_clipped_vector_path)
    temp_clipped_layer = (
        temp_clipped_vector.CreateLayer(
            os.path.splitext(temp_clipped_vector_path)[0],
            base_spatial_reference, ogr.wkbPolygon))
    temp_clipped_defn = temp_clipped_layer.GetLayerDefn()

    # clip global polygon to global clipping box
    for feature_id in landmass_vector_rtree.intersection(base_clipping_box):
        landmass_feature = landmass_layer.GetFeature(feature_id)
        landmass_shapely = shapely.wkb.loads(
            landmass_feature.GetGeometryRef().ExportToWkb())
        intersection_shapely = base_clipping_shapely.intersection(
            landmass_shapely)
        try:
            clipped_geometry = ogr.CreateGeometryFromWkt(
                intersection_shapely.wkt)
            clipped_feature = ogr.Feature(temp_clipped_defn)
            clipped_feature.SetGeometry(clipped_geometry)
            temp_clipped_layer.CreateFeature(clipped_feature)
            clipped_feature = None
        except Exception:
            LOGGER.warn(
                "Couldn't process this intersection %s", intersection_shapely)
    temp_clipped_layer.SyncToDisk()
    temp_clipped_layer = None
    temp_clipped_vector = None

    # project global clipped polygons to UTM
    LOGGER.debug("reprojecting grid %s", base_shore_point_vector_path)
    pygeoprocessing.reproject_vector(
        temp_clipped_vector_path, utm_spatial_reference.ExportToWkt(),
        temp_utm_clipped_vector_path)

    clipped_geometry_shapely_list = []
    temp_utm_clipped_vector = ogr.Open(temp_utm_clipped_vector_path)
    temp_utm_clipped_layer = temp_utm_clipped_vector.GetLayer()
    for tmp_utm_feature in temp_utm_clipped_layer:
        tmp_utm_geometry = tmp_utm_feature.GetGeometryRef()
        clipped_geometry_shapely_list.append(
            shapely.wkb.loads(tmp_utm_geometry.ExportToWkb()))
        tmp_utm_geometry = None
    temp_utm_clipped_layer = None
    temp_utm_clipped_vector = None
    landmass_shapely = shapely.ops.cascaded_union(
        clipped_geometry_shapely_list)
    clipped_geometry_shapely_list = None

    # load land geometry into shapely object
    landmass_shapely_prep = shapely.prepared.prep(landmass_shapely)

    # create fetch rays
    temp_fetch_rays_vector = esri_shapefile_driver.CreateDataSource(
        temp_fetch_rays_path)
    temp_fetch_rays_layer = (
        temp_fetch_rays_vector.CreateLayer(
            os.path.splitext(temp_clipped_vector_path)[0],
            utm_spatial_reference, ogr.wkbLineString))
    temp_fetch_rays_defn = temp_fetch_rays_layer.GetLayerDefn()
    temp_fetch_rays_layer.CreateField(ogr.FieldDefn(
        'fetch_dist', ogr.OFTReal))

    target_shore_point_vector = ogr.Open(target_fetch_point_vector_path, 1)
    target_shore_point_layer = target_shore_point_vector.GetLayer()
    target_shore_point_layer.CreateField(ogr.FieldDefn('REI', ogr.OFTReal))

    n_fetch_rays = 16
    shore_point_logger = _make_logger_callback("Shore point %.2f%% complete.")
    for shore_point_feature in target_shore_point_layer:
        shore_point_logger(
            float(shore_point_feature.GetFID()) /
            target_shore_point_layer.GetFeatureCount())
        rei_value = 0.0
        for sample_index in xrange(n_fetch_rays):
            compass_theta = float(sample_index) / n_fetch_rays * 360
            rei_pct = shore_point_feature.GetField(
                'REI_PCT%d' % int(compass_theta))
            rei_v = shore_point_feature.GetField(
                'REI_V%d' % int(compass_theta))
            cartesian_theta = -(compass_theta - 90)
            delta_x = math.cos(cartesian_theta * math.pi / 180)
            delta_y = math.sin(cartesian_theta * math.pi / 180)

            shore_point_geometry = shore_point_feature.GetGeometryRef()
            point_a_x = (
                shore_point_geometry.GetX() + delta_x * smallest_feature_size)
            point_a_y = (
                shore_point_geometry.GetY() + delta_y * smallest_feature_size)
            point_b_x = point_a_x + delta_x * (
                max_fetch_distance - smallest_feature_size)
            point_b_y = point_a_y + delta_y * (
                max_fetch_distance - smallest_feature_size)

            shore_point_geometry = None

            ray = ogr.Geometry(ogr.wkbLineString)
            ray.AddPoint(point_a_x, point_a_y)
            ray.AddPoint(point_b_x, point_b_y)

            # TEST RAY AGAINST LANDMASS FOR POSSIBLE CLIP
            ray_shapely = shapely.wkb.loads(ray.ExportToWkb())
            ray_length = 0.0
            if not landmass_shapely_prep.intersects(ray_shapely):
                ray_feature = ogr.Feature(temp_fetch_rays_defn)
                ray_feature.SetGeometry(ray)
                ray_feature.SetField(
                    'fetch_dist', ray_shapely.length)
                temp_fetch_rays_layer.CreateFeature(ray_feature)
            else:
                intersection_ray = ray_shapely.difference(landmass_shapely)
                # if there's a difference it's be a line or a multiline
                # anything else will be an empty set
                if intersection_ray.geom_type == "LineString":
                    # ensure they have the same starting points, otherwise
                    # it's probably further down the line ray
                    if intersection_ray.coords[0] == ray_shapely.coords[0]:
                        ray_feature = ogr.Feature(temp_fetch_rays_defn)
                        ray_feature.SetGeometry(
                            ogr.CreateGeometryFromWkt(
                                intersection_ray.wkt))
                        ray_length = intersection_ray.length
                        ray_feature.SetField(
                            'fetch_dist', ray_length)
                        temp_fetch_rays_layer.CreateFeature(ray_feature)
                elif intersection_ray.geom_type == "MultiLineString":
                    # the first line segment is the originating ray
                    if (intersection_ray.geoms[0].coords[0] ==
                            ray_shapely.coords[0]):
                        ray_feature = ogr.Feature(temp_fetch_rays_defn)
                        ray_feature.SetGeometry(
                            ogr.CreateGeometryFromWkt(
                                intersection_ray.geoms[0].wkt))
                        ray_length = intersection_ray.geoms[0].length
                        ray_feature.SetField('fetch_dist', ray_length)
                        temp_fetch_rays_layer.CreateFeature(ray_feature)
            ray_feature = None
            rei_value += ray_length * rei_pct * rei_v
        shore_point_feature.SetField('REI', rei_value)
        target_shore_point_layer.SetFeature(shore_point_feature)

    target_shore_point_layer.SyncToDisk()
    target_shore_point_layer = None
    target_shore_point_vector = None
    temp_fetch_rays_layer.SyncToDisk()
    temp_fetch_rays_layer = None
    temp_fetch_rays_vector = None
    with open(work_complete_token_path, 'w') as done_file:
        done_file.write('Done.\n')


def _create_shore_points(
        sample_grid_vector_path, grid_id, landmass_bounding_rtree_path,
        landmass_vector_path, wwiii_vector_path, smallest_feature_size,
        temp_workspace, target_shore_point_vector_path,
        work_complete_token_path):
    """Create points that lie on the coast line of the landmass.

    Parameters:
        sample_grid_vector_path (string): path to vector containing grids
            that are used for discrete global sampling of the landmass
            polygon.
        grid_id (integer): feature ID in `sample_grid_vector_path`'s layer to
            operate on.
        landmass_bounding_rtree_path (string): path to an rtree index that has
            bounding box indexes of the polygons in `landmass_vector_path`.
        landmass_vector_path (string): path to polygon vector representing
            landmass.
        wwiii_vector_path (string): path to point shapefile representing
            the Wave Watch III data.
        smallest_feature_size (float): smallest feature size to grid a shore
            point on.
        temp_workspace (string): path to a directory that can be created
            during run to hold temporary files.  Will be deleted on successful
            function completion.
        target_shore_point_vector_path (string): path to a point vector that
            will be created and contain points on the shore of the landmass.

    Returns:
        None.
    """
    # create the spatial reference from the base vector
    landmass_spatial_reference = osr.SpatialReference()
    landmass_spatial_reference.ImportFromWkt(
        pygeoprocessing.get_vector_info(landmass_vector_path)['projection'])

    if os.path.exists(temp_workspace):
        shutil.rmtree(temp_workspace)
    os.makedirs(temp_workspace)

    temp_clipped_vector_path = os.path.join(
        temp_workspace, 'clipped_geometry_vector.shp')
    temp_grid_raster_path = os.path.join(temp_workspace, 'grid.tif')
    temp_convolution_raster_path = os.path.join(
        temp_workspace, 'convolution.tif')
    temp_utm_clipped_vector_path = os.path.join(
        temp_workspace, 'clipped_geometry_utm.shp')
    temp_shore_kernel_path = os.path.join(
        temp_workspace, 'kernel.tif')
    temp_shore_raster_path = os.path.join(
        temp_workspace, 'shore.tif')

    for path in [target_shore_point_vector_path,
                 temp_clipped_vector_path,
                 temp_grid_raster_path]:
        if os.path.exists(path):
            os.remove(path)

    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")

    # this will hold the clipped landmass geometry
    temp_clipped_vector = esri_shapefile_driver.CreateDataSource(
        temp_clipped_vector_path)
    temp_clipped_layer = (
        temp_clipped_vector.CreateLayer(
            os.path.splitext(temp_clipped_vector_path)[0],
            landmass_spatial_reference, ogr.wkbPolygon))
    temp_clipped_defn = temp_clipped_layer.GetLayerDefn()

    # this will hold the output sample points on the shore
    target_shore_point_vector = esri_shapefile_driver.CreateDataSource(
        target_shore_point_vector_path)
    target_shore_point_layer = target_shore_point_vector.CreateLayer(
        os.path.splitext(target_shore_point_vector_path)[0],
        landmass_spatial_reference, ogr.wkbPoint)

    wwiii_vector = ogr.Open(wwiii_vector_path)
    wwiii_layer = wwiii_vector.GetLayer()
    wwiii_defn = wwiii_layer.GetLayerDefn()
    field_names = []
    for field_index in xrange(wwiii_defn.GetFieldCount()):
        field_defn = wwiii_defn.GetFieldDefn(field_index)
        field_name = field_defn.GetName()
        if field_name in ['I', 'J']:
            continue
        field_names.append(field_name)
        target_shore_point_layer.CreateField(field_defn)
    target_shore_point_defn = target_shore_point_layer.GetLayerDefn()

    landmass_vector = ogr.Open(landmass_vector_path)
    landmass_layer = landmass_vector.GetLayer()

    grid_vector = ogr.Open(sample_grid_vector_path)
    grid_layer = grid_vector.GetLayer()
    grid_feature = grid_layer.GetFeature(grid_id)
    grid_shapely = shapely.wkb.loads(
        grid_feature.GetGeometryRef().ExportToWkb())

    landmass_vector_rtree = rtree.index.Index(
        os.path.splitext(landmass_bounding_rtree_path)[0])

    # project global polygon clip to UTM
    # transform lat/lng box to utm -> local box
    utm_spatial_reference = _get_utm_spatial_reference(grid_shapely.bounds)
    utm_bounding_box = pygeoprocessing.transform_bounding_box(
        grid_shapely.bounds, landmass_spatial_reference.ExportToWkt(),
        utm_spatial_reference.ExportToWkt(), edge_samples=11)

    # transform local box back to lat/lng -> global clipping box
    utm_clipping_box = pygeoprocessing.transform_bounding_box(
        utm_bounding_box, utm_spatial_reference.ExportToWkt(),
        landmass_spatial_reference.ExportToWkt(), edge_samples=11)
    utm_clipping_shapely = shapely.geometry.box(*utm_clipping_box)

    # clip global polygon to global clipping box
    for feature_id in landmass_vector_rtree.intersection(utm_clipping_box):
        base_feature = landmass_layer.GetFeature(feature_id)
        base_geometry = base_feature.GetGeometryRef()
        base_shapely = shapely.wkb.loads(
            base_geometry.ExportToWkb())
        base_geometry = None
        intersection_shapely = utm_clipping_shapely.intersection(base_shapely)
        try:
            target_geometry = ogr.CreateGeometryFromWkt(
                intersection_shapely.wkt)
            target_feature = ogr.Feature(temp_clipped_defn)
            target_feature.SetGeometry(target_geometry)
            temp_clipped_layer.CreateFeature(target_feature)
            target_feature = None
            target_geometry = None
        except Exception:
            LOGGER.warn(
                "Couldn't process this intersection %s",
                intersection_shapely)
    temp_clipped_layer.SyncToDisk()
    temp_clipped_layer = None
    temp_clipped_vector = None

    # create grid for underlying local utm box
    pygeoprocessing.reproject_vector(
        temp_clipped_vector_path, utm_spatial_reference.ExportToWkt(),
        temp_utm_clipped_vector_path)

    byte_nodata = 255
    pygeoprocessing.create_raster_from_vector_extents(
        temp_utm_clipped_vector_path,
        temp_grid_raster_path, (
            smallest_feature_size / 2.0, -smallest_feature_size / 2.0),
        gdal.GDT_Byte, byte_nodata, fill_value=0)

    # rasterize utm global clip to grid
    pygeoprocessing.rasterize(
        temp_utm_clipped_vector_path, temp_grid_raster_path, [1], None)

    # grid shoreline from raster
    _make_shore_kernel(temp_shore_kernel_path)
    pygeoprocessing.convolve_2d(
        (temp_grid_raster_path, 1), (temp_shore_kernel_path, 1),
        temp_convolution_raster_path, target_datatype=gdal.GDT_Byte)

    temp_grid_nodata = pygeoprocessing.get_raster_info(
        temp_grid_raster_path)['nodata'][0]

    def _shore_mask_op(shore_convolution):
        """Mask values on land that border water."""
        result = numpy.empty(shore_convolution.shape, dtype=numpy.uint8)
        result[:] = byte_nodata
        valid_mask = shore_convolution != temp_grid_nodata
        # If a pixel is on land, it gets at least a 9, but if it's all on
        # land it gets an 17 (8 neighboring pixels), so we search between 9
        # and 17 to determine a shore pixel
        result[valid_mask] = numpy.where(
            (shore_convolution[valid_mask] >= 9) &
            (shore_convolution[valid_mask] < 17), 1, byte_nodata)
        return result

    pygeoprocessing.raster_calculator(
        [(temp_convolution_raster_path, 1)], _shore_mask_op,
        temp_shore_raster_path, gdal.GDT_Byte, byte_nodata)

    shore_geotransform = pygeoprocessing.get_raster_info(
        temp_shore_raster_path)['geotransform']

    utm_to_base_transform = osr.CoordinateTransformation(
        utm_spatial_reference, landmass_spatial_reference)

    wwiii_rtree = rtree.index.Index()

    for wwiii_feature in wwiii_layer:
        wwiii_geometry = wwiii_feature.GetGeometryRef()
        wwiii_x = wwiii_geometry.GetX()
        wwiii_y = wwiii_geometry.GetY()
        wwiii_rtree.insert(
            wwiii_feature.GetFID(), (wwiii_x, wwiii_y, wwiii_x, wwiii_y))

    for offset_info, data_block in pygeoprocessing.iterblocks(
            temp_shore_raster_path):
        row_indexes, col_indexes = numpy.mgrid[
            offset_info['yoff']:offset_info['yoff']+offset_info['win_ysize'],
            offset_info['xoff']:offset_info['xoff']+offset_info['win_xsize']]
        valid_mask = data_block == 1
        x_coordinates = (
            shore_geotransform[0] +
            shore_geotransform[1] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[2] * (row_indexes[valid_mask] + 0.5))
        y_coordinates = (
            shore_geotransform[3] +
            shore_geotransform[4] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[5] * (row_indexes[valid_mask] + 0.5))
        for x_coord, y_coord in zip(x_coordinates, y_coordinates):
            shore_point_geometry = ogr.Geometry(ogr.wkbPoint)
            shore_point_geometry.AddPoint(x_coord, y_coord)
            shore_point_geometry.Transform(utm_to_base_transform)
            if grid_shapely.intersects(shapely.geometry.Point(
                    shore_point_geometry.GetX(), shore_point_geometry.GetY())):
                shore_point_feature = ogr.Feature(target_shore_point_defn)
                shore_point_feature.SetGeometry(shore_point_geometry)

                nearest_points = list(wwiii_rtree.nearest(
                    (shore_point_geometry.GetX(),
                     shore_point_geometry.GetY(),
                     shore_point_geometry.GetX(),
                     shore_point_geometry.GetY()), 3))[0:3]
                shore_point_shapely = shapely.geometry.Point(
                    (shore_point_geometry.GetX(),
                     shore_point_geometry.GetY()))
                for field_name in field_names:
                    value = 0.0
                    total_distance = 0.0
                    for fid in nearest_points:
                        wwiii_feature = wwiii_layer.GetFeature(fid)
                        wwiii_shapely = shapely.wkb.loads(
                            wwiii_feature.GetGeometryRef().ExportToWkb())
                        distance = shore_point_shapely.distance(
                            wwiii_shapely)
                        field_value = wwiii_feature.GetField(field_name)
                        value += (float(field_value) * distance)
                        total_distance += distance
                    shore_point_feature.SetField(
                        field_name, value / total_distance)

                target_shore_point_layer.CreateFeature(shore_point_feature)
                shore_point_feature = None
    with open(work_complete_token_path, 'w') as work_complete_token_file:
        work_complete_token_file.write("Complete!\n")
    #shutil.rmtree(temp_workspace)


def _grid_edges_of_vector(
        base_bounding_box, base_vector_path,
        base_feature_bounding_box_rtree_path, target_grid_vector_path,
        target_grid_size, done_token_path):
    """Build a sparse grid covering the edges of the base polygons.

    Parameters:
        base_bounding_box (list/tuple): format [minx, miny, maxx, maxy]
            represents the bounding box of the underlying grid cells to be
            created.
        base_vector_path (string): path to shapefile of polygon features,
            a grid cell will be created if it contains the *edge* of any
            feature.
        base_feature_bounding_box_rtree_path (string): path to an rtree that
            indexes the bounding boxes of the polygon features in
            `base_vector_path`.
        target_grid_vector_path (string): path to shapefile grid that will be
            created by this function.
        target_grid_size (float): length of side of the grid cell to be
            created in the target_grid_vector.
        done_token_path (string): path to a file to create when the function
            is complete.
    """
    LOGGER.info("Building global grid.")
    n_rows = int((
        base_bounding_box[3]-base_bounding_box[1]) / float(
            target_grid_size))
    n_cols = int((
        base_bounding_box[2]-base_bounding_box[0]) / float(
            target_grid_size))

    if os.path.exists(target_grid_vector_path):
        os.remove(target_grid_vector_path)

    # create the target vector as an ESRI shapefile
    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
    target_grid_vector = esri_shapefile_driver.CreateDataSource(
        target_grid_vector_path)

    # create the spatial reference from the base vector
    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromWkt(
        pygeoprocessing.get_vector_info(base_vector_path)['projection'])

    target_grid_layer = target_grid_vector.CreateLayer(
        os.path.splitext(target_grid_vector_path)[0],
        spatial_reference, ogr.wkbPolygon)

    target_grid_defn = target_grid_layer.GetLayerDefn()

    base_vector = ogr.Open(base_vector_path)
    base_layer = base_vector.GetLayer()

    # the input path has a .dat extension, but the `rtree` package only uses
    # the basename.  It's a quirk of the library, so we'll deal with it by
    # cutting off the extension.
    target_feature_rtree_index = rtree.index.Index(
        os.path.splitext(base_feature_bounding_box_rtree_path)[0])

    logger_callback = _make_logger_callback(
        'Cell coverage %.2f%% complete')

    prepared_geometry = {}
    for cell_index in xrange(n_rows * n_cols):
        logger_callback(float(cell_index) / (n_rows * n_cols))
        row_index = cell_index / n_cols
        col_index = cell_index % n_cols
        # format of bounding box is  [xmin, ymin, xmax, ymax]
        cell_bounding_box = [
            base_bounding_box[0]+col_index*target_grid_size,
            base_bounding_box[1]+row_index*target_grid_size,
            base_bounding_box[0]+(col_index+1)*target_grid_size,
            base_bounding_box[1]+(row_index+1)*target_grid_size]

        intersections = list(
            target_feature_rtree_index.intersection(cell_bounding_box))
        if len(intersections) == 0:
            # skip this cell if no intersections with the bounding boxes occur
            continue

        # construct cell geometry both in OGR and Shapely formats
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for i, j in [(0, 1), (2, 1), (2, 3), (0, 3), (0, 1)]:
            ring.AddPoint(cell_bounding_box[i], cell_bounding_box[j])
        cell_geometry = ogr.Geometry(ogr.wkbPolygon)
        cell_geometry.AddGeometry(ring)
        cell_geometry_shapely = shapely.wkb.loads(
            cell_geometry.ExportToWkb())

        for fid in intersections:
            if fid not in prepared_geometry:
                base_feature = base_layer.GetFeature(fid)
                base_feature_geometry = base_feature.GetGeometryRef()
                prepared_geometry[fid] = shapely.wkb.loads(
                    base_feature_geometry.ExportToWkb())
                base_feature_geometry = None

            if (prepared_geometry[fid].intersects(
                    cell_geometry_shapely) and
                not prepared_geometry[fid].contains(
                    cell_geometry_shapely)):
                # add cell to target layer if it intersects the edge of a
                # base polygon feature
                target_feature = ogr.Feature(target_grid_defn)
                target_feature.SetGeometry(cell_geometry)
                target_grid_layer.CreateFeature(target_feature)
                target_feature = None
                # no need to test the rest if one intersects
                break
    target_grid_layer.SyncToDisk()
    target_grid_layer = None
    target_grid_vector = None

    with open(done_token_path, 'w') as done_token_file:
        done_token_file.write("Done!\n")


def _get_utm_spatial_reference(bounding_box):
    """Determine UTM spatial reference given lat/lng bounding box.

    Parameter:
        bounding_box (list/tuple): UTM84 coordinate bounding box in the
            format [min_lng, min_lat, max_lng, max_lat]

    Returns:
        An osr.SpatialReference that corresponds to the UTM zone in the
            median point of the bounding box.
    """
    # project lulc_map to UTM zone median
    mean_long = (bounding_box[0] + bounding_box[2]) / 2
    mean_lat = (bounding_box[1] + bounding_box[3]) / 2
    utm_code = (math.floor((mean_long + 180)/6) % 60) + 1

    # Determine if north or sourth
    lat_code = 6 if mean_lat > 0 else 7

    # and this is the format of an EPSG UTM code:
    epsg_code = int('32%d%02d' % (lat_code, utm_code))

    utm_sr = osr.SpatialReference()
    utm_sr.ImportFromEPSG(epsg_code)
    return utm_sr


def _build_feature_bounding_box_rtree(
        vector_path, target_rtree_path, done_token_path):
    """Builds an r-tree index of the global feature envelopes.

    Parameter:
        vector_path (string): path to vector to build bounding box index for
        target_rtree_path (string): path to ".dat" file to store the saved
            r-tree.
        done_token_path (string): path to a file to create when function is
            complete.

    Returns:
        None.
    """
    LOGGER.info("Building rtree index.")
    # the input path has a .dat extension, but the `rtree` package only uses
    # the basename.  It's a quirk of the library, so we'll deal with it by
    # cutting off the extension.
    global_feature_index_base = os.path.splitext(
        target_rtree_path)[0]
    global_feature_index = rtree.index.Index(global_feature_index_base)

    global_vector = ogr.Open(vector_path)
    global_layer = global_vector.GetLayer()
    n_features = global_layer.GetFeatureCount()

    logger_callback = _make_logger_callback(
        'rTree construction %.2f%% complete')

    for feature_index, global_feature in enumerate(global_layer):
        feature_geometry = global_feature.GetGeometryRef()
        # format of envelope is [minx, maxx, miny, maxy]
        feature_envelope = feature_geometry.GetEnvelope()
        # format of tree bounding box is [minx, miny, maxx, maxy]
        global_feature_index.insert(
            global_feature.GetFID(), (
                feature_envelope[0], feature_envelope[2],
                feature_envelope[1], feature_envelope[3]))
        logger_callback(float(feature_index) / n_features)
    global_feature_index.close()

    with open(done_token_path, 'w') as done_token_file:
        done_token_file.write("Done!\n")


def _make_shore_kernel(kernel_path):
    """Make a 3x3 raster with a 9 in the middle and 1s on the outside."""
    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_path.encode('utf-8'), 3, 3, 1,
        gdal.GDT_Byte)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([0, 1, 0, 0, 0, -1])
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(127)
    kernel_band.WriteArray(numpy.array([[1, 1, 1], [1, 9, 1], [1, 1, 1]]))


def _make_logger_callback(message):
    """Build a timed logger callback that prints `message` replaced.

    Parameters:
        message (string): a string that expects a %f replacement variable for
            `proportion_complete`.

    Returns:
        Function with signature:
            logger_callback(proportion_complete, psz_message, p_progress_arg)
    """
    def logger_callback(proportion_complete):
        """The argument names come from the GDAL API for callbacks."""
        try:
            current_time = time.time()
            if ((current_time - logger_callback.last_time) > 5.0 or
                    (proportion_complete == 1.0 and
                     logger_callback.total_time >= 5.0)):
                LOGGER.info(message, proportion_complete * 100)
                logger_callback.last_time = current_time
                logger_callback.total_time += current_time
        except AttributeError:
            logger_callback.last_time = time.time()
            logger_callback.total_time = 0.0

    return logger_callback


if __name__ == '__main__':
    main()
