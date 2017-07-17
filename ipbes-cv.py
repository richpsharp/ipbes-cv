"""IPBES global coastal vulnerability calculation."""
import traceback
import shutil
import time
import os
import math
import logging
import multiprocessing

import numpy
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import rtree
import shapely
import shapely.wkb
import shapely.ops
import shapely.speedups
import pygeoprocessing

import Task

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

_TARGET_WORKSPACE = "ipbes_cv_workspace"

_GLOBAL_POLYGON_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\Land\global_polygon.shp"

_GLOBAL_WWIII_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\CoastalProtection\Input\WaveWatchIII.shp"

_GLOBAL_DEM_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\DEMs\global_dem"

# layer name, (layer path, layer rank, protection distance)
_GLOBAL_HABITAT_LAYER_PATHS = {
    'mangrove': (r"C:\Users\rpsharp\Dropbox\ipbes-data\cv\habitat\DataPack-14_001_WCMC010_MangrovesUSGS2011_v1_3\01_Data\14_001_WCMC010_MangroveUSGS2011_v1_3.shp", 1, 1000.0),
    'saltmarsh': (r"C:\Users\rpsharp\Dropbox\ipbes-data\cv\habitat\Datapack-14_001_WCMC027_Saltmarsh2017_v4\01_Data\14_001_WCMC027_Saltmarsh_py_v4.shp", 2, 1000.0),
    'coralreef': (r"C:\Users\rpsharp\Dropbox\ipbes-data\cv\habitat\DataPack-14_001_WCMC008_CoralReef2010_v1_3\01_Data\14_001_WCMC008_CoralReef2010_v1_3.shp", 1, 2000.0),
    'seagrass': (r"C:\Users\rpsharp\Dropbox\ipbes-data\cv\habitat\DataPack-14_001_WCMC013_014_SeagrassPtPy_v4\01_Data\WCMC_013_014_SeagrassesPy_v4.shp", 4, 500.0),
}

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

# Wave Watch III data does not cover the planet.  Make sure we don't deal
# with a point that's not in range of said point.  I'm picking 160km since
# that's double the diagonal distance between two WWIII points
_MAX_WWIII_DISTANCE = 160000.0

_N_FETCH_RAYS = 16

_GLOBAL_GRID_VECTOR_FILE_PATTERN = 'global_grid.shp'
_LANDMASS_BOUNDING_RTREE_FILE_PATTERN = 'global_feature_index.dat'
_GLOBAL_WWIII_RTREE_FILE_PATTERN = 'wwiii_rtree.dat'
_GRID_POINT_FILE_PATTERN = 'grid_points_%d.shp'
_WIND_EXPOSURE_POINT_FILE_PATTERN = 'rei_points_%d.shp'
_WAVE_EXPOSURE_POINT_FILE_PATTERN = 'wave_points_%d.shp'
_HABITAT_PROTECTION_POINT_FILE_PATTERN = 'habitat_protection_points_%d.shp'
_RELIEF_POINT_FILE_PATTERN = 'relief_%d.shp'
_GLOBAL_REI_POINT_FILE_PATTERN = 'global_rei_points.shp'
_GLOBAL_WAVE_POINT_FILE_PATTERN = 'global_wave_points.shp'
_GLOBAL_HABITAT_PROTECTION_FILE_PATTERN = (
    'global_habitat_protection_points.shp')
_GLOBAL_FETCH_RAY_FILE_PATTERN = 'global_fetch_rays.shp'
_WORK_COMPLETE_TOKEN_PATH = os.path.join(
    _TARGET_WORKSPACE, 'work_tokens')
_WIND_EXPOSURE_WORKSPACES = os.path.join(
    _TARGET_WORKSPACE, 'wind_exposure_workspaces')
_WAVE_EXPOSURE_WORKSPACES = os.path.join(
    _TARGET_WORKSPACE, 'wave_exposure_workspaces')
_HABITAT_PROTECTION_WORKSPACES = os.path.join(
    _TARGET_WORKSPACE, 'habitat_protection_workspaces')
_RELIEF_WORKSPACES = os.path.join(
    _TARGET_WORKSPACE, 'relief_workspaces')
_GRID_WORKSPACES = os.path.join(
    _TARGET_WORKSPACE, 'grid_workspaces')

_SMALLEST_FEATURE_SIZE = 2000
_MAX_FETCH_DISTANCE = 60000


def main():
    """Entry point."""
    logger = logging.getLogger('ipbes-cv.main')
    if not os.path.exists(_TARGET_WORKSPACE):
        os.makedirs(_TARGET_WORKSPACE)

    task_graph = Task.TaskGraph(
        _WORK_COMPLETE_TOKEN_PATH, multiprocessing.cpu_count())

    wwiii_rtree_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_WWIII_RTREE_FILE_PATTERN)

    build_wwiii_task = task_graph.add_task(
        target=build_wwiii_rtree, args=(
            _GLOBAL_WWIII_PATH, wwiii_rtree_path))

    simplified_vector_path = os.path.join(
        _TARGET_WORKSPACE, 'simplified_geometry.shp')
    # make an approximation of smallest feature size in degrees
    smallest_feature_size_degrees = 1. / 111000 * _SMALLEST_FEATURE_SIZE / 2.0
    simplify_geometry_task = task_graph.add_task(
        target=simplify_geometry, args=(
            _GLOBAL_POLYGON_PATH, smallest_feature_size_degrees,
            simplified_vector_path))

    simplify_habitat_task_list = []
    simplified_habitat_vector_lookup = {}
    for habitat_id, (habitat_path, habitat_rank, habitat_dist) in (
            _GLOBAL_HABITAT_LAYER_PATHS.iteritems()):
        smallest_feature_size_degrees = 1. / 111000 * habitat_dist / 2.0
        simplified_habitat_vector_lookup[habitat_id] = (
            os.path.join(_TARGET_WORKSPACE, '%s.shp' % habitat_id),
            habitat_rank, habitat_dist)
        simplify_habitat_task = task_graph.add_task(
            target=simplify_geometry, args=(
                habitat_path, smallest_feature_size_degrees,
                simplified_habitat_vector_lookup[habitat_id][0]))
        simplify_habitat_task_list.append(simplify_habitat_task)

    landmass_bounding_rtree_path = os.path.join(
        _TARGET_WORKSPACE, _LANDMASS_BOUNDING_RTREE_FILE_PATTERN)

    build_rtree_task = task_graph.add_task(
        target=build_feature_bounding_box_rtree,
        args=(simplified_vector_path, landmass_bounding_rtree_path),
        dependent_task_list=[simplify_geometry_task])

    global_grid_vector_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_GRID_VECTOR_FILE_PATTERN)

    grid_edges_of_vector_task = task_graph.add_task(
        target=grid_edges_of_vector, args=(
            _GLOBAL_BOUNDING_BOX_WGS84, simplified_vector_path,
            landmass_bounding_rtree_path, global_grid_vector_path,
            _WGS84_GRID_SIZE), dependent_task_list=[build_rtree_task])

    task_graph.join()

    global_grid_vector = ogr.Open(global_grid_vector_path)
    global_grid_layer = global_grid_vector.GetLayer()
    grid_count = global_grid_layer.GetFeatureCount()
    global_grid_layer = None
    global_grid_vector = None

    local_rei_point_path_list = []
    wind_exposure_task_list = []
    wave_exposure_task_list = []
    local_fetch_ray_path_list = []
    local_wave_point_path_list = []
    habitat_protection_task_list = []
    local_habitat_protection_path_list = []
    relief_task_list = []
    local_relief_path_list = []
    #for grid_id in xrange(grid_count):
    for grid_id in xrange(grid_count):
        logger.info("Calculating grid %d of %d", grid_id, grid_count)

        shore_points_workspace = os.path.join(
            _GRID_WORKSPACES, 'grid_%d' % grid_id)
        grid_point_path = os.path.join(
            shore_points_workspace, _GRID_POINT_FILE_PATTERN % (grid_id))

        create_shore_points_task = task_graph.add_task(
            target=create_shore_points, args=(
                global_grid_vector_path, grid_id, landmass_bounding_rtree_path,
                simplified_vector_path, _GLOBAL_WWIII_PATH, wwiii_rtree_path,
                _SMALLEST_FEATURE_SIZE, shore_points_workspace,
                grid_point_path),
            dependent_task_list=[grid_edges_of_vector_task, build_wwiii_task])

        wind_exposure_workspace = os.path.join(
            _WIND_EXPOSURE_WORKSPACES, 'wind_exposure_%d' % grid_id)
        target_wind_exposure_point_path = os.path.join(
            wind_exposure_workspace,
            _WIND_EXPOSURE_POINT_FILE_PATTERN % grid_id)
        wind_exposure_task = task_graph.add_task(
            target=calculate_wind_exposure, args=(
                grid_point_path, landmass_bounding_rtree_path,
                simplified_vector_path, wind_exposure_workspace,
                _SMALLEST_FEATURE_SIZE, _MAX_FETCH_DISTANCE,
                target_wind_exposure_point_path),
            dependent_task_list=[create_shore_points_task])
        wind_exposure_task_list.append(
            wind_exposure_task)
        local_rei_point_path_list.append(
            target_wind_exposure_point_path)
        local_fetch_ray_path_list.append(
            os.path.join(wind_exposure_workspace, 'fetch_rays.shp'))

        wave_exposure_workspace = os.path.join(
            _WAVE_EXPOSURE_WORKSPACES, 'wave_exposure_%d' % grid_id)
        target_wave_exposure_point_path = os.path.join(
            wave_exposure_workspace,
            _WAVE_EXPOSURE_POINT_FILE_PATTERN % grid_id)
        wave_exposure_task = task_graph.add_task(
            target=calculate_wave_exposure, args=(
                target_wind_exposure_point_path, _MAX_FETCH_DISTANCE,
                wave_exposure_workspace,
                target_wave_exposure_point_path),
            dependent_task_list=[wind_exposure_task])
        wave_exposure_task_list.append(
            wave_exposure_task)
        local_wave_point_path_list.append(
            target_wave_exposure_point_path)

        habitat_protection_workspace = os.path.join(
            _HABITAT_PROTECTION_WORKSPACES, 'habitat_protection_%d' % grid_id)
        target_habitat_protection_point_path = os.path.join(
            habitat_protection_workspace,
            _HABITAT_PROTECTION_POINT_FILE_PATTERN % grid_id)
        habitat_protection_task = task_graph.add_task(
            target=calculate_habitat_protection, args=(
                grid_point_path,
                simplified_habitat_vector_lookup,
                habitat_protection_workspace,
                target_habitat_protection_point_path),
            dependent_task_list=[create_shore_points_task])
        habitat_protection_task_list.append(habitat_protection_task)
        local_habitat_protection_path_list.append(
            target_habitat_protection_point_path)

        relief_workspace = os.path.join(
            _RELIEF_WORKSPACES, 'relief_%d' % grid_id)
        target_relief_point_vector_path = os.path.join(
            relief_workspace, _RELIEF_POINT_FILE_PATTERN % grid_id)
        relief_task = task_graph.add_task(
            target=calculate_relief, args=(
                grid_point_path, _GLOBAL_DEM_PATH,
                relief_workspace,
                target_relief_point_vector_path),
            dependent_task_list=[create_shore_points_task])
        relief_task_list.append(relief_task)
        local_relief_path_list.append(
            target_relief_point_vector_path)

    target_merged_rei_points_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_REI_POINT_FILE_PATTERN)
    target_spatial_reference_wkt = pygeoprocessing.get_vector_info(
        _GLOBAL_POLYGON_PATH)['projection']
    _ = task_graph.add_task(
        target=merge_vectors, args=(
            local_rei_point_path_list,
            target_spatial_reference_wkt,
            target_merged_rei_points_path,
            ['REI']),
        dependent_task_list=wind_exposure_task_list)

    target_merged_wave_points_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_WAVE_POINT_FILE_PATTERN)
    _ = task_graph.add_task(
        target=merge_vectors, args=(
            local_wave_point_path_list,
            target_spatial_reference_wkt,
            target_merged_wave_points_path,
            ['Ew']),
        dependent_task_list=wave_exposure_task_list)

    target_merged_relief_points_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_RELIEF_POINT_FILE_PATTERN)
    _ = task_graph.add_task(
        target=merge_vectors, args=(
            local_relief_path_list,
            target_spatial_reference_wkt,
            target_merged_relief_points_path,
            ['relief']),
        dependent_task_list=relief_task_list)

    target_habitat_protection_points_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_HABITAT_PROTECTION_FILE_PATTERN)
    _ = task_graph.add_task(
        target=merge_vectors, args=(
            local_habitat_protection_path_list,
            target_spatial_reference_wkt,
            target_habitat_protection_points_path,
            ['Rhab']),
        dependent_task_list=habitat_protection_task_list)

    target_merged_fetch_rays_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_FETCH_RAY_FILE_PATTERN)
    target_spatial_reference_wkt = pygeoprocessing.get_vector_info(
        _GLOBAL_POLYGON_PATH)['projection']
    _ = task_graph.add_task(
        target=merge_vectors, args=(
            local_fetch_ray_path_list,
            target_spatial_reference_wkt,
            target_merged_fetch_rays_path,
            []),
        dependent_task_list=wind_exposure_task_list)

    task_graph.join()


def simplify_geometry(
        base_vector_path, tolerance, target_simplified_vector_path):
    """Simplify all the geometry in the vector.

    Parameters:
        base_vector_path (string): path to base vector.
        tolerance (float): all new vertices in the geometry will be within
            this distance (in units of the vector's projection).
        target_simplified_vector_path (string): path to desired simplified
            vector.

    Returns:
        None
    """
    logger = logging.getLogger('ipbes-cv.simplify_geometry')
    base_vector = ogr.Open(base_vector_path)
    base_layer = base_vector.GetLayer()

    if os.path.exists(target_simplified_vector_path):
        os.remove(target_simplified_vector_path)

    esri_driver = ogr.GetDriverByName('ESRI Shapefile')

    target_simplified_vector = esri_driver.CreateDataSource(
        target_simplified_vector_path)
    target_simplified_layer = target_simplified_vector.CreateLayer(
        os.path.splitext(os.path.basename(target_simplified_vector_path))[0],
        base_layer.GetSpatialRef(), ogr.wkbPolygon)

    for feature in base_layer:
        target_feature = ogr.Feature(target_simplified_layer.GetLayerDefn())
        feature_geometry = feature.GetGeometryRef()
        simplified_geometry = feature_geometry.Simplify(tolerance)
        feature_geometry = None
        if (simplified_geometry is not None and
                simplified_geometry.GetArea() > 0):
            target_feature.SetGeometry(simplified_geometry)
            target_simplified_layer.CreateFeature(target_feature)

    target_simplified_layer.SyncToDisk()
    target_simplified_layer = None
    target_simplified_vector = None


def calculate_habitat_protection(
        base_shore_point_vector_path,
        habitat_layer_lookup, workspace_dir,
        target_habitat_protection_point_vector_path):
    """Calculate habitat protection at a set of points.

    Parameters:
        base_shore_point_vector_path (string):  path to a point shapefile to
            analyze habitat protection at.
        habitat_layer_lookup: a dictionary mapping habitat id to a
            (path, rank, distance) tuple
        workspace_dir (string): path to a directory to make local calculations
            in
        target_habitat_protection_point_vector_path (string): path to desired
            output vector.  after completion will have a rank for each
            habitat ID, and a field called Rhab with a value from 1-5
            indicating relative level of protection of that point.
    Returns:
        None.
    """
    try:
        logger = logging.getLogger('ipbes-cv.calculate_habitat_protection')
        if not os.path.exists(os.path.dirname(
                target_habitat_protection_point_vector_path)):
            os.makedirs(
                os.path.dirname(target_habitat_protection_point_vector_path))
        if os.path.exists(target_habitat_protection_point_vector_path):
            os.remove(target_habitat_protection_point_vector_path)

        base_ref_wkt = pygeoprocessing.get_vector_info(
            base_shore_point_vector_path)['projection']
        base_spatial_reference = osr.SpatialReference()
        base_spatial_reference.ImportFromWkt(base_ref_wkt)

        # reproject base_shore_point_vector_path to utm coordinates
        base_shore_info = pygeoprocessing.get_vector_info(
            base_shore_point_vector_path)
        base_shore_bounding_box = base_shore_info['bounding_box']

        utm_spatial_reference = get_utm_spatial_reference(
            base_shore_info['bounding_box'])
        base_spatial_reference = osr.SpatialReference()
        base_spatial_reference.ImportFromWkt(base_shore_info['projection'])

        pygeoprocessing.reproject_vector(
            base_shore_point_vector_path, utm_spatial_reference.ExportToWkt(),
            target_habitat_protection_point_vector_path)

        utm_bounding_box = pygeoprocessing.get_vector_info(
            target_habitat_protection_point_vector_path)['bounding_box']

        # extend bounding box for max fetch distance
        utm_bounding_box = [
            utm_bounding_box[0] - _MAX_FETCH_DISTANCE,
            utm_bounding_box[1] - _MAX_FETCH_DISTANCE,
            utm_bounding_box[2] + _MAX_FETCH_DISTANCE,
            utm_bounding_box[3] + _MAX_FETCH_DISTANCE]

        # get lat/lng bounding box of utm projected coordinates

        # get global polygon clip of that utm box
        # transform local box back to lat/lng -> global clipping box
        lat_lng_clipping_box = pygeoprocessing.transform_bounding_box(
            utm_bounding_box, utm_spatial_reference.ExportToWkt(),
            base_spatial_reference.ExportToWkt(), edge_samples=11)
        if (base_shore_bounding_box[0] > 0 and
                lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[2] += 360
        elif (base_shore_bounding_box[0] < 0 and
              lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[0] -= 360
        elif base_shore_bounding_box == [0, 0, 0, 0]:
            # this case guards for an empty shore point in case there are
            # very tiny islands or such
            lat_lng_clipping_box = (0, 0, 0, 0)
        lat_lng_clipping_shapely = shapely.geometry.box(*lat_lng_clipping_box)
        habitat_shapely_lookup = {}

        for habitat_id in habitat_layer_lookup:
            habitat_vector = ogr.Open(habitat_layer_lookup[habitat_id][0])
            habitat_layer = habitat_vector.GetLayer()

            # this will hold the clipped landmass geometry
            esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
            temp_clipped_vector_path = os.path.join(
                workspace_dir, 'clipped_habitat_%s.shp' % habitat_id)
            utm_clipped_vector_path = os.path.join(
                workspace_dir, 'utm_clipped_habitat_%s.shp' % habitat_id)
            for path in [temp_clipped_vector_path, utm_clipped_vector_path]:
                if os.path.exists(path):
                    os.remove(path)
            temp_clipped_vector = esri_shapefile_driver.CreateDataSource(
                temp_clipped_vector_path)
            temp_clipped_layer = (
                temp_clipped_vector.CreateLayer(
                    os.path.splitext(temp_clipped_vector_path)[0],
                    base_spatial_reference, ogr.wkbPolygon))
            temp_clipped_defn = temp_clipped_layer.GetLayerDefn()

            # clip global polygon to global clipping box
            for habitat_feature in habitat_layer:
                habitat_shapely = shapely.wkb.loads(
                    habitat_feature.GetGeometryRef().ExportToWkb())
                intersection_shapely = lat_lng_clipping_shapely.intersection(
                    habitat_shapely)
                try:
                    clipped_geometry = ogr.CreateGeometryFromWkt(
                        intersection_shapely.wkt)
                    clipped_feature = ogr.Feature(temp_clipped_defn)
                    clipped_feature.SetGeometry(clipped_geometry)
                    temp_clipped_layer.CreateFeature(clipped_feature)
                    clipped_feature = None
                except Exception:
                    logger.warn(
                        "Couldn't process this intersection %s",
                        intersection_shapely)
            temp_clipped_layer.SyncToDisk()
            temp_clipped_layer = None
            temp_clipped_vector = None

            pygeoprocessing.reproject_vector(
                temp_clipped_vector_path, utm_spatial_reference.ExportToWkt(),
                utm_clipped_vector_path)

            clipped_geometry_shapely_list = []
            temp_utm_clipped_vector = ogr.Open(utm_clipped_vector_path)
            temp_utm_clipped_layer = temp_utm_clipped_vector.GetLayer()
            for tmp_utm_feature in temp_utm_clipped_layer:
                tmp_utm_geometry = tmp_utm_feature.GetGeometryRef()
                shapely_geometry = shapely.wkb.loads(
                    tmp_utm_geometry.ExportToWkb())
                if shapely_geometry.is_valid:
                    clipped_geometry_shapely_list.append(shapely_geometry)
                tmp_utm_geometry = None
            temp_utm_clipped_layer = None
            temp_utm_clipped_vector = None
            habitat_shapely_lookup[habitat_id] = shapely.ops.cascaded_union(
                clipped_geometry_shapely_list)

        esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
        os.path.dirname(base_shore_point_vector_path)
        target_habitat_protection_point_vector = ogr.Open(
            target_habitat_protection_point_vector_path, 1)
        target_habitat_protection_point_layer = (
            target_habitat_protection_point_vector.GetLayer())
        target_habitat_protection_point_layer.CreateField(ogr.FieldDefn(
            'Rhab', ogr.OFTReal))
        for habitat_id in habitat_layer_lookup:
            target_habitat_protection_point_layer.CreateField(ogr.FieldDefn(
                habitat_id, ogr.OFTReal))

        for target_feature in target_habitat_protection_point_layer:
            target_feature_geometry = (
                target_feature.GetGeometryRef().Clone())

            target_feature_shapely = shapely.wkb.loads(
                target_feature_geometry.ExportToWkb())
            min_rank = 5
            sum_sq_rank = 0.0
            for habitat_id, (_, rank, protection_distance) in (
                    habitat_layer_lookup.iteritems()):
                if (target_feature_shapely.distance(
                        habitat_shapely_lookup[habitat_id]) >
                        protection_distance):
                    # 5 is the worst rank
                    rank = 5
                if rank < min_rank:
                    min_rank = rank
                target_feature.SetField(habitat_id, rank)
                sum_sq_rank += (5 - rank) ** 2

            # Equation 4
            target_feature.SetField(
                'Rhab', 4.8 - 0.5 * (
                    (1.5 * (5-min_rank))**2 +
                    (sum_sq_rank - (5-min_rank))**2)**0.5)
            target_feature.SetGeometry(target_feature_geometry)
            target_habitat_protection_point_layer.SetFeature(
                target_feature)
        target_habitat_protection_point_layer.SyncToDisk()
        target_habitat_protection_point_layer = None
        target_habitat_protection_point_vector = None
    except Exception as e:
        traceback.print_exc()
        raise


def calculate_wave_exposure(
        base_fetch_point_vector_path, max_fetch_distance, workspace_dir,
        target_wave_exposure_point_vector_path):
    """Calculate the wave exposure index.

    Parameters:
        base_fetch_point_vector_path (string): path to a point shapefile that
            contains 16 'WavP_[direction]' fields, 'WavPPCT[direction]'
            fields, 'fdist_[direction]' fields, a single H, and a single T
            field.
        max_fetch_distance (float): max fetch distance before a wind fetch ray
            is terminated
        target_wave_exposure_point_vector_path (string): path to an output
            shapefile that will contain a field called 'Ew' which is the
            maximum of ocean or wind waves occurring at that point.

    Returns:
        None
    """
    # this will hold the clipped landmass geometry
    try:
        if not os.path.exists(os.path.dirname(
                target_wave_exposure_point_vector_path)):
            os.makedirs(
                os.path.dirname(target_wave_exposure_point_vector_path))
        if os.path.exists(target_wave_exposure_point_vector_path):
            os.remove(target_wave_exposure_point_vector_path)

        base_ref_wkt = pygeoprocessing.get_vector_info(
            base_fetch_point_vector_path)['projection']
        base_spatial_reference = osr.SpatialReference()
        base_spatial_reference.ImportFromWkt(base_ref_wkt)

        esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
        os.path.dirname(base_fetch_point_vector_path)
        target_wave_exposure_point_vector = (
            esri_shapefile_driver.CreateDataSource(
                target_wave_exposure_point_vector_path))
        target_wave_exposure_point_layer = (
            target_wave_exposure_point_vector.CreateLayer(
                os.path.splitext(target_wave_exposure_point_vector_path)[0],
                base_spatial_reference, ogr.wkbPoint))
        target_wave_exposure_point_layer.CreateField(ogr.FieldDefn(
            'Ew', ogr.OFTReal))
        target_wave_exposure_point_defn = (
            target_wave_exposure_point_layer.GetLayerDefn())

        base_fetch_point_vector = ogr.Open(base_fetch_point_vector_path)
        base_fetch_point_layer = base_fetch_point_vector.GetLayer()
        for base_fetch_point_feature in base_fetch_point_layer:
            target_feature = ogr.Feature(target_wave_exposure_point_defn)
            target_feature.SetGeometry(
                base_fetch_point_feature.GetGeometryRef().Clone())
            e_local_wedge = (
                0.5 *
                float(base_fetch_point_feature.GetField('H_10PCT'))**2 *
                float(base_fetch_point_feature.GetField('H_10PCT'))) / float(
                    _N_FETCH_RAYS)
            e_ocean = 0.0
            e_local = 0.0
            for sample_index in xrange(_N_FETCH_RAYS):
                compass_degree = int(sample_index * 360 / 16.)
                fdist = base_fetch_point_feature.GetField(
                    'fdist_%d' % compass_degree)
                if numpy.isclose(fdist, max_fetch_distance):
                    e_ocean += (
                        base_fetch_point_feature.GetField(
                            'WavP_%d' % compass_degree) *
                        base_fetch_point_feature.GetField(
                            'WavPPCT%d' % compass_degree))
                elif fdist > 0.0:
                    e_local += e_local_wedge
            target_feature.SetField('Ew', max(e_ocean, e_local))
            target_wave_exposure_point_layer.CreateFeature(target_feature)
    except Exception as e:
        traceback.print_exc()
        raise


def calculate_wind_exposure(
        base_shore_point_vector_path,
        landmass_bounding_rtree_path, landmass_vector_path, workspace_dir,
        smallest_feature_size, max_fetch_distance,
        target_fetch_point_vector_path):
    """Calculate wind exposure for each shore point.

    Parameters:
        base_shore_point_vector_path (string): path to a point shapefile
            representing shore points that should be sampled for wind
            exposure.
        landmass_bounding_rtree_path (string): path to an rtree bounding box
            for the landmass polygons.
        landmass_vector_path (string): path to landmass polygon vetor.
        workspace_dir (string): path to a directory that can be created for
            temporary workspace files
        smallest_feature_size (float): smallest feature size to detect in
            meters.
        max_fetch_distance (float): maximum fetch distance for a ray in
            meters.
        target_fetch_point_vector_path (string): path to target point file,
            will be a copy of `base_shore_point_vector_path`'s geometry with
            an 'REI' (relative exposure index) field added.

    Returns:
        None
    """
    try:
        logger = logging.getLogger('ipbes-cv.calculate_wind_exposure')
        if os.path.exists(workspace_dir):
            shutil.rmtree(workspace_dir)
        os.makedirs(workspace_dir)

        utm_clipped_vector_path = os.path.join(
            workspace_dir, 'utm_clipped_landmass.shp')
        temp_fetch_rays_path = os.path.join(
            workspace_dir, 'fetch_rays.shp')

        # reproject base_shore_point_vector_path to utm coordinates
        base_shore_info = pygeoprocessing.get_vector_info(
            base_shore_point_vector_path)
        base_shore_bounding_box = base_shore_info['bounding_box']

        utm_spatial_reference = get_utm_spatial_reference(
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
        lat_lng_clipping_box = pygeoprocessing.transform_bounding_box(
            utm_bounding_box, utm_spatial_reference.ExportToWkt(),
            base_spatial_reference.ExportToWkt(), edge_samples=11)
        if (base_shore_bounding_box[0] > 0 and
                lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[2] += 360
        elif (base_shore_bounding_box[0] < 0 and
              lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[0] -= 360
        elif base_shore_bounding_box == [0, 0, 0, 0]:
            # this case guards for an empty shore point in case there are
            # very tiny islands or such
            lat_lng_clipping_box = (0, 0, 0, 0)
        lat_lng_clipping_shapely = shapely.geometry.box(*lat_lng_clipping_box)

        landmass_vector_rtree = rtree.index.Index(
            os.path.splitext(landmass_bounding_rtree_path)[0])

        landmass_vector = ogr.Open(landmass_vector_path)
        landmass_layer = landmass_vector.GetLayer()

        # this will hold the clipped landmass geometry
        esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
        temp_clipped_vector_path = os.path.join(
            workspace_dir, 'clipped_geometry_vector.shp')
        temp_clipped_vector = esri_shapefile_driver.CreateDataSource(
            temp_clipped_vector_path)
        temp_clipped_layer = (
            temp_clipped_vector.CreateLayer(
                os.path.splitext(temp_clipped_vector_path)[0],
                base_spatial_reference, ogr.wkbPolygon))
        temp_clipped_defn = temp_clipped_layer.GetLayerDefn()

        # clip global polygon to global clipping box
        for feature_id in landmass_vector_rtree.intersection(
                lat_lng_clipping_box):
            landmass_feature = landmass_layer.GetFeature(feature_id)
            landmass_shapely = shapely.wkb.loads(
                landmass_feature.GetGeometryRef().ExportToWkb())
            intersection_shapely = lat_lng_clipping_shapely.intersection(
                landmass_shapely)
            try:
                clipped_geometry = ogr.CreateGeometryFromWkt(
                    intersection_shapely.wkt)
                clipped_feature = ogr.Feature(temp_clipped_defn)
                clipped_feature.SetGeometry(clipped_geometry)
                temp_clipped_layer.CreateFeature(clipped_feature)
                clipped_feature = None
            except Exception:
                logger.warn(
                    "Couldn't process this intersection %s", intersection_shapely)
        temp_clipped_layer.SyncToDisk()
        temp_clipped_layer = None
        temp_clipped_vector = None

        # project global clipped polygons to UTM
        logger.info("reprojecting grid %s", base_shore_point_vector_path)
        pygeoprocessing.reproject_vector(
            temp_clipped_vector_path, utm_spatial_reference.ExportToWkt(),
            utm_clipped_vector_path)

        clipped_geometry_shapely_list = []
        temp_utm_clipped_vector = ogr.Open(utm_clipped_vector_path)
        temp_utm_clipped_layer = temp_utm_clipped_vector.GetLayer()
        for tmp_utm_feature in temp_utm_clipped_layer:
            tmp_utm_geometry = tmp_utm_feature.GetGeometryRef()
            shapely_geometry = shapely.wkb.loads(
                tmp_utm_geometry.ExportToWkb())
            if shapely_geometry.is_valid:
                clipped_geometry_shapely_list.append(shapely_geometry)
            tmp_utm_geometry = None
        temp_utm_clipped_layer = None
        temp_utm_clipped_vector = None
        landmass_shapely = shapely.ops.cascaded_union(
            clipped_geometry_shapely_list)
        clipped_geometry_shapely_list = None

        # load land geometry into shapely object
        landmass_shapely_prep = shapely.prepared.prep(landmass_shapely)

        # explode landmass into lines for easy intersection
        temp_polygon_segements_path = os.path.join(
            workspace_dir, 'polygon_segments.shp')
        temp_polygon_segments_vector = esri_shapefile_driver.CreateDataSource(
            temp_polygon_segements_path)
        temp_polygon_segments_layer = (
            temp_polygon_segments_vector.CreateLayer(
                os.path.splitext(temp_clipped_vector_path)[0],
                utm_spatial_reference, ogr.wkbLineString))
        temp_polygon_segments_defn = temp_polygon_segments_layer.GetLayerDefn()

        polygon_line_rtree = rtree.index.Index()
        polygon_line_index = []
        shapely_line_index = []
        line_id = 0
        for line in geometry_to_lines(landmass_shapely):
            segment_feature = ogr.Feature(temp_polygon_segments_defn)
            segement_geometry = ogr.Geometry(ogr.wkbLineString)
            segement_geometry.AddPoint(*line.coords[0])
            segement_geometry.AddPoint(*line.coords[1])
            segment_feature.SetGeometry(segement_geometry)
            temp_polygon_segments_layer.CreateFeature(segment_feature)

            if (line.bounds[0] == line.bounds[2] and
                    line.bounds[1] == line.bounds[3]):
                continue
            polygon_line_rtree.insert(line_id, line.bounds)
            line_id += 1
            polygon_line_index.append(segement_geometry)
            shapely_line_index.append(shapely.wkb.loads(
                segement_geometry.ExportToWkb()))

        temp_polygon_segments_layer.SyncToDisk()
        temp_polygon_segments_layer = None
        temp_polygon_segments_vector = None

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

        target_shore_point_vector = ogr.Open(
            target_fetch_point_vector_path, 1)
        target_shore_point_layer = target_shore_point_vector.GetLayer()
        target_shore_point_layer.CreateField(
            ogr.FieldDefn('REI', ogr.OFTReal))
        for ray_index in xrange(_N_FETCH_RAYS):
            compass_degree = int(ray_index * 360 / 16.)
            target_shore_point_layer.CreateField(
                ogr.FieldDefn('fdist_%d' % compass_degree, ogr.OFTReal))

        shore_point_logger = _make_logger_callback(
            "Wind exposure %.2f%% complete.", logger)
        # Iterate over every shore point
        for shore_point_feature in target_shore_point_layer:
            shore_point_logger(
                float(shore_point_feature.GetFID()) /
                target_shore_point_layer.GetFeatureCount())
            rei_value = 0.0
            # Iterate over every ray direction
            for sample_index in xrange(_N_FETCH_RAYS):
                compass_degree = int(sample_index * 360 / 16.)
                compass_theta = float(sample_index) / _N_FETCH_RAYS * 360
                rei_pct = shore_point_feature.GetField(
                    'REI_PCT%d' % int(compass_theta))
                rei_v = shore_point_feature.GetField(
                    'REI_V%d' % int(compass_theta))
                cartesian_theta = -(compass_theta - 90)

                # Determine the direction the ray will point
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

                # build ray geometry so we can intersect it later
                ray_geometry = ogr.Geometry(ogr.wkbLineString)
                ray_geometry.AddPoint(point_a_x, point_a_y)
                ray_geometry.AddPoint(point_b_x, point_b_y)

                # keep a shapely version of the ray so we can do fast intersection
                # with it and the entire landmass
                ray_point_origin_shapely = shapely.geometry.Point(
                    point_a_x, point_a_y)

                ray_length = 0.0
                if not landmass_shapely_prep.intersects(
                        ray_point_origin_shapely):
                    # the origin is in ocean

                    # This algorithm searches for intersections, if one is found
                    # the ray updates and a smaller intersection set is determined
                    # by experimentation I've found this is significant, but not
                    # an order of magnitude, faster than looping through all
                    # original possible intersections.  Since this algorithm
                    # will be run for a long time, it's worth the additional
                    # complexity
                    tested_indexes = set()
                    while True:
                        intersection = False
                        ray_envelope = ray_geometry.GetEnvelope()
                        for poly_line_index in polygon_line_rtree.intersection(
                                [ray_envelope[i] for i in [0, 2, 1, 3]]):
                            if poly_line_index in tested_indexes:
                                continue
                            tested_indexes.add(poly_line_index)
                            line_segment = (
                                polygon_line_index[poly_line_index])
                            if ray_geometry.Intersects(line_segment):
                                # if the ray intersects the poly line, test if
                                # the intersection is closer than any known
                                # intersection so far
                                intersection_point = ray_geometry.Intersection(
                                    line_segment)
                                # offset the dist with smallest_feature_size
                                # update the endpoint of the ray
                                ray_geometry = ogr.Geometry(ogr.wkbLineString)
                                ray_geometry.AddPoint(point_a_x, point_a_y)
                                ray_geometry.AddPoint(
                                    intersection_point.GetX(),
                                    intersection_point.GetY())
                                intersection = True
                                break
                        if not intersection:
                            break
                    # when we get here `min_point` and `ray_length` are the
                    # minimum intersection points for the ray and the landmass
                    ray_feature = ogr.Feature(temp_fetch_rays_defn)
                    ray_length = ray_geometry.Length()
                    ray_feature.SetField('fetch_dist', ray_length)
                    ray_feature.SetGeometry(ray_geometry)
                    temp_fetch_rays_layer.CreateFeature(ray_feature)
                shore_point_feature.SetField(
                    'fdist_%d' % compass_degree, ray_length)
                ray_feature = None
                ray_geometry = None
                rei_value += ray_length * rei_pct * rei_v
            shore_point_feature.SetField('REI', rei_value)
            target_shore_point_layer.SetFeature(shore_point_feature)

        target_shore_point_layer.SyncToDisk()
        target_shore_point_layer = None
        target_shore_point_vector = None
        temp_fetch_rays_layer.SyncToDisk()
        temp_fetch_rays_layer = None
        temp_fetch_rays_vector = None
    except Exception as e:
        traceback.print_exc()
        raise


def calculate_relief(
        base_shore_point_vector_path, global_dem_path, workspace_dir,
        target_relief_point_vector_path):
    """Calculate DEM relief as average coastal land area within 5km.

    Parameters:
    TODO fill in

    Returns:
        None.
    """
    try:
        logger = logging.getLogger('ipbes-cv.calculate_relief')
        if not os.path.exists(os.path.dirname(
                target_relief_point_vector_path)):
            os.makedirs(
                os.path.dirname(target_relief_point_vector_path))
        if not os.path.exists(workspace_dir):
            os.makedirs(workspace_dir)
        if os.path.exists(target_relief_point_vector_path):
            os.remove(target_relief_point_vector_path)

        base_ref_wkt = pygeoprocessing.get_vector_info(
            base_shore_point_vector_path)['projection']
        base_spatial_reference = osr.SpatialReference()
        base_spatial_reference.ImportFromWkt(base_ref_wkt)

        # reproject base_shore_point_vector_path to utm coordinates
        base_shore_info = pygeoprocessing.get_vector_info(
            base_shore_point_vector_path)
        base_shore_bounding_box = base_shore_info['bounding_box']

        utm_spatial_reference = get_utm_spatial_reference(
            base_shore_info['bounding_box'])
        base_spatial_reference = osr.SpatialReference()
        base_spatial_reference.ImportFromWkt(base_shore_info['projection'])

        pygeoprocessing.reproject_vector(
            base_shore_point_vector_path, utm_spatial_reference.ExportToWkt(),
            target_relief_point_vector_path)

        utm_bounding_box = pygeoprocessing.get_vector_info(
            target_relief_point_vector_path)['bounding_box']

        # extend bounding box for max fetch distance
        utm_bounding_box = [
            utm_bounding_box[0] - _MAX_FETCH_DISTANCE,
            utm_bounding_box[1] - _MAX_FETCH_DISTANCE,
            utm_bounding_box[2] + _MAX_FETCH_DISTANCE,
            utm_bounding_box[3] + _MAX_FETCH_DISTANCE]

        # get lat/lng bounding box of utm projected coordinates

        # get global polygon clip of that utm box
        # transform local box back to lat/lng -> global clipping box
        lat_lng_clipping_box = pygeoprocessing.transform_bounding_box(
            utm_bounding_box, utm_spatial_reference.ExportToWkt(),
            base_spatial_reference.ExportToWkt(), edge_samples=11)
        if (base_shore_bounding_box[0] > 0 and
                lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[2] += 360
        elif (base_shore_bounding_box[0] < 0 and
              lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
            lat_lng_clipping_box[0] -= 360
        elif base_shore_bounding_box == [0, 0, 0, 0]:
            # this case guards for an empty shore point in case there are
            # very tiny islands or such
            lat_lng_clipping_box = (0, 0, 0, 0)

        clipped_lat_lng_dem_path = os.path.join(
            workspace_dir, 'clipped_lat_lng_dem.tif')

        target_pixel_size = pygeoprocessing.get_raster_info(
            global_dem_path)['pixel_size']
        pygeoprocessing.warp_raster(
            global_dem_path, target_pixel_size, clipped_lat_lng_dem_path,
            'bilinear', target_bb=base_shore_bounding_box)
        clipped_utm_dem_path = os.path.join(
            workspace_dir, 'clipped_utm_dem.tif')
        pygeoprocessing.warp_raster(
            clipped_lat_lng_dem_path, target_pixel_size, clipped_utm_dem_path,
            'bilinear', target_sr_wkt=utm_spatial_reference.ExportToWkt())
        # mask out all DEM < 0 to 0
        # convolve over a 5km radius
        # can i aggregate by point? if so, that's the relief
    except Exception as e:
        traceback.print_exc()
        raise




def create_shore_points(
        sample_grid_vector_path, grid_id, landmass_bounding_rtree_path,
        landmass_vector_path, wwiii_vector_path, wwiii_rtree_path,
        smallest_feature_size,
        workspace_dir, target_shore_point_vector_path):
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
        wwiii_rtree_path (string): path to an rtree index that has
            the points of `wwiii_vector_path` indexed.
        smallest_feature_size (float): smallest feature size to grid a shore
            point on.
        workspace_dir (string): path to a directory that can be created
            during run to hold temporary files.  Will be deleted on successful
            function completion.
        target_shore_point_vector_path (string): path to a point vector that
            will be created and contain points on the shore of the landmass.

    Returns:
        None.
    """
    logger = logging.getLogger('ipbes-cv.create_shore_points')
    logger.info("Creating shore points for grid %s", grid_id)
    # create the spatial reference from the base vector
    landmass_spatial_reference = osr.SpatialReference()
    landmass_spatial_reference.ImportFromWkt(
        pygeoprocessing.get_vector_info(landmass_vector_path)['projection'])

    if os.path.exists(workspace_dir):
        shutil.rmtree(workspace_dir)
    os.makedirs(workspace_dir)

    lat_lng_clipped_vector_path = os.path.join(
        workspace_dir, 'clipped_geometry_lat_lng.shp')
    grid_raster_path = os.path.join(workspace_dir, 'grid.tif')
    convolution_raster_path = os.path.join(
        workspace_dir, 'convolution.tif')
    utm_clipped_vector_path = os.path.join(
        workspace_dir, 'clipped_geometry_utm.shp')
    shore_kernel_path = os.path.join(
        workspace_dir, 'shore_kernel.tif')
    shore_raster_path = os.path.join(
        workspace_dir, 'shore_raster.tif')

    for path in [target_shore_point_vector_path,
                 lat_lng_clipped_vector_path,
                 grid_raster_path]:
        if os.path.exists(path):
            os.remove(path)

    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")

    # this will hold the clipped landmass geometry
    lat_lng_clipped_vector = esri_shapefile_driver.CreateDataSource(
        lat_lng_clipped_vector_path)
    lat_lng_clipped_layer = (
        lat_lng_clipped_vector.CreateLayer(
            os.path.splitext(lat_lng_clipped_vector_path)[0],
            landmass_spatial_reference, ogr.wkbPolygon))
    lat_lng_clipped_defn = lat_lng_clipped_layer.GetLayerDefn()

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
    grid_geometry_ref = grid_feature.GetGeometryRef()
    grid_shapely = shapely.wkb.loads(grid_geometry_ref.ExportToWkb())

    landmass_vector_rtree = rtree.index.Index(
        os.path.splitext(landmass_bounding_rtree_path)[0])

    # project global polygon clip to UTM
    # transform lat/lng box to utm -> local box
    utm_spatial_reference = get_utm_spatial_reference(grid_shapely.bounds)
    utm_bounding_box = pygeoprocessing.transform_bounding_box(
        grid_shapely.bounds, landmass_spatial_reference.ExportToWkt(),
        utm_spatial_reference.ExportToWkt(), edge_samples=11)

    # add a pixel buffer so we clip land that's a little outside the grid
    pixel_buffer = 1
    utm_bounding_box[0] -= pixel_buffer * smallest_feature_size
    utm_bounding_box[1] -= pixel_buffer * smallest_feature_size
    utm_bounding_box[2] += pixel_buffer * smallest_feature_size
    utm_bounding_box[3] += pixel_buffer * smallest_feature_size

    # transform local box back to lat/lng -> global clipping box
    lat_lng_clipping_box = pygeoprocessing.transform_bounding_box(
        utm_bounding_box, utm_spatial_reference.ExportToWkt(),
        landmass_spatial_reference.ExportToWkt(), edge_samples=11)
    # see if we're wrapped on the dateline
    if (grid_shapely.bounds[0] > 0 and
            lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
        lat_lng_clipping_box[2] += 360
    elif (grid_shapely.bounds[0] < 0 and
          lat_lng_clipping_box[0] > lat_lng_clipping_box[2]):
        lat_lng_clipping_box[0] -= 360
    lat_lng_clipping_shapely = shapely.geometry.box(*lat_lng_clipping_box)

    # clip global polygon to utm clipping box
    logger.info(
        "clip global polygon to utm clipping box for grid %s", grid_id)
    for feature_id in landmass_vector_rtree.intersection(
            lat_lng_clipping_box):
        base_feature = landmass_layer.GetFeature(feature_id)
        base_geometry = base_feature.GetGeometryRef()
        base_shapely = shapely.wkb.loads(base_geometry.ExportToWkb())
        base_geometry = None
        intersection_shapely = lat_lng_clipping_shapely.intersection(
            base_shapely)
        try:
            target_geometry = ogr.CreateGeometryFromWkt(
                intersection_shapely.wkt)
            target_feature = ogr.Feature(lat_lng_clipped_defn)
            target_feature.SetGeometry(target_geometry)
            lat_lng_clipped_layer.CreateFeature(target_feature)
            target_feature = None
            target_geometry = None
        except Exception:
            logger.warn(
                "Couldn't process this intersection %s",
                intersection_shapely)
    lat_lng_clipped_layer.SyncToDisk()
    lat_lng_clipped_layer = None
    lat_lng_clipped_vector = None

    # create grid for underlying local utm box
    pygeoprocessing.reproject_vector(
        lat_lng_clipped_vector_path, utm_spatial_reference.ExportToWkt(),
        utm_clipped_vector_path)

    byte_nodata = 255

    pygeoprocessing.create_raster_from_vector_extents(
        utm_clipped_vector_path,
        grid_raster_path, (
            smallest_feature_size / 2.0, -smallest_feature_size / 2.0),
        gdal.GDT_Byte, byte_nodata, fill_value=0)

    # rasterize utm global clip to grid
    pygeoprocessing.rasterize(
        utm_clipped_vector_path, grid_raster_path, [1], None)

    # grid shoreline from raster
    make_shore_kernel(shore_kernel_path)
    pygeoprocessing.convolve_2d(
        (grid_raster_path, 1), (shore_kernel_path, 1),
        convolution_raster_path, target_datatype=gdal.GDT_Byte)

    temp_grid_nodata = pygeoprocessing.get_raster_info(
        grid_raster_path)['nodata'][0]

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
        [(convolution_raster_path, 1)], _shore_mask_op,
        shore_raster_path, gdal.GDT_Byte, byte_nodata)

    shore_geotransform = pygeoprocessing.get_raster_info(
        shore_raster_path)['geotransform']

    utm_to_base_transform = osr.CoordinateTransformation(
        utm_spatial_reference, landmass_spatial_reference)

    # rtree index loads without the extension
    wwiii_rtree_base_path = os.path.splitext(
        wwiii_rtree_path)[0]
    wwiii_rtree = rtree.index.Index(wwiii_rtree_base_path)
    wwiii_field_lookup = {}

    logger.info(
        "Interpolating shore points with Wave Watch III data for grid %s",
        grid_id)
    feature_lookup = {}
    for offset_info, data_block in pygeoprocessing.iterblocks(
            shore_raster_path):
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
            # make sure shore point is within the bounding box of the gri
            if grid_geometry_ref.Contains(shore_point_geometry):
                shore_point_feature = ogr.Feature(target_shore_point_defn)
                shore_point_feature.SetGeometry(shore_point_geometry)

                # get the nearest wave watch III points from the shore point
                nearest_points = list(wwiii_rtree.nearest(
                    (shore_point_geometry.GetX(),
                     shore_point_geometry.GetY(),
                     shore_point_geometry.GetX(),
                     shore_point_geometry.GetY()), 3))[0:3]

                # create placeholders for point geometry and field values
                wwiii_points = numpy.empty((3, 2))
                wwiii_values = numpy.empty((3, len(field_names)))
                for fid_index, fid in enumerate(nearest_points):
                    wwiii_feature = wwiii_layer.GetFeature(fid)
                    wwiii_geometry = wwiii_feature.GetGeometryRef()
                    wwiii_points[fid_index] = numpy.array(
                        [wwiii_geometry.GetX(), wwiii_geometry.GetY()])
                    try:
                        wwiii_values[fid_index] = wwiii_field_lookup[fid]
                    except KeyError:
                        wwiii_field_lookup[fid] = numpy.array(
                            [float(wwiii_feature.GetField(field_name))
                             for field_name in field_names])
                        wwiii_values[fid_index] = wwiii_field_lookup[fid]

                distance = numpy.linalg.norm(
                    wwiii_points - numpy.array(
                        (shore_point_geometry.GetX(),
                         shore_point_geometry.GetY())))

                # make sure we're within a valid data distance
                if distance > _MAX_WWIII_DISTANCE:
                    continue

                wwiii_values *= distance
                wwiii_values = numpy.mean(wwiii_values, axis=0)
                wwiii_values /= numpy.sum(distance)

                for field_name_index, field_name in enumerate(field_names):
                    shore_point_feature.SetField(
                        field_name, wwiii_values[field_name_index])

                target_shore_point_layer.CreateFeature(shore_point_feature)
                shore_point_feature = None
    del feature_lookup
    logger.info("All done with shore points for grid %s", grid_id)


def grid_edges_of_vector(
        base_bounding_box, base_vector_path,
        base_feature_bounding_box_rtree_path, target_grid_vector_path,
        target_grid_size):
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
    logger = logging.getLogger('ipbes-cv.grid_edges_of_vector')
    logger.info("Building global grid.")
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
        'Cell coverage %.2f%% complete', logger)

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


def get_utm_spatial_reference(bounding_box):
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


def build_feature_bounding_box_rtree(vector_path, target_rtree_path):
    """Builds an r-tree index of the global feature envelopes.

    Parameter:
        vector_path (string): path to vector to build bounding box index for
        target_rtree_path (string): path to ".dat" file to store the saved
            r-tree.  A ValueError is raised if this file already exists

    Returns:
        None.
    """
    # the input path has a .dat extension, but the `rtree` package only uses
    # the basename.  It's a quirk of the library, so we'll deal with it by
    # cutting off the extension.
    logger = logging.getLogger('ipbes-cv.build_feature_bounding_box_rtree')
    global_feature_index_base = os.path.splitext(
        target_rtree_path)[0]
    logger.info("Building rtree index at %s", global_feature_index_base)
    if os.path.exists(target_rtree_path):
        raise ValueError("rtree storage path %s already exists.")
    global_feature_index = rtree.index.Index(global_feature_index_base)

    global_vector = ogr.Open(vector_path)
    global_layer = global_vector.GetLayer()
    n_features = global_layer.GetFeatureCount()

    logger_callback = _make_logger_callback(
        'rTree construction %.2f%% complete', logger)

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


def make_shore_kernel(kernel_path):
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


def _make_logger_callback(message, logger):
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
                logger.info(message, proportion_complete * 100)
                logger_callback.last_time = current_time
                logger_callback.total_time += current_time
        except AttributeError:
            logger_callback.last_time = time.time()
            logger_callback.total_time = 0.0

    return logger_callback


def build_wwiii_rtree(wwiii_vector_path, wwiii_rtree_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    wwiii_rtree = rtree.index.Index(os.path.splitext(wwiii_rtree_path)[0])

    wwiii_vector = ogr.Open(wwiii_vector_path)
    wwiii_layer = wwiii_vector.GetLayer()
    for wwiii_feature in wwiii_layer:
        wwiii_geometry = wwiii_feature.GetGeometryRef()
        wwiii_x = wwiii_geometry.GetX()
        wwiii_y = wwiii_geometry.GetY()
        wwiii_rtree.insert(
            wwiii_feature.GetFID(), (wwiii_x, wwiii_y, wwiii_x, wwiii_y))
    wwiii_layer = None
    wwiii_vector = None


def merge_vectors(
        base_vector_path_list, target_spatial_reference_wkt,
        target_merged_vector_path, field_list_to_copy):
    """Merge all the vectors in the `base_vector_path_list`.

    Parameters:
        base_vector_path_list (list): a list of OGR DataSources.  Should
            all be single layer identical feature definitions.
        target_merged_vector_path (string): path to desired output vector.
            will contain a single layer with all merged features from the
            base vectors.
        field_list_to_copy (string): list of field names to copy from the base
            vectors to the target vectors

    Returns:
        None
    """
    target_spatial_reference = osr.SpatialReference()
    target_spatial_reference.ImportFromWkt(target_spatial_reference_wkt)

    esri_driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(target_merged_vector_path):
        esri_driver.DeleteDataSource(target_merged_vector_path)
    base_vector = ogr.Open(base_vector_path_list[0])
    base_layer = base_vector.GetLayer()
    base_layer_defn = base_layer.GetLayerDefn()

    target_vector = esri_driver.CreateDataSource(target_merged_vector_path)
    target_layer = target_vector.CreateLayer(
        target_merged_vector_path, srs=target_spatial_reference,
        geom_type=base_layer.GetGeomType())
    for field_name in field_list_to_copy:
        target_layer.CreateField(
            base_layer_defn.GetFieldDefn(
                base_layer_defn.GetFieldIndex(field_name)))

    base_layer = None
    base_vector = None

    for base_vector_path in base_vector_path_list:
        base_vector = ogr.Open(base_vector_path)
        base_layer = base_vector.GetLayer()

        base_spatial_reference = base_layer.GetSpatialRef()
        base_to_target_transform = osr.CoordinateTransformation(
            base_spatial_reference, target_spatial_reference)

        for feature in base_layer:
            target_feat = ogr.Feature(target_layer.GetLayerDefn())
            target_geometry = feature.GetGeometryRef().Clone()
            target_geometry.Transform(base_to_target_transform)
            target_feat.SetGeometry(target_geometry)
            for field_name in field_list_to_copy:
                target_feat.SetField(field_name, feature.GetField(field_name))
            target_layer.CreateFeature(target_feat)
            target_feat = None
    target_layer.SyncToDisk()
    target_layer = None
    target_vector = None


def geometry_to_lines(geometry):
    if geometry.type == 'Polygon':
        return polygon_to_lines(geometry)
    elif geometry.type == 'MultiPolygon':
        line_list = []
        for geom in geometry.geoms:
            line_list.extend(geometry_to_lines(geom))
        return line_list
    else:
        return []


def polygon_to_lines(geometry):
    """Return a list of shapely lines given higher order shapely geometry."""
    line_list = []
    last_point = geometry.exterior.coords[0]
    for point in geometry.exterior.coords[1::]:
        if point == last_point:
            continue
        line_list.append(shapely.geometry.LineString([last_point, point]))
        last_point = point
    line_list.append(shapely.geometry.LineString([
        last_point, geometry.exterior.coords[0]]))
    for interior in geometry.interiors:
        last_point = interior.coords[0]
        for point in interior.coords[1::]:
            if point == last_point:
                continue
            line_list.append(shapely.geometry.LineString([last_point, point]))
            last_point = point
        line_list.append(shapely.geometry.LineString([
            last_point, interior.coords[0]]))
    return line_list

if __name__ == '__main__':
    main()
