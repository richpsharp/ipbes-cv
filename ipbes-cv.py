"""IPBES global coastal vulnerability calculation."""
import shutil
import time
import os
import math
import logging

import numpy
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import rtree
import shapely
import shapely.wkb
import pygeoprocessing

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')

_TARGET_WORKSPACE = "ipbes_cv_workspace"

_GLOBAL_POLYGON_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\Land\global_polygon.shp"

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


def _make_task(func, args, expected_output_path_list, dependant_task_list):
    """Make a task.

    Parameters:
        func (function): a function that takes the argument list `args`
        args (tuple): a list of arguments to pass to `func`.
        expected_output_path_list (list): a list of strings representing
            expected file path outputs.
        dependant_task_list (list of function): a list of other functions
            generated by `_make_task` that contain an .execute() method.
            These are all invoked before `func(args)` is invoked.

    Returns:
        a function that contains an `execute` method.  When invoked it will
        first invoke the `execute` methods in `dependant_task_list`.  Once
        those are complete, a check to see if the file paths in
        `expected_output_path_list` exist on disk.  If not, func(args) is
        executed.
    """
    def execute():
        """Invoke this method when ready to execute task."""
        if all(os.path.exists(p) for p in expected_output_path_list):
            LOGGER.info(
                "All expected files exist for %s so not executing",
                func.__name__)
            return

        for task in dependant_task_list:
            task()
        func(*args)

    return execute


def main():
    """Entry point."""
    if not os.path.exists(_TARGET_WORKSPACE):
        os.makedirs(_TARGET_WORKSPACE)

    global_rtree_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_FEATURE_INDEX_FILE_PATTERN)

    build_feature_bounding_box_rtree_task = _make_task(
        _build_feature_bounding_box_rtree, (
            _GLOBAL_POLYGON_PATH, global_rtree_path),
        [global_rtree_path], [])

    global_grid_vector_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_GRID_VECTOR_FILE_PATTERN)

    grid_edges_of_vector_task = _make_task(
        _grid_edges_of_vector, (
            _GLOBAL_BOUNDING_BOX_WGS84, _GLOBAL_POLYGON_PATH,
            global_rtree_path, global_grid_vector_path,
            _WGS84_GRID_SIZE),
        [global_grid_vector_path], [build_feature_bounding_box_rtree_task])

    grid_edges_of_vector_task()

    LOGGER.info("Find shore points in clipped polygon.")
    grid_id = 122
    grid_point_path = os.path.join(
        _TARGET_WORKSPACE, _GRID_POINT_FILE_PATTERN % (grid_id))
    smallest_feature_size = [250, -250]
    temp_workspace = os.path.join(
        _TARGET_WORKSPACE, 'grid_%d' % grid_id)
    _create_shore_points(
        global_grid_vector_path, grid_id, global_rtree_path,
        _GLOBAL_POLYGON_PATH, smallest_feature_size, temp_workspace,
        grid_point_path)


def _create_shore_points(
        grid_vector_path, grid_id, rtree_path, base_vector_path,
        smallest_feature_size, temp_workspace,
        target_sample_point_vector_path):
    """Create points that lie on the coast line.

    Parameters:
        grid_vector_path (string): path to vector containing grids
        grid_id (integer): feature ID in `grid_vector_path`'s layer to operate
            on.
        rtree_path (string): path to an rtree index that has bounding box
            indexes of polygons that might intersect a grid in question.
        base_vector_path (string): path to polygon vector that we're analyzing
            over.
        smallest_feature_size (float): smallest feature size to grid a shore
            point on.
        temp_workspace (string): path to a directory that can be created
            during run to hold temporary files.  Will be deleted on successful
            function completion.
        target_sample_point_vector_path (string): path to a point vector that
            samples the edges of the polygon

    Returns:
        None.
    """
    # create the spatial reference from the base vector
    base_spatial_reference = osr.SpatialReference()
    base_spatial_reference.ImportFromWkt(
        pygeoprocessing.get_vector_info(base_vector_path)['projection'])

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

    for path in [target_sample_point_vector_path,
                 temp_clipped_vector_path,
                 temp_grid_raster_path]:
        if os.path.exists(path):
            os.remove(path)

    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")

    temp_clipped_vector = esri_shapefile_driver.CreateDataSource(
        temp_clipped_vector_path)
    temp_clipped_layer = (
        temp_clipped_vector.CreateLayer(
            os.path.splitext(temp_clipped_vector_path)[0],
            base_spatial_reference, ogr.wkbPolygon))

    target_sample_point_vector = esri_shapefile_driver.CreateDataSource(
        target_sample_point_vector_path)
    target_sample_point_layer = target_sample_point_vector.CreateLayer(
        os.path.splitext(target_sample_point_vector_path)[0],
        base_spatial_reference, ogr.wkbPoint)

    target_sample_point_defn = target_sample_point_layer.GetLayerDefn()

    base_vector = ogr.Open(base_vector_path)
    base_layer = base_vector.GetLayer()

    grid_vector = ogr.Open(grid_vector_path)
    grid_layer = grid_vector.GetLayer()
    grid_feature = grid_layer.GetFeature(grid_id)

    base_vector_rtree = rtree.index.Index(
        os.path.splitext(rtree_path)[0])

    grid_shapely = shapely.wkb.loads(
        grid_feature.GetGeometryRef().ExportToWkb())

    # project global polygon clip to UTM
    # transform lat/lng box to utm -> local box
    utm_spatial_reference = _get_utm_spatial_reference(grid_shapely.bounds)
    utm_bounding_box = pygeoprocessing.transform_bounding_box(
        grid_shapely.bounds, base_spatial_reference.ExportToWkt(),
        utm_spatial_reference.ExportToWkt(), edge_samples=11)

    # transform local box back to lat/lng -> global clipping box
    base_clipping_box = pygeoprocessing.transform_bounding_box(
        utm_bounding_box, utm_spatial_reference.ExportToWkt(),
        base_spatial_reference.ExportToWkt(), edge_samples=11)
    base_clipping_shapely = shapely.geometry.box(*base_clipping_box)

    # clip global polygon to global clipping box
    for feature_id in base_vector_rtree.intersection(base_clipping_box):
        base_feature = base_layer.GetFeature(feature_id)
        base_shapely = shapely.wkb.loads(
            base_feature.GetGeometryRef().ExportToWkb())
        intersection_shapely = base_clipping_shapely.intersection(base_shapely)
        try:
            target_geometry = ogr.CreateGeometryFromWkt(
                intersection_shapely.wkt)
            target_feature = ogr.Feature(target_sample_point_defn)
            target_feature.SetGeometry(target_geometry)
            temp_clipped_layer.CreateFeature(target_feature)
            target_feature = None
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
        temp_grid_raster_path, [i / 2.0 for i in smallest_feature_size],
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

    LOGGER.warn('TODO: delete temporary files')


def _grid_edges_of_vector(
        base_bounding_box, base_vector_path,
        base_feature_bounding_box_rtree_path, target_grid_vector_path,
        target_grid_size):
    """Builds a sparse grid covering the edges of the base polygons.

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
        vector_path, target_rtree_path):
    """Builds an r-tree index of the global feature envelopes.

    Parameter:
        vector_path (string): path to vector to build bounding box index for
        target_rtree_path (string): path to ".dat" file to store the saved
            r-tree.

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
