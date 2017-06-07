"""IPBES global coastal vulnerability calculation."""
import time
import os
import math
import logging

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
_WGS84_GRID_SIZE = 10.0

# Once global grid is cut, it is reprojected into UTM with this underlying
# square grid cell size
_UTM_GRID_SIZE = 250

_GLOBAL_GRID_VECTOR_FILE_PATTERN = 'global_grid.shp'
_GLOBAL_FEATURE_INDEX_FILE_PATTERN = 'global_feature_index'


def main():
    """Entry point."""
    if not os.path.exists(_TARGET_WORKSPACE):
        os.makedirs(_TARGET_WORKSPACE)

    global_feature_index_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_FEATURE_INDEX_FILE_PATTERN)
    LOGGER.info("Building rtree index.")
    _build_feature_bounding_box_rtree(
        _GLOBAL_POLYGON_PATH, global_feature_index_path)

    global_grid_vector_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_GRID_VECTOR_FILE_PATTERN)
    LOGGER.info("Building global grid.")
    _grid_edges_of_vector(
        _GLOBAL_BOUNDING_BOX_WGS84, _GLOBAL_POLYGON_PATH,
        global_feature_index_path, global_grid_vector_path,
        _WGS84_GRID_SIZE)


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

    target_feature_rtree_index = rtree.index.Index(
        base_feature_bounding_box_rtree_path)

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
