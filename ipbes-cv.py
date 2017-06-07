"""IPBES global coastal vulnerability calculation."""
import os
import math

from osgeo import osr
from osgeo import ogr
import rtree
import shapely
import shapely.wkt

_TARGET_WORKSPACE = "ipbes_cv_workspace"

_GLOBAL_POLYGON_PATH = r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\Land\global_polygon.shp"

# The global bounding box to do the entire analysis
# This range was roughly picked to avoid the poles
# [min_lat, min_lng, max_lat, max_lng]
_GLOBAL_BOUNDING_BOX_WGS84 = [-180, -60, 180, 60]

# This is the lat/lng grid size to slice the runs into, annoying since it's
# lat/lng, but if you have a better idea lets hear it.
_WGS84_GRID_SIZE = 0.5

# Once global grid is cut, it is reprojected into UTM with this underlying
# square grid cell size
_UTM_GRID_SIZE = 250

_GLOBAL_GRID_VECTOR_FILE_PATTERN = 'global_grid.shp'
_GLOBAL_FEATURE_INDEX_FILE_PATTERN = 'global_feature_index'

def main():
    """Entry point."""
    if not os.path.exists(_TARGET_WORKSPACE):
        os.makedirs(_TARGET_WORKSPACE)

    global_vector = ogr.Open(_GLOBAL_POLYGON_PATH)
    global_layer = global_vector.GetLayer()

    global_feature_index_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_FEATURE_INDEX_FILE_PATTERN)
    index_pre_built = os.path.exists(global_feature_index_path + '.dat')
    global_feature_index = rtree.index.Index(global_feature_index_path)

    if not index_pre_built:
        print 'building global index'
        for global_feature in global_layer:
            feature_geometry = global_feature.GetGeometryRef()
            # minx maxx miny maxy
            feature_envelope = feature_geometry.GetEnvelope()
            # left, bottom, right, top
            global_feature_index.insert(
                global_feature.GetFID(), (
                    feature_envelope[0], feature_envelope[2],
                    feature_envelope[1], feature_envelope[3]))
    else:
        print 'global index pre-built'

    print global_feature_index.bounds

    n_rows = int((
        _GLOBAL_BOUNDING_BOX_WGS84[3]-_GLOBAL_BOUNDING_BOX_WGS84[1]) / float(
            _WGS84_GRID_SIZE))
    n_cols = int((
        _GLOBAL_BOUNDING_BOX_WGS84[2]-_GLOBAL_BOUNDING_BOX_WGS84[0]) / float(
            _WGS84_GRID_SIZE))

    global_grid_vector_path = os.path.join(
        _TARGET_WORKSPACE, _GLOBAL_GRID_VECTOR_FILE_PATTERN)

    if os.path.exists(global_grid_vector_path):
        os.remove(global_grid_vector_path)


    # set up the shapefile driver
    esri_shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")

    # create the data source
    global_grid_vector = esri_shapefile_driver.CreateDataSource(
        global_grid_vector_path)

    # create the spatial reference, WGS84
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)

    global_grid_layer = global_grid_vector.CreateLayer(
        os.path.splitext(_GLOBAL_GRID_VECTOR_FILE_PATTERN)[0],
        wgs84_srs, ogr.wkbPolygon)

    global_grid_defn = global_grid_layer.GetLayerDefn()

    prepared_geometry = {}
    for cell_index in xrange(n_rows * n_cols):
        row_index = cell_index / n_cols
        col_index = cell_index % n_cols
        # xmin, ymin xmax, ymax
        local_bounding_box_wgs84 = [
            _GLOBAL_BOUNDING_BOX_WGS84[0]+col_index*_WGS84_GRID_SIZE,
            _GLOBAL_BOUNDING_BOX_WGS84[1]+row_index*_WGS84_GRID_SIZE,
            _GLOBAL_BOUNDING_BOX_WGS84[0]+(col_index+1)*_WGS84_GRID_SIZE,
            _GLOBAL_BOUNDING_BOX_WGS84[1]+(row_index+1)*_WGS84_GRID_SIZE]

        intersections = list(
            global_feature_index.intersection(local_bounding_box_wgs84))

        if len(intersections) == 0:
            continue
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(
            local_bounding_box_wgs84[0], local_bounding_box_wgs84[1])
        ring.AddPoint(
            local_bounding_box_wgs84[2], local_bounding_box_wgs84[1])
        ring.AddPoint(
            local_bounding_box_wgs84[2], local_bounding_box_wgs84[3])
        ring.AddPoint(
            local_bounding_box_wgs84[0], local_bounding_box_wgs84[3])
        ring.AddPoint(
            local_bounding_box_wgs84[0], local_bounding_box_wgs84[1])
        cell_geometry = ogr.Geometry(ogr.wkbPolygon)
        cell_geometry.AddGeometry(ring)
        cell_geometry_shapely = shapely.wkt.loads(
            cell_geometry.ExportToWkt())

        print 'testing %s' % intersections
        for fid in intersections:
            global_feature = global_layer.GetFeature(fid)
            feature_geometry = global_feature.GetGeometryRef()
            if fid not in prepared_geometry:
                prepared_geometry[fid] = shapely.wkt.loads(
                    feature_geometry.ExportToWkt())
            prepared_shapely_feature = prepared_geometry[fid]
            if prepared_shapely_feature.intersects(cell_geometry_shapely):
                # add new geom to layer
                grid_feature = ogr.Feature(global_grid_defn)
                grid_feature.SetGeometry(cell_geometry)
                global_grid_layer.CreateFeature(grid_feature)
                grid_feature = None
                break

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

if __name__ == '__main__':
    main()
