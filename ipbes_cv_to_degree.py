"""IPBES global sea level rise post SLR analytics https://docs.google.com/document/d/1yX7aIOehd04v3CZt9FDgRsKsUjR8DSlkpv_XhKDXpC8/edit#heading=h.6fz9zq50ghb3."""
import pickle
import os
import logging
import taskgraph
import numpy
import collections

from osgeo import osr
from osgeo import ogr
from osgeo import gdal

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')
LOGGER.setLevel(logging.DEBUG)

POSSIBLE_DROPBOX_LOCATIONS = [
    r'D:\Dropbox',
    r'C:\Users\Rich\Dropbox',
    r'C:\Users\rpsharp\Dropbox',
    r'E:\Dropbox']

LOGGER.info("checking dropbox locations")
for path in POSSIBLE_DROPBOX_LOCATIONS:
    print path
    if os.path.exists(path):
        BASE_DROPBOX_DIR = path
        break
LOGGER.info("found %s", BASE_DROPBOX_DIR)

GLOBAL_CV_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_cv_results', 'global_cv_rt_ssp_recalc.shp')

TARGET_SUMMARY_DEGREE_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_cv_results', 'sum_to_degree')

SLR_LIST_PICKLE_PATH = 'slr.pickle'


def sort_points(summary_field_list, pickle_path):
    if os.path.exists(pickle_path):
        return
    vector = gdal.OpenEx(GLOBAL_CV_VECTOR_PATH, gdal.OF_VECTOR)
    layer = vector.GetLayer()

    summary_grid_index_path = collections.defaultdict(list)

    while True:
        feature = layer.GetNextFeature()
        if not feature:
            break
        centroid = feature.GetGeometryRef().Centroid()
        grid_index = (int(centroid.GetX()), int(centroid.GetY()))
        summary_fields = dict([
            (x, feature.GetField(x)) for x in summary_field_list])
        summary_grid_index_path[grid_index].append(summary_fields)

    with open(pickle_path, 'wb') as pickle_file:
        pickle.dump(summary_grid_index_path, pickle_file)


def main():
    """Entry point."""
    try:
        os.makedirs(TARGET_SUMMARY_DEGREE_DIR)
    except OSError:
        pass

    cur_ssp_list = ['cur', 'ssp1', 'ssp3', 'ssp5']

    summary_field_list = (
        ['Rt', 'rnom_cur'] +
        ['rt_ssp%d' % x for x in [1, 3, 5]] +
        ['crt_ssp%d' % x for x in [1, 3, 5]] +
        ['rpop_%s' % x for x in cur_ssp_list] +
        ['rpov_%s' % x for x in cur_ssp_list] +
        ['rage_%s' % x for x in cur_ssp_list] +
        ['crnom_ssp%s' % x for x in [1, 3, 5]] +
        ['rpop_ssp%d' % x for x in [1, 3, 5]] +
        ['rpov_ssp%d' % x for x in [1, 3, 5]] +
        ['rage_ssp%d' % x for x in [1, 3, 5]])

    task_graph = taskgraph.TaskGraph(
        os.path.join(TARGET_SUMMARY_DEGREE_DIR, 'taskgraph_cache'), -1)

    summary_grid_index_path_path = os.path.join(
        TARGET_SUMMARY_DEGREE_DIR, 'summary_grid_index_path.pickle')

    grid_index_task = task_graph.add_task(
        func=sort_points,
        args=(summary_field_list, summary_grid_index_path_path,),
        target_path_list=[summary_grid_index_path_path])

    print 'join grid index'
    grid_index_task.join()

    for summary_field in summary_field_list:
        print 'summary %s' % summary_field
        summary_raster_path = os.path.join(
            TARGET_SUMMARY_DEGREE_DIR, '%s.tif' % summary_field)
        wgs84_sr = osr.SpatialReference()
        wgs84_sr.ImportFromEPSG(4326)
        driver = gdal.GetDriverByName('GTiff')
        print 'create raster'
        summary_raster = driver.Create(
            summary_raster_path, 361, 181, 1, gdal.GDT_Float32)
        summary_raster.SetProjection(wgs84_sr.ExportToWkt())
        wgs84_gt = [-180.0, 1.0, 0, 90., 0, -1]
        summary_raster.SetGeoTransform(wgs84_gt)
        summary_band = summary_raster.GetRasterBand(1)
        nodata = -1
        summary_band.SetNoDataValue(nodata)
        base_array = numpy.empty((181, 361), dtype=numpy.float32)
        base_array[:] = nodata
        inv_gt = summary_band.InvGeoTransform()

        print 'load summary pickle'
        with open(summary_grid_index_path_path, 'r') as pickle_file:
            summary_grid_index_map = pickle.load(summary_grid_index_path_path)

        for grid_x, grid_y in summary_grid_index_map:
            i_x = inv_gt[0] + grid_x * inv_gt[1]
            i_y = inv_gt[3] + grid_y * inv_gt[5]
            print grid_x, grid_y, i_x, i_y

    print 'done'

if __name__ == '__main__':
    main()
