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
    "ipbes-cv-post-processing-workspace",
    "POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp")

TARGET_SUMMARY_DEGREE_DIR = os.path.join(
    "ipbes-cv-post-processing-workspace", 'sum_to_degree')

SLR_LIST_PICKLE_PATH = 'slr.pickle'

WORKING_DIR = 'ipbes_cv_to_degree_workspace'


def sort_points(summary_field_list, pickle_dir):
    """Read all features from GLOBAL_CV_VECTOR_PATH and sort to grid.

    Parameters:
        summary_field_list (string): field ID to record per grid.
        pickle_dir (string): path to a directoty that will create pickle
            files named field_name.pickle. Each pickle file will be a
            dictionary mapping
            (long, lat) -> [feature0[summary_field], feature1[summary_field]...]

    Returns:
        None.
    """
    try:
        os.makedirs(pickle_dir)
    except OSError:
        pass

    vector = gdal.OpenEx(GLOBAL_CV_VECTOR_PATH, gdal.OF_VECTOR)
    print GLOBAL_CV_VECTOR_PATH
    layer = vector.GetLayer()

    field_list_grid_index_map = collections.defaultdict(
        lambda: collections.defaultdict(list))

    while True:
        feature = layer.GetNextFeature()
        if not feature:
            break
        centroid = feature.GetGeometryRef().Centroid()
        grid_index = (int(centroid.GetX()), int(centroid.GetY()))
        for field_id in summary_field_list:
            field_list_grid_index_map[field_id][grid_index].append(
                feature.GetField(field_id))

    for field_id in summary_field_list:
        pickle_path = os.path.join(pickle_dir, '%s.pickle' % field_id)
        with open(pickle_path, 'wb') as pickle_file:
            pickle.dump(field_list_grid_index_map[field_id], pickle_file)


def main():
    """Entry point."""
    for dir_path in [WORKING_DIR, TARGET_SUMMARY_DEGREE_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    cur_ssp_list = ['cur', 'ssp1', 'ssp3', 'ssp5']

    summary_field_list = (
        ['Rt_%s' % x for x in cur_ssp_list] +
        ['cRtssp%d' % x for x in [1, 3, 5]] +
        ['Serv_%s' % x for x in cur_ssp_list] +
        ['cServssp%d' % x for x in [1, 3, 5]] +
        ['pRisk_%s' % x for x in cur_ssp_list] +
        ['nRisk_%s' % x for x in cur_ssp_list] +
        ['pServ_%s' % x for x in cur_ssp_list] +
        ['nServ_%s' % x for x in cur_ssp_list] +
        ['cpRiskssp%d' % x for x in [1, 3, 5]] +
        ['cpServssp%d' % x for x in [1, 3, 5]] +
        ['cnRiskssp%d' % x for x in [1, 3, 5]] +
        ['cnServssp%d' % x for x in [1, 3, 5]])

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKING_DIR, 'taskgraph_cache'), -1)

    summary_grid_pickle_dir = os.path.join(WORKING_DIR, 'pickles')
    task = task_graph.add_task(
        func=sort_points,
        args=(summary_field_list, summary_grid_pickle_dir,),
        target_path_list=[summary_grid_pickle_dir])
    """
    sort_points(summary_field_list, summary_grid_pickle_dir)
    """

    task_graph.join()

    for summary_field in summary_field_list:
        print 'summary %s' % summary_field
        summary_raster_path = os.path.join(
            WORKING_DIR, '%s.tif' % summary_field)
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
        summary_band.Fill(nodata)
        base_array = numpy.empty((181, 361), dtype=numpy.float32)
        base_array[:] = nodata
        inv_gt = gdal.InvGeoTransform(wgs84_gt)

        print 'load summary pickle'
        pickle_path = os.path.join(
            summary_grid_pickle_dir, '%s.pickle' % summary_field)
        with open(pickle_path, 'r') as pickle_file:
            summary_grid_index_map = pickle.load(pickle_file)

        for grid_x, grid_y in summary_grid_index_map:
            i_x = int(inv_gt[0] + grid_x * inv_gt[1])
            i_y = int(inv_gt[3] + grid_y * inv_gt[5])
            value = sorted(
                summary_grid_index_map[(grid_x, grid_y)])[
                    int(.9*len(summary_grid_index_map[(grid_x, grid_y)]))]
            summary_band.WriteArray(numpy.array([[value]]), i_x, i_y)

    print 'done'

if __name__ == '__main__':
    main()
