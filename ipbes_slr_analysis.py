"""IPBES global sea level rise risk calculation."""
import pickle
import os
import logging
import taskgraph
import numpy

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
    'ipbes_cv_results', 'global_cv_rt_ssp.shp')

TARGET_GLOBAL_CV_SLR_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_cv_results', 'global_cv_rt_ssp_slr_risk.shp')


SLR_LIST_PICKLE_PATH = 'slr.pickle'


def calculate_slr_list(sst_vector_path, slr_pickle_path):
    global_cv_sst_vector = ogr.Open(sst_vector_path)
    global_cv_sst_layer = global_cv_sst_vector.GetLayer()

    slr_list = []
    while True:
        feature = global_cv_sst_layer.GetNextFeature()
        if not feature:
            break
        for slr_field_id in ['slr_rcp26', 'slr_rcp60', 'slr_rcp85']:
            slr_list.append(feature.GetField(slr_field_id))

    with open(slr_pickle_path, 'wb') as slr_p_file:
        pickle.dump(slr_list, slr_p_file)


def main():
    """Entry point."""

    task_graph = taskgraph.TaskGraph('taskgraph_cache', -1, 5.0)

    task_graph.add_task(
        func=calculate_slr_list,
        args=(GLOBAL_CV_VECTOR_PATH, SLR_LIST_PICKLE_PATH),
        target_path_list=[SLR_LIST_PICKLE_PATH],
        task_name='calculate_slr_list')

    task_graph.join()

    with open(SLR_LIST_PICKLE_PATH, 'r') as slr_file:
        slr_list = pickle.load(slr_file)

    slr_count, slr_thresholds = numpy.histogram(sorted(slr_list), bins=5)
    LOGGER.debug(slr_thresholds)
    LOGGER.debug(numpy.count_nonzero(slr_list > slr_thresholds[-2]) / float(
        len(slr_list)))

    esri_driver = gdal.GetDriverByName('ESRI Shapefile')
    if os.path.exists(TARGET_GLOBAL_CV_SLR_VECTOR_PATH):
        esri_driver.Delete(TARGET_GLOBAL_CV_SLR_VECTOR_PATH)

    base_vector = gdal.OpenEx(GLOBAL_CV_VECTOR_PATH)
    target_vector = esri_driver.CreateCopy(
        TARGET_GLOBAL_CV_SLR_VECTOR_PATH, base_vector)

    target_layer = target_vector.GetLayer()
    for slr_risk_field_id, slr_id in [
            ('slr_ssp1_r', 'slr_rcp26'),
            ('slr_ssp3_r', 'slr_rcp60'),
            ('slr_ssp5_r', 'slr_rcp85')]:
        target_layer.CreateField(
            ogr.FieldDefn(slr_risk_field_id, ogr.OFTInteger))

    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break
        for slr_risk_field_id, slr_id in [
                ('slr_ssp1_r', 'slr_rcp26'),
                ('slr_ssp3_r', 'slr_rcp60'),
                ('slr_ssp5_r', 'slr_rcp85')]:
            slr_val = feature.GetField(slr_id)
            slr_risk = numpy.searchsorted(
                slr_thresholds, slr_val, side='right')
            feature.SetField(slr_risk_field_id, slr_risk)
        target_layer.SetFeature(feature)


if __name__ == '__main__':
    main()
