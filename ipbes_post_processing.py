"""This script is meant to be run on the data produced by `ipbes-cv.py`."""
import os
import pickle
import logging
import numpy

import taskgraph
from osgeo import ogr
from osgeo import gdal

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-post-processing')
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

WORKSPACE_DIR = 'ipbes-cv-post-processing-workspace'

BASE_CV_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'cv_results',
    'global_cv_risk_and_aggregate_analysis.shp')

TARGET_CV_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'cv_results',
    'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

SLRISE_LIST_PICKLE_PATH = os.path.join(WORKSPACE_DIR, 'slrise.pickle')
SLRATE_LIST_PICKLE_PATH = os.path.join(WORKSPACE_DIR, 'slrate.pickle')


def main():
    """Entry point."""

    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass
    task_graph = taskgraph.TaskGraph(
        'post_processing_taskgraph_cache', -1)

    task_graph.add_task(
        func=calculate_slr_list,
        args=(BASE_CV_VECTOR_PATH, 'rise', SLRISE_LIST_PICKLE_PATH),
        target_path_list=[SLRISE_LIST_PICKLE_PATH],
        task_name='calculate_slrise_list')

    task_graph.join()

    with open(SLRISE_LIST_PICKLE_PATH, 'r') as slr_file:
        slr_list = pickle.load(slr_file)

    slr_count, slr_thresholds = numpy.histogram(sorted(slr_list), bins=5)

    esri_driver = gdal.GetDriverByName('ESRI Shapefile')
    if os.path.exists(TARGET_CV_VECTOR_PATH):
        esri_driver.Delete(TARGET_CV_VECTOR_PATH)

    base_vector = gdal.OpenEx(BASE_CV_VECTOR_PATH)
    target_vector = esri_driver.CreateCopy(
        TARGET_CV_VECTOR_PATH, base_vector)

    target_layer = target_vector.GetLayer()

    # 1) Sea level rise
    for slr_risk_field_id in [
            'Rslr_cur', 'Rslr_ssp1', 'Rslr_ssp3', 'Rslr_ssp5']:
        target_layer.CreateField(
            ogr.FieldDefn(slr_risk_field_id, ogr.OFTInteger))

    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break
        for slr_risk_field_id, slr_id in [
                ('Rslr_cur', 'SLRrise_c'),
                ('Rslr_ssp1', 'SLRrise_1'),
                ('Rslr_ssp3', 'SLRrise_3'),
                ('Rslr_ssp5', 'SLRrise_5')]:
            # retrieve slr value and find what percentile it's in for risk
            slr_val = feature.GetField(slr_id)
            slr_risk = numpy.searchsorted(
                slr_thresholds, slr_val, side='right')
            feature.SetField(slr_risk_field_id, slr_risk)


def calculate_slr_list(sst_vector_path, slr_type, slr_pickle_path):
    global_cv_sst_vector = ogr.Open(sst_vector_path)
    global_cv_sst_layer = global_cv_sst_vector.GetLayer()

    slr_list = []
    while True:
        feature = global_cv_sst_layer.GetNextFeature()
        if not feature:
            break
        for slr_field_id in [
                'SLR%s_%s' % (slr_type, c) for c in ['c', '1', '3', '5']]:
            slr_list.append(feature.GetField(slr_field_id))

    with open(slr_pickle_path, 'wb') as slr_p_file:
        pickle.dump(slr_list, slr_p_file)


if __name__ == '__main__':
    main()
