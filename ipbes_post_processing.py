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
    WORKSPACE_DIR, 'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

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

    LOGGER.info('pickling slr list')
    task_graph.add_task(
        func=calculate_slr_list,
        args=(BASE_CV_VECTOR_PATH, 'rise', SLRISE_LIST_PICKLE_PATH),
        target_path_list=[SLRISE_LIST_PICKLE_PATH],
        task_name='calculate_slrise_list')

    task_graph.join()
    LOGGER.info('loading pickle file')

    with open(SLRISE_LIST_PICKLE_PATH, 'r') as slr_file:
        slr_list = pickle.load(slr_file)

    slr_count, slr_thresholds = numpy.histogram(sorted(slr_list), bins=5)

    LOGGER.info('copying CV vector')
    esri_driver = gdal.GetDriverByName('ESRI Shapefile')
    base_vector = gdal.OpenEx(BASE_CV_VECTOR_PATH)
    if not os.path.exists(TARGET_CV_VECTOR_PATH):
        #esri_driver.Delete(TARGET_CV_VECTOR_PATH)
        LOGGER.info('creating vector')
        target_vector = esri_driver.CreateCopy(
            TARGET_CV_VECTOR_PATH, base_vector)
    else:
        LOGGER.info('loading existing vector')
        target_vector = gdal.OpenEx(
            TARGET_CV_VECTOR_PATH, gdal.OF_VECTOR | gdal.OF_UPDATE)

    target_layer = target_vector.GetLayer()

    # Add new fields that we're making in this file
    for new_field_id in [
            'Rslr_cur', 'Rslr_ssp1', 'Rslr_ssp3', 'Rslr_ssp5',
            'Serv_cur', 'Serv_ssp1', 'Serv_ssp3', 'Serv_ssp5',
            'dSLRrt_sp1', 'dSLRrt_sp3', 'dSLRrt_sp5', 'Rt_cur_nh',
            'Rt_ssp1', 'Rt_ssp3', 'Rt_ssp5',
            'Rt_ssp1_nh', 'Rt_ssp3_nh', 'Rt_ssp5_nh',
            'pdnrc_ssp1', 'pdnrc_ssp3', 'pdnrc_ssp5']:
        if target_layer.FindFieldIndex(new_field_id, 1) == -1:
            target_layer.CreateField(
                ogr.FieldDefn(new_field_id, ogr.OFTReal))

    risk_ids = ['Rwind', 'Rwave', 'Rrelief', 'Rsurge']

    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break
        for slr_risk_field_id, rhab_id, slr_id, rt_hab_id, rt_nohab_id, serv_id in [
                ('Rslr_cur', 'Rhab_cur', 'SLRrise_c', 'Rt_cur', 'Rt_cur_nh', 'Serv_cur'),
                ('Rslr_ssp1', 'Rhab_ssp1', 'SLRrise_1', 'Rt_ssp1', 'Rt_ssp1_nh', 'Serv_ssp1'),
                ('Rslr_ssp3', 'Rhab_ssp3', 'SLRrise_3', 'Rt_ssp3', 'Rt_ssp3_nh', 'Serv_ssp3'),
                ('Rslr_ssp5', 'Rhab_ssp5', 'SLRrise_5', 'Rt_ssp5', 'Rt_ssp5_nh', 'Serv_ssp5')]:
            # retrieve slr value and find what percentile it's in for risk
            slr_val = feature.GetField(slr_id)
            slr_risk = numpy.searchsorted(
                slr_thresholds, slr_val, side='right')
            feature.SetField(slr_risk_field_id, slr_risk)

            base_risk_list = [
                float(feature.GetField(x)) for x in risk_ids + [slr_id]]

            #LOGGER.debug(rhab_id)
            hab_risk = float(feature.GetField(rhab_id))
            #LOGGER.debug([hab_risk] + base_risk_list)
            # calculate geometric mean with hab
            rt_slr_risk = numpy.prod(
                [hab_risk] + base_risk_list) ** (1. / (len(base_risk_list) + 1))
            feature.SetField(rt_hab_id, rt_slr_risk)

            # calculate geometric mean with no hab
            rt_slr_nohab_risk = numpy.prod(
                [5.] + base_risk_list) ** (
                    1. / (len(base_risk_list) + 1))
            feature.SetField(rt_nohab_id, rt_slr_nohab_risk)

            feature.SetField(
                serv_id, feature.GetField(rt_nohab_id) -
                feature.GetField(rt_hab_id))

        for ssp_id in [1, 3, 5]:
            feature.SetField(
                'dSLRrt_sp%d' % ssp_id,
                feature.GetField('SLRrise_%d' % ssp_id) -
                feature.GetField('SLRrise_c'))

            if feature.GetField('pdn_2010') > 0:
                feature.SetField(
                    'pdnrc_ssp%d' % ssp_id,
                    feature.GetField('pdn_ssp%d' % ssp_id) /
                    feature.GetField('pdn_2010') * feature.GetField('pdn_gpw'))
            else:
                feature.SetField('pdnrc_ssp%d' % ssp_id, 0.0)

        target_layer.SetFeature(feature)



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
