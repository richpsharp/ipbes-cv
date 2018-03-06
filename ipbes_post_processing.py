"""This script is meant to be run on the data produced by `ipbes-cv.py`."""
import os
import pickle
import logging
import numpy

import taskgraph
from osgeo import ogr
from osgeo import gdal

gdal.UseExceptions()

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

WORKSPACE_DIR = r'C:\fast_dir\ipbes-cv-post-processing-workspace'

BASE_CV_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'cv_results',
    'global_cv_risk_and_aggregate_analysis.shp')

TARGET_CV_VECTOR_PATH = os.path.join(
    WORKSPACE_DIR, 'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

SLRISE_LIST_PICKLE_PATH = os.path.join(WORKSPACE_DIR, 'slrise.pickle')

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
    if os.path.exists(TARGET_CV_VECTOR_PATH):
        esri_driver.Delete(TARGET_CV_VECTOR_PATH)
    LOGGER.info('creating vector')
    target_vector = esri_driver.CreateCopy(
        TARGET_CV_VECTOR_PATH, base_vector)
    target_layer = target_vector.GetLayer()

    # Add new fields that we're making in this file
    for new_field_id in [
            'Rslr_cur', 'Rslr_ssp1', 'Rslr_ssp3', 'Rslr_ssp5',
            'Serv_cur', 'Serv_ssp1', 'Serv_ssp3', 'Serv_ssp5',
            'dSLRrt_sp1', 'dSLRrt_sp3', 'dSLRrt_sp5', 'Rt_cur_nh',
            'Rt_ssp1', 'Rt_ssp3', 'Rt_ssp5',
            'Rt_ssp1_nh', 'Rt_ssp3_nh', 'Rt_ssp5_nh',
            'pdnrc_ssp1', 'pdnrc_ssp3', 'pdnrc_ssp5', 'pRisk_cur',
            'pRisk_ssp1', 'pRisk_ssp3', 'pRisk_ssp5',
            'pServ_cur', 'pServ_ssp1', 'pServ_ssp3', 'pServ_ssp5',
            'aRisk_cur', 'aRisk_ssp1', 'aRisk_ssp3', 'aRisk_ssp5',
            'aServ_cur', 'aServ_ssp1', 'aServ_ssp3', 'aServ_ssp5',
            'vRisk_cur', 'vRisk_ssp1', 'vRisk_ssp3', 'vRisk_ssp5',
            'vServ_cur', 'vServ_ssp1', 'vServ_ssp3', 'vServ_ssp5',
            'nRisk_cur', 'nRisk_ssp1', 'nRisk_ssp3', 'nRisk_ssp5',
            'nServ_cur', 'nServ_ssp1', 'nServ_ssp3', 'nServ_ssp5',
            'cRtssp1', 'cRtssp3', 'cRtssp5',
            'cServssp1', 'cServssp3', 'cServssp5',
            'cpRiskssp1', 'cpRiskssp3', 'cpRiskssp5',
            'cpServssp1', 'cpServssp3', 'cpServssp5',
            'cnRiskssp1', 'cnRiskssp3', 'cnRiskssp5',
            'cnServssp1', 'cnServssp3', 'cnServssp5',
            'SvRt_cur', 'SvRt_ssp1', 'SvRt_ssp3', 'SvRt_ssp5',
            'pSvRt_cur', 'pSvRt_ssp1', 'pSvRt_ssp3', 'pSvRt_ssp5',
            'aSvRt_cur', 'aSvRt_ssp1', 'aSvRt_ssp3', 'aSvRt_ssp5',
            'vSvRt_cur', 'vSvRt_ssp1', 'vSvRt_ssp3', 'vSvRt_ssp5',
            'cSvRt_ssp1', 'cSvRt_ssp3', 'cSvRt_ssp5',
            'nSvRt_cur', 'nSvRt_ssp1', 'nSvRt_ssp3', 'nSvRt_ssp5',
            'cnSvRtssp1', 'cnSvRtssp3', 'cnSvRtssp5',
            'cpSvRtssp1', 'cpSvRtssp3', 'cpSvRtssp5',
            'logpop_cur', 'logpop_s1', 'logpop_s3', 'logpop_s5'
            ]:
        target_layer.CreateField(
            ogr.FieldDefn(new_field_id, ogr.OFTReal))

    risk_ids = ['Rwind', 'Rwave', 'Rrelief', 'Rsurge']

    min_max_id = {}

    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break
        #logpop_cur = log(pdn_gpw)
        pdn_gpw = feature.GetField('pdn_gpw')
        if pdn_gpw > 0:
            feature.SetField('logpop_cur', numpy.log(pdn_gpw))
        else:
            feature.SetField('logpop_cur', 0.0)

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

            # logpop_s[1|3|5] = log(pdnrc_ssp[1|3|5])
            if feature.GetField('pdnrc_ssp%d' % ssp_id) > 0:
                feature.SetField(
                    'logpop_s%d' % ssp_id, numpy.log(
                        feature.GetField('pdnrc_ssp%d' % ssp_id)))
            else:
                feature.SetField('logpop_s%d' % ssp_id, 0.0)

            pRisk = (
                feature.GetField('pdnrc_ssp%d' % ssp_id) *
                feature.GetField('Rt_ssp%d' % ssp_id))
            if 'pRisk' not in min_max_id:
                min_max_id['pRisk'] = [pRisk, pRisk]
            else:
                min_max_id['pRisk'][0] = min(pRisk, min_max_id['pRisk'][0])
                min_max_id['pRisk'][1] = max(pRisk, min_max_id['pRisk'][1])
            feature.SetField('pRisk_ssp%d' % ssp_id, pRisk)

            pServ = (
                feature.GetField('pdnrc_ssp%d' % ssp_id) *
                feature.GetField('Serv_ssp%d' % ssp_id))
            if 'pServ' not in min_max_id:
                min_max_id['pServ'] = [pServ, pServ]
            else:
                min_max_id['pServ'][0] = min(pServ, min_max_id['pServ'][0])
                min_max_id['pServ'][1] = max(pServ, min_max_id['pServ'][1])
            feature.SetField('pServ_ssp%d' % ssp_id, pServ)

            aRisk = feature.GetField('Rt_ssp%d' % ssp_id) * (
                feature.GetField('14bt_pop') +
                feature.GetField('65plus_pop'))
            if 'aRisk' not in min_max_id:
                min_max_id['aRisk'] = [aRisk, aRisk]
            else:
                min_max_id['aRisk'][0] = min(aRisk, min_max_id['aRisk'][0])
                min_max_id['aRisk'][1] = max(aRisk, min_max_id['aRisk'][1])
            feature.SetField('aRisk_ssp%d' % ssp_id, aRisk)

            aServ = feature.GetField('Serv_ssp%d' % ssp_id) * (
                feature.GetField('14bt_pop') +
                feature.GetField('65plus_pop'))
            if 'aServ' not in min_max_id:
                min_max_id['aServ'] = [aServ, aServ]
            else:
                min_max_id['aServ'][0] = min(aServ, min_max_id['aServ'][0])
                min_max_id['aServ'][1] = max(aServ, min_max_id['aServ'][1])
            feature.SetField('aServ_ssp%d' % ssp_id, aServ)

            vRisk = (
                feature.GetField('poverty_p') *
                feature.GetField('Rt_ssp%d' % ssp_id))
            if 'vRisk' not in min_max_id:
                min_max_id['vRisk'] = [vRisk, vRisk]
            else:
                min_max_id['vRisk'][0] = min(vRisk, min_max_id['vRisk'][0])
                min_max_id['vRisk'][1] = max(vRisk, min_max_id['vRisk'][1])
            feature.SetField('vRisk_ssp%d' % ssp_id, vRisk)

            vServ = (
                feature.GetField('poverty_p') *
                feature.GetField('Serv_ssp%d' % ssp_id))
            if 'vServ' not in min_max_id:
                min_max_id['vServ'] = [vServ, vServ]
            else:
                min_max_id['vServ'][0] = min(vServ, min_max_id['vServ'][0])
                min_max_id['vServ'][1] = max(vServ, min_max_id['vServ'][1])
            feature.SetField('vServ_ssp%d' % ssp_id, vServ)

            if feature.GetField('Rt_ssp%d' % ssp_id) != 0:
                SvRt = (
                    feature.GetField('Serv_ssp%d' % ssp_id) /
                    feature.GetField('Rt_ssp%d' % ssp_id))
            else:
                SvRt = 0
            if 'SvRt' not in min_max_id:
                min_max_id['SvRt'] = [SvRt, SvRt]
            else:
                min_max_id['SvRt'][0] = min(SvRt, min_max_id['SvRt'][0])
                min_max_id['SvRt'][1] = max(SvRt, min_max_id['SvRt'][1])
            feature.SetField('SvRt_ssp%d' % ssp_id, SvRt)

            pSvRt = (
                feature.GetField('pdnrc_ssp%d' % ssp_id) *
                feature.GetField('SvRt_ssp%d' % ssp_id))
            if 'pSvRt' not in min_max_id:
                min_max_id['pSvRt'] = [pSvRt, pSvRt]
            else:
                min_max_id['pSvRt'][0] = min(pSvRt, min_max_id['pSvRt'][0])
                min_max_id['pSvRt'][1] = max(pSvRt, min_max_id['pSvRt'][1])
            feature.SetField('pSvRt_ssp%d' % ssp_id, pSvRt)

            aSvRt = feature.GetField('SvRt_ssp%d' % ssp_id) * (
                feature.GetField('14bt_pop') +
                feature.GetField('65plus_pop'))

            if 'aSvRt' not in min_max_id:
                min_max_id['aSvRt'] = [aSvRt, aSvRt]
            else:
                min_max_id['aSvRt'][0] = min(aSvRt, min_max_id['aSvRt'][0])
                min_max_id['aSvRt'][1] = max(aSvRt, min_max_id['aSvRt'][1])
            feature.SetField('aSvRt_ssp%d' % ssp_id, aSvRt)

            vSvRt = feature.GetField('SvRt_ssp%d' % ssp_id) * (
                feature.GetField('poverty_p'))

            if 'vSvRt' not in min_max_id:
                min_max_id['vSvRt'] = [vSvRt, vSvRt]
            else:
                min_max_id['vSvRt'][0] = min(vSvRt, min_max_id['vSvRt'][0])
                min_max_id['vSvRt'][1] = max(vSvRt, min_max_id['vSvRt'][1])
            feature.SetField('vSvRt_ssp%d' % ssp_id, vSvRt)

        pRisk = feature.GetField('pdn_gpw') * feature.GetField('Rt_cur')
        if 'pRisk' not in min_max_id:
            min_max_id['pRisk'] = [pRisk, pRisk]
        else:
            min_max_id['pRisk'][0] = min(pRisk, min_max_id['pRisk'][0])
            min_max_id['pRisk'][1] = max(pRisk, min_max_id['pRisk'][1])
        feature.SetField('pRisk_cur', pRisk)

        pServ = feature.GetField('pdn_gpw') * feature.GetField('Serv_cur')
        if 'pServ' not in min_max_id:
            min_max_id['pServ'] = [pServ, pServ]
        else:
            min_max_id['pServ'][0] = min(pServ, min_max_id['pServ'][0])
            min_max_id['pServ'][1] = max(pServ, min_max_id['pServ'][1])
        feature.SetField('pServ_cur', pServ)

        aRisk = feature.GetField('Rt_cur') * (
            feature.GetField('14bt_pop') + feature.GetField('65plus_pop'))
        if 'aRisk' not in min_max_id:
            min_max_id['aRisk'] = [aRisk, aRisk]
        else:
            min_max_id['aRisk'][0] = min(aRisk, min_max_id['aRisk'][0])
            min_max_id['aRisk'][1] = max(aRisk, min_max_id['aRisk'][1])
        feature.SetField('aRisk_cur', aRisk)

        aServ = feature.GetField('Serv_cur') * (
            feature.GetField('14bt_pop') + feature.GetField('65plus_pop'))
        if 'aServ' not in min_max_id:
            min_max_id['aServ'] = [aServ, aServ]
        else:
            min_max_id['aServ'][0] = min(aServ, min_max_id['aServ'][0])
            min_max_id['aServ'][1] = max(aServ, min_max_id['aServ'][1])
        feature.SetField('aServ_cur', aServ)

        vServ = (
            feature.GetField('poverty_p') * feature.GetField('Serv_cur'))
        if 'vServ' not in min_max_id:
            min_max_id['vServ'] = [vServ, vServ]
        else:
            min_max_id['vServ'][0] = min(vServ, min_max_id['vServ'][0])
            min_max_id['vServ'][1] = max(vServ, min_max_id['vServ'][1])
        feature.SetField('vServ_cur', vServ)

        vRisk = (
            feature.GetField('poverty_p') * feature.GetField('Rt_cur'))
        if 'vRisk' not in min_max_id:
            min_max_id['vRisk'] = [vRisk, vRisk]
        else:
            min_max_id['vRisk'][0] = min(vRisk, min_max_id['vRisk'][0])
            min_max_id['vRisk'][1] = max(vRisk, min_max_id['vRisk'][1])
        feature.SetField('vRisk_cur', vRisk)

        if feature.GetField('Rt_cur') != 0:
            SvRt = (
                feature.GetField('Serv_cur') /
                feature.GetField('Rt_cur'))
        else:
            SvRt = 0
        if 'SvRt' not in min_max_id:
            min_max_id['SvRt'] = [SvRt, SvRt]
        else:
            min_max_id['SvRt'][0] = min(SvRt, min_max_id['SvRt'][0])
            min_max_id['SvRt'][1] = max(SvRt, min_max_id['SvRt'][1])
        feature.SetField('SvRt_cur', SvRt)

        pSvRt = (
            feature.GetField('pdn_gpw') *
            feature.GetField('SvRt_cur'))
        if 'pSvRt' not in min_max_id:
            min_max_id['pSvRt'] = [pSvRt, pSvRt]
        else:
            min_max_id['pSvRt'][0] = min(pSvRt, min_max_id['pSvRt'][0])
            min_max_id['pSvRt'][1] = max(pSvRt, min_max_id['pSvRt'][1])
        feature.SetField('pSvRt_cur', pSvRt)

        aSvRt = feature.GetField('SvRt_cur') * (
            feature.GetField('14bt_pop') +
            feature.GetField('65plus_pop'))

        if 'aSvRt' not in min_max_id:
            min_max_id['aSvRt'] = [aSvRt, aSvRt]
        else:
            min_max_id['aSvRt'][0] = min(aSvRt, min_max_id['aSvRt'][0])
            min_max_id['aSvRt'][1] = max(aSvRt, min_max_id['aSvRt'][1])
        feature.SetField('aSvRt_cur', aSvRt)

        vSvRt = feature.GetField('SvRt_cur') * (
            feature.GetField('poverty_p'))

        if 'vSvRt' not in min_max_id:
            min_max_id['vSvRt'] = [vSvRt, vSvRt]
        else:
            min_max_id['vSvRt'][0] = min(vSvRt, min_max_id['vSvRt'][0])
            min_max_id['vSvRt'][1] = max(vSvRt, min_max_id['vSvRt'][1])
        feature.SetField('vSvRt_cur', vSvRt)

        for ssp_id in [1, 3, 5]:
            if feature.GetField('Rt_cur') != 0:
                feature.SetField(
                    'cRtssp%d' % ssp_id, (
                        feature.GetField('Rt_ssp%d' % ssp_id) -
                        feature.GetField('Rt_cur')) /
                    feature.GetField('Rt_cur'))
            else:
                feature.SetField('cRtssp%d' % ssp_id, 0)

            if feature.GetField('Serv_cur') != 0:
                feature.SetField(
                    'cServssp%d' % ssp_id, (
                        feature.GetField('Serv_ssp%d' % ssp_id) -
                        feature.GetField('Serv_cur')) /
                    feature.GetField('Serv_cur'))
            else:
                feature.SetField('cServssp%d' % ssp_id, 0)

            if feature.GetField('pRisk_cur') != 0:
                feature.SetField(
                    'cpRiskssp%d' % ssp_id, (
                        feature.GetField('pRisk_ssp%d' % ssp_id) -
                        feature.GetField('pRisk_cur')) /
                    feature.GetField('pRisk_cur'))
            else:
                feature.SetField('cpRiskssp%d' % ssp_id, 0)

            if feature.GetField('pServ_cur') != 0:
                feature.SetField(
                    'cpServssp%d' % ssp_id, (
                        feature.GetField('pServ_ssp%d' % ssp_id) -
                        feature.GetField('pServ_cur')) /
                    feature.GetField('pServ_cur'))
            else:
                feature.SetField('cpServssp%d' % ssp_id, 0)

            if feature.GetField('SvRt_cur') != 0:
                feature.SetField(
                    'cSvRt_ssp%d' % ssp_id, (
                        feature.GetField('SvRt_ssp%d' % ssp_id) -
                        feature.GetField('SvRt_cur')) /
                    feature.GetField('SvRt_cur'))
            else:
                feature.SetField('cSvRt_ssp%d' % ssp_id, 0)

            if feature.GetField('pSvRt_cur') != 0:
                feature.SetField(
                    'cpSvRtssp%d' % ssp_id, (
                        feature.GetField('pSvRt_ssp%d' % ssp_id) -
                        feature.GetField('pSvRt_cur')) /
                    feature.GetField('pSvRt_cur'))
            else:
                feature.SetField('cpSvRt%d' % ssp_id, 0)

        target_layer.SetFeature(feature)

    target_layer.ResetReading()
    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break

        feature.SetField(
            'nRisk_cur', numpy.mean([
                (feature.GetField('pRisk_cur') - min_max_id['pRisk'][0]) / (min_max_id['pRisk'][1]-min_max_id['pRisk'][0]),
                (feature.GetField('vRisk_cur') - min_max_id['vRisk'][0]) / (min_max_id['vRisk'][1]-min_max_id['vRisk'][0]),
                (feature.GetField('aRisk_cur') - min_max_id['aRisk'][0]) / (min_max_id['aRisk'][1]-min_max_id['aRisk'][0])]))

        feature.SetField(
            'nServ_cur', numpy.mean([
                (feature.GetField('pServ_cur') - min_max_id['pServ'][0]) / (min_max_id['pServ'][1]-min_max_id['pServ'][0]),
                (feature.GetField('vServ_cur') - min_max_id['vServ'][0]) / (min_max_id['vServ'][1]-min_max_id['vServ'][0]),
                (feature.GetField('aServ_cur') - min_max_id['aServ'][0]) / (min_max_id['aServ'][1]-min_max_id['aServ'][0])]))

        feature.SetField(
            'nSvRt_cur', numpy.mean([
                (feature.GetField('pSvRt_cur') - min_max_id['pSvRt'][0]) / (min_max_id['pSvRt'][1]-min_max_id['pSvRt'][0]),
                (feature.GetField('vSvRt_cur') - min_max_id['vSvRt'][0]) / (min_max_id['vSvRt'][1]-min_max_id['vSvRt'][0]),
                (feature.GetField('aSvRt_cur') - min_max_id['aSvRt'][0]) / (min_max_id['aSvRt'][1]-min_max_id['aSvRt'][0])]))

        for ssp_id in [1, 3, 5]:
            feature.SetField(
                'nRisk_ssp%d' % ssp_id, numpy.mean([
                    (feature.GetField('pRisk_ssp%d' % ssp_id) - min_max_id['pRisk'][0]) / (min_max_id['pRisk'][1]-min_max_id['pRisk'][0]),
                    (feature.GetField('vRisk_ssp%d' % ssp_id) - min_max_id['vRisk'][0]) / (min_max_id['vRisk'][1]-min_max_id['vRisk'][0]),
                    (feature.GetField('aRisk_ssp%d' % ssp_id) - min_max_id['aRisk'][0]) / (min_max_id['aRisk'][1]-min_max_id['aRisk'][0])]))

            feature.SetField(
                'nServ_ssp%d' % ssp_id, numpy.mean([
                    (feature.GetField('pServ_ssp%d' % ssp_id) - min_max_id['pServ'][0]) / (min_max_id['pServ'][1]-min_max_id['pServ'][0]),
                    (feature.GetField('vServ_ssp%d' % ssp_id) - min_max_id['vServ'][0]) / (min_max_id['vServ'][1]-min_max_id['vServ'][0]),
                    (feature.GetField('aServ_ssp%d' % ssp_id) - min_max_id['aServ'][0]) / (min_max_id['aServ'][1]-min_max_id['aServ'][0])]))

            feature.SetField(
                'nSvRt_ssp%d' % ssp_id, numpy.mean([
                    (feature.GetField('pSvRt_ssp%d' % ssp_id) - min_max_id['pSvRt'][0]) / (min_max_id['pSvRt'][1]-min_max_id['pSvRt'][0]),
                    (feature.GetField('vSvRt_ssp%d' % ssp_id) - min_max_id['vSvRt'][0]) / (min_max_id['vSvRt'][1]-min_max_id['vSvRt'][0]),
                    (feature.GetField('aSvRt_ssp%d' % ssp_id) - min_max_id['aSvRt'][0]) / (min_max_id['aSvRt'][1]-min_max_id['aSvRt'][0])]))

            if feature.GetField('nRisk_cur') != 0:
                feature.SetField(
                    'cnRiskssp%d' % ssp_id, (
                        feature.GetField('nRisk_ssp%d' % ssp_id) -
                        feature.GetField('nRisk_cur')) /
                    feature.GetField('nRisk_cur'))
            else:
                feature.SetField('cpRiskssp%d' % ssp_id, 0)

            if feature.GetField('nServ_cur') != 0:
                feature.SetField(
                    'cnServssp%d' % ssp_id, (
                        feature.GetField('nServ_ssp%d' % ssp_id) -
                        feature.GetField('nServ_cur')) /
                    feature.GetField('nServ_cur'))
            else:
                feature.SetField('cnServssp%d' % ssp_id, 0)

            if feature.GetField('nServ_cur') != 0:
                feature.SetField(
                    'cnSvRtssp%d' % ssp_id, (
                        feature.GetField('nSvRt_ssp%d' % ssp_id) -
                        feature.GetField('nSvRt_cur')) /
                    feature.GetField('nSvRt_cur'))
            else:
                feature.SetField('cnSvRtssp%d' % ssp_id, 0)

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
