"""IPBES global sea level rise post SLR analytics https://docs.google.com/document/d/1yX7aIOehd04v3CZt9FDgRsKsUjR8DSlkpv_XhKDXpC8/edit#heading=h.6fz9zq50ghb3."""
import pickle
import os
import logging
import taskgraph
import numpy
import collections

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
    'ipbes_cv_results', 'without slr incorporated into RT',
    'global_cv_rt_ssp.shp')

TARGET_GLOBAL_CV_SLR_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_cv_results', 'global_cv_rt_ssp_recalc.shp')

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

    #1) Sea level rise
    # Set R_slr_cur to 1 for whole globe
    target_layer.CreateField(ogr.FieldDefn('R_slr_cur', ogr.OFTInteger))

    target_layer.CreateField(ogr.FieldDefn('rnom_cur', ogr.OFTReal))

    for slr_risk_field_id, slr_id in [
            ('slr_ssp1_r', 'slr_rcp26'),
            ('slr_ssp3_r', 'slr_rcp60'),
            ('slr_ssp5_r', 'slr_rcp85')]:
        target_layer.CreateField(
            ogr.FieldDefn(slr_risk_field_id, ogr.OFTInteger))

    #Population - recalibrate scenarios population to current gpw pdn (see Beneficiaries Dilemma design doc)
    #New columns: pdnrc_ssp[1|3|5]
    for ssp_id in [1, 3, 5]:
        target_layer.CreateField(
            ogr.FieldDefn('pdnrc_ssp%d' % ssp_id, ogr.OFTReal))
        target_layer.CreateField(
            ogr.FieldDefn('crnom_ssp%d' % ssp_id, ogr.OFTReal))

    rt_slr_field_id_ssp_id_map = {}
    rt_nohab_slr_field_id_ssp_id_map = {}
    for ssp_id in [1, 3, 5]:
        # "Rtslrsp[n]" means total risk with sea level rise for SSP[n]
        rt_slr_field_id = 'Rtslrsp%d' % ssp_id
        # "Rtslrsp[n]nh" means total risk with sea level rise but no habitat for SSP[n]
        rt_nohab_slr_field_id = 'Rtslrsp%dnh' % ssp_id
        rt_slr_field_id_ssp_id_map[ssp_id] = rt_slr_field_id
        rt_nohab_slr_field_id_ssp_id_map[ssp_id] = rt_nohab_slr_field_id
        target_layer.CreateField(ogr.FieldDefn(rt_slr_field_id, ogr.OFTReal))
        target_layer.CreateField(
            ogr.FieldDefn(rt_nohab_slr_field_id, ogr.OFTReal))

    risk_ids = ['Rwind', 'Rwave', 'Rrelief', 'Rsurge']
    hab_id = 'Rhab'

    crpop_ssp_list = collections.defaultdict(list)
    crage_ssp_list = collections.defaultdict(list)
    crpov_ssp_list = collections.defaultdict(list)
    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break
        feature.SetField('R_slr_cur', 1)

        # Rt/Rt_nohab = adding slr into geometric mean (with SLR = 1 everywhere)
        # load up the existing risks for the feature with a SLR risk of 1.0
        hab_risk = feature.GetField(hab_id)
        base_risk_list = [
            feature.GetField(x) for x in risk_ids] + [1.0]

        # calculate geometric mean with hab
        rt_slr_risk = numpy.prod(
            [hab_risk] + base_risk_list) ** (
                1. / (len(base_risk_list) + 1))

        # calculate geometric mean with no hab
        rt_slr_nohab_risk = numpy.prod(
            [5.] + base_risk_list) ** (
                1. / (len(base_risk_list) + 1))

        feature.SetField('Rt', rt_slr_risk)
        feature.SetField('Rtnohab', rt_slr_nohab_risk)

        # rpop_cur = Rt * pdn_gpw
        rpop_cur = feature.GetField('Rt') * feature.GetField('pdn_gpw')
        feature.SetField('rpop_cur', rpop_cur)
        crpop_ssp_list['cur'].append(rpop_cur)

        # rage_cur  =  rt_cur * poverty_p (ssp versions below)
        rage_cur = feature.GetField('Rt') * feature.GetField('poverty_p')
        feature.SetField('rage_cur', rage_cur)
        crage_ssp_list['cur'].append(rage_cur)

        # rpov_cur|ssp[1|3|5]  =  rt_cur|ssp[1|3|5] * (14bt_pop+65plus_pop)
        rpov_cur = feature.GetField('Rt') * (
                feature.GetField('14bt_pop') + feature.GetField('65plus_pop'))
        feature.SetField('rpov_cur', rpov_cur)
        crpov_ssp_list['cur'].append(rpov_cur)

        for ssp_id in [1, 3, 5]:
            # Population - recalibrate scenarios population to current gpw pdn (see Beneficiaries Dilemma design doc)
            # New columns: pdnrc_ssp[1|3|5]
            # pdn_sspn / pdn_2010 * pdn_gpw
            try:
                pdnrc_ssp_val = feature.GetField('pdn_gpw') * (
                    feature.GetField('pdn_ssp%d' % ssp_id) / feature.GetField('pdn_2010'))
            except ZeroDivisionError:
                pdnrc_ssp_val = 0.0
            feature.SetField('pdnrc_ssp%d' % ssp_id, pdnrc_ssp_val)

        # Histogram values in slr_rcp[26 | 60 | 85] columns across all rcp scenarios; then rank 1-5 (with 5 being highest risk)
        for slr_risk_field_id, slr_id, ssp_id in [
                ('slr_ssp1_r', 'slr_rcp26', 1),
                ('slr_ssp3_r', 'slr_rcp60', 3),
                ('slr_ssp5_r', 'slr_rcp85', 5)]:
            # retrieve slr value and find what percentile it's in for risk
            slr_val = feature.GetField(slr_id)
            slr_risk = numpy.searchsorted(
                slr_thresholds, slr_val, side='right')
            feature.SetField(slr_risk_field_id, slr_risk)

            # load up the existing risks for the feature
            hab_risk = feature.GetField(hab_id)
            base_risk_list = [
                feature.GetField(x) for x in risk_ids] + [slr_risk]

            # calculate geometric mean with hab
            rt_slr_risk = numpy.prod(
                [hab_risk] + base_risk_list) ** (
                    1. / (len(base_risk_list) + 1))

            # calculate geometric mean with no hab
            rt_slr_nohab_risk = numpy.prod(
                [5.] + base_risk_list) ** (
                    1. / (len(base_risk_list) + 1))

            feature.SetField(
                rt_slr_field_id_ssp_id_map[ssp_id], rt_slr_risk)
            feature.SetField(
                rt_nohab_slr_field_id_ssp_id_map[ssp_id], rt_slr_nohab_risk)

            # rt_ssp[1|3|5] = Rtslrsp[1|3|5] when urbp_ssp[1|3|5] < 0.2; else Rtslrsp[1|3|5]nh
            if feature.GetField('urbp_ssp%d' % ssp_id) < 0.2:
                feature.SetField('rt_ssp%d' % ssp_id, rt_slr_risk)
            else:
                feature.SetField('rt_ssp%d' % ssp_id, rt_slr_nohab_risk)

            # rpop_ssp[1|3|5] = rt_ssp[1|3|5] *  pdnrc_ssp[1|3|5]
            feature.SetField(
                'rpop_ssp%s' % ssp_id,
                feature.GetField('rt_ssp%d' % ssp_id) * feature.GetField('pdnrc_ssp%d' % ssp_id))

            # rage_ssp[1|3|5]  =  rt|ssp[1|3|5] * poverty_p (current versions above)
            feature.SetField(
                'rage_ssp%s' % ssp_id,
                feature.GetField('rt_ssp%d' % ssp_id) * feature.GetField('poverty_p'))

            # rpov_cur|ssp[1|3|5]  =  rt_cur|ssp[1|3|5] * (14bt_pop+65plus_pop)
            feature.SetField(
                'rpov_ssp%s' % ssp_id,
                feature.GetField('rt_ssp%d' % ssp_id) * (
                    feature.GetField('14bt_pop') + feature.GetField('65plus_pop')))

            # crt_ssp[1|3|5] = (rt_ssp[1|3|5] - rt_cur)/rt_cur
            feature.SetField(
                'crt_ssp%s' % ssp_id, (
                    feature.GetField('rt_ssp%d' % ssp_id) -
                    feature.GetField('Rt')) / feature.GetField('Rt'))

            # crpop_ssp[1|3|5] = (rpop_ssp[1|3|5] - rpop_cur)/ rpop_cur
            try:
                crpop_ssp = (
                    feature.GetField('rpop_ssp%d' % ssp_id) -
                    feature.GetField('rpop_cur')) / feature.GetField('rpop_cur')
            except ZeroDivisionError:
                crpop_ssp = 0.0

            feature.SetField('crpop_ssp%s' % ssp_id, crpop_ssp)
            crpop_ssp_list[ssp_id].append(crpop_ssp)

            # crage_ssp[1|3|5] = (rage_ssp[1|3|5] - rage_cur) / rage_cur
            try:
                crage_ssp = (
                    feature.GetField('rage_ssp%d' % ssp_id) -
                    feature.GetField('rage_cur')) / feature.GetField('rage_cur')
            except ZeroDivisionError:
                crage_ssp = 0.0

            feature.SetField(
                'crage_ssp%s' % ssp_id, crage_ssp)
            crage_ssp_list[ssp_id].append(crage_ssp)

            # crpov_ssp[1|3|5] = (rpov_ssp[1|3|5] -rpov_cur) / rpov_cur
            try:
                crpov_ssp = (
                    feature.GetField('rpov_ssp%d' % ssp_id) -
                    feature.GetField('rpov_cur')) / feature.GetField('rpov_cur')
            except ZeroDivisionError:
                crpov_ssp = 0.0
            feature.SetField('crpov_ssp%s' % ssp_id, crpov_ssp)
            crpov_ssp_list[ssp_id].append(crpov_ssp)

        target_layer.SetFeature(feature)

    # loop through one more time to set the normalized results
    crpop_ssp_90th = {
        'cur': sorted(crpop_ssp_list['cur'])[int(len(crpop_ssp_list['cur']) * .9)],
        1: sorted(crpop_ssp_list[1])[int(len(crpop_ssp_list[1]) * .9)],
        3: sorted(crpop_ssp_list[3])[int(len(crpop_ssp_list[3]) * .9)],
        5: sorted(crpop_ssp_list[5])[int(len(crpop_ssp_list[5]) * .9)],
    }

    crage_ssp_90th = {
        'cur': sorted(crage_ssp_list['cur'])[int(len(crage_ssp_list['cur']) * .9)],
        1: sorted(crage_ssp_list[1])[int(len(crage_ssp_list[1]) * .9)],
        3: sorted(crage_ssp_list[3])[int(len(crage_ssp_list[3]) * .9)],
        5: sorted(crage_ssp_list[5])[int(len(crage_ssp_list[5]) * .9)],
    }

    crpov_ssp_90th = {
        'cur': sorted(crpov_ssp_list['cur'])[int(len(crpov_ssp_list['cur']) * .9)],
        1: sorted(crpov_ssp_list[1])[int(len(crpov_ssp_list[1]) * .9)],
        3: sorted(crpov_ssp_list[3])[int(len(crpov_ssp_list[3]) * .9)],
        5: sorted(crpov_ssp_list[5])[int(len(crpov_ssp_list[5]) * .9)],
    }

    target_layer.ResetReading()

    while True:
        feature = target_layer.GetNextFeature()
        if not feature:
            break

        rnom_cur = (
            feature.GetField('rpop_cur') / crpop_ssp_90th['cur'] +
            feature.GetField('rage_cur') / crage_ssp_90th['cur'] +
            feature.GetField('rpov_cur') / crpov_ssp_90th['cur']) / 3.
        feature.SetField('rnom_cur', rnom_cur)

        for ssp_id in [1, 3, 5]:
            crnom_ssp = (
                feature.GetField('crpop_ssp%d' % ssp_id) / crpop_ssp_90th[ssp_id] +
                feature.GetField('crage_ssp%d' % ssp_id) / crage_ssp_90th[ssp_id] +
                feature.GetField('crpov_ssp%d' % ssp_id) / crpov_ssp_90th[ssp_id]) / 3.
            feature.SetField(
                'crnom_ssp%d' % ssp_id, crnom_ssp)

        target_layer.SetFeature(feature)

if __name__ == '__main__':
    main()