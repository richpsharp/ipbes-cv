"""Create a factors and a result shapefile with only the desired fields."""
import logging
import os

import taskgraph
from osgeo import gdal
from osgeo import ogr

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('clean-vector-files')
LOGGER.setLevel(logging.DEBUG)

BASE_CV_VECTOR_PATH = os.path.join(
    'C:', 'fast_dir', 'ipbes-cv-post-processing-workspace',
    'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

WORKSPACE_DIR = os.path.join("C:", 'fast_dir', 'clean_vector_dir')

TARGET_FACTOR_VECTOR_PATH = os.path.join(WORKSPACE_DIR, 'CV_factors.shp')
TARGET_OUTPUTS_VECTOR_PATH = os.path.join(WORKSPACE_DIR, 'CV_outputs.shp')
TARGET_SVRT_VECTOR_PATH = os.path.join(WORKSPACE_DIR, 'CV_SvRt.shp')
TARGET_CV_POPOUTPUTS_PATH = os.path.join(WORKSPACE_DIR, 'CV_popoutputs.shp')


def clean_cv_vector(base_path, target_path, field_set):
    """Copy base and remove the fields in `field_set`."""
    esri_driver = gdal.GetDriverByName('ESRI Shapefile')
    if os.path.exists(target_path):
        esri_driver.Delete(target_path)

    base_vector = gdal.OpenEx(base_path, gdal.OF_VECTOR)
    base_layer = base_vector.GetLayer()
    field_definitions = base_layer.GetLayerDefn()

    LOGGER.info('copying %s to memory', base_path)
    target_vector = esri_driver.CreateCopy(target_path, base_vector)

    target_layer = target_vector.GetLayer()
    target_defn = target_layer.GetLayerDefn()

    base_field_names = set(
        [field_definitions.GetFieldDefn(i).GetName()
         for i in xrange(field_definitions.GetFieldCount())])
    field_definitions = None

    unknown_field_names = field_set.subtract(base_field_names)
    if unknown_field_names:
        raise ValueError(
            "The following field names were unknown: %s",
            unknown_field_names)

    for i in reversed(xrange(target_defn.GetFieldCount())):
        if target_defn.GetFieldDefn(i).GetName() not in field_set:
            LOGGER.info(
                'deleting %s from %s',
                target_defn.GetFieldDefn(i).GetName(), base_path)
            target_layer.DeleteField(i)

    LOGGER.info('done saving %s', target_path)


def main():
    """Entry point."""
    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    factor_field_set = set(
        ['Rhab_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['Rslra_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['Rslri_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['urbp_%s' % x for x in ['2015', 'ssp1', 'ssp3', 'ssp5']] +
        ['pdn_%s' % x for x in ['2010', 'ssp1', 'ssp3', 'ssp5']] +
        ['SLRrate_%s' % x for x in ['c', '1', '3', '5']] +
        ['SLRrise_%s' % x for x in ['1', '3', '5']] +
        ['curpb_ssp%s' % x for x in ['1', '3', '5']] +
        ['cpdn_ssp%s' % x for x in ['1', '3', '5']] +
        ['pdnrc_ssp%s' % x for x in ['1', '3', '5']] +
        ['slr_rcp%s' % x for x in ['26', '60', '85']] +
        ['dSLRrt_sp%s' % x for x in ['1', '3', '5']] +
        ['dSLRrs_sp%s' % x for x in ['1', '3', '5']] +
        ['pdn_gpw'] +
        ['14bt_pop'] +
        ['65plus_pop'] +
        ['poverty_p'] +
        ['Rwind'] +
        ['Rwave'] +
        ['Rrelief'] +
        ['Rsurge'])

    outputs_field_set = set(
        ['Rt_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['Rt_%s_nh' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['Serv_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['cRtssp%s' % x for x in ['1', '3', '5']] +
        ['cServssp%s' % x for x in ['1', '3', '5']] +
        ['cpRiskssp%s' % x for x in ['1', '3', '5']] +
        ['cpServssp%s' % x for x in ['1', '3', '5']] +
        ['cnRiskssp%s' % x for x in ['1', '3', '5']] +
        ['cnServssp%s' % x for x in ['1', '3', '5']] +
        ['pRisk_ssp%s' % x for x in ['1', '3', '5']] +
        ['pServ_ssp%s' % x for x in ['1', '3', '5']] +
        ['aRisk_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['vRisk_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['aServ_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['vServ_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['nRisk_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['nServ_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['pServ_cur'] +
        ['pRisk_cur'])

    cv_svrt_field_set = set(
        ['SvRt_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['aSvRt_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['vSvRt_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['nSvRt_%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['cSvRt%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['cpSvRt%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['cnSvRt%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']] +
        ['pSvRt_ssp%s' % x for x in ['cur', 'ssp1', 'ssp3', 'ssp5']])

    cv_pop_outputs_set = set(
        ['SvRt_cur'] +
        ['cSvRt_ssp%d' % x for x in [1,3,5]] +
        ['logpop_cur'] +
        ['logpop_s%d' % x for x in [1,3,5]] +
        ['c_logpop_s%d' % x for x in [1,3,5]] +
        ['logage'] +
        ['logpRt_cur'] +
        ['logpServ_cur'] +
        ['logpSvRt_cur'] +
        ['logaRt_cur'] +
        ['logaServ_cur'] +
        ['logaSvRt_cur'] +
        ['pcRt_s%d' % x for x in [1,3,5]] +
        ['pcServ_s%d' % x for x in [1,3,5]] +
        ['pcSvRt_s%d' % x for x in [1,3,5]] +
        ['acRt_s%d' % x for x in [1,3,5]] +
        ['acServ_s%d' % x for x in [1,3,5]] +
        ['acSvRt_s%d' % x for x in [1,3,5]])

    task_graph = taskgraph.TaskGraph('clean_vector_taskgraph_dir', 4)
    for path, field_set in [
            (TARGET_FACTOR_VECTOR_PATH, factor_field_set),
            (TARGET_OUTPUTS_VECTOR_PATH, outputs_field_set),
            (TARGET_SVRT_VECTOR_PATH, cv_svrt_field_set),
            (TARGET_CV_POPOUTPUTS_PATH, cv_pop_outputs_set)]:
        task_graph.add_task(
            func=clean_cv_vector,
            args=(BASE_CV_VECTOR_PATH, path, field_set),
            task_name=os.path.basename(path))

    task_graph.close()
    task_graph.join()

if __name__ == '__main__':
    main()
