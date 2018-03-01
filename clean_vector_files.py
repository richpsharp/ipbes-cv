"""Create a factors and a result shapefile with only the desired fields."""
import logging
import os

from osgeo import gdal
from osgeo import ogr

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('clean-vector-files')
LOGGER.setLevel(logging.DEBUG)

BASE_CV_VECTOR_PATH = os.path.join(
    'ipbes-cv-post-processing-workspace',
    'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

WORKSPACE_DIR = 'clean_vector_dir'

TARGET_FACTOR_VECTOR_PATH = os.path.join(WORKSPACE_DIR, 'CV_factors.shp')
TARGET_OUTPUTS_VECTOR_PATH = os.path.join(WORKSPACE_DIR, 'CV_outputs.shp')


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
        ['SLRrise_%s' % x for x in |['1', '3', '5']] +
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
        ['pServ_cur'] +
        ['nRisk'] +
        ['pRisk_cur'] +
        ['nServ'])

    esri_driver = gdal.GetDriverByName('ESRI Shapefile')

    base_vector = gdal.OpenEx(BASE_CV_VECTOR_PATH, gdal.OF_VECTOR)

    for path, field_set in [
            (TARGET_FACTOR_VECTOR_PATH, factor_field_set),
            (TARGET_OUTPUTS_VECTOR_PATH, outputs_field_set)]:
        if os.path.exists(path):
            esri_driver.Delete(path)

        target_vector = esri_driver.CreateCopy(path, base_vector)
        target_layer = target_vector.GetLayer()

        for i in reversed(xrange(target_layer.GetFieldCount())):
            if self.GetFieldDefnRef(i).GetName() not in field_set:
                LOGGER.info(
                    'deleting %s from %s', self.GetFieldDefnRef(i).GetName(),
                    path)
                target_layer.DeleteField(i)

if __name__ == '__main__':
    main()
