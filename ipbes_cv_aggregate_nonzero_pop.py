"""
This script aggregates a simplified point shapefile for CV based off the

'IPBES_CV.shp' file described in https://docs.google.com/document/d/1yX7aIOehd04v3CZt9FDgRsKsUjR8DSlkpv_XhKDXpC8/edit#heading=h.1mo37dt7cnav
"""
import logging
import os
import shapely
import shapely.wkb
import shapely.prepared

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes_cv_aggregate_nonzero_pop')
LOGGER.setLevel(logging.DEBUG)

POSSIBLE_DROPBOX_LOCATIONS = [
    r'D:\Dropbox',
    r'C:\Users\Rich\Dropbox',
    r'C:\Users\rpsharp\Dropbox',
    r'E:\Dropbox']

LOGGER.info("checking dropbox locations")
for path in POSSIBLE_DROPBOX_LOCATIONS:
    if os.path.exists(path):
        BASE_DROPBOX_DIR = path
        LOGGER.info("found %s", BASE_DROPBOX_DIR)
        break

ROOT_CV_SHAPEFILE_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_cv_results',
    'cv_results_3_6_2018',
    'POST_PROCESS_global_cv_risk_and_aggregate_analysis.shp')

EEZ_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff', 'summary table shapefile',
    'EEZv8_WVS_DIS_V3_ALL_final_v7disIPBES',
    'EEZv8_WVS_DIS_V3_ALL_final_v7disIPBES.shp')

WORKSPACE_DIR = 'aggregate_nzpop_workspace'

TARGET_PATH = os.path.join(WORKSPACE_DIR, 'IPBES_CV.shp')


def main():
    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    eez_vector = gdal.OpenEx(EEZ_PATH, gdal.OF_VECTOR)
    eez_layer = eez_vector.GetLayer()

    LOGGER.info('load IPBES regions')
    prep_poly_region_map = {}
    while True:
        eez_feature = eez_layer.GetNextFeature()
        if eez_feature is None:
            break
        eez_geom = eez_feature.GetGeometryRef()
        shapely_eez_geom_prep = shapely.prepared.prep(shapely.wkb.loads(
            eez_geom.ExportToWkb()))
        eez_geom = None
        prep_poly_region_map[
            (eez_feature.GetField('IPBES_regi'),
             eez_feature.GetField('IPBES_sub'))] = shapely_eez_geom_prep
        shapely_eez_geom_prep = None

    LOGGER.info('load create target file')
    root_vector = gdal.OpenEx(ROOT_CV_SHAPEFILE_PATH, gdal.OF_VECTOR)
    root_layer = root_vector.GetLayer()
    LOGGER.info(root_layer.GetFeatureCount())

    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    esri_driver = gdal.GetDriverByName("ESRI Shapefile")
    if os.path.exists(TARGET_PATH):
        esri_driver.Delete(TARGET_PATH)
    target_vector = esri_driver.Create(TARGET_PATH, 0, 0, 0, gdal.GDT_Unknown)
    target_layer = target_vector.CreateLayer(
        'ipbes_cv', wgs84_srs, ogr.wkbPoint)

    field = ogr.FieldDefn("IPBES_regi", ogr.OFTString)
    field.SetWidth(24)
    target_layer.CreateField(field)
    field = ogr.FieldDefn("IPBES_sub", ogr.OFTString)
    field.SetWidth(24)
    target_layer.CreateField(field)
    field = ogr.FieldDefn("country", ogr.OFTString)
    field.SetWidth(24)
    target_layer.CreateField(field)

    for field_name in [
            'cRtssp1', 'cRtssp3', 'cRtssp5', 'cSvRt_ssp1', 'cSvRt_ssp3',
            'cSvRt_ssp5', 'pdnrc_ssp1', 'pdnrc_ssp3', 'pdnrc_ssp5',
            '14bt_pop', '65plus_pop', 'poverty_p']:
        target_layer.CreateField(ogr.FieldDefn(field_name, ogr.OFTReal))

    LOGGER.info('process file')
    nonzero_count = 0
    while True:
        feature = root_layer.GetNextFeature()
        if feature is None:
            break
        if feature.GetField('pdnrc_ssp1') == 0.0:
            nonzero_count += 1
        target_feature = ogr.Feature(target_layer.GetLayerDefn())
        geom = feature.GetGeometryRef()
        target_feature.SetGeometry(geom.Clone())
        shapely_point = shapely.wkb.loads(geom.ExportToWkb())
        geom = None

        for (region, subregion), prepped_poly in prep_poly_region_map.iteritems():
            if prepped_poly.contains(shapely_point):
                target_feature.SetField('IPBES_regi', region)
                target_feature.SetField('IPBES_sub', subregion)
                break

        for field_name in [
                'cRtssp1', 'cRtssp3', 'cRtssp5', 'cSvRt_ssp1', 'cSvRt_ssp3',
                'cSvRt_ssp5', 'pdnrc_ssp1', 'pdnrc_ssp3', 'pdnrc_ssp5',
                '14bt_pop', '65plus_pop', 'poverty_p']:
            target_feature.SetField(
                field_name, feature.GetField(field_name))
        target_layer.CreateFeature(target_feature)
        target_feature = None
    target_vector = None

if __name__ == '__main__':
    main()

