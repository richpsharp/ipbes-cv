"""Sort CV fields into alphabetical order."""
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

CV_VECTOR_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_cv_results', 'global_cv_rt_ssp_recalc - Copy.shp')


def main():
    vector = gdal.OpenEx(CV_VECTOR_PATH, gdal.OF_VECTOR | gdal.OF_UPDATE)
    layer = vector.GetLayer()

    feature = layer.GetNextFeature()
    for last_index in xrange(feature.GetFieldCount()):
        layer.ResetReading()
        feature = layer.GetNextFeature()
        field_name_list = []
        for i in xrange(last_index, feature.GetFieldCount()):
            field_name_list.append(feature.GetFieldDefnRef(i).GetName())
        min_field_name = min(field_name_list)
        print min_field_name, last_index, feature.GetFieldCount()
        layer.ReorderField(
            field_name_list.index(
                min_field_name, key=lambda x: x.lower())+last_index,
            last_index)

    print 'syncing'
    layer.SyncToDisk()
    print 'done'

if __name__ == '__main__':
    main()
