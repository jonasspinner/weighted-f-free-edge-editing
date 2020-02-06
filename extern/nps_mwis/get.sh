#!/bin/bash
cd "$(dirname "$0")" || exit

ILS_MWIS_DIR=.
ILS_MWIS_DOWNLOAD_DIR=$ILS_MWIS_DIR/download

mkdir -p $ILS_MWIS_DOWNLOAD_DIR

curl -L -o $ILS_MWIS_DOWNLOAD_DIR/ils-mwis.tar.gz https://sites.google.com/site/nogueirabruno/software/ils-mwis.tar.gz
tar -C $ILS_MWIS_DIR -xf $ILS_MWIS_DOWNLOAD_DIR/ils-mwis.tar.gz
rm $ILS_MWIS_DOWNLOAD_DIR/ils-mwis.tar.gz

for file in $ILS_MWIS_DIR/patches/*.patch; do
    patch -p1 -i $file
done