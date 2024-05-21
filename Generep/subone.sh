#!/bin/bash

BIN=/blue/dtran/share/alberto/TCGA/bin
CONF=$1
DONE=${CONF}.done
submit -o --mem=5G,--qos=icbrbi,--time=120:00:00 -p $CONF -done $DONE generic.qsub $BIN/runone.sh $CONF
