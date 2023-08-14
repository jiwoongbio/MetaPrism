#!/usr/bin/env bash
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

codePath=`dirname $0`
dataPath="$codePath/data"

time perl $codePath/MetaPrism_data.gene.pl -p 32 -b -r -t 2,2157,4751 -i 100
time perl $codePath/MetaPrism_data.taxon.pl -t 2,2157,4751
time perl $codePath/MetaPrism_data.gene.base.pl -p 32
