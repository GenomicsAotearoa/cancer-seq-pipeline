#!/bin/bash

#ml BEDTools/2.26.0-gimkl-2017a


echo $1
echo $2

echo "bedtools multicov -bams $1 $2  -bed ../data/hg19.3kb.bed > ../data/prosper01_3kb.count"
