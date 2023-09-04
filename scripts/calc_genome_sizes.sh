#!/bin/bash
bins=$1
outfile=$2
ext=$3
echo "ID, length" > ${outfile}
for file in ${bins}/*.${ext}; do filename=$(basename $file) ; echo "$(echo ${filename##*mes/} | sed -e "s/\.${ext}//"),$(cat $file | grep -v ">" | wc -c)" >> ${outfile}; done