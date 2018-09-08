#!/bin/sh

for i in `find ./ -type f|grep vcf|grep "-"`
do
    echo $i
    subtractBed -a $i -b NOG.bed > tmp
    mv tmp $i
done