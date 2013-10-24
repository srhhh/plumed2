#!/bin/bash
echo -n "unset key; plot " >check.gplt
n=`ls -1d CV_* | wc -l | awk '{print $1}' `
for i in ` ls CV_* | cut -d_ -f2 | sort -n `
do 
echo -n " \"CV_$i\" u 2:3 w p " >>check.gplt
if [ "$i" -lt "$n" ]; then 
echo -n "," >>check.gplt
fi
done
gnuplot -persist check.gplt
