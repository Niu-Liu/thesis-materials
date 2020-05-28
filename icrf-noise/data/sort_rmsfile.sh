#########################################################################
# File Name: sort_rmsfile.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Tue Jul 10 15:22:42 2018
#########################################################################
#!/bin/bash

rmsfile=$1
sortfile=${rmsfile}.sort

head -n 3 ${rmsfile} >${sortfile}

# For data among 1979.0 ~ 1999.0
for year in {79..99}
do
    grep "\$${year}" $rmsfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done

# for data among 2000 ~ 2009
for year in {0..9}
do
    grep "\$0${year}" $rmsfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done

for year in {10..18}
do
    grep "\$${year}" $rmsfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done
