#########################################################################
# File Name: sort_arcfile.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Tue Jul 10 15:22:42 2018
#########################################################################
#!/bin/bash

arcfile=$1
sortfile=${arcfile}.sort

rm ${sortfile}
#cp ${arcfile} ${arcfile}.bk
sed -i "" "/^*/d" ${arcfile}

# For data among 1979.0 ~ 1999.0
for year in {79..99}
do
    grep "\$${year}" $arcfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done

# for data among 2000 ~ 2009
for year in {0..9}
do
    grep "\$0${year}" $arcfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done

for year in {10..18}
do
    grep "\$${year}" $arcfile | sort -b -M -k 2.4,2.6 >> ${sortfile}
done

awk 'NR % 2 == 1' ${sortfile} >deci-test1c.arc
awk 'NR % 2 == 0' ${sortfile} >deci-test2c.arc
