#########################################################################
# File Name: get_result.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Sat May  5 15:16:07 2018
#########################################################################
#!/bin/bash


cat=('GA15' 'GA15-noGrad' 'GA00' 'GA00-GR' 'noGA')

vsh01="opa-sx-180425_Gaiadr2_vsh01.log"
ga="opa-sx-180425_Gaiadr2_GA.log"

if [ -e $vsh01 ]
then
	rm $vsh01
fi
if [ -e $ga ]
then
	rm $ga
fi

for cati in ${cat[@]};
do
	echo "#${cati}" >>${vsh01} 
 	grep "used_for_grep_all" opa-sx-180425-${cati}_GaiaDR2_vsh.log >>${vsh01}
	echo "#${cati}" >>${ga} 
 	grep "used_for_grep_GA" opa-sx-180425-${cati}_GaiaDR2_vsh.log >>${ga} 
done
