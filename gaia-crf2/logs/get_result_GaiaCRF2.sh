#########################################################################
# File Name: get_result.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Sat May  5 15:16:07 2018
#########################################################################
#!/bin/bash


cat=('GaiaCRF2')

for cati in ${cat[@]};
do
	echo $cati
	echo "wrt. GaiaDR2"
 	sed -n "12p" ${cati}_GaiaDR2_vsh_param.tex
	echo "wrt. GaiaDR1"
	sed -n "12p" ${cati}_GaiaDR1_vsh_param.tex
	echo "wrt. ICRF2"
 	sed -n "12p" ${cati}_icrf2_vsh_param.tex
done
