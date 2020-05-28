#########################################################################
# File Name: fig_cp.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Sat Jun 30 22:56:59 2018
#########################################################################
#!/bin/bash

figs=("x-eema.eps" \
	  "DSM.eps" \
	  "SBL.eps" \
	  "sess_num_sf.eps" \
	  "dec_nf.eps" \
	  "dec_sf.eps" \
	  "catalog_comparison.eps")

for fig in ${figs[@]}
do
#	cp ../plots/${fig} ../notes/20180530/figures/
# cp ../plots/${fig} ../notes/20180821/figures/
    cp ../plots/${fig} ~/Desktop/TBD/thesis-note/chapter4/figs/noise/
done
