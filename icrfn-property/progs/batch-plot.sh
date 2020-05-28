#!/bin/bash
#########################################################################
# File Name: batch-plot.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Sat May 18 10:37:55 2019
#########################################################################

# Plot figures
ipython median-error.py
echo "Plot for median-error: done!"

ipython icrf-formal-error-bgt.py
echo "Plot for icrf-formal-error-bgt: done!"

ipython ra-dependent-error.py
echo "Plot for ra-dependent-error: done!"

ipython dec-dependent-error.py
echo "Plot for dec-dependent-error: done!"

ipython pos-diff-vs-dec.py
echo "Plot for pos-diff-vs-dec: done!"

# ipython gaia-vlbi-separation.py
# echo "Plot for gaia-vlbi-separation: done!"

ipython gaia-vlbi-separation-bgt.py
echo "Plot for gaia-vlbi-separation-bgt: done!"

ipython vsh_param.py
echo "Plot for vsh_param: done!"

# ipython nor-sep-dist.py
# echo "Plot for nor-sep-dist: done!"

# ipython effect-of-g-limiting-magnitude-plot.py
# echo "Plot for effect-of-g-limiting-magnitude-plot: done!"

ipython nor-sep-dist-bgt.py
echo "Plot for nor-sep-dist-bgt: done!"

ipython angular-separation-icrf-vs-icrf3.py
echo "Plot for angular-separation-icrf-vs-icrf3: done!"


# Copy figures
ipython fig_copy.py
