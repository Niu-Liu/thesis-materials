#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: fig_copy.py
"""
Created on Mon May 13 00:37:57 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os

# -----------------------------  FUNCTIONS -----------------------------
figfiles = [
    "gdr2-G-sky-dist.eps",
    "gdr2-sigma-vs-G.eps",
    "gdr2-nosou-vs-G.eps",
    "gdr2-rotation-vs-G.eps",
    "gdr2-glide-vs-G.eps",
    "gdr2-EM20-vs-G.eps"]

dir1 = "../plots/"
dir2 = "/Users/Neo/Desktop/TBD/thesis-note/chapter5/figs/"

print("Input directory: ", dir1)
print("Output directory: ", dir2)

# Remove old figures
os.system("rm {}*".format(dir2))

for fig in figfiles:
    # os.system("cp {}{} {}".format(dir1, fig, dir2))
    os.system("cp {}{} {}".format(dir1, fig, dir2))
    print("Copy {}: done!".format(fig))

# --------------------------------- END --------------------------------
