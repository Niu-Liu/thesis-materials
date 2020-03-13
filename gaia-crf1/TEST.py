##!/usr/bin/env python3
## -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:45:58 2017

@author: Neo
"""

OUTFILE = open('test.log', 'w')
print('a b c', file = OUTFILE)
OUTFILE.close()
