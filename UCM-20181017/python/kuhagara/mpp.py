#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 06:04:17 2018

@author: vdoan
"""

import numpy as np
kme = 100
z_top = 2000.
dz1 = 2
dz2 = z_top/20.

z_t = np.zeros(kme)
z_t[0] = -dz1/2.
z_t[1] =  dz1/2
dz_t = np.zeros(kme)
print dz1, dz2
dz_t[0] = 2.

a = (2000. - dz1/2 - 2* (kme - 2)) * 2. / (kme-2)/(kme-3)
a = 2.*2.*z_top / (kme-2)/(kme-3)

for i in range(1,kme-1):
    #z_t[i] = z_t[i-1] + a*2*(i-2)
    dz_t[i] = dz_t[i-1] + i**1.0001
    z_t[i+1] = z_t[i] + dz_t[i]
    print dz_t[i], z_t[i]
