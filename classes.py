# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""

class Sphere:
    def __init__(self, center, radius, ref_ind):
        self.center = center
        self.radius = radius
        self.ref_ind = ref_ind
        
material = {
'water': 1.333,
'ethanol': 1.36,
'ice': 1.309,
'diamond': 2.42
}