# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Particle:
    def __init__(self,refractive_index):
        self.ind = refractive_index
    
    def intersection_distance(self,ray): pass

class Sphere(Particle):
    def __init__(self, radius, center=np.array([0,0,0]), refractive_index=1.333):
        Particle.__init__(self,refractive_index)
        self.r = radius
        self.c = center
        
    def normal(self,point):
        
        return
        
    def is_surface_point(self,point):
        return True
    
    def intersection_distance(self,ray):
        a = np.dot(ray.d, ray.d)
        b = 2*np.dot(ray.d - self.c, ray.o)
        c = np.dot(ray.o - self.c, ray.o - self.c) - self.r**2
        
        disc = b*b - 4*a*c
        
        if disc<0:
            return -1
        
        dist_sqrt = np.sqrt(disc)
        if b<0:
            q = (-b - dist_sqrt)/2
        else:
            q = (-b + dist_sqrt)/2
        
        t0 = q / a
        t1 = c / q
        
        if t0>t1:
            tmp = t0
            t0 = t1
            t1 = tmp
        
        if (t1<0):
            return -1
        
        if t0<0:
            return t1
            
        return t0
        
material = {
'water': 1.333,
'ethanol': 1.36,
'ice': 1.309,
'diamond': 2.42
}

class Ray:
    def __init__(self, origin, direction):
        self.o = origin
        self.d = direction
        self.stokes = np.array([1,0,0,0])
    
    def intersection_point(self,distance):
        return self.o + distance*self.d

def main():
    sph = Sphere(10)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for x in range(200):
        ray = Ray(np.array([rnd.uniform(-10,10),rnd.uniform(-10,10),-15]), np.array([0,0,1]))
        pnt = ray.intersection_point(sph.intersection_distance(ray))
        plot_point(ax,pnt)
        

def plot_point(ax, p):
    ax.scatter(p[0],p[1],p[2])

if __name__ == '__main__':
    main()