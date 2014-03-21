# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import numpy as np

class Particle:
    def __init__(self,refractive_index):
        self.ind = refractive_index
    
    def intersection_distance(self,ray): pass

class Sphere(Particle):
    def __init__(self, radius, center=np.array([0,0,0]), refractive_index=1.333):
        Particle.__init__(self,refractive_index)
        self.r = radius
        self.c = center
    
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
    
    def intersection_point(self,distance):
        return self.o + distance*self.d

def reflect(ray,surf,point):
    return

def refract(ray,surf,point):
    return

def main():
    sph = Sphere(4)
    ray = Ray(np.array([1,2,-5]), np.array([0,0,1]))

if __name__ == '__main__':
    main()
