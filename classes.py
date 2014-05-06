# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from sys import float_info
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d

class Particle:
    def __init__(self,refractive_index):
        self.m = refractive_index
    
    def intersection_distance(self,ray): pass

    def has_point(self,point): pass

    def unit_normal(self,point): pass

class Sphere(Particle):
    def __init__(self, radius, center=np.array([0,0,0]), refractive_index=1.333):
        Particle.__init__(self,refractive_index)
        self.r = radius
        self.c = center
        
    def unit_normal(self,point):
        return Ray.unit(point-self.c)
        
    def has_point(self,point):
        d_center = np.linalg.norm(point-self.c)
        return d_center + float_info.epsilon < self.r
    
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
        
    def plot(self,ax):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x=self.r*np.cos(u)*np.sin(v)
        y=self.r*np.sin(u)*np.sin(v)
        z=self.r*np.cos(v)
        ax.plot_wireframe(x, y, z, color="r")
        
m = {
'water': 1.333,
'ethanol': 1.36,
'ice': 1.309,
'diamond': 2.42
}

class Ray:
    def __init__(self, origin, direction=np.array([0,0,1])):
        self.o = origin
        self.d = Ray.unit(direction)
        self.stokes = np.array([1,0,0,0])
        self.l = -1
    
    @staticmethod
    def unit(vector):
        return vector/np.linalg.norm(vector)
    
    def intersection_point(self,particle):
        distance = particle.intersection_distance(self)
        if distance < 0:
            return np.array([])
        return self.o + distance*self.d
        
    def intersects(self,particle):
        return particle.intersection_distance(self) > 0
        
    def reflection(self, particle, point):
        if point.size == 0:
            return
        dnsign = 1 if particle.has_point(self.o) else -1
        di = -self.d
        dn = dnsign*particle.unit_normal(point)
        ds = 2*(np.dot(dn,di))*dn-di
        self.l = self.distance_from_origin(point)
        return Ray(point,ds)
    
    def distance_from_origin(self,p):
        return np.linalg.norm(p-self.o)
    
    def plot(self,ax):
        Ray.plot_line(ax,self.o,self.endpoint())

    def endpoint(self):
        length = 10 if self.l<0 else self.l
        return self.o + length*self.d
        
    def plot_origin(self, ax):
        ax.scatter(self.o[0],self.o[1],self.o[2],color='#6633ff')
    
    @staticmethod    
    def plot_line(ax, p1, p2):
        l = art3d.Line3D([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]])
        ax.add_line(l)

def main():
    sph = Sphere(10)
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', projection='3d')
    for x in range(100):
        ray = Ray(np.array([rnd.uniform(-10,10),rnd.uniform(-10,10),-15]))
        pnt = ray.intersection_point(sph)
        if pnt.size == 0:
            continue
        ray_r = ray.reflection(sph,pnt)
        ray.plot(ax)
        ray_r.plot(ax)
        ray_r.plot_origin(ax)
    sph.plot(ax)
    plt.show()

if __name__ == '__main__':
    main()