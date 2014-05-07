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
    EPS = 100*float_info.epsilon
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
        return d_center < self.r + self.EPS 
        
    def has_surfpoint(self,point):
        d_center = np.linalg.norm(point-self.c)
        return abs(d_center - self.r) < self.EPS
    
    def intersection_distance(self,ray):
        ray_originates_surf = self.has_surfpoint(ray.o)
        if ray_originates_surf:
            ray_displacement = 100
        else:        
            ray_displacement = 0
        ray_o = ray.o - ray_displacement*ray.d
        a = np.dot(ray.d, ray.d)
        b = 2*np.dot(ray_o - self.c, ray.d)
        c = np.dot(ray_o - self.c, ray_o - self.c) - self.r**2
        
        disc = b**2 - 4*a*c
        
        if disc<0:
            return -1
        
        dist_sqrt = np.sqrt(disc)
        if b<0:
            q = (-b - dist_sqrt)/2
        else:
            q = (-b + dist_sqrt)/2
        
        t0 = q / a - ray_displacement
        t1 = c / q - ray_displacement
        
        if t0>t1:
            tmp = t0
            t0 = t1
            t1 = tmp
        
        if t1<0:
            return -1
        
        if t0<0 or ray_originates_surf:
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
        
    def __str__(self):
        return "ray origin: %s, direction: %s" % (self.o, self.d)
    
    @staticmethod
    def unit(vector):
        return vector/np.linalg.norm(vector)
    
    def intersection_point(self,particle):
        distance = particle.intersection_distance(self)
        if distance < 0:
            return np.array([])
        self.l = distance
        return self.o + distance*self.d
        
    def intersects(self,particle):
        return particle.intersection_distance(self) > 0
        
    def ref_ction(self, particle):
        point = self.intersection_point(particle)
        if point.size == 0:
            return (None, None)
        ray_is_inside = particle.has_point(self.o)
        if ray_is_inside:
            dnsign = -1
            r = particle.m
        else:
            dnsign = 1
            r = 1/particle.m
        di = self.d
        dn = dnsign*particle.unit_normal(point)
        drefl = di-2*(np.dot(dn,di))*dn
        c = np.dot(-dn,di)
        drefr = r*di+(r*c-np.sqrt(1-r**2*(1-c**2)))*dn
        reflection = Ray(point,drefl)
        refraction = Ray(point,drefr)
        return (reflection,refraction)
        
    def refadefa(self,particle):
        point = self.intersection_point(particle)
        if point.size == 0:
            return (None, None)
        ray_is_inside = particle.has_point(self.o)
        if ray_is_inside:
            dnsign = -1
            r = particle.m
        else:
            dnsign = 1
            r = 1/particle.m
        di = self.d
        dn = dnsign*particle.unit_normal(point)
        ci = dn*np.dot(-di,dn)
        si = ci + di
        dr = ci + si
        st = r*si
        ct = -dn*np.sqrt(1-np.dot(st,st))
        dt = ct + st
        reflection = Ray(point,dr)
        refraction = Ray(point,dt)
        return (reflection,refraction)
    
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
    for x in range(200):
        ray = Ray(np.array([rnd.uniform(-10,10),rnd.uniform(-10,10),-15]))
        ray_refl, ray_refr = ray.ref_ction(sph)
        if ray_refl is None:
            continue
        ray_refl2, ray_refr2 = ray_refr.refadefa(sph)
        ray.plot(ax)
        ray_refl.plot(ax)
        ray_refr.plot(ax)
        ray_refl.plot_origin(ax)
        ray_refr2.plot_origin(ax)
        ray_refl2.plot(ax)
        ray_refr2.plot(ax)
    sph.plot(ax)
    plt.show()

if __name__ == '__main__':
    main()