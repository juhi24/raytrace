# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from sys import float_info
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d, axes3d

class Particle:
    EPS = 100*float_info.epsilon
    def __init__(self,refractive_index):
        self.m = refractive_index
    
    def intersection_distance(self,ray): pass

    def has_point(self,point): pass

    def unit_normal(self,point): pass

class Sphere(Particle):
    def __init__(self, radius, center=np.array([0,0,0]), refractive_index=complex(1.333,1e-3)):
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

class Stokes:
    def __init__(self,s0,s1,s2,s3):
        self.s0 = s0
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        
    def psi(self):
        return 0.5*np.arctan(self.s2/self.s1)
        
    def chi(self):
        return 0.5*np.arctan(self.s3/np.sqrt(self.s1**2 + self.s2**2))
        
    def p(self):
        return np.sqrt(self.s1**2 + self.s2**2 + self.s3**2)/self.s0

class Ray:
    def __init__(self, origin, direction=np.array([0,0,1]), stokes=np.matrix([[1],[0],[0],[0]]), wavelength=0.5):
        self.o = origin
        self.d = Ray.unit(direction)
        self.stokes = stokes
        self.l = -1
        self.lam = wavelength
        
    def __str__(self):
        return "ray origin: %s, direction: %s, intensity: %s" % (self.o, self.d, self.stokes[0])
    
    @staticmethod
    def unit(vector):
        return vector/np.linalg.norm(vector)
    
    @staticmethod    
    def angle(v1,v2):
        return np.arccos(np.dot(Ray.unit(v1),Ray.unit(v2)))
    
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
            r = particle.m.real
            self.attenuate(particle)
        else:
            dnsign = 1
            r = 1/particle.m.real
        di = self.d
        dn = dnsign*particle.unit_normal(point)
        dr = di-2*(np.dot(dn,di))*dn
        c = np.dot(-dn,di)
        dt = r*di+(r*c-np.sqrt(1-r**2*(1-c**2)))*dn
        angle_i = Ray.angle(-1*dn,di)
        angle_t = Ray.angle(-1*dn,dt)
        #print('i: %s, t: %s' % (angle_i,angle_t))
        ang_plus = angle_i + angle_t
        ang_minus = angle_i - angle_t
        cos2_minus = np.cos(ang_minus)**2
        cos2_plus = np.cos(ang_plus)**2
        r_11_22 = cos2_minus + cos2_plus
        r_12_21 = cos2_minus - cos2_plus
        r_33_44 = -2*np.cos(ang_plus)*np.cos(ang_minus)
        mueller_r = 0.5*np.tan(ang_minus/ang_plus)**2*np.matrix([
            [r_11_22, r_12_21, 0, 0],
            [r_12_21, r_11_22, 0, 0],
            [0, 0, r_33_44, 0],
            [0, 0, 0, r_33_44]
        ])
        t_33_44 = 2*np.cos(ang_minus)
        mueller_t = np.sin(2*angle_i)*np.sin(2*angle_t)/2*(np.sin(np.sin(ang_plus)*np.cos(ang_minus)))**2*np.matrix([
            [cos2_minus+1, cos2_minus-1, 0, 0],
            [cos2_minus-1, cos2_minus+1, 0, 0],
            [0, 0, t_33_44, 0],
            [0, 0, 0, t_33_44]        
        ])
        stokes_r = mueller_r*self.stokes
        stokes_t = mueller_t*self.stokes
        reflection = Ray(point, dr, stokes=stokes_r, wavelength=self.lam)
        refraction = Ray(point, dt, stokes=stokes_t, wavelength=self.lam)
        return (reflection,refraction)
        
    def attenuate(self,particle):
        if self.l < 0:
            raise Exception('Cannot attenuate: Ray length not defined.')
        c_abs = 4*np.pi*particle.m.imag/self.lam
        self.stokes[0] = self.stokes[0]*np.exp(-c_abs*self.l)
    
    def distance_from_origin(self,p):
        return np.linalg.norm(p-self.o)
    
    def plot(self,ax):
        Ray.plot_line(ax,self.o,self.endpoint(),a=np.log(1+10*self.stokes[0])/np.log(11))

    def endpoint(self):
        length = 2000 if self.l<0 else self.l
        return self.o + length*self.d
        
    def plot_origin(self, ax):
        ax.scatter(self.o[0],self.o[1],self.o[2],color='#6633ff')
    
    @staticmethod    
    def plot_line(ax, p1, p2, a=1):
        l = art3d.Line3D([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],alpha=a)
        ax.add_line(l)
        
def trace(ray,particle,ax=None,draw_rays=True):
    intensity = ray.stokes[0]
    if intensity < 1e-5 or np.isnan(intensity):
        return
    ray_r, ray_t = ray.ref_ction(particle)
    if (ray_r is not None or intensity < 1) and draw_rays:
        ray.plot(ax)
        #ray.plot_origin(ax)
    #print(ray)
    if ray_r is None: # ray misses particle
        #print('Missed!')
        return
    trace(ray_t,particle,ax,draw_rays)
    trace(ray_r,particle,ax,draw_rays)
    
def particle_radius(x,wavelength):
    return 0.5*x*wavelength/np.pi
    
def main():
    draw3d = True
    lam = 0.5 #um
    x = 1000
    
    sph_r = particle_radius(x,lam)
    sph = Sphere(sph_r)
    if draw3d:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal', projection='3d')
        sph.plot(ax)
    else:
        ax = None
    for x in range(500): # initial number of rays
        ray = Ray(np.array([rnd.uniform(-sph_r,sph_r),rnd.uniform(-sph_r,sph_r),-sph_r-1]),wavelength=lam)
        trace(ray,sph,ax,draw3d) # recursive
    if draw3d:
        plt.show()

if __name__ == '__main__':
    main()