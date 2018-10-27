import numpy as np
from math import gamma
from vec4d import *

def floorToZero(a,N=0):
    return a
    return int(a*10**N)*10**-N

class Rambo(object):
    def __init__(self, nout, Ecms):
        self.nout = nout
        self.Ecms = Ecms

        self.ps_volume = (np.pi/2.)**(nout-1) * Ecms**(2*nout-4)/gamma(nout)/gamma(nout-1)

    def generate_weight(self):
        return self.ps_volume
    
    def generate_fourvector(self):
        ran = np.random.rand(4)
        c = 2.*ran[0]-1.
        phi = 2.*np.pi*ran[1]
        
        q = Mom4D()
        q[0] = -np.log(ran[2]*ran[3])
        q[1] = floorToZero(q[0]*np.sqrt(1-c**2)*np.cos(phi),12)
        q[2] = floorToZero(q[0]*np.sqrt(1-c**2)*np.sin(phi),12)
        q[3] = floorToZero(q[0]*c,12)
        
        return q
    
    def generate_point(self):
        q = []
        Q = Mom4D()
        for i in range(self.nout):
            q.append(self.generate_fourvector())
            Q += q[i]
        
        M = Q.m
        b = -Q.mom3d/M
        x = self.Ecms/M
        gamma = Q.E/M
        a = 1./(1.+gamma)
        
        p = []
        for i in range(self.nout):
            p.append(Mom4D())
            bdotq = b.dot(q[i].mom3d)
            p[i].E = x*(gamma*q[i].E + bdotq)    
            p[i].mom3d = x*(q[i].mom3d + b*q[i].E + a*bdotq*b)
    
        return p

    def generate_massive_point(self,pmasses,FixedMasses=False):
        p = self.generate_point()

        masses = pmasses.copy()
        w = 0
        nmassiv = 0 
        if not FixedMasses:
            for i, m in enumerate(masses):
                if m == 1:
                    nmassiv += 1
                    masses[i] = np.random.random()*p[i].E
                    w += 1/p[i].E

            w /= nmassiv
        if w == 0:
            w = 1


        Ki_abs = np.array([np.sqrt(p[i].E**2-masses[i]**2) for i in range(self.nout)])
        
        Xi = sum(Ki_abs)
        Xi /= self.Ecms

        k = []
        for i  in range(self.nout):
            k.append(Mom4D())
            k[i].E = np.sqrt(masses[i]**2+Xi**2*p[i].E**2)
            k[i].mom3d = Xi * p[i].mom3d 

        W = self.ps_volume
        W *= Xi**(2*self.nout-3)    
        temp = 0
        for i in range(self.nout):
            W *= Ki_abs[i]/k[i].E
            temp += Ki_abs[i]**2/k[i].E
        W /= temp
        W *=self.Ecms

        #W *= w

        return k ,W
