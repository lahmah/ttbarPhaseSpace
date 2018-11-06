import numpy as np
from math import gamma
from functools import partial
from scipy.optimize import newton
from vec4d import *

def floorToZero(a,N=0):
    return a
    return int(a*10**N)*10**-N

class Rambo(object):
    def __init__(self, nout, Ecms):
        self.nout = nout
        self.Ecms = Ecms

        self.ps_volume = (np.pi/2.)**(nout-1) \
                * Ecms**(2*nout-4) \
                /gamma(nout)/gamma(nout-1)# \
                #* (2 * np.pi) **(4-3*nout)

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

    def fxi(self,xi):
        val = 0
        for i in range(self.nout):
            val += np.sqrt(self.masses[i]**2+xi**2*self.p[i].E**2)
        return val -self.Ecms

    def dfxi(self,xi):
        val = 0
        for i in range(self.nout):
            denom = np.sqrt(self.masses[i]**2+xi**2*self.p[i].E**2)
            val += xi * self.p[i].E**2 / denom
        return val


    def generate_massive_point(self,masses):
        p = self.generate_point() 
        self.masses = masses
        self.p = p


        mass_sum = np.sum(masses)

        f = partial(self.fxi)
        df = partial(self.dfxi)

        xi_0 = np.sqrt(1-(mass_sum / self.Ecms)**2)
        #print(mass_sum)
        #print(xi_0)

        #print(self.fxi(xi_0))
        #print(self.dfxi(xi_0))


        Xi = newton(f,xi_0,df)

        #print(Xi)




        #Ki_abs = np.array([np.sqrt(p[i].E**2-masses[i]**2) for i in range(self.nout)])
        
        #Xi = sum(Ki_abs)
        #Xi /= self.Ecms

        k = []
        term1 = 0
        term2 = 0
        term3 = 1

        for i  in range(self.nout):
            k.append(Mom4D())
            k[i].E = np.sqrt(masses[i]**2+Xi**2*p[i].E**2)
            k[i].mom3d = Xi * p[i].mom3d 

            modulus = np.sqrt( np.sum((k[i].mom3d)**2))
            
            term1 += modulus / self.Ecms
            term2 += modulus**2 / k[i].E
            term3 *= modulus / k[i].E

        term1 = term1**(2*self.nout-3)
        term2 = term2**-1
        WW = self.ps_volume * term1*term2*term3*self.Ecms

        """
        W = self.ps_volume
        W *= Xi**(2*self.nout-3)    
        temp = 0
        for i in range(self.nout):
            W *= Ki_abs[i]/k[i].E
            temp += Ki_abs[i]**2/k[i].E
        W /= temp
        W *=self.Ecms
        """
        #print(W,WW)

        

        return k ,WW
