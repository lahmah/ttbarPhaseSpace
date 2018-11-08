import numpy as np
from math import gamma
from functools import partial
from scipy.optimize import newton
from vec4d import *

def floorToZero(a,N=0):
    return a
    return int(a*10**N)*10**-N

class topGenerator (object):
    def __init__(self, nout, Ecms):
        self.nout = nout
        self.Ecms = Ecms

        self.ps_volume = (np.pi/2.)**(nout-1) \
                * Ecms**(2*nout-4) \
                /gamma(nout)/gamma(nout-1)# \

    def generate_weight(self,pout):
        term1 = 0
        term2 = 0
        term3 = 1

        for i  in range(self.nout): 
            modulus = np.sqrt( np.sum((pout[:,i].mom3d)**2))
            
            term1 += modulus / self.Ecms
            term2 += modulus**2 / pout[:,i].E
            term3 *= modulus / pout[:,i].E

        term1 = term1**(2*self.nout-3)
        term2 = term2**-1
        return self.ps_volume * term1*term2*term3*self.Ecms
    
    def generate_fourvectors(self,n):
        ran = np.random.rand(4,n)
        c = 2.*ran[0]-1.
        phi = 2.*np.pi*ran[1]
        
        q = np.empty((4,n))
        q[0] = -np.log(ran[2]*ran[3])
        q[1] = q[0]*np.sqrt(1-c**2)*np.cos(phi)
        q[2] = q[0]*np.sqrt(1-c**2)*np.sin(phi)
        q[3] = q[0]*c
      
        return Mom4D(q)
    
    def generate_p(self,n):
        q = generate_fourvectors(n) 
        Q = Mom4D(np.sum(q._arr,1))
       
        M = Q.m
        b = -Q.mom3d/M
        x = self.Ecms/M
        gamma = Q.E/M
        a = 1./(1.+gamma)
        
        p = np.empty(4,n) 

        bdotq = np.sum(b * qq.mom3d.T,axis = 1) 
        p[0] = x* (gamma*q.E+bdotq)
        p[1:] = x * (q.mom3d + b[:,None]*q.E + a*bdotq*b[:,None])
   
        return Mom4D(p)

    def fxi(self,p,masses,xi):
        nout = len(masses)
        val = 0
        for i in range(nout):
            val += np.sqrt(masses[i]**2+xi**2*p[:,i].E**2)
        return val -self.Ecms

    def dfxi(self,p,masses,xi):
        nout = len(masses)
        val = 0
        for i in range(nout):
            denom = np.sqrt(masses[i]**2+xi**2*p[:,i].E**2)
            val += xi * p[:,i].E**2 / denom
        return val


    def generate_massiv_point(self,masses):
        n = len(masses)

        p = self.generate_p(n) 
        
        mass_sum = np.sum(masses)

        f = partial(self.fxi,p,masses)
        df = partial(self.dfxi,p,masses)

        xi_0 = np.sqrt(1-(mass_sum / self.Ecms)**2)

        Xi = newton(f,xi_0,df)

        k = np.empty((4,n))
        k[0] = np.sqrt(P.E**2*Xi**2+masses**2)
        
        return Mom4D(k)


    def generate_ttbar(self): 
        if self.nout == 2:
            particles=("6OS","-6OS")
        else:
            particles=("6","-6")

         
        masses,BW_weight = self.generate_Masses(self.Ecms,*particles)

        ttbar = generate_massiv_point(masses)

        return ttbar , BW_weight
