import numpy as np
from numpy import pi, sin, cos, arctan, tan, log, sqrt, array
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

        self.ps_volume = (pi/2.)**(nout-1) \
                * Ecms**(2*nout-4) \
                /gamma(nout)/gamma(nout-1)# \

        # Particle masses and widths
        MW = 80.38
        GW = 2.09
        MT = 173
        GT = 1.41
        rWmin = arctan(-(pow(MW,2))/(MW*GW))
        rTmin = arctan(-(pow(MT,2))/(MT*GT))
        
        if nout == 3:
            rTmin =  arctan(MW-(pow(MT,2))/(MT*GT))
       
        self.particles={             
        5  : (0,),#4.18,        # b-quark 
        6  : (MT,GT,rTmin),   # t-quark
        11 : (0,),              # elektron
        12 : (0,),              # e-eutrino
        24 : (MW,GW,rWmin)    # wboson
        }


    def generate_weight(self,pout):
        term1 = 0
        term2 = 0
        term3 = 1

        for i  in range(self.nout): 
            modulus = sqrt( np.sum((pout[:,i].mom3d)**2))
            
            term1 += modulus / self.Ecms
            term2 += modulus**2 / pout[:,i].E
            term3 *= modulus / pout[:,i].E

        term1 = term1**(2*self.nout-3)
        term2 = term2**-1
        return self.ps_volume * term1*term2*term3*self.Ecms
    
    def generate_Masses(self,Ecms,*argsv):
        masses = []             
        remainingE = Ecms
        
        W = 1
        for arg in argsv:
            if arg.endswith("OS"):
                id = abs(int(arg[:-2]))
                m = self.particles.get(id)
                if isinstance(m,tuple): 
                    m = m[0]
                remainingE -= m
        for arg in argsv:        
            if arg.endswith("OS"):
                id = abs(int(arg[:-2]))
                m = self.particles.get(id)
                if isinstance(m,tuple):
                    m = m[0]
            else:    
                id = abs(int(arg)) 
                m = self.particles.get(id)
                if isinstance(m,tuple):
                    if id == 6:
                        remainingE = 500
                    rmax = arctan((remainingE**2-m[0]**2)/(m[0]*m[1])) 
                    r = m[2] + np.random.random()*(rmax-m[2])
                    s = m[0]*m[1]*tan(r)+m[0]*m[0]
                    
                    W *= (rmax-m[2])*((s-m[0]**2)**2+(m[0]*m[1])**2)/(m[0]*m[1])  # inverse breit wigner weight

                    W *=  1./(2.*remainingE*128.*pow(pi,3))*(1.-s/remainingE**2) 

                    m = sqrt(s)        

                remainingE -= m
            masses.append(m)
        return array(masses) , W



    def generate_fourvectors(self,n):
        ran = np.random.rand(4,n)
        c = 2.*ran[0]-1.
        phi = 2.*pi*ran[1]
        
        q = np.empty((4,n))
        q[0] = -log(ran[2]*ran[3])
        q[1] = q[0]*sqrt(1-c**2)*cos(phi)
        q[2] = q[0]*sqrt(1-c**2)*sin(phi)
        q[3] = q[0]*c
      
        return Mom4D(q)
    
    def generate_p(self,n):
        q = self.generate_fourvectors(n) 
        Q = Mom4D(np.sum(q._arr,1))
       
        M = Q.m
        b = -Q.mom3d/M
        x = self.Ecms/M
        gamma = Q.E/M
        a = 1./(1.+gamma)
        
        p = np.empty((4,n)) 

        bdotq = np.sum(b * q.mom3d.T,axis = 1) 
        p[0] = x* (gamma*q.E+bdotq)
        p[1:] = x * (q.mom3d + b[:,None]*q.E + a*bdotq*b[:,None])
   
        return Mom4D(p)

    def fxi(self,p,masses,xi):
        nout = len(masses)
        val = 0
        for i in range(nout):
            val += sqrt(masses[i]**2+xi**2*p[:,i].E**2)
        return val -self.Ecms

    def dfxi(self,p,masses,xi):
        nout = len(masses)
        val = 0
        for i in range(nout):
            denom = sqrt(masses[i]**2+xi**2*p[:,i].E**2)
            val += xi * p[:,i].E**2 / denom
        return val


    def generate_k(self,masses):
        n = len(masses)

        p = self.generate_p(n) 
        
        mass_sum = np.sum(masses)

        f = partial(self.fxi,p,masses)
        df = partial(self.dfxi,p,masses)

        xi_0 = sqrt(1-(mass_sum / self.Ecms)**2)

        Xi = newton(f,xi_0,df)

        k = np.empty((4,n))
        k[0] = sqrt(p.E**2*Xi**2+masses**2)
        k[1:] = Xi*p.mom3d
        
        return Mom4D(k)


    def generate_ttbar(self): 
        if self.nout == 2:
            particles=("6OS","-6OS")
        else:
            particles=("6","-6")

         
        masses,BW_weight = self.generate_Masses(self.Ecms,*particles)

        ttbar = self.generate_k(masses)

        return ttbar , BW_weight


