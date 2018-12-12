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
    def __init__(self, nout, Ecms,cexp, debug = False):
        self.nout = nout
        self.Ecms = Ecms
        self.debug = debug

        self.ps_volume = (pi/2.)**(nout-1) \
                * Ecms**(2*nout-4) \
                /gamma(nout)/gamma(nout-1)

        #self.ps_volume = ((2.*pi)**(4.-3*nout))
        self.cexp = cexp


        # Particle masses and widths
        self.MW = 80.385
        self.GW = 2.085
        self.MT = 173.21
        self.GT = 2
        self.Mb = 0#4.8
        self.Me = 0#0.000511
        self.Mmu = 0#0.105

        # minkowski metrik
        self.mk=np.array([1,-1,-1,-1])

    def generate_point(self,size=1):
        pin = np.array([[self.Ecms],[0],[0],[0]])
        weight = self.ps_volume

        if self.nout == 2:
            masses = np.full((2,size),self.MT)
            ttbar, Decay_Weight= self.decay(pin,masses,self.Ecms,True)
            weight *= Decay_Weight

            return ttbar, weight

        if self.nout == 3:
            masses = np.empty((2,size))
            masses[0] = self.MT
            masses[1], BW_Weight = self.generate_mass(self.MT,self.GT,self.Ecms-self.MT,self.MW+self.Mb,size)
            ttbar, Decay_Weight = self.decay(pin,masses,self.Ecms,True)
            weight *=BW_Weight * Decay_Weight
            masses[0] = self.MW
            masses[1] = self.Mb
            wb, Decay_Weight = self.decay(ttbar[:,1],masses,self.MT)
            weight *= Decay_Weight

            pout = np.append(np.append(wb[:,0,None],ttbar[:,0,None],axis=1),wb[:,1,None],axis=1)
            return pout , weight

        if self.nout == 4:
            masses = np.empty((2,size))
            masses[0] = self.MT
            masses[1], BW_Weight = self.generate_mass(self.MT,self.GT,self.Ecms-self.MT,self.Me+self.Mb,size)
            ttbar, Decay_Weight = self.decay(pin,masses,self.Ecms,True)
            weight *=BW_Weight * Decay_Weight
            masses[0] = self.generate_mass(self.MW,self.GW,masses[1]-self.Mb,self.Me,size)
            masses[1] = self.Mb
            wb, Decay_Weight = self.decay(ttbar[:,1],masses,self.MT)
            weight *= Decay_Weight
            masses[0] = self.Me
            masses[1] = 0
            enu, Decay_Weight = self.decay(wb[:,0],masses,self.MW)
            weight *= Decay_Weight
            
            pout = np.append(np.append(np.append(enu[:,1,None],ttbar[:,0,None],axis=1),enu[:,0,None],axis=1),wb[:,1,None],axis=1)
            return pout, weight

 







    def generate_mass(self,mass,gamma,mmax,mmin=0,size=1):
        if hasattr(mmax, "__len__"):
            size = len(mmax)

        mass2 = mass**2
        mw = mass*gamma
        smin = mmin**2
        smax = mmax**2
        ymax=arctan((smin-mass2)/mw)
        ymin=arctan((smax-mass2)/mw)
        
        s = mass2+mw*tan(ymin + np.random.rand(size)*(ymax-ymin))

        wt= (ymin-ymax)*((s-mass2)**2+mw**2)/mw
        wt/=(smax-smin)
        #wt /= np.sqrt(s)

        return np.sqrt(s), wt
    
    def mass(self,p):
        mk = np.reshape(self.mk,(4,*tuple(np.ones(len(np.shape(p))-1,dtype=int))))
        m2 = np.sum(p**2*mk,axis=0)
        m = np.sign(m2)*np.sqrt(np.abs(m2))
        return m
       
    def sqlambda(self,s,s1,s2):
        return np.sqrt(((s-s1-s2)**2-4*s1*s2))/s

    def Lambda(self,a,b,c):
        L = (a-b-c)**2 -4*b*c
        return sqrt(L)/2/np.sqrt(a)

    def DecayMassWeight(self,s,sp,b,c):
        return self.Lambda(sp,b,c)/self.Lambda(s,b,c) *s/sp


    def decay(self,pint, masses,m_pin,pin_const =False):
        Ecms = self.mass(pint)
        if pin_const:
            size = len(masses[0])
        else:
            size = len(Ecms)
        

        s = self.mass(pint)**2  
        s1 = masses[0]**2
        s2 = masses[1]**2 

        p = np.empty((4,2,size))
        p[0,0] = (s+s1-s2)/Ecms/2
        p1m = Ecms * self.sqlambda(s,s1,s2)/2

        ct  = 2*np.random.random(size)-1
        st  = sqrt(1-ct**2)
        phi   = 2*pi*np.random.random(size)

        p[1:,0]=p1m * array([st*sin(phi),st*cos(phi),ct])
#        for i in range(size): 
#            rot = self.rotat(np.array([1,0,0,1]),pint[:,i])
 #           p[:,0,i] = self.rotat_inv(p[:,0,i],rot)

        p[:,0] = self.boost(pint,p[:,0])
        p[:,1] = pint - p[:,0]

        #Decay_Weight = pi*self.sqlambda(s,s1,s2)/2
        Decay_Weight = self.DecayMassWeight(m_pin**2,Ecms**2,masses[0]**2,masses[1]**2) 

        return p, Decay_Weight

    def PeakedDist(self,a, cn, cxm, cxp,k,size):
        ce = 1-cn
        res = np.random.rand(size)

        CTZ = np.isclose(0,ce)
        NCTZ = np.logical_not(CTZ)

        res[CTZ]  = k *( (a+k*cxm)*( (a+k*cxp)/(a+k*cxm))**res[CTZ]  - a)
        res[NCTZ] = k * ((res[NCTZ]*(a+k*cxp)**ce+(1.-res[NCTZ])*(a+k*cxm)**ce)**(1/ce)-a)
        
        weight = np.empty(size)
        weight[CTZ] = log((a+k*cxp)/(a+k*cxm))/k
        weight[NCTZ] = ((a+k*cxp)**ce-(a+k*cxm)**ce)/(k*ce)

        weight *= (a+res)**cn
        return res, weight

    def AnisotropicDecay(self, pint, masses, ctexp=.8, ctmin=-1, ctmax=1):
        Ecms = self.mass(pint)
        s = Ecms**2
        s1 = masses[0]**2
        s2 = masses[1]**2
        
        ctexp=self.cexp

        size = len(Ecms)
        
        p = np.empty((4,2,size))
        p[0,0] = (s+s1-s2)/Ecms/2.

        p1m = Ecms*self.sqlambda(s,s1,s2)/2
        pim = sqrt(pint[0]**2-s)
        a   = pint[0]*p[0,0]/pim/p1m
        a[np.logical_and(1.>=a,a>=0.)] = 1.0000000001
        ct, PeakedWeight     = self.PeakedDist(a,ctexp,ctmin,ctmax,1,size)
        st = sqrt(1.-ct**2)
        phi = 2.*pi*np.random.rand(size)
        p[1:,0]=p1m * array([st*sin(phi),st*cos(phi),ct])
        print(self.mass(p[:,0]))
  
        pref = np.empty((4,size))
        pref[0] = pint[0]
        pref[1:3] = 0
        pref[3] = pim
        
        ref = np.array([1,0,0,1])
        p[:,0] = self.Poincare(ref[:,None],pint,p[:,0])
        p[:,0] = self.boost(pint,p[:,0])
        #p[:,0] = self.boost(pref,p[:,0])
        #self.Poincare(pref,pint,p[:,0])
        
        p[:,1] = pint - p[:,0]
        print(self.mass(p[:,0]))
        print(self.mass(p[:,1]))

        Decay_Weight =  (pi*self.sqlambda(s,s1,s2)/2* PeakedWeight)

        return p, Decay_Weight





            
    def boost(self, q, ph):
    #
    # Boost of a 4-vector ( relative speed q/q(0) ):
    #
    # ph is the 4-vector in the rest frame of q
    # p is the corresponding 4-vector in the lab frame
    #
    # INPUT     OUTPUT
    # q, ph     p

        p = np.empty(np.shape(ph))
    
        rsq = self.mass(q)
    
        p[0]  = (q[0]*ph[0]+np.sum(q[1:]*ph[1:],axis=0)) / rsq
        c1 = (ph[0]+p[0]) / (rsq+q[0])
        p[1:] = ph[1:] + c1*q[1:]
    
        return p
    
    def boost_inv(self, q, p):
    #
    # Boost of a 4-vector ( relative speed q/q(0) ):
    #
    # ph is the 4-vector in the rest frame of q
    # p is the corresponding 4-vector in the lab frame
    #
    # INPUT     OUTPUT
    # q, p      ph
    
        ph = np.empty(np.shape(q))
        rsq = self.mass(q)
    
        ph[0] = q*p/rsq
        c1 = (p[0]+ph[0]) / (rsq+q[0])
        ph[1:] = p[1:] - c1*q[1:]
    
        return ph



    def Poincare(self,v1,v2,v3): 
        b = v1[1:]  
        a = v2[1:] 

        ct=np.sum(a*b,axis=0)/np.sqrt(np.sum(a**2,axis=0) * np.sum(b**2,axis=0));
        n  = np.cross(a,b,axis=0) 

        nsq=np.sum(n**2,axis=0)
        nabs=sqrt(nsq)
        ct[ct>1]=1
        ct[ct<-1]=-1

        st=-sqrt(1.0-ct**2) 

        c=v3[1:]
        at = n*np.sum(n*c,axis=0)/nsq
        ap = c-at
        return np.append(v3[0,None],at+ct*ap+st/nabs*np.cross(n,ap,axis=0),axis=0)


    def rotat(self, p1, p2):
    # Rotation of a 4-vector:
    #
    #            p1 = rot*p2
    #
    # INPUT     OUTPUT
    #
    # p1, p2    rot

        rot = np.empty((3, 3))

        r = np.empty((2, 3, 3))
        pm = np.empty(2)
        sp = np.empty(2)
        cp = np.empty(2)
        st = np.empty(2)
        ct = np.empty(2)
        pp = np.empty((2,4))

        pm[0] = sqrt(p1[1:].dot(p1[1:]))
        pm[1] = sqrt(p2[1:].dot(p2[1:]))
        pp[0] = (1./pm[0])*p1
        pp[1] = (1./pm[1])*p2

        for i in range(2):
            ct[i] = pp[i][3]
            st[i] = sqrt(1.-ct[i]*ct[i])
            if np.isclose(abs(ct[i]), 1.):
                cp[i] = 1.
                sp[i] = 0.
            else:
                cp[i] = pp[i][2] / st[i]
                sp[i] = pp[i][1] / st[i]

            r[i, 0, 0] = cp[i]
            r[i, 0, 1] = sp[i]*ct[i]
            r[i, 0, 2] = st[i]*sp[i]
            r[i, 1, 0] = -sp[i]
            r[i, 1, 1] = ct[i]*cp[i]
            r[i, 1, 2] = cp[i]*st[i]
            r[i, 2, 0] = 0.
            r[i, 2, 1] = -st[i]
            r[i, 2, 2] = ct[i]


        rot = np.matmul(r[0],r[1].T) 
        return rot

    
    
    def rotat_inv(self, p2, rot):
    # Rotation of a 4-vector:
    #
    #            p1 = rot*p2
    #
    # INPUT     OUTPUT
    #
    # p2, rot   p1

        p1 = np.zeros(4)
    
        p1[0] = p2[0]
        p1[1:] = np.matmul(rot,p2[1:])
    
        return p1

