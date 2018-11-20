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
    def __init__(self, nout, Ecms, debug = False):
        self.nout = nout
        self.Ecms = Ecms
        self.debug = debug

        self.ps_volume = (pi/2.)**(nout-1) \
                * Ecms**(2*nout-4) \
                /gamma(nout)/gamma(nout-1)# \


        self.labsystem = Mom4D([0,0,0,1])

        # Particle masses and widths
        self.MW = 80.385
        self.GW = 2.085
        self.MT = 173.21
        self.GT = 2
        self.Mb = 0#4.8
        self.Me = 0#0.000511
        self.Mmu = 0#0.105

    def generate_weight(self,pout):
        term1 = 0
        term2 = 0
        term3 = 1

        for i  in range(self.nout): 
            modulus = sqrt( np.sum((pout[i].mom3d)**2))
            
            term1 += modulus / self.Ecms
            term2 += modulus**2 / pout[i].E
            term3 *= modulus / pout[i].E

        term1 = term1**(2*self.nout-3)
        term2 = term2**-1
        return self.ps_volume * term1*term2*term3*self.Ecms
    
    def generate_Mass(self,mass,gamma,mmax,mmin = 0):
        
        Smin = mmin**2
        Smax = mmax**2
        M2 = mass**2
        MG = mass*gamma
        rmin = arctan((Smin-M2)/MG)
        rmax = arctan((Smax-M2)/MG)

        # generate r
        r = rmin+np.random.random()*(rmax-rmin)
        # generate s                                   
        s = MG*tan(r)+M2     
        
        weight = (rmax-rmin)*((s-M2)**2+MG**2)/MG  
        #weight *= np.pi/(rmax-rmin)

        weight /= MG*(tan(rmax) - tan(rmin))

        if np.sqrt(s) == 0:
            print("s = 0 - ",mass,gamma)

        return np.sqrt(s) , weight

    def generate_point(self):
        Weight = 1
        ttbar,tt_BW_Weight = self.generate_ttbar()

        Weight *= tt_BW_Weight
        
        if self.nout == 2:
            pout = ttbar
        elif self.nout == 3: # decay anti top - to W b
            masses = array([self.MW,self.Mb])
            w1 , b1, wb_D_weight = self.decay(ttbar[1],masses)
            #Weight *= wb_D_weight
            pout = [w1, ttbar[0], b1]
        elif self.nout == 4:
            mw1, w1_BW_Weight = self.generate_Mass(self.MW,self.GW,ttbar[1].m-self.Mb,self.Me)
            Weight *= w1_BW_Weight

            masses = array([mw1,self.Mb])
            w1, b1, wb_D_weight = self.decay(ttbar[1],masses)
            masses = array([self.Me,0])
            e, enu, enu_D_weight = self.decay(w1,masses)

            #Weight *= wb_D_weight*enu_D_weight
            pout = [enu, ttbar[0], e, b1]
        elif self.nout == 6:
            mw1, w1_BW_Weight = self.generate_Mass(self.MW,self.GW,ttbar[1].m-self.Mb,self.Me)
            Weight *= w1_BW_Weight

            masses = array([mw1,self.Mb])
            w1, b1, wb1_D_weight = self.decay(ttbar[1],masses)
            masses = array([self.Me,0])
            e, enu, enu_D_weight = self.decay(w1,masses)

            mw2, w2_BW_Weight = self.generate_Mass(self.MW,self.GW,ttbar[0].m-self.Mb,self.Mmu)
            Weight *= w2_BW_Weight

            masses = array([mw2,self.Mb])
            w2, b2, wb2_D_weight = self.decay(ttbar[0],masses)
            masses = array([self.Mmu,0])
            mu, munu, munu_D_weight = self.decay(w2,masses)
            #Weight *= wb1_D_weight*wb2_D_weight*enu_D_weight*munu_D_weight




            #Weight *= wb_D_weight*enu_D_weight
            pout = [munu,b2,enu,mu,e,b1]            
        
        Weight *= self.generate_weight(pout)

        if self.debug:
            print("Output")
            for i in range(len(pout)):
                print(pout[i], " - ", pout[i].m)
            for i in range(len(pout)):
                for j in range(i+1,len(pout)):
                    print(i, " - ",j ,": ",angle(pout[i],pout[j]))
            print(Weight)
            print("----------------")

        return pout, Weight

    def generate_q(self,n):
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
        q = self.generate_q(n) 
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

        try:
            Xi = newton(f,xi_0,df)
        except:
            Xi =xi_0
            print("Faild Newton with: ",p)

        k = np.empty((4,n))
        k[0] = sqrt(p.E**2*Xi**2+masses**2)
        k[1:] = Xi*p.mom3d

        
        return [Mom4D(k[:,0]) ,Mom4D(k[:,1])]


    def generate_ttbar(self): 
        
        if self.nout == 2:
            m1 = self.MT
            m2 = m1
            BW_weight = 1
        elif self.nout in [3,4]:  
            Mmax = self.Ecms
            if self.nout == 3:
                Mmin = self.MW+self.Mb
            else:
                Mmin = self.Mb+self.Me
            m1 = self.MT
            m2, BWw2 = self.generate_Mass(self.MT,self.GT,Mmax-m1,Mmin)
            BW_weight = BWw2 

        else: 
            Mmax = self.Ecms
            Mmin = self.Mb
            m1, BWw1 = self.generate_Mass(self.MT,self.GT,Mmax,Mmin+self.Mmu)
            m2, BWw2 = self.generate_Mass(self.MT,self.GT,Mmax-m1,Mmin+self.Me)
            BW_weight = BWw1*BWw2 

        masses=np.array([m1,m2])

        ttbar = self.generate_k(masses)

        return ttbar , BW_weight

    def decay(self, pint, masses):
        Ecms = pint.m

        p = [Mom4D(),Mom4D()] 

        momentum = Ecms/2 * (1-sum(array(masses)**2)/Ecms**2)

        s = masses[0]**2-masses[1]**2 
        p[0].E = (1+s/Ecms**2) * Ecms/2
        p[1].E = (1-s/Ecms**2) * Ecms/2

        cosT = 2*np.random.random()-1
        theta = np.arccos(cosT)
        phi = np.random.random()
        phi   *= 2*pi

        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)

        momentum *= array([x,y,z])
        p[0].mom3d = +momentum
        p[1].mom3d = -momentum
        
        if self.debug:
            print("befor boost:")
            print(p[0], " - " , p[0].m)
            print(p[1], " - " , p[1].m)


        m_rot = self.rotat(pint,self.labsystem)
        p[0] = self.rotat_inv(p[0],m_rot)
        p[1] = self.rotat_inv(p[1],m_rot)
        p[0] = self.boost(pint,p[0])
        p[1] = self.boost(pint,p[1])

                
        


        if self.debug:
            print("after boost:") 
            print(p[0], " - " , p[0].m)
            print(p[1], " - " , p[1].m)

        def lamda(P,p1,p2):
            S = P.m**2
            s1 = p1.m**2
            s2 = p2.m**2
            return (S-s1-s2)**2-4*s1*s2

        Decay_Weight = np.pi * np.sqrt(lamda(pint,p[0],p[1]))/2/pint.m**2
        return p[0], p[1], Decay_Weight



    def boost(self, q, ph):
    #                                      _
    # Boost of a 4-vector ( relative speed q/q(0) ):
    #
    # ph is the 4-vector in the rest frame of q
    # p is the corresponding 4-vector in the lab frame
    #
    # INPUT     OUTPUT
    # q, ph     p

        p = Mom4D()
    
        rsq = q.m
    
        p.E = (q.E*ph.E+np.dot(q.mom3d, ph.mom3d)) / rsq
        c1 = (ph.E+p.E) / (rsq+q.E)
        p.mom3d = ph.mom3d + c1*q.mom3d
    
        return p
    
    def boost_inv(self, q, p):
    #                                      _
    # Boost of a 4-vector ( relative speed q/q(0) ):
    #
    # ph is the 4-vector in the rest frame of q
    # p is the corresponding 4-vector in the lab frame
    #
    # INPUT     OUTPUT
    # q, p      ph
    
        ph = Mom4D()
        rsq = q.m
    
        ph.E = q*p/rsq
        c1 = (p.E+ph.E) / (rsq+q.E)
        ph.mom3d = p.mom3d - c1*q.mom3d
    
        return ph

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
        pp = [Mom4D(), Mom4D]
    
        pm[0] = sqrt(p1.mom3d.dot(p1.mom3d))
        pm[1] = sqrt(p2.mom3d.dot(p2.mom3d))
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

        p1 = Mom4D()
    
        p1.E = p2.E
        p1.mom3d = np.matmul(rot,p2.mom3d)
    
        return p1

