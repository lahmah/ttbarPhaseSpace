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
        MW = 80 #80.38
        GW = 2# 2.09
        MT = 175 #173
        GT = 1.5 #1.41
        rWmin = arctan(-(pow(MW,2))/(MW*GW))
        rTmin = arctan(-(pow(MT,2))/(MT*GT))
        
        if nout == 3:
            rTmin =  arctan(MW-(pow(MT,2))/(MT*GT))
       
        self.particles={             
        5  : (0,),#4.18,        # b-quark 
        6  : (MT,GT,rTmin),   # t-quark
        11 : (0,),              # elektron
        12 : (0,),              # e-neutrino
        13 : (0,),              # muon
        14 : (0,),              # mu-neutrino

        24 : (MW,GW,rWmin)    # wboson
        }


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
    
    def generate_Masses(self,Ecms,*argsv,FixedE=False):
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
                    if FixedE:
                        remainingE = Ecms
                    #mmax = m[0]+30*m[1]
                    #remainingE = mmax
                    mmin = 0 #m[0]-30*m[1]
                    rmax = arctan((remainingE**2-m[0]**2)/(m[0]*m[1])) 
                    rmin = arctan(((mmin)**2-m[0]**2)/(m[0]*m[1])) 
                    m = (m[0],m[1],rmin)

                    r = m[2] + np.random.random()*(rmax-m[2])
                    s = m[0]*m[1]*tan(r)+m[0]*m[0]
                    
                    W *= (rmax-m[2])*((s-m[0]**2)**2+(m[0]*m[1])**2)/(m[0]*m[1])  # inverse breit wigner weight

                    W /= m[0]*m[1]*tan(rmax)+m[0]**2 - (m[0]*m[1]*tan(rmin)+m[0]**2)


                    m = sqrt(s)        

                remainingE -= m
            masses.append(m)
        return array(masses) , W


    def generate_point(self):
        ttbar, BW_tt_weight   = self.generate_ttbar()

        Weight = BW_tt_weight

        if self.nout == 3 or True:
            OS = "OS"
        else:
            OS = ""
        if self.nout > 2:
            wb1,   BW_wb1_weight  = self.decay(ttbar[:,0],"24"+OS,"5OS")
            Weight *= BW_wb1_weight
        if self.nout > 3:
            enu,_  = self.decay(wb1[:,0],"11OS","12OS")
        if self.nout == 6:
            wb2,   BW_wb2_weight  = self.decay(ttbar[:,1],"24","5OS")
            munu,_ = self.decay(wb2[:,0],"11OS","12OS")
            Weight *=BW_wb2_weight

        if self.nout == 2:
            pout = [ttbar[:,0],ttbar[:,1]]
        elif self.nout == 3:
            pout = [wb1[:,0],  ttbar[:,1], wb1[:,1]] 
        elif self.nout == 4:
            pout = [wb1[:,1],ttbar[:,0],enu[:,0],enu[:,1]]
        elif self.nout == 6:
            pout = [wb1[:,1],wb2[:,1],munu[:,0],munu[:,1],enu[:,0],enu[:,1]]

        Weight *= self.generate_weight(pout)
        return pout, Weight, BW_tt_weight, BW_wb1_weight, self.generate_weight(pout)





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

         
        masses,BW_weight = self.generate_Masses(self.Ecms/2,*particles,FixedE=True)

        if self.debug:
            print("ttbar mass: ",masses)
        ttbar = self.generate_k(masses)

        return ttbar , BW_weight

    def decay(self, pint, *pout):
        Ecms = pint.m

        masses,BW_weight = self.generate_Masses(Ecms,*pout) 

        p = Mom4D(np.empty((4,2)))

        momentum = Ecms/2 * (1-sum(array(masses)**2)/Ecms**2)
        s = masses[0]**2-masses[1]**2 
        p[0,0] = (1+s/Ecms**2) * Ecms/2
        p[0,1] = (1-s/Ecms**2) * Ecms/2

        theta, phi = np.random.random(2)
        theta *= pi
        phi   += 2*pi

        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)

        momentum *= array([x,y,z])
        p[1:,0] = +momentum
        p[1:,1] = -momentum
        
        if self.debug:
            print("w,b befor boost:")
            print(p, " - " , p.m)




        m_rot = self.rotat(pint,self.labsystem)
        p[:,0] = self.rotat_inv(p[:,0],m_rot)._arr
        p[:,1] = self.rotat_inv(p[:,1],m_rot)._arr
        p[:,0] = self.boost(pint,p[:,0])._arr
        p[:,1] = self.boost(pint,p[:,1])._arr


        if self.debug:
            print("w,b after boost:")
            print(p, " - " , p.m)


        return p,BW_weight



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

