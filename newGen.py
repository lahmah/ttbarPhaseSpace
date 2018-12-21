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


    def generate_point(self,size):
        Weight = 1
        ttbar,tt_BW_Weight = self.generate_ttbar(size)

        Weight*=tt_BW_Weight

        if self.nout == 2:
            pout = ttbar
        elif self.nout == 3:
            masses = np.empty((2,size))
            masses[0] = self.MW
            masses[1] = self.Mb

            wb, wb_D_weight = self.decay(ttbar[:,1],masses,self.MT) 
            #wb, wb_D_weight = self.AnisotropicDecay(ttbar[:,1],masses)
            Weight *= wb_D_weight

            pout = np.append(np.append(wb[:,0,None],ttbar[:,0,None],axis=1),wb[:,1,None],axis=1)

        elif self.nout == 4:
            masses = np.empty((2,size))
            masses[0], w1_BW_Weight = self.generate_mass(self.MW,self.GW,self.mass(ttbar[:,1])-self.Mb,self.Me)
            masses[1] = self.Mb
            Weight *= w1_BW_Weight

            wb, wb_D_weight = self.decay(ttbar[:,1],masses,self.MT)
            masses[0] = self.Me
            masses[1] = 0
            enu, enu_D_weight= self.decay(wb[:,0],masses,self.MW)

            Weight *= wb_D_weight*enu_D_weight
            pout = np.append(np.append(np.append(enu[:,1,None],ttbar[:,0,None],axis=1),enu[:,0,None],axis=1),wb[:,1,None],axis=1)
            #Weight /= np.prod(pout[0],axis=0)
        elif self.nout == 6:
            masses = np.empty((2,size))
            masses[0], w1_BW_Weight = self.generate_mass(self.MW,self.GW,self.mass(ttbar[:,1])-self.Mb,self.Me)
            masses[1] = self.Mb
            Weight *= w1_BW_Weight
            
            
            wb1, wb1_D_weight = self.decay(ttbar[:,1],masses,self.MT)
            masses[0] = self.Me
            masses[1] = 0
            enu, enu_D_weight = self.decay(wb1[:,0],masses,self.MW)

            masses[0], w2_BW_Weight = self.generate_mass(self.MW,self.GW,self.mass(ttbar[:,0])-self.Mb,self.Mmu)
            masses[1] = self.Mb
            Weight *= w2_BW_Weight
            
            wb2, wb2_D_weight = self.decay(ttbar[:,0],masses,self.MT)
            masses[0] = self.Mmu
            masses[1] = 0
            munu, munu_D_weight = self.decay(wb2[:,0],masses,self.MW)
            Weight *= wb1_D_weight*wb2_D_weight*enu_D_weight*munu_D_weight

            pout = np.append(munu[:,1,None],wb2[:,1,None],axis=1)
            pout = np.append(pout,enu[:,1,None],axis=1)
            pout = np.append(pout,munu[:,0,None],axis=1)
            pout = np.append(pout,enu[:,0,None],axis=1)
            pout = np.append(pout,wb1[:,1,None],axis=1)

        
        Weight *= self.generate_weight(pout)
        

        return pout, Weight    
    
    def generate_weight(self,pout):
        if self.nout == 2:
            modulus = sqrt(np.sum(pout[1:]**2,axis=0))
            term1 = np.sum(modulus/self.Ecms,axis=0)**(2*self.nout-3)
            term2 = 1/np.sum(modulus**2/pout[0],axis=0)
            term3 = np.prod(modulus / pout[0],axis=0)
            return self.ps_volume * term1*term2*term3*self.Ecms
        return self.ps_volume #* term1*term2*term3*self.Ecms 


    def pol3(self,s,m,smax):
        x = sqrt(s)
        w1 =(2*np.pi)**-2 *(1-(x/m))**2+1
        w1 = 1/ w1
        norm1 = np.sum(w1)/len(x)
        w2 = sqrt(1-s/smax)
        norm2 = np.sum(w2)/len(x)
        w = w1*w2
        norm = norm1 * norm2
        print("pol norm: ",norm)
        return w/norm

    def nop_generate_mass(self,mass,gamma,mmax,mmin = 0,size=1):
        if hasattr(mmax, "__len__"):
            size = len(mmax)
        Smin = mmin**2
        Smax = mmax**2
        M2 = mass**2
        MG = mass*gamma
        rmin = arctan((Smin-M2)/MG)
        rmax = arctan((Smax-M2)/MG)

        # generate r
        r = rmin+np.random.random(size)*(rmax-rmin)
        # generate s                                   
        s = MG*tan(r)+M2     
        
        #weight = (rmax-rmin)*((s-M2)**2+MG**2)/MG  
        ##weight *= np.pi/(rmax-rmin)

        #weight /= MG*(tan(rmax) - tan(rmin))
        ##weight /= np.sqrt(s)
        ##weight *= self.pol(np.sqrt(s)) 

        upper  = (Smax-M2)/MG      
        lower  = (Smin-M2)/MG        
        weight=((s-M2)**2+MG**2)/MG
        weight*=arctan(upper)-arctan(lower)             


        #print( min(weight), " - " ,max(weight))
        #print( min(np.sqrt(s)), " - " ,max(np.sqrt(s)))
        #print(Smin)
        #print(Smax)

        return np.sqrt(s) , weight

    def generate_mass(self,mass,gamma,mmax,mmin=0,size=1):
        if hasattr(mmax, "__len__"):
            size = len(mmax)

        mass2 = mass**2
        mw = mass*gamma
        smin = mmin**2
        smax = mmax**2
        ymin=arctan((smin-mass2)/mw)
        ymax=arctan((smax-mass2)/mw)  
        
        r = ymin+np.random.random(size)*(ymax-ymin)
    
        s = mass2+mw*tan(r)

        wt= (ymax-ymin)*((s-mass2)**2+mw**2)/mw
        wt/=(smax-smin)
        #wt /= np.sqrt(s)

        print("BW: Norm: ",np.sum(wt)/size)

        if self.nout == 3 or mass == self.MT:
            wt*=self.pol3(s,mass,smax)

        if mass == self.MT:
            wt *= s/mmax**2
        #if mass == self.MW:
        #    wt *= np.sqrt(1-s/mmax**2)

        return np.sqrt(s), wt#/norm

    def nop_generate_mass(self,mass,gamma,mmax,mmin=0,size=1):
        if hasattr(mmax, "__len__"):
            size = len(mmax)

        smax = mmax**2
        smin = mmin**2
        s = np.random.uniform(smin,smax,size)
        weight = 1/(smax-smin)

        #m = np.random.uniform(mmin,mmax,size)
        #s = m**2
        #weight = 1/(mmax-mmin)

        return np.sqrt(s), weight 







    
    def mass(self,p):
        mk = np.reshape(self.mk,(4,*tuple(np.ones(len(np.shape(p))-1,dtype=int))))
        m2 = np.sum(p**2*mk,axis=0)
        m = np.sign(m2)*np.sqrt(np.abs(m2))
        return m
    
    def generate_q(self,n):
        ran = np.random.rand(4,*n)
        c = 2.*ran[0]-1.
        phi = 2.*pi*ran[1]
        
        q = np.empty((4,*n))
        q[0] = -log(ran[2]*ran[3])
        q[1] = q[0]*sqrt(1-c**2)*cos(phi)
        q[2] = q[0]*sqrt(1-c**2)*sin(phi)
        q[3] = q[0]*c
      
        return q
    
        
    def generate_p(self,n):
        q = self.generate_q(n) 
        Q = np.sum(q,1)
       
        M = self.mass(Q)
        b = -Q[1:]/M
        x = self.Ecms/M
        gamma = Q[0]/M
        a = 1./(1.+gamma)
        bdotq = np.sum(b[:,None]*q[1:],axis=0)
        
        p = np.empty((4,*n)) 
        
        p[0] = x*(q[0]*gamma+bdotq)
        p[1:] = x * (q[1:]+ b[:,None]*q[0] + a*bdotq*b[:,None])
   
        return p
    
    def fxi(self,p,masses,xi):
        return np.sum(sqrt(masses**2+xi**2*p[0]**2))-self.Ecms

    def dfxi(self,p,masses,xi):
        denom = sqrt(masses**2+xi**2*p[0]**2)
        return np.sum(xi * p[0]**2 / denom)
    
    def generate_k(self,masses):
        n = np.shape(masses)
        p = self.generate_p(n) 
        mass_sum = np.sum(masses,axis=0)
        xi_0 = sqrt(1-(mass_sum / 1000)**2)

        xi = np.empty(n[-1])
        for i in range(len(xi)):
            f = partial(self.fxi,p[:,:,i],masses[:,i])
            df = partial(self.dfxi,p[:,:,i],masses[:,i])
            try:
                xi[i] = newton(f,xi_0[i],df)
            except:
                xi[i] =xi_0[i]

        k = np.empty((4,*n))
        k[0] = sqrt(p[0]**2*xi**2+masses**2)
        k[1:] = xi*p[1:]

        return k
    
    def generate_ttbar(self,size=1): 
        
        if self.nout == 2:
            masses = np.full((2,size),self.MT)
            BW_weight = 1
        elif self.nout in [3,4]:
            
            masses = np.full((2,size),self.MT)
            Mmax = self.Ecms
            if self.nout == 3:
                Mmin = self.MW+self.Mb
            else:
                Mmin = self.Mb+self.Me

            
            masses[1], BWw2 = self.generate_mass(self.MT,self.GT,Mmax-self.MT,Mmin,size=size)
            BW_weight = BWw2 
            #BW_weight *= self.pol(masses[1])


        else: 
            Mmax = self.Ecms
            Mmin = self.Mb
            masses = np.empty((2,size))
            masses[0], BWw1 = self.generate_mass(self.MT,self.GT,Mmax,Mmin+self.Mmu,size=size)
            masses[1], BWw2 = self.generate_mass(self.MT,self.GT,Mmax-masses[0],Mmin+self.Me,size=size)
            BW_weight = BWw1*BWw2 
            
            #BW_weight *= self.pol(masses[0])
            #BW_weight *= self.pol(masses[1])




        ttbar = self.generate_k(masses) 
        
        return ttbar , BW_weight



    def lamda(self,P,p1,p2):
            S  = self.mass(P)**2
            s1 = self.mass(p1)**2
            s2 = self.mass(p2)**2
            return (S-s1-s2)**2-4*s1*s2

    def nop_decay(self, pint, masses):
        Ecms = self.mass(pint)
        size = len(Ecms)
        momentum = Ecms/2 * (1-np.sum(array(masses)**2,axis=0)/Ecms**2)
        
        s = masses[0]**2-masses[1]**2
        p = np.empty((4,2,size))
        p[0,0] = (1+s/Ecms**2) * Ecms/2
        p[0,1] = (1-s/Ecms**2) * Ecms/2

        cosT  = 2*np.random.random(size)-1
        theta = np.arccos(cosT)
        phi   = np.random.random(size)
        phi  *= 2*pi
        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)
        momentum = momentum * array([x,y,z])
        p[1:,0] = momentum
        p[1:,1] = -momentum

        for i in range(size):
            #rot = self.rotat(pint[:,i],np.array([1,0,0,1]))
            rot = self.rotat(np.array([1,0,0,1]),pint[:,i])

            p[:,0,i] = self.rotat_inv(p[:,0,i],rot)
            p[:,1,i] = self.rotat_inv(p[:,1,i],rot)
        p[:,0] = self.boost(pint,p[:,0])
        p[:,1] = self.boost(pint,p[:,1])

        Decay_Weight = np.pi * np.sqrt(self.lamda(pint,p[:,0],p[:,1]))/2/self.mass(pint)**2
        return p, Decay_Weight
       
    def sqlamda(self,s,s1,s2):
        return np.sqrt(((s-s1-s2)**2-4*s1*s2))/s

    def Lambda(self,a,b,c):
        L = (a-b-c)**2 -4*b*c
        return sqrt(L)/2/sqrt(a)

    def DecayMassWeight(self,s,sp,b,c):
        return self.Lambda(sp,b,c)/self.Lambda(s,b,c) *s/sp


    def decay(self,pint, masses,m_pin):
        Ecms = self.mass(pint)
        size = len(Ecms)

        s = self.mass(pint)**2  
        s1 = masses[0]**2
        s2 = masses[1]**2 

        p = np.empty((4,2,size))
        p[0,0] = (s+s1-s2)/Ecms/2
        p1m = Ecms * self.sqlamda(s,s1,s2)/2  

        ct  = 2*np.random.random(size)-1
        st  = sqrt(1-ct**2) 
        phi   = 2*pi*np.random.random(size)



        p[1:,0]=p1m * array([st*cos(phi),st*sin(phi),ct]) 
        #for i in range(size): 
        #    #rot = self.rotat(np.array([1,0,0,1]),pint[:,i])
        #    rot = self.rotat(pint[:,i],np.array([1,0,0,1]))
         #   p[:,0,i] = self.rotat_inv(p[:,0,i],rot)

        #ref = np.array([1,0,0,1])
        #p[:,0] = self.Poincare(ref[:,None],pint,p[:,0])
        p[:,0] = self.boost(pint,p[:,0])
        p[:,1] = pint - p[:,0]

        Decay_Weight = pi*self.sqlamda(s,s1,s2)/2
        #Decay_Weight = self.DecayMassWeight(m_pin**2,Ecms**2,masses[0]**2,masses[1]**2) 
   
        norm = np.sum(Decay_Weight)/size
        print("Decay Norm: ",norm)

        return p, Decay_Weight/norm

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
        return res, weight

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

    def AnisotropicDecay(self, pint, masses, ctexp=.8, ctmin=-1, ctmax=1):
        Ecms = self.mass(pint)
        s = Ecms**2
        s1 = masses[0]**2
        s2 = masses[1]**2
        
        ctexp=self.cexp

        size = len(Ecms)
        
        p = np.empty((4,2,size))
        p[0,0] = (s+s1-s2)/Ecms/2.

        p1m = Ecms*self.sqlamda(s,s1,s2)/2
        pim = sqrt(pint[0]**2-s)
        a   = pint[0]*p[0,0]/pim/p1m
        a[np.logical_and(1.>=a,a>=0.)] = 1.0000000001
        ct, PeakedWeight     = self.PeakedDist(a,ctexp,ctmin,ctmax,1,size)
        st = sqrt(1.-ct**2)
        phi = 2.*pi*np.random.rand(size)
        p[1:,0]=p1m * array([st*sin(phi),st*cos(phi),ct])
        
  
        pref = np.empty((4,size))
        pref[0] = pint[0]
        pref[1:3] = 0
        pref[3] = pim
        
        #ref = np.array([1,0,0,1])
        #p[:,0] = self.Poincare(ref[:,None],pint,p[:,0])
        p[:,0] = self.boost(pint,p[:,0])
        #p[:,0] = self.boost(pref,p[:,0])
        #self.Poincare(pref,pint,p[:,0])
        
        p[:,1] = pint - p[:,0]
      
        

        Decay_Weight =  (pi*self.sqlamda(s,s1,s2)/4*(a+ct)**ctexp * PeakedWeight)

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

