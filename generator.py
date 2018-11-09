import numpy as np
from numpy import sqrt, pi, log ,exp, tan, arctan, sin, cos, array
from vec4d import *
from math import gamma

def floorToZero(a,N=8):
    return int(a*10**N)*10**-N

class topGenerator(object):
    def __init__(self, nin, nout, Ecms, debug = False):
        self.nin = nin
        self.nout = nout
        self.Ecms = Ecms
        self.debug = debug

        self.labsystem = Mom4D([0,0,0,1])


        MW = 80.38
        GW = 2.09
        MT = 173
        GT = 1.41
        rWmin = arctan(-(pow(MW,2))/(MW*GW))
        rTmin = arctan(-(pow(MT,2))/(MT*GT))
        
        if nout == 3:
            rTmin =  arctan(MW-(pow(MT,2))/(MT*GT))
       
        self.particles={             
        5  : 0,#4.18,        # b-quark 
        6  : (MT,GT,rTmin),   # t-quark
        11 : 0,              # elektron
        12 : 0,              # e-eutrino
        24 : (MW,GW,rWmin)    # wboson
        }
                       

    def generate_weight(self,pout):
        Ecms = self.Ecms
        nout = self.nout

        ps_volume = (np.pi/2.)**(nout-1) * Ecms**(2*nout-4)/gamma(nout)/gamma(nout-1)

        return ps_volume

    def generate_point(self):
        out = []

        weight = 1
        ttbar,w1 = self.generate_ttbar()
        weight *=w1

        if self.nout ==2:
            out = ttbar
        else:
            wb1,w2 = self.top_decay(ttbar[1])
            weight*=w2
       
        if self.nout == 3:
            out.append(wb1[0])
            out.append(ttbar[0])
            out.append(wb1[1])             
            
            p = [wb1[0],ttbar[0]]
            masses = [wb1[0].m,wb1[1].m]
            Ki_abs = array([sqrt(p[i].E**2-masses[i]**2) for i in range(2)])
            Xi = sum(Ki_abs)/self.Ecms
 
            # weight for mass shift of top quark
            nout = 2
            weight *= Xi**(2*nout-3)    
            temp = 0
            for i in range(nout):
                weight *= Ki_abs[i]/p[i].E
                temp += Ki_abs[i]**2/p[i].E
            weight /= temp
            weight *=self.Ecms

            
        elif self.nout == 4:
            enu = self.w_decay(wb1[0])
            out.append(enu[0])
            out.append(enu[1])
            out.append(ttbar[0])
            out.append(wb1[1]) 


        
        elif self.nout == 6: 
            wb2,w3= self.top_decay(ttbar[0])
            weight*=w3
            enu = self.w_decay(wb1[0])
            munu = self.w_decay(wb2[0])
            
            #out.append(wb2[1])
            #out.append(munu[1])
            #out.append(enu[1])
            #out.append(munu[0])
            #out.append(enu[0])
            #out.append(wb1[1])
            out.append(wb2[1])
            out.append(wb1[1])
            out.append(munu[0])
            out.append(munu[1])
            out.append(enu[0])
            out.append(enu[1])
                                                 

        weight *= self.generate_weight(out)

        if self.nout ==2 or self.nout == 3 or self.nout == 4:
            w3 = 1
            if self.nout == 2:
                w2 = 1

        return out ,weight, w1,w2, w3
    
    def generate_ttbar(self): 
        q = []
        Q = Mom4D()
        for i in range(2):
            q.append(self.generate_fourvector())
            Q += q[i]
        
        M = Q.m
        b = -Q.mom3d/M
        x = self.Ecms/M
        gamma = Q.E/M
        a = 1./(1.+gamma)
        
        p = []
        for i in range(2):
            p.append(Mom4D())
            bdotq = b.dot(q[i].mom3d)
            p[i].E = x*(gamma*q[i].E + bdotq)    
            p[i].mom3d = x*(q[i].mom3d + b*q[i].E + a*bdotq*b)

        if self.nout == 2:
            masses,weight = self.generate_Masses(self.Ecms,"6OS","-6OS")

        if self.nout == 3 or self.nout == 4:
            masses,weight = self.generate_Masses(self.Ecms,"6","-6")
        if self.nout == 6:
            masses,weight = self.generate_Masses(self.Ecms,"6","-6")
            



        Ki_abs = array([sqrt(p[i].E**2-masses[i]**2) for i in range(2)])
 
        Xi = sum(Ki_abs)/self.Ecms
        
        for i in range(2):
            p[i].mom3d *= Xi
            p[i].E = np.sqrt(masses[i]**2+Xi**2*p[i].E**2)

         
        if self.nout == 2:
            # weight for mass shift of top quark
            nout = 2
            weight *= Xi**(2*nout-3)    
            temp = 0
            for i in range(nout):
                weight *= Ki_abs[i]/p[i].E
                temp += Ki_abs[i]**2/p[i].E
            weight /= temp
            weight *=self.Ecms
        
        

        if self.debug:
            print("ttbar:")
            print("masses: ", masses) 
            print("t: ",p[0]," - ",p[0].m)
            print("tb: ",p[1]," - ",p[1].m)
        return p, weight

    def top_decay(self, pint):
        E_CM = pint.m

        if self.nout == 3:
            masses,weight = self.generate_Masses(E_CM,"24OS","5OS") 
        if self.nout == 4 or self.nout == 6:
             masses,weight = self.generate_Masses(E_CM,"24","5OS") 

        if self.debug:
            print("top decay:")
            print("mass in: ",E_CM)
            print("masses out: ", masses) 

        p = [Mom4D(), Mom4D()]

        momentum = E_CM/2 * (1-sum(array(masses)**2)/E_CM**2)
        s = masses[0]**2-masses[1]**2 
        p[0].E = (1+s/E_CM**2) * E_CM/2
        p[1].E = (1-s/E_CM**2) * E_CM/2

        theta, phi = np.random.random(2)
        theta *= pi
        phi   += 2*pi

        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)

        momentum *= array([x,y,z])
        p[0].mom3d = +momentum
        p[1].mom3d = -momentum
        
        if self.debug:
            print("w,b befor boost:")
            print("w: ",p[0], " - " , p[0].m)
            print("b: ",p[1], " - " , p[1].m)



        m_rot = self.rotat(pint,self.labsystem)
        p[0] = self.rotat_inv(p[0],m_rot)
        p[1] = self.rotat_inv(p[1],m_rot)
        p[0] = self.boost(pint,p[0])
        p[1] = self.boost(pint,p[1])


        if self.debug:
            print("w,b after boost:")
            print("w: ",p[0], " - " , p[0].m)
            print("b: ",p[1], " - " , p[1].m)

        return p,weight


    def w_decay(self, pint):
        E_CM = pint.m
        
        p = [Mom4D(),Mom4D()]

        p[0].E = E_CM/2
        p[1].E = E_CM/2


        theta, phi = np.random.random(2)
        theta *= pi
        phi   += 2*pi

        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)
        momentum = E_CM/2 * array([x,y,z])

        p[0].mom3d = momentum
        p[1].mom3d = -momentum

        m_rot = self.rotat(pint,self.labsystem)
        p[0] = self.rotat_inv(p[0],m_rot)
        p[1] = self.rotat_inv(p[1],m_rot)
        p[0] = self.boost(pint,p[0])
        p[1] = self.boost(pint,p[1])


        if self.debug:
            print("W decay")

            print("e / nu: ",p[1], " - " , p[1].m)
            print("nu / nu: ", p[0], " - " , p[0].m)


        return p


    # generation of masses according to Breit-Wigner
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
                    remainingE = m[0]+5
                    rmax = arctan((remainingE**2-m[0]**2)/(m[0]*m[1])) 
                    rmin = arctan(((m[0]-5)**2-m[0]**2)/(m[0]*m[1])) 
                    m = (m[0],m[1],rmin)

                    r = m[2] + np.random.random()*(rmax-m[2])
                    s = m[0]*m[1]*tan(r)+m[0]*m[0]
                    
                    W *= (rmax-m[2])*((s-m[0]**2)**2+(m[0]*m[1])**2)/(m[0]*m[1])  # inverse breit wigner weight

                    #W *=  1./(2.*remainingE*128.*pow(pi,3))*(1.-s/remainingE**2) 

                    W /= m[0]*m[1]*tan(rmax)+m[0]**2 - (m[0]*m[1]*tan(rmin)+m[0]**2)


                    m = sqrt(s)        

                remainingE -= m
            masses.append(m)
        return masses , W

    def generate_fourvector(self):
        ran = np.random.rand(4)
        c = 2.*ran[0]-1.
        phi = 2.*pi*ran[1]
        
        q = Mom4D()
        q[0] = -log(ran[2]*ran[3])
        q[1] = q[0]*sqrt(1-c**2)*cos(phi)
        q[2] = q[0]*sqrt(1-c**2)*sin(phi)
        q[3] = q[0]*c
        
        return q


    
    def polytope(self, m):
    # Produces a uniform random distribution inside a polytope with
    # | x_k | < 1, | x_k-x_l | < 1, see physics/0003078
        
        x = np.empty(m+1)
    
        # number of negative values
        k = int((m+1)*np.random.random())
        x[0] = 0.
    
        if k == 0:
            for i in range(1, m+1):
                x[i] = np.random.random()
        
        elif k == m:
            for i in range(1, m+1):
                x[i] = -np.random.random()
    
        else:
            prod = 1.
            for i in range(1, k+1):
                prod *= np.random.random()
    
            v1 = -log(prod)
            prod = 1.
            for i in range(1, m-k+2):
                prod *= np.random.random()
    
            v2 = -log(prod)
            y1 = v1/(v1+v2)
            x[1] = -y1
            x[m] = (1-y1) * random.random()**(1./(m-k))
            for i in range(2, k+1):
                x[i] = x[1]*random.random()
    
            for i in range(k+1, m):
                x[i] = x[m]*np.random.random()
    
        return np.random.permutation(x)
    
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
    
            for i in range(3):
                for l in range(3):
                    rot[i, l] = 0.
                    for k in range(3):
                        rot[i, l] += r[0, i, k] * r[1, l, k]
    
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
        for i in range(3):
            p1[i+1] = 0.
            for j in range(3):
                p1[i+1] += rot[i, j] * p2[j+1]
    
        return p1

