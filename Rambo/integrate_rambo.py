#!/usr/bin/env python3
import Sherpa
from rambo import *
from vec4d import *
import numpy as np
from itertools import combinations

### run parameters ############################
seed = 1234
N = 100000 # number of events
E_CM = 1000.
nin = 2
nout = 3
pT_min = 5.


FixMass = True 
TopMass = 173
WMass   = 80
###############################################

if FixMass:
    if nout == 2:
        masses = [TopMass,TopMass]
    if nout == 3:
        masses = [TopMass,WMass,0]
    if nout == 4:
        masses = [TopMass,0,0,0]
    if nout == 6:
        masses = [0,0,0,0,0,0]
else: # 1 = massiv particle ; 0 = massles particle
    if nout == 2:
        masses = [1,1]
    if nout == 3:
        masses = [1,1,0]
    if nout == 4:
        masses = [1,0,0,0]
    if nout == 6:
        masses = [0,0,0,0]

################################################

np.random.seed(seed)

conversion = 0.389379*1e9 # convert to pb

## Initialize Sherpa
Generator=Sherpa.Sherpa()
Generator.InitializeTheRun(4, [''.encode('ascii'), ('RUNDATA=../runcards/ttbar'+str(nout)+'.dat').encode('ascii'), 'INIT_ONLY=2'.encode('ascii'), 'OUTPUT=0'.encode('ascii')])
Process=Sherpa.MEProcess(Generator)

# Incoming flavors must be added first!
Process.Initialize();

# First argument corresponds to particle index:
# index 0 corresponds to particle added first, index 1 is the particle added second, and so on...
Process.SetMomentum(0, E_CM/2, 0., 0., E_CM/2)
Process.SetMomentum(1, E_CM/2, 0., 0., -E_CM/2)

## Initialize Rambo
PSGenerator = Rambo(nout, E_CM)


sum = 0
sum2 = 0.
N_gen = 0
N_nan = 0
N_acc = 0
last_print = 0
weights = []
while N_acc<N:
    N_gen +=1

    

    momenta,weight = PSGenerator.generate_massive_point(masses,FixMass)

    

    for i, momentum in enumerate(momenta):
        Process.SetMomentum(i+2, momentum[0], momentum[1], momentum[2], momentum[3])

    
    #Process.SetMomentum(4, momenta[0][0], momenta[0][1], momenta[0][2], momenta[0][3])
    #Process.SetMomentum(2, momenta[1][0], momenta[1][1], momenta[1][2], momenta[0][3])
    #Process.SetMomentum(3, momenta[2][0], momenta[2][1], momenta[2][2], momenta[0][3])
    
    me = Process.CSMatrixElement()


    cosPhi = [angle(p1, p2) for p1, p2 in combinations(momenta, 2)][0]
    p = [p for p in momenta]
    p1 = [p[0].E,p[0].m, p[0].mom3d]
    p2 = [p[0].E,p[0].m, p[0].mom3d]
    if np.isnan(me): 
        N_gen -= 1
        N_nan += 1
        
        continue
    N_acc += 1
    
    sum += me*weight
  
    sum2 += (me*weight)**2

    weights.append(weight)

    if N_acc%1000 == 0 and N_acc > last_print:
        print('Event', N_acc,end="\r")
        last_print = N_acc


x0 = sum / N_acc
x00 = sum2 / N_acc
sigma2 = x00-x0**2
error = np.sqrt(sigma2/N_acc)

xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) *x0 



sigma_xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) * error 

print('xs: ', xs, "pb  +/- " , sigma_xs , "pb = " , sigma_xs/xs*100,' %')


print("generated events:", N_gen)
print("accepted events", N_acc)
print('acceptance rate: ', N_acc/N_gen)
print('nan detection: ', N_nan)


