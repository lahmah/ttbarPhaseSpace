#!/usr/bin/env python3
import Sherpa
from newGen import *
from vec4d import *
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from export_hepmc import export_hepmc


import sys
if "debug" in sys.argv:
    DEBUG = True
    BREAK = True
else:
    DEBUG = False

if "nobreak" in sys.argv:
    BREAK = False 

if "noweights" in sys.argv:
    NOWEIGHTS = True
else:
    NOWEIGHTS = False

if "plot" in sys.argv:
    PLOT =True
else:
    PLOT = False

### run parameters ############################
seed = 1234
N = 100000 # number of events
E_CM = 1000.
nin = 2
nout = 3 
###############################################

np.random.seed(seed)

conversion = 0.389379*1e9 # convert to pb

## Initialize Sherpa
Generator=Sherpa.Sherpa()

Generator.InitializeTheRun(4, [''.encode('ascii'), ('RUNDATA=runcards/ttbar'+str(nout)+'.dat').encode('ascii'), 'INIT_ONLY=2'.encode('ascii'), 'OUTPUT=0'.encode('ascii')])
Process=Sherpa.MEProcess(Generator)

# Incoming flavors must be added first!
Process.Initialize();


# First argument corresponds to particle index:
# index 0 corresponds to particle added first, index 1 is the particle added second, and so on...
Process.SetMomentum(0, E_CM/2, 0., 0., E_CM/2)
Process.SetMomentum(1, E_CM/2, 0., 0., -E_CM/2)

## Initialize

PSGenerator = topGenerator(nout,E_CM,debug=DEBUG)
pin1 = Mom4D([E_CM/2, 0., 0., E_CM/2])
pin2 = Mom4D([E_CM/2, 0., 0., -E_CM/2])
pin = [pin1, pin2]


sum = 0.
sum2 = 0.
N_gen = 0
N_nan = 0
N_acc = 0
last_print = 0

WAll = 0
W1 = 0
W2 = 0
W3 = 0
W12 = 0

ALLP = [] 
ALLW = []
while N_acc<N:
    N_gen +=1

    if N_acc == -1:
        PSGenerator.debug=True
        DEBUG= True
        BREAK = True
    momenta,weight,w1,w2,w3 = PSGenerator.generate_point() 
    WAll += weight
    W1 += w1
    W2 += w2
    W12 += w1*w2
    W3 += w3

    if DEBUG:
        print("sum of all Momenta:")
        print(np.sum(momenta))
        print("------")
        if BREAK:
            input()

    temp = []
    for i, momentum in enumerate(momenta):
        if DEBUG:
            print(momentum, " - ", momentum.m)
        Process.SetMomentum(i+2, momentum[0], momentum[1], momentum[2], momentum[3])
        temp.append(momentum._arr)
    ALLP.append(temp)

        
    me = Process.CSMatrixElement()
    

    if np.isnan(me):
        N_gen -= 1
        N_nan += 1
        #print(N_nan,end="\r")
        continue


    N_acc += 1
    sum += weight*me
    sum2 += (weight*me)**2
    ALLW.append(weight*me)
        
    if N_acc%1000 == 0 and N_acc > last_print:
        print('Event', N_acc,end ="\r")
        last_print = N_acc


x0 = sum / N_acc
x00 = sum2 / N_acc
sigma2 = x00-x0**2
error = np.sqrt(sigma2/N_acc)


xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) *x0
sigma_xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) *error 

    
print('xs: ', xs, "pb  +/- " , sigma_xs , "pb = " , sigma_xs/xs*100,' %')
print("generated events:", N_gen)
print("accepted events", N_acc)
print('acceptance rate: ', N_acc/N_gen)
print('generated nan: ', N_nan)

print(W1/N_acc)
print(W2/N_acc)
print(W12/N_acc)
print(W3/N_acc)
print(WAll/N_acc)
export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), ALLW, "./Rambo/newGen.hepmc")

