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
N = 50000 # number of events
E_CM = 1000.
nin = 2
nout = 4 
filename="own"
n_per_run = 50000
###############################################
if "nout" in sys.argv:
    index = sys.argv.index("nout")
    nout = int(sys.argv[index+1])
if "file" in sys.argv:
    index = sys.argv.index("file")
    filename = sys.argv[index+1]
if "N" in sys.argv:
    index = sys.argv.index("N")
    N = int(sys.argv[index+1])

if "npr" in sys.argv:
    index = sys.argv.index("npr")
    n_per_run= int(sys.argv[index+1])

if "c" in sys.argv:
    index = sys.argv.index("c")
    cexp = float(sys.argv[index+1])
else:
    cexp = -1

n_per_run = min(N,n_per_run)


mask = np.arange(nout)

if "mask" in sys.argv:
    index = sys.argv.index("mask")
    for i in range(nout):
        mask[i] = int(sys.argv[index+1+i])
        

if nout == 2:
    pids = [6,-6]
if nout == 3:
    pids = [-24,6,-5]
if nout == 4:
    pids = np.array([-12,6,11,-5])
    pids = pids[mask]
    print(pids)
if nout == 6:
    pids = np.array([14,5,-12,-13,11,-5])

if "seed" in sys.argv:
    index = sys.argv.index("seed")
    seed = int(sys.argv[index+1])

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

PSGenerator = topGenerator(nout,E_CM,cexp,debug=DEBUG)
pin1 = Mom4D([E_CM/2, 0., 0., E_CM/2])
pin2 = Mom4D([E_CM/2, 0., 0., -E_CM/2])
pin = [pin1, pin2]


sum = 0.
sum2 = 0.
sumW = 0
sumW2 = 0
N_gen = 0
N_nan = 0
N_acc = 0
N_run = -1
last_print = 0

ALLP = [] 
ALLW = []
ALLWNoM = []

while N_acc<N:
    if N-N_acc < n_per_run:
        n_per_run=N-N_acc

    N_gen +=n_per_run
    N_run += 1

    momenta,weight = PSGenerator.generate_point(n_per_run) 
    ALLW = np.append(ALLW,weight)
    ALLWNoM = np.append(ALLW,weight)




    for i in range(n_per_run):
        temp = []
        momentum = momenta[:,:,i]
        for j in range(nout):
            Process.SetMomentum(j+2, momentum[0,j], momentum[1,j], momentum[2,j], momentum[3,j])
            temp.append(momentum[:,j])

                    
        me = Process.CSMatrixElement()
        #    print(temp)
        #    print(weight)
        #    print(me)
        #    input()


        if np.isnan(me):
            N_gen -= 1
            N_nan += 1
            #print(N_nan,end="\r")
            continue



        N_acc += 1
        sum += weight[i]*me
        sum2 += (weight[i]*me)**2
        sumW += weight[i]
        sumW2 += weight[i]**2
        ALLW[N_run*n_per_run+i]*=me
        ALLP.append(temp)

        
    if N_acc%100 == 0 and N_acc > last_print:
        print('Event', N_acc,end ="\r")
        last_print = N_acc


x0 = sum / N_acc
x00 = sum2 / N_acc
sigma2 = x00-x0**2
error = np.sqrt(sigma2/N_acc)


xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) *x0
sigma_xs = conversion * (2*np.pi)**(4-3.*nout)/(2*E_CM**2) *error 

ps = sumW/ N_acc
ps2 = sumW2 / N_acc
sigma2W = ps2-ps**2
errorW = np.sqrt(sigma2W/N_acc)
    
print('xs: %e' % xs, "pb  +/- " , sigma_xs , "pb = " , sigma_xs/xs*100,' %')
print('ps: ', ps, " +/- ", errorW)
print("generated events:", N_gen)
print("accepted events", N_acc)
print('acceptance rate: ', N_acc/N_gen)
print('generated nan: ', N_nan)

export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), ALLW,pids, "./"+filename+".hepmc")
export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), np.ones(len(ALLP)),pids, "./"+filename+"NW.hepmc")
export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), ALLWNoM,pids, "./"+filename+"NM.hepmc")

#np.save("wbtr",wbtr)

