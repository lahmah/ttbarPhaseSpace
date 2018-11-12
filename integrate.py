#!/usr/bin/env python3
import Sherpa
from generator import *
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
pT_min = 5
angle_min = np.cos(0.3)
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

PSGenerator = topGenerator(nin, nout,E_CM,debug=DEBUG)
pin1 = Mom4D([E_CM/2, 0., 0., E_CM/2])
pin2 = Mom4D([E_CM/2, 0., 0., -E_CM/2])
pin = [pin1, pin2]
s0 = 2*pT_min**2*min(1-np.cos(angle_min), 1./(1+np.sqrt(1-pT_min**2/E_CM**2))) # cut on inv. mass

sum = 0.
sum2 = 0.
N_gen = 0
N_nan = 0
N_acc = 0
last_print = 0

EW = []
ET = []
EB = []
EB2 = []

EE = []
ENU = []
EMU = []
MUNU = []

MT = []

W = []
W1 = []
W2 = []
W3 = []


ALLP = []

while N_acc<N:
    N_gen +=1

    if N_acc == -1:
        PSGenerator.debug=True
        DEBUG= True
        BREAK = True
    momenta,weight, w1,w2,w3 = PSGenerator.generate_point() 
    W1.append(w1)
    W2.append(w2)
    W3.append(w3)



    if DEBUG:
        print("sum of all Momenta:")
        print(np.sum(momenta))
        print("------")
        if BREAK:
            input()


    temp = []
    for i, momentum in enumerate(momenta):
        temp.append([momentum._arr])
        Process.SetMomentum(i+2, momentum[0], momentum[1], momentum[2], momentum[3])
    ALLP.append(temp)
        
    me = Process.CSMatrixElement()
    

    if np.isnan(me):
        N_gen -= 1
        N_nan += 1
        #print(N_nan,end="\r")
        continue


    if nout == 3:
        EW.append(momenta[0].E)

    if nout == 4:
        EE.append(momenta[0].E)
        ENU.append(momenta[1].E)

    if nout == 3 or nout == 4:
        ET.append(momenta[-2].E)
        EB.append(momenta[-1].E)
        MT.append(momenta[-2].m)

    if nout == 6:
       EE.append(momenta[0].E) 
       EMU.append(momenta[1].E)
       MUNU.append(momenta[2].E) 
       ENU.append(momenta[3].E)
       EB2.append(momenta[4].E)
       EB.append(momenta[5].E)




    
    if NOWEIGHTS:
        W.append(1)
    else:
        W.append(weight*me)




    N_acc += 1
    sum += weight*me
    sum2 += (weight*me)**2

        
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

W1 = np.array(W1)
W2 = np.array(W2)
W3 = np.array(W3)

print(np.sum(W1)/N)
print(np.sum(W2)/N)
print(np.sum(W3)/N)

print(np.sum(W1*W2)/N)
print(np.sum(W1*W3)/N)
print(np.sum(W2*W3)/N)
print(np.sum(W1*W2*W3)/N)

print(np.sum(W)/N)




export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), W, "./oldGen.hepmc")
export_hepmc(E_CM, np.array(ALLP).reshape(N,nout*4), np.zeros(len(ALLP)), "./oldGenNoWeights.hepmc")




if PLOT:
    gs = gridspec.GridSpec(2, 4)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])
    ax5 = plt.subplot(gs[4])
    ax6 = plt.subplot(gs[5])
    ax7 = plt.subplot(gs[6])
    ax8 = plt.subplot(gs[7])

    ax1.hist(EB, label ="E b",density=True,weights=W)
    ax1.set_title("E b")
    if nout == 6:
        ax2.hist(EB2,density=True,weights=W)
        ax2.set_title("E b2")
    else:
        ax2.hist(ET, label ="E T",density=True,weights=W)
        ax2.set_title("E T")
        ax3.hist(MT, label="M T",density=True,weights=W)
        ax3.set_title("M T")
    if nout == 3:
        ax4.hist(EW, label="E W",density=True,weights=W)
        ax4.set_title("E W")
    if nout == 4 or nout == 6:
        ax5.hist(EE, label ="E E",density=True,weights=W)
        ax5.set_title("E E")
        ax6.hist(ENU, label="E Nu",density=True,weights=W)
        ax6.set_title("E Nu")
    if  nout == 6:
        ax7.hist(EMU, label ="E Mu",density=True,weights=W)
        ax7.set_title("E Mu")
        ax8.hist(MUNU, label="E MuNu",density=True,weights=W)
        ax8.set_title("E MuNu")




    plt.legend(loc="best")
    plt.show()

    gs = gridspec.GridSpec(1, 3)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.hist(W1)
    ax1.set_title("weight ttbar")
    ax2.hist(W2,range=[0,10000])
    ax2.set_title("w1")
    ax3.hist(W3,range=[0,10000])
    ax3.set_title("w2")
    plt.show()
    



