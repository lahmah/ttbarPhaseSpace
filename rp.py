from os import system
from  multiprocessing import Process
from numpy import log

command = "python newInt.py nout 4 N 10000000 npr 1000000"
filename = " file is4polSFT-"

processes = []
for i in range(5,15):
    seed = int(log((33*(i+1))**6)**3)
    seed = " seed "+str(seed)
    logf = " > log"+str(i)
    print(seed)
    t = Process(target=system, args=(command+filename+str(i)+seed+logf,))
    processes.append(t)
    t.start()

for p in processes:
    p.join()

print("Done")
