from os import system
from  multiprocessing import Process
from numpy import log

command = "python newInt.py nout 4 N 2000000 npr 1000000"
filename = " file is4polSFTPT-"

processes = []
for i in range(5):
    seed = " seed"+str(log((33*i)**3))
    log = " > log"+str(i)
    print(seed)
    t = Process(target=system, args=(command+filename+str(i)+seed+log,))
    processes.append(t)
    t.start()

for p in processes:
    p.join()

print("Done")
