import sys
from numpy import pi, sqrt


nout = float(sys.argv[1])
E_CM = float(sys.argv[2])

filenames = sys.argv[3:]

conversion = 0.389379*1e9 # convert to pb


sum = 0 
sum2 = 0
N_gen = 0

for filename in filenames:
    with open(filename, "r") as file:
        for line in file:
            if line[0] == "E":
                split = line.split()
                N_gen += 1
                sum += float(split[-1])
                sum2 +=  float(split[-1])**2

                x0 = sum / N_gen
                x00 = sum2 / N_gen
                sigma2 = x00-x0**2
                error = sqrt(sigma2/N_gen)

                xs = conversion * (2*pi)**(4-3.*nout)/(2*E_CM**2) *x0 
                sigma_xs = conversion * (2*pi)**(4-3.*nout)/(2*E_CM**2) *error 

print(N_gen,": ", xs, "\t", sigma_xs)
