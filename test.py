from newGen import *
gen = topGenerator(6,1000)
A = gen.generate_fourvectors(2)
print(A)
print(A.m)
print(np.sum(A.E))
B = gen.generate_p(2)
print(B)
print(B.m)
print(np.sum(B.E))
B = gen.generate_p(2)
masses = np.array([175,175])
C = gen.generate_k(masses)
print(C)
print(C.m)
print(np.sum(C.E))
ttbar, BW_weight = gen.generate_ttbar()
print(ttbar)
print(ttbar.m)
print(sum(ttbar.E))
print(BW_weight)

W = 0
for i in range(1,10000000):
    ttbar ,BW_weight = gen.generate_ttbar()
    W += BW_weight
    print(W/i, " - " ,i,end="\r")

print()
