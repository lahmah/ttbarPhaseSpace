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

tt = 0
ww = 0
ttww = 0
tt2 = 0
wW2 = 0
ttww2 = 0
precision = "{:.4f}"
print("\n\n\n")
for i in range(1,10000000):
    ttbar ,ttw = gen.generate_ttbar()
    W, ww1 = gen.generate_Masses(175,"24")
    W, ww2 = gen.generate_Masses(175,"24")

    tt += ttw
    ww += ww1*ww2
    ttww += ttw*ww1*ww2
    tt2 += ttw**2
    wW2 += (ww1*ww2)**2
    ttww2 += (ttw*ww1*ww2)**2

    
    if i % 1000 == 0:
        print(end="\033[F\033[F\033[F\033[F")
        print(i)
        print("tt: ",precision.format(tt/i), "\t-\t", precision.format(np.sqrt((tt2/i-(tt/i)**2)/i)))
        print("ww: ",precision.format(ww/i), "\t-\t", precision.format(np.sqrt((wW2/i-(ww/i)**2)/i)))
        print("tw: ",precision.format(ttww/i), "\t-\t", precision.format(np.sqrt((ttww2/i-(ttww/i)**2)/i)))


print("\n\n\n\n")
