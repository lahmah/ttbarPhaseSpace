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


def testBWWeights():
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


wb, _ = gen.decay(ttbar[:,0],"24","5OS")
print(wb)
print(wb.m)

munu,_ = gen.decay(wb[:,0],"13OS","14OS")
print(munu)
print(munu.m)

gen.nout = 2
p,w =gen.generate_point()
print(p)
print(p[0].m, " - " ,p[1].m)
gen.nout = 3
p,w =gen.generate_point()
print(p)
print(p[0].m, " - " ,p[1].m, " - " ,p[2].m)
gen.nout = 4 
p,w =gen.generate_point()
print(p)
print(p[0].m, " - " ,p[1].m , " - ", p[2].m, " - ",p[3].m)
gen.nout = 6 
p,w =gen.generate_point()
print(p)
print(p[0].m, " - " ,p[1].m , " - ", p[2].m, " - ",p[3].m, " - " ,p[4].m , " - ", p[5].m)
#print("\n\n\n\n")
#testBWWeights()

