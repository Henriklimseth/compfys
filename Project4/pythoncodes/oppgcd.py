import matplotlib.pylab as plt
import numpy as np




infile = open("20x20energies_highT_2.txt")
E = []

for line in infile:
	E.append(float(line.split()[0]))


plt.hist(E, bins=121, normed=True)
plt.xlabel(r"$E$")
plt.ylabel(r"$P(E)$")
plt.title(r"$T = 2.4$")

plt.savefig("P_E_highT_2.png")
plt.show()


"""
infile1 = open("20x20run_ordered_lowT.txt")
infile2 = open("20x20run_disordered_lowT.txt")

infile3 = open("20x20run_ordered_highT.txt")
infile4 = open("20x20run_disordered_highT.txt")



E1 = []; M1 = []; C1 = []

E2 = []; M2 = []; C2 = []

E3 = []; M3 = []; C3 = []

E4 = []; M4 = []; C4 = []



for line in infile1:
	E1.append(float(line.split()[0]))
	M1.append(float(line.split()[1]))
	C1.append(float(line.split()[2]))

for line in infile2:
	E2.append(float(line.split()[0]))
	M2.append(float(line.split()[1]))
	C2.append(float(line.split()[2]))


for line in infile3:
	E3.append(float(line.split()[0]))
	M3.append(float(line.split()[1]))
	C3.append(float(line.split()[2]))

for line in infile4:
	E4.append(float(line.split()[0]))
	M4.append(float(line.split()[1]))
	C4.append(float(line.split()[2]))

plt.plot(M1[0:2000])
plt.hold('on')
plt.plot(M2[0:2000])
plt.axis([0,2000,0,1.1])
plt.xlabel("Number of MC cycles")
plt.ylabel(r"$\langle|M|\rangle/L^2$")
plt.title(r"$T=1.0$")
plt.legend(["Ordered", "Disordered"], loc="lower right")
plt.savefig("meanmagT10.png")
plt.show()


plt.subplot(2,1,1)
plt.plot(C1[0:2000])
plt.hold('on')
plt.plot(C3[0:2000])
plt.title("Ordered")
#plt.xlabel("Number of MC cycles")
plt.ylabel(r"$\frac{Accepted\ configurations}{MC cycles}$")
plt.legend([r"$T=1.0$", r"$T=2.4$"], loc="lower right")

plt.subplot(2,1,2)
plt.plot(C2[0:2000])
plt.hold('on')
plt.plot(C4[0:2000])
plt.title("Disordered")
plt.xlabel("Number of MC cycles")
plt.ylabel(r"$\frac{Accepted\ configurations}{MC cycles}$")
plt.legend([r"$T=1.0$", r"$T=2.4$"], loc="lower right")

plt.savefig("configurations.png")
#plt.axis([0,2000,0,1.1])
plt.show()
"""
