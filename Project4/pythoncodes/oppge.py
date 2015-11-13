import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import *

infile80 = open("total80x80.txt")

infile60 = open("total60x60.txt")

infile40 = open("total40x40.txt")

infile20 = open("total20x20.txt")

T=[]
E80=[]; M80=[]; C_V80=[]; X80=[]
E60=[]; M60=[]; C_V60=[]; X60=[]
E40=[]; M40=[]; C_V40=[]; X40=[]
E20=[]; M20=[]; C_V20=[]; X20=[]

for line in infile80:
	T.append(float(line.split()[0]))
	E80.append(float(line.split()[1]))
	M80.append(float(line.split()[2]))
	C_V80.append(float(line.split()[3]))
	X80.append(float(line.split()[4]))


for line in infile60:
	E60.append(float(line.split()[1]))
	M60.append(float(line.split()[2]))
	C_V60.append(float(line.split()[3]))
	X60.append(float(line.split()[4]))


for line in infile40:
	E40.append(float(line.split()[1]))
	M40.append(float(line.split()[2]))
	C_V40.append(float(line.split()[3]))
	X40.append(float(line.split()[4]))

for line in infile20:
	E20.append(float(line.split()[1]))
	M20.append(float(line.split()[2]))
	C_V20.append(float(line.split()[3]))
	X20.append(float(line.split()[4]))



for i in range(len(T)-1):
	if X80[i+1]>X80[i]:
		n80 = i+1
	if X60[i+1]>X60[i]:
		n60 = i+1
	if X40[i+1]>X40[i]:
		n40 = i+1
	if X20[i+1]>X20[i]:
		n20 = i+1

func = lambda p,x: p[0] + p[1]*x
L_inv = np.array([1./80,1./60, 1./40, 1./20])
T_C = np.array([T[n80],T[n60], T[n40], T[n20]])
x0 = np.array([0.0,0.0])
Errorfunc = lambda p,x,y: func(p,x)-y
T_C_coeffs = leastsq(Errorfunc, x0, args=(L_inv, T_C))[0]
T_C_reg = lambda coeffs, x: coeffs[0] + coeffs[1]*x
L_inv_array = np.linspace(0,0.5,100)
plt.plot(L_inv_array, T_C_reg(T_C_coeffs, L_inv_array))
plt.xlabel(r"$1/L$")
plt.ylabel(r"$T_C(1/L)$")
plt.savefig("criticaltemp.png")
plt.show()


"""
plt.plot(T, M20)
plt.hold('on')
plt.plot(T, M40)
plt.plot(T, M60)
plt.plot(T, M80)
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle |M|\rangle$")
plt.legend([r"$20\times20$",r"$40\times40$",r"$60\times60$",r"$80\times80$"], loc="Upper left")
plt.savefig("magnetization.png")
plt.show()

"""



