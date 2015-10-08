from matplotlib.pylab import *
import numpy as np

infile1 = open('Eigenvectors_nonint_omega_r0.010000.txt')
infile2 = open('Eigenvectors_nonint_omega_r0.500000.txt')
infile3 = open('Eigenvectors_nonint_omega_r1.000000.txt')
infile4 = open('Eigenvectors_nonint_omega_r5.000000.txt')

rho_list1 = []
arm_list1 = []
for line in infile1:
	x = line.split()
	rho_list1.append(float(x[0]))
	arm_list1.append(float(x[1])**2)

rho_list2 = []
arm_list2 = []
for line in infile2:
	x = line.split()
	rho_list2.append(float(x[0]))
	arm_list2.append(float(x[1])**2)

rho_list3 = []
arm_list3 = []
for line in infile3:
	x = line.split()
	rho_list3.append(float(x[0]))
	arm_list3.append(float(x[1])**2)

rho_list4 = []
arm_list4 = []
for line in infile4:
	x = line.split()
	rho_list4.append(float(x[0]))
	arm_list4.append(float(x[1])**2)



arm1 = np.array(arm_list1)
rho1 = np.array(rho_list1)

arm2 = np.array(arm_list2)
rho2 = np.array(rho_list2)


arm3 = np.array(arm_list3)
rho3 = np.array(rho_list3)

arm4 = np.array(arm_list4)
rho4 = np.array(rho_list4)

"""
def u_sq1(x):
	return x**2*exp(-x**2/4.)*(1+x/2.)**2

def u_sq2(x):
	return x**2*exp(-x**2/20.)*(1+x/2.+x**2/20)**2
"""
integral_arm1 = 0.
integral_arm2 = 0.
integral_arm3 = 0.
integral_arm4 = 0.

h1 = rho1[1]-rho1[0]
h2 = rho2[1]-rho2[0]
h3 = rho3[1]-rho3[0]
h4 = rho4[1]-rho4[0]



for i in range(len(arm1)):
	integral_arm1 += h1*arm1[i]
	integral_arm2 += h2*arm2[i]
	integral_arm3 += h3*arm3[i]
	integral_arm4 += h4*arm4[i]




arm_norm1 = arm1/integral_arm1
arm_norm2 = arm2/integral_arm2
arm_norm3 = arm3/integral_arm3
arm_norm4 = arm4/integral_arm4



rc('text', usetex='True')

subplot(2,1,1)
plot(rho1, arm_norm1, 'b')
legend([r'$\omega_r = 0.01$'])
ylabel(r'$|u(\rho)|^2$')
xlabel(r'$\rho$')

subplot(2,1,2)
plot(rho2, arm_norm2, 'b')
legend([r'$\omega_r = 0.5$'])
xlabel(r'$\rho$')
ylabel(r'$|u(\rho)|^2$')

savefig('ikkevvlitenomega.png')
show()
