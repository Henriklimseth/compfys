from matplotlib.pylab import *
import numpy as np

infile1 = open('Eigenvectors_n501.txt')
rho = []
jacobi1 = []
jacobi2 = []
jacobi3 = []
for line in infile1:
	x = line.split()
	rho.append(float(x[0]))
	jacobi1.append(float(x[1])**2)
	jacobi2.append(float(x[2])**2)
	jacobi3.append(float(x[3])**2)

print sum(jacobi1)

infile2 = open('Eigenvectors_armadillo_n501.txt')
arma1 = []
arma2 = []
arma3 = []

for line in infile2:
	x = line.split()
	arma1.append(float(x[0])**2)
	arma2.append(float(x[1])**2)
	arma3.append(float(x[2])**2)

# Analytical solutions:
def psisq0(x):
	return 4./np.sqrt(np.pi)*x**2*np.exp(-x**2)

def psisq1(x):
	return 8./(3*np.sqrt(np.pi))*x**2*(3./2-x**2)**2*np.exp(-x**2)

def psisq2(x):
	return 8./(15*np.sqrt(np.pi))*x**2*(x**4-5*x**2+15./4.)**2*np.exp(-x**2)

rhoarray = array(rho)
#Normalize the eigenfunctions:
int1 = 0
int2 = 0
int3 = 0

h = rho[1]-rho[0]
for i in range(len(jacobi1)):
	int1 += jacobi1[i]*h
	int2 += jacobi2[i]*h
	int3 += jacobi3[i]*h
	


ja1 = array(jacobi1)/int1
ja2 = array(jacobi2)/int2
ja3 = array(jacobi3)/int3
aa1 = array(arma1)/int1
aa2 = array(arma2)/int2
aa3 = array(arma3)/int3


psi1 = psisq0(rhoarray)
psi2 = psisq1(rhoarray)
psi3 = psisq2(rhoarray)

rc('text', usetex='True')

subplot(3,1,1)
plot(rho, ja1, 'r')
hold('on')
plot(rho, aa1, 'b')
hold('on')
plot(rhoarray, psi1, 'g')
legend(['Jacobi', 'Armadillo', 'Analytisk'])
ylabel(r'$|u_0(\rho)|^2$')

subplot(3,1,2)
plot(rho, ja2, 'r')
hold('on')
plot(rho, aa2, 'b')
hold('on')
plot(rhoarray, psi2, 'g')
legend(['Jacobi', 'Armadillo', 'Analytisk'])
ylabel(r'$|u_1(\rho)|^2$')

subplot(3,1,3)
plot(rho, ja3, 'r')
hold('on')
plot(rho, aa3, 'b')
hold('on')
plot(rhoarray, psi3, 'g')
legend(['Jacobi', 'Armadillo', 'Analytisk'])
xlabel(r'$\rho$')
ylabel(r'$|u_2(\rho)|^2$')
savefig('ettelektron.png')
show()





"""
rho_test = [1., 2., 3., 4., 4.5, 5., 6., 7., 8.]
It_rho_n100 = [17114., 16755., 16588., 16405., 16270., 16217., 16021., 15835., 15696.]

n_test = [100., 200., 300., 400., 500]
it_n_rho5 = [16475., 66828., 150798., 269021., 422320.]
t_n_rho5 = [0.19, 2.94, 14.4, 44.66, 108.05]

plt.rc('text', usetex='True')

plt.subplot(2,1,1)
plt.plot(rho_test, It_rho_n100)
plt.xlabel(r'$\rho_{max}$')
plt.ylabel('Antall transformasjoner')
plt.subplot(2,1,2)
plt.plot(n_test, it_n_rho5)
plt.ylabel('Antall transformasjoner')
plt.xlabel('n')
plt.savefig('transformasjoner.png')

plt.subplot(2,1,1)
plt.plot(n_test, t_n_rho5)
plt.xlabel(r'$n$')
plt.ylabel(r'$t$')
plt.subplot(2,1,2)
plt.plot(t_n_rho5,it_n_rho5) 
plt.xlabel(r'$t$')
plt.ylabel(r'Antall transformasjoner')
plt.savefig('tidsomfunkavn.png')
plt.show()
"""
