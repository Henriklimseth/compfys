from matplotlib.pylab import *
import numpy as np

data = open('plot_values_n10.txt')

temp_x = []
y = []
error=[]
analytic = []


for line in data:
	x = line.split()
	temp_x.append(float(x[0]))
	y.append(float(x[1]))
	analytic.append(float(x[2]))
	error.append(float(x[3]))

del error[0], error[len(error)-1], error[len(error)-2]


errorx = temp_x[1:-2]

x = np.array(temp_x)



plot(x,y,'r')
hold('on')
plot(x,analytic,'b')
legend(["Numeric", "Analytic"])
plt.figure()
plot(errorx, error)

show()
