from matplotlib.pylab import *
import numpy as np

data = open('plot_values_n10.txt')

lu_data = open('ludecomp_values_n10.txt')
lu_x = []
lu_y = []
lu_error = []
temp_values = []
for line in lu_data:
	x = line.split()
	if len(x) >=1:
		temp_values.append(float(x[0]))


lu_x = temp_values[:len(temp_values)/3]
lu_y = temp_values[len(temp_values)/3: len(temp_values)*2/3]
lu_error = temp_values[len(temp_values)*2/3:]

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

del error[0], error[len(error)-1]


errorx = temp_x[1:-1]

x = np.array(temp_x)



plot(x,y,'r')
hold('on')
plot(x,analytic,'b')
hold('on')
plot(lu_x, lu_y, 'g')
legend(["Numeric", "Analytic", "LU"])
plt.figure()
plot(errorx, error, 'r')
hold('on')
plot(lu_x, lu_error)
show()
