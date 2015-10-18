from matplotlib.pylab import *
import numpy as np


def f(x):
	return np.exp(-4*(x))

x = np.linspace(0,2,100)

plot(x, f(x))
hold('on')
plot(2.3, f(2.3), 'o')
xlabel('x')
ylabel('exp(-4x)')
savefig('integrandplot.png')
show()
