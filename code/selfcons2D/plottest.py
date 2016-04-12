 #import corresponding module
from dolfin import *
from meshCreator import *
from schrodinger import schrodingerEq
from electronDensity import electronOccupationState, electronDensityFunction
from scipy import constants as pc
import numpy as np
import matplotlib.pyplot as plt 
from poisson import poissonEq
import time
import sys




x = numpy.arange(0,1,0.05)
y = numpy.power(x, 2)

fig = plt.figure()
ax = fig.gca()
ax.set_xticks(numpy.arange(0,1,0.1))
ax.set_yticks(numpy.arange(0,1.,0.1))
plt.scatter(x,y)
plt.grid()
plt.show()
