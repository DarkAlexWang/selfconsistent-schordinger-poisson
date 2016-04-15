 #import corresponding module
from scipy import constants as pc
import numpy as np
import matplotlib.pyplot as plt 


x = np.arange(0,1,0.05)
y = np.power(x, 2)

fig = plt.figure()
ax = fig.gca()
ax.set_xticks(np.arange(0,1,0.1))
ax.set_yticks(np.arange(0,1.,0.1))
plt.scatter(x,y)
plt.grid()
plt.show()
