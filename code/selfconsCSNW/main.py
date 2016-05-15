'''
2D Poisson-Schrodinger equation self-consistant solver for a core-shell 
hexagonal nanowire, the mesh is produced by mshr.
Created by Zhihuan Wang
Date 04/14/2016  Modified 04/14/2016

'''

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

#start timing
start = time.clock()

#constants
#band_offset = 0.24 #eV
evA2tocm3 = 7.53858e30
q = 1.602e-19

#system input fermi energy with python main.py 0.0
fermiEnergy = float(sys.argv[1])

meshArray = meshFunc()
V0 = meshArray[0]
mesh = meshArray[5]
band_offset = meshArray[6]
V = meshArray[7]
band_offset = interpolate(band_offset, V)
# iPython debugger
#import ipdb
#ipdb.set_trace()
#ipdb.pm()
#ipdb.run('x{0} = 3')
#result =  ipdb.runcall(function, arg0, arg1, kwarg='foo')
#result = ipdb.runeval('f(1,2) - 3')

potential = interpolate(Expression("1"), V)
n = interpolate(Expression("1"), V)

#Start of Iteration
k = 0
while k < 10:
  print str(k+1) + ' iteration'
  # get the eigenvalues and eigenfunctions
  eigenvalues, eigenvectors, u, rx_list = schrodingerEq(potential, meshArray, fermiEnergy)
  # find the electron density 
  nk = electronOccupationState(eigenvalues, fermiEnergy)
  new_n = electronDensityFunction(eigenvectors, nk)
  # solve the poisson equation
  phi = poissonEq(new_n,meshArray)



  #get the new potential V(x) = -q phi(x) + band_offset
  new_potential = Function(V)
  new_potential.vector()[:] = np.array([-q*j for j in phi.vector()[:]])
 
  print len(new_potential.vector()) 
  #for i in range(len(new_potential.vector())):
   # import pdb
   # pdb.set_trace()
   # new_potential.vector()[i] = new_potential.vector()[i] + band_offset.vector()[i]
  #print len(new_potential.vector()[i])

  print '\n'

  #convergence 
  # this loop will run until we reach a convergent solution, or we reach our
  # maximum number of iterations, given by
  print 'convergence'
  v_error = []
  n_error = []

  for i in range(len(potential.vector())):
    v_error.append((new_potential.vector()[i]-potential.vector()[i]))
    n_error.append((n.vector()[i]-new_n[i]))
  
  v_error = abs(sum(v_error)/len(v_error))
  n_error = abs(sum(n_error)/len(n_error))
  n_average = sum(new_n)/len(new_n)
  e_average = sum(eigenvalues)/len(eigenvalues)

  print 'v_error: ' + str(v_error)
  print 'n_error: ' + str(n_error) + '\n'
  print 'fermiEnergy (4 K) = '  + str(fermiEnergy)
  print 'average electron density = ' + str(n_average)
  print 'average eigenvalue = ' + str(e_average)


  #update potential, electron density and iteration
  potential = new_potential
  n.vector()[:] = np.array([j for j in new_n])
  k = k+1
  
  if(n_error < 10e-5 and v_error < 10e-5):
    print 'Convergence occured at k: ' + str(k)
    break
  
#***End of Iteration***

#end timing
end = time.clock()

#print timing information
print 'time elapsed: ' + str(end - start) + ' seconds'

#with open("fermi.txt", "a") as myfile:
#    myfile.write(str(fermiEnergy) + ', ' + str(n_average) + ', ' + str(end - start) + ', ' +  str(n_error) + ', ' + str(v_error) + ', ' + str(n_error < 10e-5 and v_error < 10e-5) + ', '  + str(k) + '\n')
#    myfile.close()

#plotting
v1 = Function(V)

e1 = Function(V)
e2 = Function(V)
e3 = Function(V)

#electron density in units of 1/cm^3
n.vector()[:] = np.array([evA2tocm3 * j for j in new_n])
v1.vector()[:] = n.vector().array()

e1.vector()[:] = np.array(eigenvectors[0])
e2.vector()[:] = np.array(eigenvectors[1])
e3.vector()[:] = np.array(eigenvectors[2])

Iv1 = interpolate(v1,V)
Ie1 = interpolate(e1,V)
Ie2 = interpolate(e2,V)
Ie3 = interpolate(e3,V)

#save the output
output = open("out.txt", "w+")

coor = mesh.coordinates()
u_array = Iv1.vector().array()
if mesh.num_vertices() == len(u_array)+1:
  for i in range(mesh.num_vertices()-1):
    output.write('n(%8g,%8g) = %g\n' % (coor[i][0], coor[i][1], u_array[i]))


plot(Iv1, axes=True, title='Electron Density n(x)')
plot(Ie1, axes=True, title='Wave Function psi1(x)')
plot(Ie2, axes=True, title='Wave Function psi2(x)')
plot(Ie3, axes=True, title='Wave Function psi3(x)')

plot(new_potential, axes=True, title='V(x)')

plot(phi,  axes=True, title='phi(x)')


interactive()


