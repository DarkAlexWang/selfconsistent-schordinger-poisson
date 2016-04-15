'''
2D Poisson-Schrodinger equation self-consistant solver for a core-shell 
hexagonal nanowire, the mesh is produced by mshr.
Created by Zhihuan Wang
Date 04/14/2016  Modified 04/14/2016

'''

mport mpmath
from scipy import inf
import numpy as np
import math
from scipy import constants as pc

#constants
hb2m = 3.81 #ev*A^2
nk_factor =  1.0/pc.pi * 1.0/2/(hb2m)#pc.m_e/(pc.pi * pc.hbar * pc.hbar * 1e40)
K = 8.6173324e-5 #ev K^-1
T = 4.0 # temperature


#function to return the electron occupation state nk
def electronOccupationState(eigenvalues, fermiEnergy):
 print 'Finding n_k...'
 nk = []
 for i in range(0,len(eigenvalues)):
  ek = float(eigenvalues[i])	
  #result = nk_factor *( limit(x-ln(1+exp(x-fermiEnergy)/pc.k/300),x,oo) - (ek - ln(1+exp(ek-fermiEnergy)/pc.k/300)))
  result = mpmath.quad(lambda x: nk_factor/(1+mpmath.exp((x-fermiEnergy)/T/K)), [ek, 2*ek + fermiEnergy])
  print float(result)
  nk.append(float(result))

 return nk

#fuction to calculate the electron density
def electronDensityFunction(eigenvectors,  nk):
 print 'Finding the electron density function n(x)'
 result = []
 for i in range(len(eigenvectors)):
   kth_term = [(j * j * nk[i]) for j in eigenvectors[i]]
   result.append(kth_term)
   
 n = np.sum(result, axis=0)
 return n
