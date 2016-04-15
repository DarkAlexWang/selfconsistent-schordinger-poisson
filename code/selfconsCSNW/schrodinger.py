'''
2D Poisson-Schrodinger equation self-consistant solver for a core-shell 
hexagonal nanowire, the mesh is produced by mshr.
Created by Zhihuan Wang
Date 04/14/2016  Modified 04/14/2016

'''


from dolfin import *
from meshCreator import *
from scipy import constants as pc
from scipy import sqrt
import numpy as np
import math

#constants
hb2m_g = Constant(3.81/0.063) #ev*Angstros^2
hb2m_a = Constant(3.81/0.0879)
band_offset_value = 0.24 #ev

#create mesh
#nx = 20;  ny = 20
mesh = RectangleMesh(Point(x0, y0), Point(xf, yf), nx, ny)

def schrodingerEq(potential, meshArray, fermiEnergy):  
 
 # Define boundary condition
 # Create classes for defining parts of the boundaries and the interior
 # of the domain
 class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], x0)
        
 class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], xf)
        
 class GaAs1(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (x0, gaas_thickness1))
        
 class AlGaAs1(SubDomain):
     def inside(self, x, on_boundary):
         return (between(x[0], (gaas_thickness1, gaas_thickness1 + algaas_thickness1)))
                
 class AlGaAs2(SubDomain):
     def inside(self, x, on_boundary):
         return (between(x[0], (gaas_thickness1 + algaas_thickness1, gaas_thickness1 + algaas_thickness1 + algaas_thickness2)))
         
 class GaAs2(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (gaas_thickness1 + algaas_thickness1 + algaas_thickness2, xf))
 
 # Sub domain for Periodic boundary condition
 class PeriodicBoundary(SubDomain):
 
   # Left boundary is "target domain" G
   def inside(self, x, on_boundary):
       return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS)
    # Map right boundary (H) to left boundary (G)
   def map(self, x, y):
         y[1] = x[1] - yf -DOLFIN_EPS
         y[0] = x[0]
 
 class BoundaryPoint(SubDomain):
   def __init__(self, bpt):
       SubDomain.__init__(self)
       self.bpt = bpt
   def inside(self,x,on_boundary):
       return near(x[0],self.bpt)
       
 # Initialize sub-domain instances
 top = Top()
 bottom = Bottom()
 gaas1 = GaAs1()
 algaas1 = AlGaAs1()
 algaas2 = AlGaAs2()
 gaas2 = GaAs2()
 
  
 # Initialize mesh function for interior domains
 domains = CellFunction("size_t", mesh)
 domains.set_all(0)
 gaas1.mark(domains, 1)
 algaas1.mark(domains, 2)
 algaas2.mark(domains, 3)
 gaas2.mark(domains, 4)
 
 # Initialize mesh function for boundary domains
 boundaries = FacetFunction("size_t", mesh)
 boundaries.set_all(0)
 #top.mark(boundaries, 1)
 #bottom.mark(boundaries, 2)
 
  
 #Define Function space
 V =FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())
 
 
 #define functions
 u = TrialFunction(V)
 v = TestFunction(V)
 Vpot = Function(V)
 
 Vpot.vector()[:] = np.array([i for i in potential.vector()])
 
 #plot(Vpot, interactive=True)
 
 #boundary conditions
 start = BoundaryPoint(0.0)
 left = BoundaryPoint(gaas_thickness1)
 right = BoundaryPoint(gaas_thickness1 + algaas_thickness1 + algaas_thickness2)
 end = BoundaryPoint(xf)
 start.mark(boundaries,0)
 left.mark(boundaries, 1)
 right.mark(boundaries,2)
 end.mark(boundaries, 3)
 bcs = [DirichletBC(V, 0.0, start), DirichletBC(V, 0.0, right)]#, DirichletBC(V, 0.0, left)]#, DirichletBC(V, 0.0, end)]
 #bcs.apply(A)
 #bcs.apply(M)
 
 # Define new measures associated with the interior domains and
 # exterior boundaries
 dx = Measure("dx")[domains]
 #ds = Measure("ds")[boundaries]
  
 #define problem
 a = (inner(hb2m_g * grad(u), grad(v)) \
      +  Vpot*u*v)*dx(1) + (inner(hb2m_a * grad(u), grad(v)) \
      +  Vpot*u*v)*dx(2) + (inner(hb2m_a * grad(u), grad(v)) \
      +  Vpot*u*v)*dx(3) + (inner(hb2m_g * grad(u), grad(v)) \
      +  Vpot*u*v)*dx(4) + (inner(hb2m_g * grad(u), grad(v)) \
      +  Vpot*u*v)*dx(0)
 m = u*v*dx(1) + u*v*dx(2) + u*v*dx(3) + u*v*dx(4) + u*v*dx(0)
 
 #assemble stiffness matrix
 A = PETScMatrix() 
 M = PETScMatrix()
 _ = PETScVector()
 L = Constant(0.)*v*dx(0) + Constant(0.)*v*dx(1) + Constant(0.)*v*dx(2) + Constant(0.)*v*dx(3) + Constant(0.)*v*dx(4)
 #assemble(a, tensor=A)
 #assemble(m, tensor=M)
 assemble_system(a, L, A_tensor=A, b_tensor=_)
 #assemble_system(m, L,bc, A_tensor=M, b_tensor=_)
 assemble_system(m, L, A_tensor=M, b_tensor=_)
 
 #create eigensolver
 eigensolver = SLEPcEigenSolver(A,M)
 eigensolver.parameters['spectrum'] = 'smallest magnitude'
 eigensolver.parameters['solver']   = 'lapack'
 eigensolver.parameters['tolerance'] = 1.e-15
 
 #solve for eigenvalues
 print 'solving Schrodinger\'s equation...'
 eigensolver.solve()
 
 
 #assign eigenvector to function
 u = Function(V)
 eigenvectors = []
 eigenvalues = []
 rx_list = []

 
 #extract first eigenpair
 r, c, rx, cx = eigensolver.get_eigenpair(0)
 
 #assign eigenvector to function
 #rx = np.array([q*xsi_factor for q in rx])
 #u.vector()[:] = rx
 #eigenvalues.append(r)
 #eigenvectors.append([j for j in u.vector()])
     
 #rx_list.append(rx) 
 #plot eigenfunction and probability density
 #print 'eigenvalue: ' + str(r)
 #plt  = plot(u, axes=True)
 #interactive()
  
 
 i = 0
 while (i < eigensolver.get_number_converged() and r < 1.8):#fermiEnergy):
     #extract next eigenpair
     r, c, rx, cx = eigensolver.get_eigenpair(i)
     #assign eigenvector to function
     #rx = np.array([q*xsi_factor for q in rx])
     #if(r > 0.24):
     #  break
     
     u.vector()[:] = rx
     eigenvalues.append(r)
     eigenvectors.append([j for j in u.vector()])
     
     rx_list.append(rx) 
     
     #plot eigenfunction
     print 'eigenvalue: ' + str(r)
     #plot(u).update(u)
     #interactive()
  
     #increment i
     i = i+1
 return [eigenvalues, eigenvectors,u,rx_list]
