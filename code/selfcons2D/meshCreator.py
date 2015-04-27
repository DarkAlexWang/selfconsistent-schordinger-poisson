from dolfin import *
import sys, numpy, math

#**********************************************************

#defining the scale of substrates

gaas_thickness1 = 150.0
algaas_thickness1 = 200.0
algaas_thickness2 = 50.0
gaas_thickness2 = 500.0

total_width = gaas_thickness1 + algaas_thickness1  + algaas_thickness2 + gaas_thickness2

#MESH																																																																								
nx = 30;  ny = 30
x0, y0 = 0.0, 0.0
xf, yf = total_width, total_width
mesh = RectangleMesh(x0, y0, xf, yf, nx, ny)

#mass of electron in MeV/c^2
m0 = 0.511 #MeV/c^2


# Define a MeshFunction over two subdomains
subdomains = MeshFunction('size_t', mesh, 2)

class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[0], (x0, gaas_thickness1)))

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
         return (between(x[0], (gaas_thickness1, gaas_thickness1 + algaas_thickness1)))
         
class Omega2(SubDomain):
    def inside(self, x, on_boundary):
         return (between(x[0], (gaas_thickness1 + algaas_thickness1, gaas_thickness1 + algaas_thickness1 + algaas_thickness2)))
        

class Omega3(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[0], (gaas_thickness1 + algaas_thickness1 + algaas_thickness2, xf)))

# Mark subdomains with numbers 0 and 1
subdomain0 = Omega0()
subdomain0.mark(subdomains, 0)
subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)
subdomain2 = Omega2()
subdomain2.mark(subdomains, 2)
subdomain3 = Omega3()
subdomain3.mark(subdomains, 3)

V0 = FunctionSpace(mesh, 'DG', 0)
#k = Function(V0)
#nd = Function(V0)
#m_eff = Function(V0)

#plot(subdomains)

# Loop over all cell numbers, find corresponding
# subdomain number and fill cell value in k
#k_values = [0.0, 0.8*0.3, 0.8*0.3, 0.0]  # values of k in the two subdomains
#nd_values = [1e16, 1e16,0.0, 0.0]  # values of nd in the two subdomains
#m_eff_values = [0.063*m0, 0.0879*m0, 0.0879*m0, 0.0]
#for cell_no in range(len(subdomains.array())):
#    subdomain_no = subdomains.array()[cell_no]
#    k.vector()[cell_no] = k_values[subdomain_no]
#    nd.vector()[cell_no] = nd_values[subdomain_no]
#    m_eff.vector()[cell_no] = m_eff_values[subdomain_no]
    
#######################################
class Step(Expression):
   def __init__(self, mesh):
      self.mesh = mesh
   def eval_cell(self, value, x, ufc_cell):
      cell = Cell(self.mesh, ufc_cell.index)
      if (between(cell.midpoint().x(),  (gaas_thickness1, gaas_thickness1 + algaas_thickness1 + algaas_thickness2))):
        value[0] = 0.0
      else:
        value[0] = (1.0-0.3)*1.4 + (.3)*2.7
      
P = Step(mesh)
band_offset = interpolate(P, V0)

plot(band_offset, title='band_offset', axes=True, interactive=True)

#plot(k, title='dielectric', axes=True)
#plot(m_eff, title='effective mass', axes=True)
#plot(nd , title='doping', axes = True)

#plot(subdomains, title='subdomains')
#***************************************************************

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] < DOLFIN_EPS or x[0] > (xf - DOLFIN_EPS)) \
                    )

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[1] = x[1] - yf -DOLFIN_EPS
        y[0] = x[0]


def meshFunc(): 
  # Create mesh and finite element
  V = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())
  
  # Create Dirichlet boundary condition
  u0 = Constant(0.0)
  dbc = DirichletBoundary()
  bc0 = DirichletBC(V, u0, dbc)

  # Collect boundary conditions
  bcs = [bc0]
  
  #return functionspace and boundary conditions
  return [V0, bcs, 0.0, 0.0, 0.0, mesh, band_offset, V]
  
