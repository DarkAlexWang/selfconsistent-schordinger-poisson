"""
2D Poisson-Schrodinger equation self-consistant solver for a core-shell 
hexagonal nanowire, the mesh is produced by mshr.
Original Creator by: Feras Aldahlawi
Contact: faldahlawi@gmail.com
Created by Zhihuan Wang
Date 04/14/2016  Modified 04/14/2016
"""

from dolfin import *
import sys, numpy, math
from mshr import *

#**********************************************************
# Define geometry constants
radius = 110 # total radius of ring
iradius = 70  # inner radius of ring

# Define 2D geometry
# Define the outter shell of nanowire
middle = Rectangle(dolfin.Point(-1/sqrt(3)*radius, -radius), \
         dolfin.Point(1/sqrt(3)*radius, radius))
left = Rectangle(dolfin.Point(-1/sqrt(3)*radius, -radius), \
       dolfin.Point(1/sqrt(3)*radius, radius))
right = Rectangle(dolfin.Point(-1/sqrt(3)*radius, -radius), \
        dolfin.Point(1/sqrt(3)*radius, radius))

left = CSGRotation(left, dolfin.Point(0, 0), -pi/3)
right = CSGRotation(right, dolfin.Point(0, 0), pi/3)

domain = middle + left + right

domain = Rectangle(dolfin.Point(-1/sqrt(3)*radius, -1*radius), \
        dolfin.Point(1/sqrt(3)*radius, 1*radius))

#domain.set_subdomain(1, Rectangle(dolfin.Point(20, 20), dolfin.Point(50, 50)))

domain.set_subdomain(1, middle+left+right)

# Define the inner core of nanowire
imiddle = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
         dolfin.Point(1/sqrt(3)*iradius, iradius))
ileft = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
       dolfin.Point(1/sqrt(3)*iradius, iradius))
iright = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
        dolfin.Point(1/sqrt(3)*iradius, iradius))

ileft = CSGRotation(ileft, dolfin.Point(0, 0), -pi/3)
iright = CSGRotation(iright, dolfin.Point(0, 0), pi/3)

innerdomain = imiddle + ileft + iright

domain.set_subdomain(2, innerdomain)

dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# Generate and plot mesh
mesh = generate_mesh(domain, 45)
print mesh
dolfin.plot(mesh, "2D mesh")

# Convert mf to mesh function for plotting
# Define a MeshFunction over two mf
mf = dolfin.MeshFunction("size_t", mesh, 2, mesh.domains())
dolfin.plot(mf, "Subdomains")

dolfin.interactive()

#defining the scale of substrates
"""
gaas_thickness1 = 150.0
algaas_thickness1 = 200.0
algaas_thickness2 = 50.0
gaas_thickness2 = 500.0

total_width = gaas_thickness1 + algaas_thickness1  + algaas_thickness2 + gaas_thickness2
"""
#MESH

nx = 30;  ny = 30
x0, y0 = -1/sqrt(3)*radius, -1*radius
xf, yf = 1/sqrt(3)*radius, 1*radius
centerx, centery = 0.0, 0.0
hori, vert = iradius, 2/sqrt(3)*iradius
#mass of electron in MeV/c^2
m0 = 0.511 #MeV/c^2

class Hexagon(SubDomain):
    def isinside(self, x, y):
        """Test a point is inside the hexagon or not"""
        
        # transform the test point locally and to quadrant 2
        q2x = math.fabs(x - centerx)
        q2y = math.fabs(y - centery)
        # bounding test ( since q2 is in quadrant 2 only, two tests are needed)
        if q2x > hori or q2y > vert * 2:
            return false
        # finally the dot product can be reduced to this due to the hexagon
        # symmetry
        return 2*vert * hori - vert * q2x - hori * q2y >= 0
"""
class Hexagon(SubDomain):
    def inside(self, x, on_boundary):
        if 
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

# Mark mf with numbers 0 and 1
subdomain0 = Omega0()
subdomain0.mark(mf, 0)
subdomain1 = Omega1()
subdomain1.mark(mf, 1)
subdomain2 = Omega2()
subdomain2.mark(mf, 2)
subdomain3 = Omega3()
subdomain3.mark(mf, 3)
"""
subdomain1 = Hexagon()
V0 = FunctionSpace(mesh, 'DG', 0)
#k = Function(V0)
#nd = Function(V0)
#m_eff = Function(V0)

#plot(mf)

# Loop over all cell numbers, find corresponding
# subdomain number and fill cell value in k
#k_values = [0.0, 0.8*0.3, 0.8*0.3, 0.0]  # values of k in the two mf
#nd_values = [1e16, 1e16,0.0, 0.0]  # values of nd in the two mf
#m_eff_values = [0.063*m0, 0.0879*m0, 0.0879*m0, 0.0]
#for cell_no in range(len(mf.array())):
#    subdomain_no = mf.array()[cell_no]
#    k.vector()[cell_no] = k_values[subdomain_no]
#    nd.vector()[cell_no] = nd_values[subdomain_no]
#    m_eff.vector()[cell_no] = m_eff_values[subdomain_no]
    
#######################################
class Step(Expression):
   def __init__(self, mesh):
      self.mesh = mesh
   def eval_cell(self, value, x, ufc_cell):
      cell = Cell(self.mesh, ufc_cell.index)
      if subdomain1.isinside(cell.midpoint().x(), cell.midpoint().y()):
        value[0] = 1.4
      else:
        value[0] = (1.0-0.3)*1.4 + (.3)*2.7
      
P = Step(mesh)
band_offset = interpolate(P, V0)

plot(band_offset, title='band_offset', axes=True, interactive=True)

#plot(k, title='dielectric', axes=True)
#plot(m_eff, title='effective mass', axes=True)
#plot(nd , title='doping', axes = True)

#plot(mf, title='mf')
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
  
