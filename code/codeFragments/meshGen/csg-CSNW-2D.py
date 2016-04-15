# Copyright (C) 2012-2016 Based on mshr demo of Benjamin Kehlet
# 
'''
Create a core-shell hexagonal nanowire mesh based on mshr
Produced by Zhihuan Wang
Date 04/14/2016  Modified 04/14/2016

'''
#
# mshr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# mshr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with mshr.  If not, see <http://www.gnu.org/licenses/>.

import dolfin
from mshr import *
from math import * 

dolfin.set_log_level(dolfin.TRACE)

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

#domain.set_subdomain(1, Rectangle(dolfin.Point(20, 20), dolfin.Point(50, 50)))


# Define the inner core of nanowire
imiddle = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
         dolfin.Point(1/sqrt(3)*iradius, iradius))
ileft = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
       dolfin.Point(1/sqrt(3)*iradius, iradius))
iright = Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), \
        dolfin.Point(1/sqrt(3)*iradius, iradius))

ileft = CSGRotation(ileft, dolfin.Point(0, 0), -pi/3)
iright = CSGRotation(iright, dolfin.Point(0, 0), pi/3)

#innerdomain = imiddle + ileft + iright
#domain.set_subdomain(2, Rectangle(dolfin.Point(-1/sqrt(3)*iradius, -iradius), dolfin.Point(1/sqrt(3)*iradius, iradius)))

# Code addaped from the demo 
#domain.set_subdomain(1, Rectangle(dolfin.Point(1., 1.), dolfin.Point(5., 3.)))
#domain.set_subdomain(2, Rectangle(dolfin.Point(2., 2.), dolfin.Point(3., 4.)))

dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
print mesh2d
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")

dolfin.interactive()
