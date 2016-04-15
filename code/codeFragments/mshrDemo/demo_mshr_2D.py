# Copyright (C) 2012-2013 Benjamin Kehlet
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2012-11-12
# Last changed: 2013-03-15
# Begin demo

from dolfin import *
from mshr import *

#if not has_cgal():
#    print "DOLFIN must be compiled with CGAL to run this demo."
#    exit(0)
dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain = Circle(dolfin.Point(1., 1.), 1.) - Rectangle(dolfin.Point(0.5, 0.5), dolfin.Point(1.5, 1.5)) 

#domain.set_subdomain(1, Circle(dolfin.Point(2., 2.), 0.5))
#domain.set_subdomain(2, Rectangle(dolfin.Point(0.5, 0.5), dolfin.Point(1.5, 1.5)))


#domain = Rectangle(dolfin.Point(-2., -2.), dolfin.Point(2., 2.))
#domain.set_subdomain(1, Rectangle(dolfin.Point(-2., -2.), dolfin.Point(2., 2.)))
#domain.set_subdomain(2, Circle(dolfin.Point(0, 0), 1))

# Test printing
#info("\nCompact output of 2D geometry:")
#info(domain)
#info("")
dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# Plot geometry
#dolfin.plot(domain, "2D Geometry (boundary)")

# Generate and plot mesh
mesh2d = generate_mesh(domain, 25)
print mesh2d
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")


dolfin.interactive()
