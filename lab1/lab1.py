#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:06:33 2021

@author: sblair

This lab exercise is intended merely to introduce the user to the basic
elements of an input file, the procedures necessary to invoke the code
and an orientation to the contents of an output file.  (all for MCNP)

The geometry is a sphere of radius 10 cm centered at the origin. The material
everywhere in the universe is void.  The sphere is contained within a larger
RPP 100 cm on an edge.  Importance inside the RPP (and inside the sphere) is set
to 1 and outside the RPP the importance is set to 0.

100,000 neutron histories are run and an f2 tally is used to get the neutron
fluence for the sphere.  MCNP output is confirmed to be equal to the theoretical
value of 1/(surface are of sphere) per source particle.
"""

import openmc

# MATERIALS

# this lab has no materials.  Let's just try to write empty materials
# to file
mf = openmc.Materials(())
mf.export_to_xml()

# GEOMETRY

# define the surfaces
sphere = openmc.Sphere(); #note length units = cm
sphere.r = 1.0;

#evidently OpenMC does not have RPP so we'll do this the hard way
box_xm = openmc.XPlane(x0=-50.0,boundary_type='vacuum');
box_xp = openmc.XPlane(x0=50.0,boundary_type='vacuum');
box_ym = openmc.YPlane(y0=-50.0,boundary_type='vacuum');
box_yp = openmc.YPlane(y0=50.0,boundary_type='vacuum');
box_zm = openmc.ZPlane(z0=-50.0,boundary_type='vacuum');
box_zp = openmc.ZPlane(z0=50.0,boundary_type='vacuum');

# define the cells
sphere_region = -sphere;
outer_region = +box_xm & -box_xp & +box_ym & -box_yp & +box_zm & -box_zp & +sphere;
#cell3 = ~cell2 & ~cell1;#<-- do I need this?

cell1 = openmc.Cell();
cell1.fill = 'void';
cell1.region = sphere_region;

cell2 = openmc.Cell()
cell2.fill = 'void';
cell2.region = outer_region;

root = openmc.Universe()
root.add_cells((cell1,cell2))

g = openmc.Geometry()
g.root_universe = root
g.export_to_xml()

# data/control section
settings = openmc.Settings();
settings.run_mode = 'fixed source';

settings.batches = 50;# even for fixed source problems you need batches
settings.particles = 10000; 

source = openmc.Source();
source.particle = 'neutron';
source.space = openmc.stats.Point(xyz=(0.,0.,0.)); # default is the origin
source.angle = openmc.stats.Isotropic();
source.energy = openmc.stats.Discrete([1.0e6],[1.0]);
settings.source = source;

settings.export_to_xml()

t = openmc.Tally(name='sphere');
s_filter = openmc.SurfaceFilter(sphere.id);
t.filters = [s_filter];
t.scores = ['flux'];

tallies = openmc.Tallies([t])
tallies.export_to_xml()


openmc.run()
