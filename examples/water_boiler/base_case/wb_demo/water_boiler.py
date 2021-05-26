#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 15:14:10 2021

@author: sblair

A module for generating an OpenMC model water boiler
"""

import openmc
import numpy as np
import scipy.interpolate

def get_UranylSulfate_solution_density(pct_us=0.299,temp=303):
    '''
    input:
    pct_us: float, weight percent UO2-SO4 in water
    temp: float, temperature in K
    
    output:
    rho_us: float, mass density (g/cc) of the solution
    
    This is derived from Table 5 of IEU-SOL-THERM-004
    
    '''
    temp = temp-273.; # convert to C
    y = [0.51, 0.399, 0.299];
    x = [15.,20.,25.,30.,35.,40.,50.];
    
    z = [[1.7959,1.7926,1.7891,1.7854,1.7816,1.777,1.7696],
        [1.5283,1.5257,1.5228,1.5199,1.5166,1.5133,1.5064],
        [1.3499,1.3479,1.3455,1.3430,1.3403,1.3376,1.3314]];
    f = scipy.interpolate.interp2d(x,y,z);
    rho_us = np.float64(f(temp,pct_us));
    return rho_us

# the atomic weight of enriched uranium is needed
def get_Aw_U(enrichment):
    assert (enrichment<=1.0), "enrichment should be entered as a percentage"
    # assumes a fixed ratio in the percentage of U234 and U235.
    # i.e. the enrichment process keeps these two isotopes in the
    # same relative abundance
    
    U234_to_U235_ratio = 0.0055/0.72;
    Aw_U235 = openmc.data.atomic_mass('U235');
    Aw_U234 = openmc.data.atomic_mass('U234');
    Aw_U238 = openmc.data.atomic_mass('U238');
    
    frac_235 = enrichment;
    frac_234 = frac_235*U234_to_U235_ratio;
    frac_238 = 1. - frac_235 - frac_234;
    
    aw = 1./(frac_235/Aw_U235 + frac_238/Aw_U238 + frac_234/Aw_U234);
    
    weight_frac = {};
    weight_frac['U234']=frac_234;
    weight_frac['U235']=frac_235;
    weight_frac['U238']=frac_238;
    
    return aw,weight_frac
    

# calculate atom densities of UO2-SO4 + H2O solution as a function
# of Uranium enrichment, water temperature and Urynal Sulfate concentration

def BoilerAtomDensities(enrich=0.1467,temp=303,conc=0.299):
    ''' 
    input:
    enrich: w/o U-235
    temp: solution temperature in K;
    conc: w/o concentration of Uranyl Sulfate in the water
    
    output:
    dictionary with the atom densities (atoms/b-cm) of all elements
    and nuclides in the solution
    
    this results in a 'not great, not terrible' agreement with the benchmark.
    must re-visit to make corrections
    
    '''
    assert (temp > 288) and (temp < 323), 'temperature not in correlation limits';
    assert (conc >= 0.299) and (conc <= 0.51), 'solution concentration not in correlation limits';
        
    
    Na = 0.60221; # Avagadro's number x10**-24
    AtomDensities = {};
    
    rho_uranyl_sulf_sol = get_UranylSulfate_solution_density(conc,temp);
    rho_water = rho_uranyl_sulf_sol*(1-conc);
    rho_uranyl = rho_uranyl_sulf_sol*conc;
    
    #Aw_234 = openmc.data.atomic_mass('U234');
    #Aw_235 = openmc.data.atomic_mass('U235');
    #Aw_238 = openmc.data.atomic_mass('U238');
    Aw_S = openmc.data.atomic_weight('S');
    Aw_H = openmc.data.atomic_weight('H');
    Aw_O = openmc.data.atomic_weight('O');
    Aw_U,U_weight_fracs = get_Aw_U(enrich);
    
    Aw_uranyl_sulf = Aw_U + 2.*Aw_O + Aw_S + 4.*Aw_O;
    Aw_h2o = 2.*Aw_H + Aw_O;
    
    mol_density_uranyl_sulf = (rho_uranyl/Aw_uranyl_sulf)*Na;
    mol_density_h2o = (rho_water/Aw_h2o)*Na;
    AtomDensities['H']=mol_density_h2o*2.;
    AtomDensities['O']=mol_density_h2o*1.;
    AtomDensities['O']+=mol_density_uranyl_sulf*4.;
    AtomDensities['S']=mol_density_uranyl_sulf*1.;
    AtomDensities['U234']=(mol_density_uranyl_sulf)*U_weight_fracs['U234']; #wrong but better
    AtomDensities['U235']=(mol_density_uranyl_sulf)*U_weight_fracs['U235'];
    AtomDensities['U238']=(mol_density_uranyl_sulf)*U_weight_fracs['U238'];
    
    
    return AtomDensities

def generate_model(sol_temp=303.,sol_conc=0.299,U_enrch=0.1467,cr_wd=0.1):
    sol_atom_densities = BoilerAtomDensities(enrich=U_enrch,temp=sol_temp,conc=sol_conc);

    sol = openmc.Material(name='sol');
    sol.add_element('H',sol_atom_densities['H'],percent_type='ao');
    sol.add_element('O',sol_atom_densities['O'],percent_type='ao');
    sol.add_element('S',sol_atom_densities['S'],percent_type='ao');
    sol.add_nuclide('U234',sol_atom_densities['U234'],percent_type='ao');
    sol.add_nuclide('U235',sol_atom_densities['U235'],percent_type='ao');
    sol.add_nuclide('U238',sol_atom_densities['U238'],percent_type='ao');
    sol.add_s_alpha_beta('c_H_in_H2O');

    ad_tot = 0.;
    for key in sol_atom_densities:
        ad_tot+=sol_atom_densities[key];
    
    sol.set_density('atom/b-cm',ad_tot);
    
    brass = openmc.Material(name='brass');
    brass.add_element('Fe',0.001002);
    brass.add_element('Cu',0.674918);
    brass.add_element('Zn',0.320956);
    brass.add_element('Sn',0.001451);
    brass.add_element('Pb',0.001673);
    brass.set_density('g/cc',8.070);
    
    cadmium = openmc.Material(name='cadmium');
    cadmium.add_element('Cd',1.0);
    cadmium.set_density('g/cc',8.65);
    
    shell = openmc.Material(name='shell');
    shell.add_element('C',0.003659);
    shell.add_element('Si',0.019559);
    shell.add_element('P',0.000798);
    shell.add_element('S',0.000514);
    shell.add_element('Cr',0.179602);
    shell.add_element('Mn',0.019998);
    shell.add_element('Fe',0.669338);
    shell.add_element('Ni',0.102952);
    shell.add_element('Nb',0.002365);
    shell.add_element('Ta',0.001214);
    shell.set_density('g/cc',8.0);
    
    beryl_ref = openmc.Material(name='beryl_ref');
    beryl_ref.add_element('O',6.6210e-2);
    beryl_ref.add_element('Be',6.6210e-2);
    beryl_ref.add_element('B',3.0637e-7);
    beryl_ref.add_element('Co',5.6202e-7);
    beryl_ref.add_element('Ag',3.0706e-8);
    beryl_ref.add_element('Cd',7.3662e-8);
    beryl_ref.add_element('In',1.4423e-8);
    beryl_ref.add_s_alpha_beta('c_Be_in_BeO');
    beryl_ref.set_density('g/cc',2.75);
    
    
    grph = openmc.Material(name='grph');
    grph.add_element('C',0.999999);
    grph.add_element('B',0.000001);
    grph.set_density('g/cc',1.7);
    grph.add_s_alpha_beta('c_Graphite');
    
    
    air = openmc.Material(name='air');
    air.add_element('C',0.000150);
    air.add_element('N',0.784431);
    air.add_element('O',0.210748);
    air.add_element('Ar',0.004671);
    air.set_density('g/cc',0.001205);
    
    materials = openmc.Materials([sol,shell,beryl_ref,grph,air,brass,cadmium]);
    
   
    rx_origin = [0.,76.3214,0.];
    ref_sphere = openmc.Sphere(y0=rx_origin[1],r=47.4210);
    tank_o = openmc.Sphere(y0=rx_origin[1],r=15.3614);
    tank_i = openmc.Sphere(y0=rx_origin[1],r=15.282);
    graph_base_cyl = openmc.YCylinder(r=47.4210);
    #fill_drain_cav = openmc.YCylinder(r=4.445/2.);
    fill_drain_o = openmc.YCylinder(r=2.06375);
    fill_drain_i = openmc.YCylinder(r=1.905);
    plate_plane = openmc.YPlane(y0=0.);
    base_plane = openmc.YPlane(y0=34.4114);
    sphere_center_plane = openmc.YPlane(y0=rx_origin[1]);
    upper_plane = openmc.YPlane(y0=118.2314);
    bbox = openmc.model.RightCircularCylinder([0.,-10.,0.],230.,60.,axis='y',boundary_type='vacuum');
    
    # surfaces for the control rod
    rod_channel_bottom = 118.2314-76.20;
    
    cr_cyl = openmc.YCylinder(x0=-18.891,z0=0.,r=1.42875);
    cr_cyl_bottom = openmc.YPlane(y0=rod_channel_bottom);
    
    # surfaces for the safety rod
    sr_right = openmc.XPlane(x0=17.3664);
    sr_left = openmc.XPlane(x0=15.4614);
    sr_front = openmc.ZPlane(z0=7.62/2.);
    sr_back = openmc.ZPlane(z0=-7.62/2.);
    
    # top plane for cr/sr
    rod_channel_top = openmc.YPlane(y0=200.);
    rod_length = 76.20;
    
    #sr_wd = 76.20;# cm, distance from fully inserted
    
    
    cr_bottom = openmc.YPlane(y0=(rod_channel_bottom+cr_wd));
    cr_top = openmc.YPlane(y0=(rod_channel_bottom+cr_wd+rod_length));
    
    #sr_bottom = openmc.YPlane(y0=(rod_channel_bottom+sr_wd));
    #sr_top = openmc.YPlane(y0=(rod_channel_bottom+sr_wd+rod_length));
    
    cr_brass_o = openmc.YCylinder(x0=-18.891,z0=0.,r=0.9525);
    cr_brass_i = openmc.YCylinder(x0=-18.891,z0=0.,r=0.7000);
    cr_cd_o = openmc.YCylinder(x0=-18.891,z0=0.0,r=1.03375);
    
    core = openmc.Cell();
    core.fill = sol;
    core.region = (-tank_i) | (-fill_drain_i) & -bbox
    
    steel_tank_and_pipe = openmc.Cell();
    steel_tank_and_pipe.fill = shell;
    #steel_tank_and_pipe.region = (+tank_i & -tank_o & ~(-fill_drain_i)) | \
    #                             (+fill_drain_i & -fill_drain_o & +tank_i) & -bbox
    steel_tank_and_pipe.region = (+tank_i & -tank_o & ~(-fill_drain_i)) | \
                                 (+fill_drain_i & -fill_drain_o & +tank_i) & -bbox;

    # make a universe for the control rod
    cr_brass = openmc.Cell();
    cr_brass.fill = brass;
    cr_brass.region = +cr_brass_i & -cr_brass_o & +cr_bottom & -cr_top;
    
    cr_cd = openmc.Cell();
    cr_cd.fill = cadmium;
    cr_cd.region = +cr_brass_o & -cr_cd_o & +cr_bottom & -cr_top;
    
    cr_air = openmc.Cell();
    cr_air.fill = air;
    cr_air.region = +cr_cyl_bottom & -rod_channel_top & -cr_cyl & ~cr_cd.region & ~cr_brass.region
    
    cr_univ = openmc.Universe();
    cr_univ.add_cells([cr_brass,cr_cd,cr_air]);
    
    cr = openmc.Cell();
    cr.fill = cr_univ;
    cr.region = -cr_cyl & +cr_cyl_bottom & -rod_channel_top;
    
    sr = openmc.Cell();
    sr.fill = air;
    sr.region = +sr_left & -sr_right & +cr_cyl_bottom & -upper_plane & -sr_front & +sr_back;
    
    ref = openmc.Cell();
    ref.fill = beryl_ref;
    ref.region = ((+tank_o & +fill_drain_o) & -ref_sphere & +base_plane & -upper_plane) & ~cr.region & ~sr.region; 
    
    outside = openmc.Cell();
    outside.fill = air;
    outside.region = -bbox & (+graph_base_cyl | (+ref_sphere & -upper_plane) | 
                             (+upper_plane & +fill_drain_o & ~cr.region ))
    
    graph_base = openmc.Cell();
    graph_base.fill = grph;
    graph_base.region = ((-graph_base_cyl & +plate_plane & -base_plane & +fill_drain_o) |
                         (-graph_base_cyl & +ref_sphere & +base_plane & -sphere_center_plane));
    root = openmc.Universe();
    root.add_cells([graph_base,ref,core,steel_tank_and_pipe,outside,cr,sr]);
    
    geometry = openmc.Geometry();
    geometry.root_universe = root;
    
    cell_filter = openmc.CellFilter(core);
    N = 1001;
    energy_bins = np.logspace(-3,7,num=N);
    energy_filter = openmc.EnergyFilter(values=energy_bins);
    
    abs_core = openmc.Tally(name='abs_core');
    abs_core.scores = ['absorption'];
    abs_core.filters = [cell_filter,energy_filter];
    
    fission = openmc.Tally(name='fission');
    fission.scores = ['fission'];
    fission.filters = [cell_filter,energy_filter];
    
    fission_by_nuclide = openmc.Tally(name='fission_by_nuclide');
    fission_by_nuclide.scores = ['fission'];
    fission_by_nuclide.nuclides = ['U234','U235','U238'];
    fission_by_nuclide.filters = [cell_filter,energy_filter];
    
    capture = openmc.Tally(name='capture');
    capture.scores = ['(n,gamma)'];
    capture.filters = [cell_filter,energy_filter];
    
    capture_by_nuclide = openmc.Tally(name='capture_by_nuclide');
    capture_by_nuclide.scores = ['(n,gamma)'];
    capture_by_nuclide.nuclides = ['U234','U238','H1','O16','S32'];
    capture_by_nuclide.filters = [cell_filter,energy_filter];
    
    
    flux = openmc.Tally(name='flux');
    flux.scores = ['flux'];
    flux.filters = [cell_filter,energy_filter];
    
    tallies = openmc.Tallies([abs_core, flux,
                              fission, capture,
                              fission_by_nuclide,
                              capture_by_nuclide]);
    
    settings = openmc.Settings();
    settings.batches = 100;
    settings.inactive = 20;
    settings.particles = 5000;
    
    R = 15.;
    y_org = 76.3214;
    
    bounds = [-R,-R+y_org,-R,R,R+y_org,R];
    uniform_dist = openmc.stats.Box(bounds[:3],bounds[3:],
                                       only_fissionable=True);
    settings.source = openmc.source.Source(space=uniform_dist);
    
    #settings.temperature['method']='interpolation';
    
    return materials, geometry, tallies, settings