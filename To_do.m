% To-do list comment out, when finished
1.1 Find different scientific papers describing litosphere structure
1.2 Describe them and formulate a hypothesis

2 Import data and compare them to literature
%2.1 Import Topographic data
%2.2 Import Gravity Data
2.3 Construct Bouguer gravity map

3 Perform Bouguer Inversion
3.1 Choose thickness D (different values, compare results)
3.2 Change density paramters
3.3 perform inversion

4 Airy model (thickness D, one formular)

5 Flexure model of crust
5.1 Perform spectrum analysis of Airy crust to Spherical Harmonics (adapt Test_GSHloop.m)
5.2 Choose appropiate boundaries and youngs-modulus and elastic thickness
5.3 Multiply SH domain with flexural response

6 Gravity response
6.1 Use GSH code: run_example_2_layer_planet
6.2 Analyze and Plot the degree variance responds of those signals with gravity observations
6.3 Discuss plots
6.4 compare all three models with literature and to each other

7 Find optimal Te value
7.1 find spectral domain that is sensitive to litosphere flexure
7.2 fitt the degree variance of flexural model with that of the observed gravity field
7.3 find optimal Te value
7.4 discuss optimal fitting crustal model and compare with literature -> tectonic state 

8 Write everything
8.1 write to Data import
8.2 Write to Bouger Model
8.3 Airy Model
8.4 Flexure Model
8.5 discussion


Use these commands so load and save data properly
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);

filename = [HOME '/Data/ggmes_100v08_sha.tab'];
load([HOME '/Results/r4_delta_g_mGal.mat'], 'delta_g_mGal')
save([HOME '/Results/crust_thickness_iterative.mat'], 'crust_thickness');