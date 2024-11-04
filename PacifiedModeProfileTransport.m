%% Initialize

clc
close all
clearvars -except modedata

addpath('.\Functions')
addpath('.\DataFiles')

if exist('modedata') == 0
    modedata = readtable("C:\Users\natef\OneDrive - University of Maryland\MATLAB\Floquet Research\DataFiles\poopcrosssectionmodes3d.txt");
end

n=[3.437208 3.408583 3.360361 3.291731 3.201465 3.087806 2.948290 2.779468 2.576498 2.332719 2.040777 1.713245];


%% Settings
ell=1.55/2/pi; %length unit

DwTotal=0.12    /ell; %um
dy=[0.1, 0.2, 0.3, 0.4, 0.6, 1, 2, 3, 4, 6, 10, 20, 30, 40, 60]/1000 /ell; %nm
if isinteger(sum(dy))
else
    fprintf('Step sizes are messed up')
end

%% Parameters
ncore=2.02;
nclad=1.46;

w=1.500/ell;
width=2*w;

ecore=ncore^2;
eclad=nclad^2;


corewidth=3.000/ell;
cladwidth=9.000/ell;

coreheight=0.400/ell;
cladheight=3.000/ell;

%% Preparations

N=round(cladwidth./dy);



pa=[ecore,eclad,ell,N,dy];

%% Running



