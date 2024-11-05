%% Initialize

clc
close all
clearvars -except rawdata

addpath('.\Functions')
addpath('.\DataFiles')

if exist('modedata','var') == 0
    rawdata = table2array(readtable("C:\Users\natef\OneDrive - University of Maryland\MATLAB\Floquet Research\DataFiles\poorcrosssectionmodes3d.txt"));
    n0=[3.437208 3.408583 3.360361 3.291731 3.201465 3.087806 2.948290 2.779468 2.576498 2.332719 2.040777 1.713245];

end



%% Settings
ell=1.55/2/pi; %length unit

dwTotal=0.100    /ell; %um

divs=(1:0.25:4); %fractions of total change to compute by

%% Parameters
ncore=2.03-0.1*1i;
nclad=1.46-0.05*1i;

corewidth=3.000     /ell;
cladwidth=9.000     /ell;

coreheight=0.400    /ell;
cladheight=3.000    /ell;

%% Preparations

ecore=ncore^2;
eclad=nclad^2;

w0=corewidth/2;
dw=dwTotal./10.^(divs);

i=1;
while rawdata(i,1)==rawdata(1,1)
    i=i+1;
end
Nx=i-1;
Nxy=size(rawdata,1);

if rem(Nxy,Nx)==0
    Ny=Nxy/Nx;
else
    fprintf('Ny \n');
    scream()
    return
end

nmodes=round((size(rawdata(1,:),2)-2)/3,1);
Nxy=Nx*Ny;

xVals=rawdata(:,1);
yVals=rawdata(:,2);

E0=zeros(Nxy,3,nmodes);
for i=1:nmodes
    E0(:,1,i)=rawdata(:,3*i);
    E0(:,2,i)=rawdata(:,3*i+1);
    E0(:,3,i)=rawdata(:,3*i+2);
end

coref=(heaviside(xVals+corewidth/2)-heaviside(xVals-corewidth/2)).*(heaviside(yVals+cladwidth/2)-heaviside(yVals-cladwidth/2));
dA=(xVals(Nx+1)-xVals(1))*(yVals(2)-yVals(1));

H0=zeros(nmodes,nmodes);
for i=1:nmodes
    H0(i,i)=sum((E0(:,:,i)')*(coref.*E0(:,:,i)),"all")*dA;
    for j=i+1:nmodes
        H0(i,j)=sum((E0(:,:,i)')*(coref.*E0(:,:,j)),"all")*dA;
        H0(j,i)=H0(i,j);
    end
end


S0=zeros(nmodes,nmodes);


imagesc(1:nmodes,1:nmodes,(abs(H0)))
colorbar
colormap('jet')
axis xy


%% Running

for d=divs
    modalbasisprofiles=zeros(nmodes,nmodes,d);
    eigs=zeros(nmodes,d);

    modalbasisprofiles(:,:,1)=eye(nmodes);
    eigs(:,1)=n0;
    
    cw=0;
    for w=w0:dWtotal/d:(w0+dwTotal)
        cw=cw+1;


    end






end


